import sys
import sqlite3
import requests
import pandas as pd
import json
import time
import os
from sqlalchemy import create_engine

__dataprepdir__ = os.path.dirname(os.path.abspath(__file__))
outfilename = os.path.join(__dataprepdir__, 'genes.json')
outdbname = 'postgap.20180324.db'

VALID_CHROMOSOMES = [*[str(chr) for chr in range(1, 23)], 'X', 'Y']


def build_db(filename):
    '''
    Build a sqlite database from a tsv file that meets the
    Open Targets format requirements.
    '''
    # creates db if db does not exist
    # conn = sqlite3.connect('postgap.db.scored')
    # conn = sqlite3.connect('postgap.20180324.v0.0.1.db')
    conn = sqlite3.connect(outdbname)
    cursor = conn.cursor()

    # # build raw table
    # build_raw(cursor, conn, filename)
    # conn.commit()
    # print('--- Built raw table from postgap pipeline results. ---')

    # # build genes table
    # build_ensembl_genes(cursor, conn)
    # conn.commit()

    # build lead variants table
    build_ensembl_lead_variants(cursor, conn)
    conn.commit()

    # build processed (merging previous three tables)
    build_processed(cursor, conn)
    conn.commit()

    # build chromosome tables (performance of locus queries)
    build_chroms(cursor, conn)
    conn.commit()

    # close the connection now we are done with it
    conn.close()


def build_raw(cursor, conn, filename):
    # create table from file (file must match OT format - see transform.py)
    castings = {
    'ld_snp_rsID': str,
    'chrom':str,
    'GRCh38_chrom':str,
    'GRCh38_gene_chrom':str,
    'gene_id': str,
    'gene_symbol': str,
    'gwas_pvalue': float,
    'ls_snp_is_gwas_snp': bool,
    'gene_id':str
    }
    df = pd.read_csv(filename, compression='gzip', sep='\t', na_values=['None'],
                     dtype=castings)
    print('--- Read file %s. Dumping to SQL... ---' % filename)
    df.to_sql('raw', conn)

    print('--- Adding indices on gene and SNP ---')
    cursor.executescript('''
    CREATE INDEX ix_raw_gene_id ON raw (gene_id);
    CREATE INDEX ix_raw_gwas_snp ON raw (gwas_snp);
    ''')

def batch(iterable, n=1):
    l = len(iterable)
    for ndx in range(0, l, n):
        yield iterable[ndx:min(ndx + n, l)]

def build_ensembl_genes(cursor, conn):
    '''queries the MySQL public ensembl database and outputs a gene lookup object
    in JSON format. It also injects into our sqlite database just so that we can
    do the processing directly there.
    '''

    #connect to Ensembl MySQL public server
    core = create_engine('mysql+mysqldb://anonymous@ensembldb.ensembl.org/homo_sapiens_core_92_38')

    q = """
    select
    et.exon_id,
    et.transcript_id,
    g.stable_id as gene_id,
    g.description,
    r.name as chr,
    g.seq_region_start as start,
    g.seq_region_end as end,
    e.seq_region_start as exon_start,
    e.seq_region_end as exon_end,
    t.seq_region_strand as fwdstrand
    from exon_transcript et, exon e, gene g, transcript t, seq_region r
    where
    g.canonical_transcript_id = et.transcript_id and
    g.seq_region_id = r.seq_region_id and
    r.coord_system_id = 4 and
    r.name NOT RLIKE 'CHR' and
    et.transcript_id = t.transcript_id and
    e.exon_id =et.exon_id
    """

    start_time = time.time()

    df = pd.read_sql_query(q, core,index_col='exon_id')
    df['exons'] = list(zip(df.exon_start, df.exon_end))
    df['fwdstrand'] = df['fwdstrand'].map({1:True,-1:False})
    df['tss'] = df.apply(lambda row: row['start'] if row['fwdstrand'] else row['end'], axis=1)
    keepcols = ['gene_id','description','tss','chr','start','end','fwdstrand']
    genes = pd.DataFrame(df.groupby(keepcols)['exons'].apply(list)).reset_index()
    genes.set_index('gene_id', inplace=True)
    print(genes['chr'].value_counts())
    genes.to_json(os.path.join(__dataprepdir__, 'genes.json'),orient='index')

    print("--- Genes table completed in %s seconds ---" % (time.time() - start_time))
    genes.loc[:,('start','end')].to_sql('gene', conn)

    # add indices
    cursor.execute('''
    CREATE INDEX ix_gene_gene_id ON gene (gene_id);
    ''')

def build_ensembl_lead_variants(cursor, conn):
    # retrieve all the unique lead variant ids from the raw table
    lead_variant_ids = [row[0] for row
                        in cursor.execute('SELECT DISTINCT gwas_snp FROM raw;')]
    print("--- Querying Ensembl Variation MySQL for %s SNPs---" % len(lead_variant_ids))
    #connect to Ensembl MySQL public server
    variation = create_engine('mysql+mysqldb://anonymous@ensembldb.ensembl.org/homo_sapiens_variation_92_38')
    start_time = time.time()

    q1 = '''
    select
    v.name as gwas_snp,
    vf.seq_region_start as position,
    r.name as seq_region_name
    from variation v, variation_feature vf, seq_region r
    where
    v.name in ('{0}')
    and v.variation_id = vf.variation_id
    and vf.seq_region_id = r.seq_region_id
    '''.format("','".join(lead_variant_ids))
    lead_variants = pd.read_sql_query(q1, variation)

    print("--- Variant table completed in %s seconds ---" % (time.time() - start_time))
    print("--- Committing Variant table to SQL ---")
    lead_variants.to_sql('lead_variant', conn)

    print("--- Adding indices for Variant table ---")
    cursor.executescript('''
    CREATE INDEX ix_lead_variants_lead_variant_id ON lead_variant (gwas_snp);
    ''')



def build_processed(cursor, conn):
    # create new table with gene and lead variant location information
    cursor.executescript('''
    CREATE TABLE processed AS
        SELECT
            raw.*,
            gene.start AS GRCh38_gene_start,
            gene.end AS GRCh38_gene_end,
            lead_variant.seq_region_name AS GRCh38_gwas_snp_chrom,
            lead_variant.position AS GRCh38_gwas_snp_pos
        FROM raw
        LEFT JOIN gene ON raw.gene_id = gene.gene_id
        LEFT JOIN lead_variant ON lead_variant.gwas_snp = raw.gwas_snp;
    ''')

    cursor.executescript('''
    CREATE INDEX ix_gene_id ON processed (gene_id);
    CREATE INDEX ix_ld_snp_rsID ON processed (ld_snp_rsID);
    CREATE INDEX ix_gwas_snp ON processed (gwas_snp);
    CREATE INDEX ix_disease_efo_id ON processed (disease_efo_id);
    CREATE INDEX ix_ld_snp_location ON processed (GRCh38_chrom, GRCh38_pos);
    CREATE INDEX ix_gene_start ON processed (GRCh38_gene_chrom, GRCh38_gene_start);
    CREATE INDEX ix_gene_end ON processed (GRCh38_gene_chrom, GRCh38_gene_end);
    CREATE INDEX ix_gwas_snp_location ON processed (GRCh38_gwas_snp_chrom, GRCh38_gwas_snp_pos);
    ''')

def build_chroms(cursor, conn):
    for chr in VALID_CHROMOSOMES:
        # setup
        table_name = 'chr_{}'.format(chr)
        create_sql = '''
        CREATE TABLE {table_name} AS
            SELECT *
            FROM processed
            WHERE
                GRCh38_gene_chrom="{chr}"
                OR GRCh38_chrom="{chr}"
                OR GRCh38_gwas_snp_chrom="{chr}";
        '''.format(chr=chr, table_name=table_name)
        indices_sql = '''
        CREATE INDEX ix{chr}_gene_start ON {table_name} (GRCh38_gene_start);
        CREATE INDEX ix{chr}_gene_end ON {table_name} (GRCh38_gene_end);
        CREATE INDEX ix{chr}_ld_snp_location ON {table_name} (GRCh38_pos);
        CREATE INDEX ix{chr}_gwas_snp_location ON {table_name} (GRCh38_gwas_snp_pos);
        '''.format(chr=chr, table_name=table_name)

        # create new table with gene and lead variant location information
        cursor.executescript(create_sql)

        # create indices
        cursor.executescript(indices_sql)


if __name__ == '__main__':
    filename = sys.argv[1]
    build_db(filename)
