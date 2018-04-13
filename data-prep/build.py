import sys
import sqlite3
import requests
import pandas as pd
import json
import time

API_INTER_CALL_PAUSE_S = 10
GENE_BATCH_SIZE = 500
LEAD_VARIANT_BATCH_SIZE = 200
GENE_TABLE_COLS = [
    'gene_id',
    'display_name',
    'description',
    'seq_region_name',
    'start',
    'end',
    'strand',
    'biotype',
    'canonical_transcript'
]
LEAD_VARIANT_TABLE_COLS = [
    'gwas_snp',
    'seq_region_name',
    'position',
    'strand',
]

def build_db(filename):
    '''
    Build a sqlite database from a tsv file that meets the
    Open Targets format requirements.
    '''
    # creates db if db does not exist
    # conn = sqlite3.connect('postgap.db.scored')
    conn = sqlite3.connect('postgap.20180324.v0.0.1.db')
    cursor = conn.cursor()

    # build raw table
    build_raw(cursor, conn)
    conn.commit()

    # build genes table
    build_ensembl_genes(cursor, conn)
    conn.commit()

    # build lead variants table
    build_ensembl_lead_variants(cursor, conn)
    conn.commit()

    # build processed (merging previous three tables)
    build_processed(cursor, conn)
    conn.commit()

    # close the connection now we are done with it
    conn.close()

def build_raw(cursor, conn):
    # create table from file (file must match OT format - see transform.py)
    df = pd.read_csv(filename, compression='gzip', sep='\t', na_values=['None'])
    df.to_sql('raw', conn)

    # add indices
    cursor.executescript('''
    CREATE INDEX ix_raw_gene_id ON raw (gene_id);
    CREATE INDEX ix_raw_gwas_snp ON raw (gwas_snp);
    ''')

def batch(iterable, n=1):
    l = len(iterable)
    for ndx in range(0, l, n):
        yield iterable[ndx:min(ndx + n, l)]

def build_ensembl_genes(cursor, conn):
    # retrieve all the unique gene ids from the raw table
    gene_ids = [row[0] for row
                in cursor.execute('SELECT DISTINCT gene_id FROM raw;')]

    # ensembl rest api supports queries of up to 1000 genes at a time
    genes = pd.DataFrame(columns=GENE_TABLE_COLS)
    counter = 1
    for gene_ids_batch in batch(gene_ids, n=GENE_BATCH_SIZE):
        genes_batch = fetch_ensembl_genes(gene_ids_batch)
        genes = pd.concat([genes, genes_batch], ignore_index=True)
        print('Called Ensembl genes ({})'.format(counter))
        counter += 1

        # delay before next call
        time.sleep(API_INTER_CALL_PAUSE_S)
        
    # create table
    genes.to_sql('gene', conn)

    # add indices
    cursor.executescript('''
    CREATE INDEX ix_gene_gene_id ON gene (gene_id);
    ''')

def fetch_ensembl_genes(gene_ids):
    print('Making request for {} genes'.format(len(gene_ids)))
    r = requests.post('https://rest.ensembl.org/lookup/id', json = {'expand': True, 'ids': gene_ids})
    r.raise_for_status()
    print('Request answered with status: {}'.format(r.status_code))
    raw_dict = r.json()
    cleaned = []
    skipped = 0
    for gene_id, gene in raw_dict.items():
        if (not gene) or ('Transcript' not in gene):
            skipped += 1
            print('No transcript for gene {}, total skipped is {}'.format(gene_id, skipped))
            continue

        # create json string for canonical transcript
        raw_transcripts = gene['Transcript']
        raw_canonical_transcript = [
            t for t in raw_transcripts
            if t['is_canonical'] == 1
        ][0]

        if ('Translation' not in raw_canonical_transcript):
            skipped += 1
            print('No translation for gene {} (but has transcript), total skipped is {}'.format(gene_id, skipped))
            continue

        tss = raw_canonical_transcript['start'] if (raw_canonical_transcript['strand'] == 1) else raw_canonical_transcript['end']
        clean_canonical_transcript = json.dumps({
            'id': raw_canonical_transcript['id'],
            'start': raw_canonical_transcript['start'],
            'end': raw_canonical_transcript['end'],
            'forwardStrand': (raw_canonical_transcript['strand'] == 1),
            'tss': tss,
            'translationStart': raw_canonical_transcript['Translation']['start'],
            'translationEnd': raw_canonical_transcript['Translation']['end'],
            'exons': [{
                'id': exon['id'],
                'start': exon['start'],
                'end': exon['end']
            } for exon in raw_canonical_transcript['Exon']]
        })

        # create tuple for insertion into db table
        clean_gene = (
            gene_id,
            gene['display_name'],
            gene['description'] if ('description' in gene) else '',
            gene['seq_region_name'],
            gene['start'],
            gene['end'],
            gene['strand'] == 1,
            gene['biotype'],
            clean_canonical_transcript
        )

        cleaned.append(clean_gene)

    print('Processed {} genes'.format(len(cleaned)))

    cleaned_df = pd.DataFrame(cleaned, columns=GENE_TABLE_COLS)
    return cleaned_df

def build_ensembl_lead_variants(cursor, conn):
    # retrieve all the unique lead variant ids from the raw table
    lead_variant_ids = [row[0] for row
                        in cursor.execute('SELECT DISTINCT gwas_snp FROM raw;')]

    # ensembl rest api supports queries of up to 200 variants at a time
    lead_variants = pd.DataFrame(columns=LEAD_VARIANT_TABLE_COLS)
    counter = 1
    for lead_variant_ids_batch in batch(lead_variant_ids, n=LEAD_VARIANT_BATCH_SIZE):
        lead_variants_batch = fetch_ensembl_variants(lead_variant_ids_batch)
        lead_variants = pd.concat([lead_variants, lead_variants_batch], ignore_index=True)
        print('Called Ensembl variants ({})'.format(counter))
        counter += 1

        # delay before next call
        time.sleep(API_INTER_CALL_PAUSE_S)
        
    # create table
    lead_variants.to_sql('lead_variant', conn)

    # add indices
    cursor.executescript('''
    CREATE INDEX ix_lead_variants_lead_variant_id ON lead_variant (gwas_snp);
    ''')

def fetch_ensembl_variants(variant_ids):
    print('Making request for {} variants'.format(len(variant_ids)))
    r = requests.post('https://rest.ensembl.org/variation/homo_sapiens', json = {'ids': variant_ids})
    r.raise_for_status()
    print('Request answered with status: {}'.format(r.status_code))
    raw_dict = r.json()
    cleaned = []
    for variant_id, variant in raw_dict.items():
        if ('mappings' not in variant) or len(variant['mappings']) == 0:
            print('No mappings for {}'.format(variant_id))
            continue

        # grab the first mapping object
        mapping = variant['mappings'][0]

        # create tuple for insertion into db table
        clean_variant = (
            variant_id,
            mapping['seq_region_name'],
            mapping['start'],
            mapping['strand']
        )

        cleaned.append(clean_variant)

    print('Processed {} lead variants'.format(len(cleaned)))

    cleaned_df = pd.DataFrame(cleaned, columns=LEAD_VARIANT_TABLE_COLS)
    return cleaned_df

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


if __name__ == '__main__':
    filename = sys.argv[1]
    build_db(filename)