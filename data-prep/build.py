import sys
import sqlite3
import requests
import pandas as pd
import json
import time

API_INTER_CALL_PAUSE_S = 5
GENE_BATCH_SIZE = 5
LEAD_VARIANT_BATCH_SIZE = 5
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
    conn = sqlite3.connect('postgap.db')
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

    # close the connection now we are done with it
    conn.close()

def build_raw(cursor, conn):
    # create table from file (file must match OT format - see transform.py)
    df = pd.read_csv(filename, compression='gzip', sep='\t', na_values=['None'])
    df.to_sql('raw', conn)

    # add indices
    cursor.executescript('''
    CREATE INDEX ix_gene_id ON raw (gene_id);
    CREATE INDEX ix_ld_snp_rsID ON raw (ld_snp_rsID);
    CREATE INDEX ix_gwas_snp ON raw (gwas_snp);
    CREATE INDEX ix_disease_efo_id ON raw (disease_efo_id);
    CREATE INDEX ix_gene_location ON raw (GRCh38_gene_chrom, GRCh38_gene_pos);
    CREATE INDEX ix_ld_snp_location ON raw (GRCh38_chrom, GRCh38_pos);
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
        if counter == 5:
            break
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
    print('Making request for {}'.format(gene_ids))
    r = requests.post('https://rest.ensembl.org/lookup/id', json = {'expand': True, 'ids': gene_ids})
    r.raise_for_status()
    print('Request answered: {}'.format(r.status_code))
    raw_dict = r.json()
    print('Request JSON: {}'.format(raw_dict))
    cleaned = []
    for gene_id, gene in raw_dict.items():
        # create json string for canonical transcript
        raw_transcripts = gene['Transcript']
        raw_canonical_transcript = [
            t for t in raw_transcripts
            if t['is_canonical'] == 1
        ][0]
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
            gene['description'],
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
        if counter == 5:
            break
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
    print('Making request for {}'.format(variant_ids))
    r = requests.post('https://rest.ensembl.org/variation/homo_sapiens', json = {'ids': variant_ids})
    r.raise_for_status()
    print('Request answered: {}'.format(r.status_code))
    raw_dict = r.json()
    print('Request JSON: {}'.format(raw_dict))
    cleaned = []
    for variant_id, variant in raw_dict.items():
        if not variant['mappings']:
            print('No mappings for {}'.format(variant_id))

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


if __name__ == '__main__':
    filename = sys.argv[1]
    build_db(filename)