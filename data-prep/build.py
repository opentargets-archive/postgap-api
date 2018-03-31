import sys
import sqlite3
import pandas as pd

def build_db(filename):
    '''
    Build a sqlite database from a tsv file that meets the
    Open Targets format requirements.
    '''
    # creates db if db does not exist
    conn = sqlite3.connect('postgap.db')
    c = conn.cursor()

    # create table from file (file must match OT format - see transform.py)
    df = pd.read_csv(filename, compression='gzip', sep='\t', na_values=['None'])
    df.to_sql('raw', conn)

    # add indices
    c.executescript('''
    CREATE INDEX ix_gene_id ON raw (gene_id);
    CREATE INDEX ix_ld_snp_rsID ON raw (ld_snp_rsID);
    CREATE INDEX ix_gwas_snp ON raw (gwas_snp);
    CREATE INDEX ix_disease_efo_id ON raw (disease_efo_id);
    CREATE INDEX ix_gene_location ON raw (GRCh38_gene_chrom, GRCh38_gene_pos);
    CREATE INDEX ix_ld_snp_location ON raw (GRCh38_chrom, GRCh38_pos);
    ''')

    # save the changes
    conn.commit()

    # close the connection now we are done with it
    conn.close()
    

if __name__ == '__main__':
    filename = sys.argv[1]
    build_db(filename)