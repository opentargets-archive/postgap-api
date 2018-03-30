import sys
import pandas as pd

VALID_CHROMOSOMES = [*[str(chr) for chr in range(23)], 'X', 'Y']
VALID_GWAS_SOURCES = ['GWAS Catalog']
MINIMUM_GTEX_VALUE = 0.999975
MAXIMUM_GENE_SNP_DISTANCE = 1000000

def open_targets_transform(filename):
    '''
    Iterate the rows in the tsv file and output a new tsv that meets the
    Open Targets format requirements.
    '''
    # load
    pg = pd.read_csv(filename, sep='\t', na_values=['None'])
    
    # filter for gwas source
    pg = pg[pg['gwas_source'].isin(VALID_GWAS_SOURCES)]

    # filter for chromosomes
    pg = pg[pg['GRCh38_chrom'].isin(VALID_CHROMOSOMES)]
    pg = pg[pg['GRCh38_gene_chrom'].isin(VALID_CHROMOSOMES)]

    # filter for gene and snp on same chromosome
    pg = pg[pg['GRCh38_chrom'] == pg['GRCh38_gene_chrom']]

    # filter for gene and snp maximum of 1MB apart
    pg = pg[abs(pg['GRCh38_pos'] - pg['GRCh38_gene_pos']) <= MAXIMUM_GENE_SNP_DISTANCE]

    # TODO: filter for gtex (carefully)
    # TODO: calculate the open targets g2v score
    # TODO: add gene start, end, canonical transcript, exons

    # write out
    pg.to_csv('{}.transformed'.format(filename), sep='\t', compression='gzip')

if __name__ == '__main__':
    filename = sys.argv[1]
    open_targets_transform(filename)
