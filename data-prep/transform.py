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
    print('{} rows (input file)'.format(pg.shape[0]))
    
    # filter for gwas source
    pg = pg[pg['gwas_source'].isin(VALID_GWAS_SOURCES)]
    print('{} rows (after filtering gwas_source)'.format(pg.shape[0]))

    # filter for chromosomes
    pg = pg[pg['GRCh38_chrom'].isin(VALID_CHROMOSOMES)]
    pg = pg[pg['GRCh38_gene_chrom'].isin(VALID_CHROMOSOMES)]
    print('{} rows (after filtering for valid chromosomes)'.format(pg.shape[0]))

    # filter for gene and snp on same chromosome
    pg = pg[pg['GRCh38_chrom'] == pg['GRCh38_gene_chrom']]
    print('{} rows (after filtering for gene-variant chromosome match)'.format(pg.shape[0]))

    # filter for gene and snp maximum of 1MB apart
    pg = pg[abs(pg['GRCh38_pos'] - pg['GRCh38_gene_pos']) <= MAXIMUM_GENE_SNP_DISTANCE]
    print('{} rows (after filtering for gene-variant distance)'.format(pg.shape[0]))

    # filter for gtex (carefully)
    valid_gtex = pg['GTEx'] > MINIMUM_GTEX_VALUE
    any_other_score = (pg['VEP'] > 0) | (pg['PCHiC'] > 0) | (pg['Fantom5'] > 0) | (pg['DHS'] > 0) | (pg['Nearest'] > 0)
    pg = pg[valid_gtex | any_other_score]
    print('{} rows (after filtering for valid gtex)'.format(pg.shape[0]))
    
    # TODO: calculate the open targets g2v score
    # TODO: add gene start, end, canonical transcript, exons

    # write out
    pg.to_csv('{}.transformed'.format(filename), sep='\t', compression='gzip')

if __name__ == '__main__':
    filename = sys.argv[1]
    open_targets_transform(filename)
