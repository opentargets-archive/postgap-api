import os
import sys
import argparse
import pandas as pd
import numpy as np

VALID_CHROMOSOMES = [*[str(chr) for chr in range(23)], 'X', 'Y']
VALID_GWAS_SOURCES = ['GWAS Catalog']
MINIMUM_GTEX_VALUE = 0.999975
MAXIMUM_GENE_SNP_DISTANCE = 1000000
MINIMUM_VEP_VALUE = 0.65
GTEX_WEIGHT = 13.0
FANTOM5_WEIGHT = 3.0
DHS_WEIGHT = 1.5
PCHIC_WEIGHT = 1.5
FUNCGEN_TOTAL_WEIGHT = GTEX_WEIGHT + FANTOM5_WEIGHT + DHS_WEIGHT + PCHIC_WEIGHT
OT_G2V_PHASE_1_MIN = 0.9
OT_G2V_PHASE_1_MAX = 1.0
OT_G2V_PHASE_2_MIN = 0.5
OT_G2V_PHASE_2_MAX = 0.9
OT_G2V_PHASE_3_VALUE = 0.5

def calculate_open_targets_score(eco_scores, vep_terms, gtex, pchic, fantom5, dhs, nearest):
    # calculate for one row the open targets g2v score

    # phase 1 (sufficient VEP term)
    vep_terms_list = [term for term in str(vep_terms).split(',') if term in eco_scores]
    vep_score = max([
        eco_scores[term] for term in vep_terms_list
    ], default=0)
    if vep_score >= MINIMUM_VEP_VALUE:
        # linearly interpolate [MINIMUM_VEP_VALUE, 1] to [OT_G2V_PHASE_1_MIN, OT_G2V_PHASE_1_MAX]
        x = vep_score - MINIMUM_VEP_VALUE
        c = OT_G2V_PHASE_1_MIN
        m = (OT_G2V_PHASE_1_MAX - OT_G2V_PHASE_1_MIN) / (1 - MINIMUM_VEP_VALUE)
        y = m * x + c
        return y

    # phase 2 (at least some functional genomics)
    funcgen_score = 0
    if (gtex > MINIMUM_GTEX_VALUE):
        funcgen_score += gtex * GTEX_WEIGHT
    if (pchic > 0):
        funcgen_score += pchic * PCHIC_WEIGHT
    if (fantom5 > 0):
        funcgen_score += fantom5 * FANTOM5_WEIGHT
    if (dhs > 0):
        funcgen_score += dhs * DHS_WEIGHT
    funcgen_score /= FUNCGEN_TOTAL_WEIGHT
    if funcgen_score > 0:
        # linearly interpolate [0, 1] to [OT_G2V_PHASE_2_MIN, OT_G2V_PHASE_2_MAX]
        x = funcgen_score
        c = OT_G2V_PHASE_2_MIN
        m = (OT_G2V_PHASE_2_MAX - OT_G2V_PHASE_2_MIN)
        y = m * x + c
        return y

    # phase 3 (nearest gene)
    if nearest > 0:
        return OT_G2V_PHASE_3_VALUE

    # phase 4 (no score - should be filtered)
    return 0


def calculate_open_targets_scores(pg):
    # load the eco_scores.tsv (which contains mappings from VEP terms to 0-1 values)
    eco_scores_filename = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'eco_scores.tsv')
    eco_scores_df = pd.read_csv(eco_scores_filename, sep='\t', na_values=['None'])
    eco_scores = pd.Series(eco_scores_df['value'].values, index=eco_scores_df['term']).to_dict()
    pg['ot_g2v_score'] = pg.apply(
        lambda row: calculate_open_targets_score(
            eco_scores, row['vep_terms'], row['GTEx'], row['PCHiC'], row['Fantom5'], row['DHS'], row['Nearest']
        ),
        axis=1
    )
    print('generated ot_g2v_score')

    # filter out zeros
    valid_ot_g2v_score = pg['ot_g2v_score'] > 0
    pg = pg[valid_ot_g2v_score]

    return pg


def open_targets_transform(filename,nrows):
    '''
    Iterate the rows in the tsv file and output a new tsv that meets the
    Open Targets format requirements.
    '''
    # load
    pg = pd.read_csv(filename, sep='\t', na_values=['None'],nrows=nrows,
                     dtype={'GRCh38_chrom':str,'GRCh38_gene_chrom':str})
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

    # calculate the open targets g2v score
    pg = calculate_open_targets_scores(pg)
    print('{} rows (after calculating/filtering for valid g2v score)'.format(pg.shape[0]))

    # write out
    pg.to_csv('{}.transformed'.format(filename), sep='\t', compression='gzip')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('filename', default='https://storage.googleapis.com/postgap-data/postgap.20180324.txt.gz')
    parser.add_argument('--sample', action='store_const', const=5000,
                        help='run on a small subsample')

    args = parser.parse_args()
    if args.sample:
        print('Running analysis on the first {} lines'.format(args.sample))
    open_targets_transform(args.filename, nrows=args.sample)
