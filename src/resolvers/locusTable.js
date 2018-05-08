const resolveLocusTable = (_, { chromosome, start, end, g2VMustHaves, g2VScore, r2, gwasPValue, selectedId, selectedType, offset, limit }, { db }) => {
    const params = {$chromosome: chromosome, $start: start, $end: end, $offset: offset, $limit: limit};
    const paramsWithoutPagination = { $chromosome: chromosome, $start: start, $end: end };
    let filtersSql = g2VMustHaves.length > 0 ? (g2VMustHaves.map(d => `AND (${d.toLowerCase()} > 0)`).join(' ')) : '';
    if (g2VScore && g2VScore.length === 2) {
        filtersSql += ` AND (ot_g2v_score >= ${g2VScore[0]}) AND (ot_g2v_score <= ${g2VScore[1]})`;
    }
    if (r2 && r2.length === 2) {
        filtersSql += ` AND (r2 >= ${r2[0]}) AND (r2 <= ${r2[1]})`;
    }
    if (gwasPValue && gwasPValue.length === 2) {
        // note: need to do reverse -log10 conversion
        filtersSql += ` AND (gwas_pvalue >= ${10 ** (-gwasPValue[1])}) AND (gwas_pvalue <= ${10 ** (-gwasPValue[0])})`;
    }

    let selectedSql = '';
    switch (selectedType) {
        case 'gene':
            selectedSql = `AND gene_id = "${selectedId}"`;
            break;
        case 'variant':
            selectedSql = `AND ld_snp_rsID = "${selectedId}"`;
            break;
        case 'leadVariant':
            selectedSql = `AND gwas_snp = "${selectedId}"`;
            break;
        case 'disease':
            selectedSql = `AND disease_efo_id = "${selectedId}"`;
            break;
        case 'geneVariant':
            const [geneId, vId] = selectedId.split('-');
            selectedSql = `AND (gene_id = "${geneId}" AND ld_snp_rsID = "${vId}")`;
            break;
        case 'variantLeadVariant':
            const [varId, lvId] = selectedId.split('-');
            selectedSql = `AND (ld_snp_rsID = "${varId}" AND gwas_snp = "${lvId}")`;
            break;
        case 'leadVariantDisease':
            const [leadVarId, diseaseId] = selectedId.split('-');
            selectedSql = `AND (gwas_snp = "${leadVarId}" AND disease_efo_id = "${diseaseId}")`;
            break;
        default:
            break;
    }

    const templateWhere = (filtersSql, selectedSql) => `
    WHERE
        (GRCh38_gene_start IS NOT NULL)
        AND (GRCh38_gene_end IS NOT NULL)
        AND (GRCh38_gwas_snp_pos IS NOT NULL)
        ${filtersSql}
        ${selectedSql}
        AND (
            (GRCh38_gene_chrom=$chromosome AND GRCh38_gene_start>=$start AND GRCh38_gene_start<=$end)
            OR (GRCh38_gene_chrom=$chromosome AND GRCh38_gene_end>=$start AND GRCh38_gene_end<=$end)
            OR (GRCh38_chrom=$chromosome AND GRCh38_pos>=$start AND GRCh38_pos<=$end)
            OR (GRCh38_gwas_snp_chrom=$chromosome AND GRCh38_gwas_snp_pos>=$start AND GRCh38_gwas_snp_pos<=$end))
    `;
    const filteredWhere = templateWhere(filtersSql, selectedSql);

    // rows
    const rowsSql = `
    SELECT
        "index",
        gene_id as geneId,
        gene_symbol as geneSymbol,
        GRCH38_gene_chrom as geneChromosome,
        GRCh38_gene_pos as geneTss,
        ld_snp_rsID as vId,
        GRCH38_chrom as variantChromosome,
        GRCh38_pos as vPos,
        gwas_snp as lvId,
        disease_efo_id as efoId,
        disease_name as efoName,
        ot_g2v_score as otG2VScore,
        VEP as vep,
        GTEx as gtex,
        PCHiC as pchic,
        Fantom5 as fantom5,
        DHS as dhs,
        Nearest as nearest,
        r2,
        gwas_pvalue as gwasPValue,
        gwas_odds_ratio as gwasOddsRatio,
        gwas_beta as gwasBeta,
        gwas_size as gwasSize,
        gwas_pmid as gwasPMId,
        gwas_study as gwasStudy
    FROM processed
    ${filteredWhere}
    LIMIT $limit
    OFFSET $offset
    `;
    const rowsQuery = db.all(rowsSql, params);

    // total rows
    const totalRowsSql = `
    SELECT COUNT(*) as total
    FROM processed
    ${filteredWhere}
    `;
    const totalRowsQuery = db.get(totalRowsSql, paramsWithoutPagination);

    // wait for all queries and return composite object
    return Promise.all([
        rowsQuery,
        totalRowsQuery,
    ]).then(([rows, { total }]) => {
        return {
            rows,
            total,
            offset,
            limit
        }
    })
};

export default resolveLocusTable;
