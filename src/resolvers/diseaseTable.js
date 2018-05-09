const resolverDiseaseTable = (_, { efoId, offset, limit }, { db }) => {
    const params = {$efoId: efoId, $offset: offset, $limit: limit};

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
        ot_g2v_score_reason as otG2VReason,
        ot_vep as vep,
        vep_terms as vepTerms,
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
    WHERE disease_efo_id=$efoId
    LIMIT $limit
    OFFSET $offset
    `;
    const rowsQuery = db.all(rowsSql, params);

    // total rows
    const totalRowsSql = `
    SELECT COUNT(*) as total
    FROM processed
    WHERE disease_efo_id=$efoId
    `;
    const totalRowsQuery = db.get(totalRowsSql, { $efoId: efoId });

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

export default resolverDiseaseTable;
