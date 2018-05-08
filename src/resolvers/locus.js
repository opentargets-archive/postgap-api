// helper methods
const getLocusFiltersSql = ({ g2VMustHaves, g2VScore, r2, gwasPValue }) => {
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
    return filtersSql
}
const templateWhere = filtersSql => `
WHERE
    (GRCh38_gene_start IS NOT NULL)
    AND (GRCh38_gene_end IS NOT NULL)
    AND (GRCh38_gwas_snp_pos IS NOT NULL)
    ${filtersSql}
    AND (
        (GRCh38_gene_start>=$start AND GRCh38_gene_start<=$end)
        OR (GRCh38_gene_end>=$start AND GRCh38_gene_end<=$end)
        OR (GRCh38_pos>=$start AND GRCh38_pos<=$end)
        OR (GRCh38_gwas_snp_pos>=$start AND GRCh38_gwas_snp_pos<=$end))
`;
const getFilteredLocusWhereSql = (args) => {
    const filtersSql = getLocusFiltersSql(args)
    return templateWhere(filtersSql);
}
const getUnfilteredLocusWhereSql = () => {
    return templateWhere('');
}
const getSelectedSql = ({ selectedId, selectedType }) => {
    let selectedSql = '';
    switch (selectedType) {
        case 'gene':
            selectedSql = `COUNT(CASE gene_id WHEN "${selectedId}" THEN 1 ELSE null END) > 0 AS selected,`;
            break;
        case 'variant':
            selectedSql = `COUNT(CASE ld_snp_rsID WHEN "${selectedId}" THEN 1 ELSE null END) > 0 AS selected,`;
            break;
        case 'leadVariant':
            selectedSql = `COUNT(CASE gwas_snp WHEN "${selectedId}" THEN 1 ELSE null END) > 0 AS selected,`;
            break;
        case 'disease':
            selectedSql = `COUNT(CASE disease_efo_id WHEN "${selectedId}" THEN 1 ELSE null END) > 0 AS selected,`;
            break;
        case 'geneVariant':
            const [geneId, vId] = selectedId.split('-');
            selectedSql = `COUNT(CASE ((gene_id = "${geneId}") AND (ld_snp_rsID = "${vId}")) WHEN 1 THEN 1 ELSE null END) > 0 AS selected,`;
            break;
        case 'variantLeadVariant':
            const [varId, lvId] = selectedId.split('-');
            selectedSql = `COUNT(CASE ((ld_snp_rsID = "${varId}") AND (gwas_snp = "${lvId}")) WHEN 1 THEN 1 ELSE null END) > 0 AS selected,`;
            break;
        case 'leadVariantDisease':
            const [leadVarId, diseaseId] = selectedId.split('-');
            selectedSql = `COUNT(CASE ((gwas_snp = "${leadVarId}") AND (disease_efo_id = "${diseaseId}")) WHEN 1 THEN 1 ELSE null END) > 0 AS selected,`;
            break;
        default:
            break;
    }
    return selectedSql;
}

const commonSetup = (args) => {
    const { chromosome, start, end, g2VMustHaves, g2VScore, r2, gwasPValue, selectedId, selectedType } = args;
    const params = { $start: start, $end: end };
    const tableName = `chr_${chromosome}`;

    const filteredWhere = getFilteredLocusWhereSql(args);
    const unfilteredWhere = getUnfilteredLocusWhereSql();

    let selectedSql = getSelectedSql(args);
    const orderBySql = selectedSql ? 'ORDER BY selected' : '';
    return {
        tableName,
        filteredWhere,
        unfilteredWhere,
        selectedSql,
        orderBySql,
        params
    };
};

const resolveGenes = ({ common }, args, { db, geneLocationsCache }) => {
    const { filteredWhere, unfilteredWhere, selectedSql, orderBySql, params, tableName } = common;
    const genesSql = `
    -- SQL_QUERY_TYPE=genes
    SELECT
        ${selectedSql}
        gene_id as id,
        gene_symbol as symbol,
        GRCH38_gene_chrom as chromosome,
        GRCh38_gene_start as start,
        GRCh38_gene_end as end
    FROM ${tableName}
    ${unfilteredWhere}
    GROUP BY gene_id
    ${orderBySql}
    `;

    const genesQuery = db.all(genesSql, params).then(genes => {
        const start = Date.now();
        const genesWithLocations = genes.map(d => {
            const geneLocation = geneLocationsCache[d.id];
            return {
                ...d,
                forwardStrand: geneLocation.fwdstrand,
                // canonicalTranscript: geneLocation.canonicalTranscript
                tss: geneLocation.tss,
                start: geneLocation.start,
                end: geneLocation.end,
                exons: geneLocation.exons
            }
        });
        const end = Date.now();
        console.log('gene-location-inject-time', 'genes', end - start);
        return genesWithLocations;
    });
    return genesQuery;
}

const resolveVariants = ({ common }, args, { db }) => {
    const { filteredWhere, unfilteredWhere, selectedSql, orderBySql, params, tableName } = common;
    const variantsSql = `
    -- SQL_QUERY_TYPE=variants
    SELECT
        ${selectedSql}
        ld_snp_rsID as id,
        GRCh38_pos as position
        FROM ${tableName}
    ${unfilteredWhere}
    GROUP BY ld_snp_rsID
    ${orderBySql}
    `;
    const variantsQuery = db.all(variantsSql, params);
    return variantsQuery;
}

const resolveLeadVariants = ({ common }, args, { db }) => {
    const { filteredWhere, unfilteredWhere, selectedSql, orderBySql, params, tableName } = common;
    const leadVariantsSql = `
    -- SQL_QUERY_TYPE=leadVariants
    SELECT
        ${selectedSql}
        gwas_snp as id,
        GRCh38_gwas_snp_pos as position
    FROM ${tableName}
    ${unfilteredWhere}
    GROUP BY gwas_snp
    ${orderBySql}
    `;
    const leadVariantsQuery = db.all(leadVariantsSql, params);
    return leadVariantsQuery;
}

const resolveDiseases = ({ common }, args, { db }) => {
    const { filteredWhere, unfilteredWhere, selectedSql, orderBySql, params, tableName } = common;
    const diseasesSql = `
    -- SQL_QUERY_TYPE=diseases
    SELECT
        ${selectedSql}
        disease_efo_id as id,
        disease_name as name
    FROM ${tableName}
    ${unfilteredWhere}
    GROUP BY disease_efo_id
    ${orderBySql}
    `;
    const diseasesQuery = db.all(diseasesSql, params);
    return diseasesQuery;
}

const resolveGeneVariants = ({ common }, args, { db, geneLocationsCache }) => {
    const { filteredWhere, unfilteredWhere, selectedSql, orderBySql, params, tableName } = common;
    const geneVariantsSql = `
    -- SQL_QUERY_TYPE=geneVariants
    SELECT
        ${selectedSql}
        (gene_id || "-" || ld_snp_rsID) AS id,
        gene_id as geneId,
        gene_symbol as geneSymbol,
        GRCH38_gene_chrom as geneChromosome,
        GRCh38_gene_pos as geneTss,
        ld_snp_rsID as vId,
        GRCh38_pos as vPos,
        ot_g2v_score as otG2VScore,
        VEP as vep,
        GTEx as gtex,
        PCHiC as pchic,
        Fantom5 as fantom5,
        DHS as dhs,
        Nearest as nearest
    FROM ${tableName}
    ${filteredWhere}
    GROUP BY gene_id, ld_snp_rsID
    ${orderBySql}
    `;
    const geneVariantsQuery = db.all(geneVariantsSql, params).then(geneVariants => {
        const start = Date.now();
        const geneVariantsWithLocations = geneVariants.map(d => {
            const geneLocation = geneLocationsCache[d.geneId];
            return {
                ...d,
                // canonicalTranscript: geneLocation.canonicalTranscript
                geneTss: geneLocation.tss
            };
        });
        const end = Date.now();
        console.log('gene-location-inject-time', 'geneVariants', end - start);
        return geneVariantsWithLocations;
    });
    return geneVariantsQuery;
}

const resolveVariantLeadVariants = ({ common }, args, { db }) => {
    const { filteredWhere, unfilteredWhere, selectedSql, orderBySql, params, tableName } = common;
    const variantLeadVariantsSql = `
    -- SQL_QUERY_TYPE=variantLeadVariants
    SELECT
        ${selectedSql}
        (ld_snp_rsID || "-" || gwas_snp) AS id,
        ld_snp_rsID as vId,
        GRCh38_pos as vPos,
        gwas_snp as lvId,
        GRCh38_gwas_snp_pos as lvPos,
        r2
    FROM ${tableName}
    ${filteredWhere}
    GROUP BY ld_snp_rsID, gwas_snp
    ${orderBySql}
    `;
    const variantLeadVariantsQuery = db.all(variantLeadVariantsSql, params);
    return variantLeadVariantsQuery;
}

const resolveLeadVariantDiseases = ({ common }, args, { db }) => {
    const { filteredWhere, unfilteredWhere, selectedSql, orderBySql, params, tableName } = common;
    const leadVariantDiseasesSql = `
    -- SQL_QUERY_TYPE=leadVariantDiseases
    SELECT
        ${selectedSql}
        (gwas_snp || "-" || disease_efo_id) AS id,
        gwas_snp as lvId,
        GRCh38_gwas_snp_pos as lvPos,
        disease_efo_id as efoId,
        disease_name as efoName,
        gwas_pvalue as gwasPValue,
        gwas_odds_ratio as gwasOddsRatio,
        gwas_beta as gwasBeta,
        gwas_study as gwasStudy,
        gwas_pmid as gwasPMId,
        gwas_size as gwasSize
    FROM ${tableName}
    ${filteredWhere}
    GROUP BY gwas_snp, disease_efo_id
    ${orderBySql}
    `;
    const leadVariantDiseasesQuery = db.all(leadVariantDiseasesSql, params);
    return leadVariantDiseasesQuery;
}

const resolveMaxGwasPValue = ({ common }, args, { db }) => {
    const { filteredWhere, unfilteredWhere, selectedSql, orderBySql, params, tableName } = common;
    const maxGwasPValueSql = `
    -- SQL_QUERY_TYPE=maxGwasPValue
    SELECT MIN(gwas_pvalue) as minGwasPValue
    FROM ${tableName}
    ${unfilteredWhere}
    `;
    const maxGwasPValueQuery = db.all(maxGwasPValueSql, params)
    .then(rows => {
        if (rows && rows.length > 0) {
            return -Math.log10(rows[0].minGwasPValue);
        } else {
            return Number.MAX_SAFE_INTEGER;
        }
    });
    return maxGwasPValueQuery;
}

export const resolveLocus = (_, args) => {
    const common = commonSetup(args);
    return { common };
};

const resolveDrawableLocus = {
    genes: resolveGenes,
    variants: resolveVariants,
    leadVariants: resolveLeadVariants,
    diseases: resolveDiseases,
    geneVariants: resolveGeneVariants,
    variantLeadVariants: resolveVariantLeadVariants,
    leadVariantDiseases: resolveLeadVariantDiseases,
    maxGwasPValue: resolveMaxGwasPValue
}

export default resolveDrawableLocus;
