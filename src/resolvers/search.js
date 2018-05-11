const resolveSearch = (_, { queryString }, { db }) => {
    // currently gene and disease are retrieved using OT API
    // search currently uses exact match on rs ids, but we could use
    // free text search in the future (see sqlite's fts3/fts4)
    const params = { queryString };
    const variantSql = `
        SELECT ld_snp_rsID AS id
        FROM processed
        WHERE ld_snp_rsID="${queryString}"
        GROUP BY ld_snp_rsID;
    `;
    const variantQuery = db.all(variantSql);
    const leadVariantSql = `
        SELECT gwas_snp AS id
        FROM processed
        WHERE gwas_snp="${queryString}"
        GROUP BY gwas_snp;
    `;
    const leadVariantQuery = db.all(leadVariantSql);
    return Promise.all([variantQuery, leadVariantQuery]).then(([variant, leadVariant]) => {
        const isVariant = variant.length > 0;
        const isLeadVariant = leadVariant.length > 0;

        if (isVariant) {
            const { id } = variant[0];
            const result = { id, name: id, type: 'variant' }
            return [result];
        } else if (isLeadVariant) {
            const { id } = leadVariant[0];
            const result = { id, name: id, type: 'variant' }
            return [result];
        } else {
            return [];
        }
    });
};

export default resolveSearch;
