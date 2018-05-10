const resolveGeneLocation = (_, { variantId }, { db }) => {
    const params = { $variantId: variantId };
    const variantSql = `
        SELECT
            ld_snp_rsID AS id,
            GRCh38_chrom AS chromosome,
            GRCh38_pos AS position
        FROM processed
        WHERE ld_snp_rsID=$variantId
        GROUP BY ld_snp_rsID;
    `;
    const variantQuery = db.all(variantSql, params);
    const leadVariantSql = `
        SELECT
            gwas_snp AS id,
            GRCh38_gwas_snp_chrom AS chromosome,
            GRCh38_gwas_snp_pos AS position
        FROM processed
        WHERE gwas_snp=$variantId
        GROUP BY gwas_snp;
    `;
    const leadVariantQuery = db.all(leadVariantSql, params);
    return Promise.all([variantQuery, leadVariantQuery]).then(([variant, leadVariant]) => {
        const isVariant = variant.length > 0;
        const isLeadVariant = leadVariant.length > 0;

        if (!isVariant && !isLeadVariant) {
            return null;
        }

        let value = { id: variantId, isVariant, isLeadVariant };
        if (isVariant) {
            const { chromosome, position } = variant[0];
            value = { ...value, chromosome, position };
        } else if (isLeadVariant) {
            const { chromosome, position } = leadVariant[0];
            value = { ...value, chromosome, position };
        }

        return value;
    });
};

export default resolveGeneLocation;
