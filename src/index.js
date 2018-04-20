import fs from 'fs';
import path from 'path';
import express from 'express';
import cors from 'cors';
import bodyParser from 'body-parser';
import { graphqlExpress, graphiqlExpress } from 'apollo-server-express';
import { buildSchema } from 'graphql';
import { makeExecutableSchema } from 'graphql-tools';
import sqlite3 from 'sqlite3';
import {promisify} from 'bluebird';

// open a connection to the database
let db = new sqlite3.Database('postgap.20180324.v0.0.1.db', sqlite3.OPEN_READONLY, err => {
    if (err) {
        console.error(err.message)
    } else {
        console.log('Connected to POSTGAP database.')
    }
})

// sqlite3 methods do not support promises out of the box
db.all = promisify(db.all);
db.get = promisify(db.get);

// total number of genes is ~20,000,
// so cache ensembl calls on an object
let geneCache = {};

// total number of lead variants is ~40,000 (which need location info),
// so cache ensembl calls on an object
let leadVariantCache = {};

// load the schema
const schemaFile = path.join(__dirname, 'schema.gql');
const typeDefs = fs.readFileSync(schemaFile, 'utf8');

// specify the resolution methods for allowed query
const resolvers = {
    Query: {
        locus: (_, { chromosome, start, end, g2VMustHaves, g2VScore, r2, gwasPValue, selectedId, selectedType }) => {
            const params = {$chromosome: chromosome, $start: start, $end: end};
            const oldWhere = `
            WHERE
                (GRCh38_gene_chrom=$chromosome AND GRCh38_gene_pos>=$start AND GRCh38_gene_pos<=$end)
                OR (GRCh38_chrom=$chromosome AND GRCh38_pos>=$start AND GRCh38_pos<=$end)
            `;
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
            const templateWhere = filtersSql => `
            WHERE
                (GRCh38_gene_start IS NOT NULL)
                AND (GRCh38_gene_end IS NOT NULL)
                AND (GRCh38_gwas_snp_pos IS NOT NULL)
                ${filtersSql}
                AND (
                    (GRCh38_gene_chrom=$chromosome AND GRCh38_gene_start>=$start AND GRCh38_gene_start<=$end)
                    OR (GRCh38_gene_chrom=$chromosome AND GRCh38_gene_end>=$start AND GRCh38_gene_end<=$end)
                    OR (GRCh38_chrom=$chromosome AND GRCh38_pos>=$start AND GRCh38_pos<=$end)
                    OR (GRCh38_gwas_snp_chrom=$chromosome AND GRCh38_gwas_snp_pos>=$start AND GRCh38_gwas_snp_pos<=$end))
            `;
            const filteredWhere = templateWhere(filtersSql);
            const unfilteredWhere = templateWhere('');
            
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
                    const [geneId, variantId] = selectedId.split('-');
                    selectedSql = `COUNT(CASE ((gene_id = "${geneId}") AND (ld_snp_rsID = "${variantId}")) WHEN 1 THEN 1 ELSE null END) > 0 AS selected,`;
                    break;
                case 'variantLeadVariant':
                    const [varId, leadVariantId] = selectedId.split('-');
                    selectedSql = `COUNT(CASE ((ld_snp_rsID = "${varId}") AND (gwas_snp = "${leadVariantId}")) WHEN 1 THEN 1 ELSE null END) > 0 AS selected,`;
                    break;
                case 'leadVariantDisease':
                    const [leadVarId, diseaseId] = selectedId.split('-');
                    selectedSql = `COUNT(CASE ((gwas_snp = "${leadVarId}") AND (disease_efo_id = "${diseaseId}")) WHEN 1 THEN 1 ELSE null END) > 0 AS selected,`;
                    break;
                default:
                    break;
            }
            const orderBySql = selectedSql ? 'ORDER BY selected' : '';

            // genes
            const genesSql = `
            SELECT
                ${selectedSql}
                gene_id as id,
                gene_symbol as symbol,
                GRCH38_gene_chrom as chromosome,
                GRCh38_gene_pos as tss,
                GRCh38_gene_start as start,
                GRCh38_gene_end as end
            FROM processed
            ${unfilteredWhere}
            GROUP BY gene_id
            ${orderBySql}
            `;
            const genesQuery = db.all(genesSql, params);

            // gene locations (from Ensembl)
            const geneLocationsQuery = genesQuery.then(genes => {
                // retrieve canonical transcript info for each gene from gene table
                const geneIds = genes.map(d => d.id);
                const geneLocationsSql = `
                SELECT DISTINCT
                    gene_id as geneId,
                    description,
                    strand,
                    canonical_transcript as canonicalTranscript
                FROM gene
                WHERE
                    gene_id IN ("${geneIds.join('","')}")
                `;
                return db.all(geneLocationsSql).then(genes => {
                    // create a lookup
                    const geneLookup = {};
                    genes.forEach(d => {
                        geneLookup[d.geneId] = {
                            geneId: d.geneId,
                            description: d.description,
                            forwardStrand: (d.strand === 1),
                            canonicalTranscript: JSON.parse(d.canonicalTranscript)
                        }
                    })
                    return geneLookup;
                })
            });

            // variants
            const variantsSql = `
            SELECT
                ${selectedSql}
                ld_snp_rsID as id,
                GRCH38_chrom as chromosome,
                GRCh38_pos as position
            FROM processed
            ${unfilteredWhere}
            GROUP BY ld_snp_rsID
            ${orderBySql}
            `;
            const variantsQuery = db.all(variantsSql, params);

            // lead variants
            const leadVariantsSql = `
            SELECT
                ${selectedSql}
                gwas_snp as id,
                GRCh38_gwas_snp_chrom as chromosome,
                GRCh38_gwas_snp_pos as position
            FROM processed
            ${unfilteredWhere}
            GROUP BY gwas_snp
            ${orderBySql}
            `;
            const leadVariantsQuery = db.all(leadVariantsSql, params);

            // diseases
            const diseasesSql = `
            SELECT
                ${selectedSql}
                disease_efo_id as id,
                disease_name as name
            FROM processed
            ${unfilteredWhere}
            GROUP BY disease_efo_id
            ${orderBySql}
            `;
            const diseasesQuery = db.all(diseasesSql, params);

            // geneVariants
            const geneVariantsSql = `
            SELECT
                ${selectedSql}
                (gene_id || "-" || ld_snp_rsID) AS id,
                gene_id as geneId,
                gene_symbol as geneSymbol,
                GRCH38_gene_chrom as geneChromosome,
                GRCh38_gene_pos as geneTss,
                ld_snp_rsID as variantId,
                GRCH38_chrom as variantChromosome,
                GRCh38_pos as variantPosition,
                ot_g2v_score as otG2VScore,
                VEP as vep,
                GTEx as gtex,
                PCHiC as pchic,
                Fantom5 as fantom5,
                DHS as dhs,
                Nearest as nearest
            FROM processed
            ${filteredWhere}
            GROUP BY gene_id, ld_snp_rsID
            ${orderBySql}
            `;
            const geneVariantsQuery = db.all(geneVariantsSql, params)

            // variantLeadVariants
            const variantLeadVariantsSql = `
            SELECT
                ${selectedSql}
                (ld_snp_rsID || "-" || gwas_snp) AS id,
                ld_snp_rsID as variantId,
                GRCH38_chrom as variantChromosome,
                GRCh38_pos as variantPosition,
                gwas_snp as leadVariantId,
                GRCh38_gwas_snp_chrom as leadVariantChromosome,
                GRCh38_gwas_snp_pos as leadVariantPosition,
                r2
            FROM processed
            ${filteredWhere}
            GROUP BY ld_snp_rsID, gwas_snp
            ${orderBySql}
            `;
            const variantLeadVariantsQuery = db.all(variantLeadVariantsSql, params)

            // leadVariantDiseases
            const leadVariantDiseasesSql = `
            SELECT
                ${selectedSql}
                (gwas_snp || "-" || disease_efo_id) AS id,
                gwas_snp as leadVariantId,
                GRCh38_gwas_snp_chrom as leadVariantChromosome,
                GRCh38_gwas_snp_pos as leadVariantPosition,
                disease_efo_id as efoId,
                disease_name as efoName,
                gwas_pvalue as gwasPValue,
                gwas_odds_ratio as gwasOddsRatio,
                gwas_beta as gwasBeta,
                gwas_study as gwasStudy,
                gwas_pmid as gwasPMId,
                gwas_size as gwasSize
            FROM processed
            ${filteredWhere}
            GROUP BY gwas_snp, disease_efo_id
            ${orderBySql}
            `;
            const leadVariantDiseasesQuery = db.all(leadVariantDiseasesSql, params)

            // max gwas p-value (unfiltered call to leadVariantDiseases)
            const maxGwasPValueSql = `
            SELECT MIN(gwas_pvalue) as minGwasPValue
            FROM processed
            ${unfilteredWhere}
            `;
            const maxGwasPValueQuery = db.all(maxGwasPValueSql, params)
            .then(rows => {
                if (rows && rows.length > 0) {
                    return -Math.log10(rows[0].minGwasPValue);
                } else {
                    return Number.MAX_SAFE_INTEGER;
                }
            })


            // wait for all queries and return composite object
            return Promise.all([
                genesQuery,
                variantsQuery,
                leadVariantsQuery,
                diseasesQuery,
                geneVariantsQuery,
                variantLeadVariantsQuery,
                leadVariantDiseasesQuery,
                geneLocationsQuery,
                maxGwasPValueQuery
            ]).then(([genes, variants, leadVariants, diseases, geneVariants, variantLeadVariants, leadVariantDiseases, geneLocations, maxGwasPValue]) => {
                const genesWithLocations = genes.map(d => {
                    const geneLocation = geneLocations[d.id];
                    return {
                        ...d,
                        forwardStrand: geneLocation.forwardStrand,
                        canonicalTranscript: geneLocation.canonicalTranscript
                    }
                })
                const geneVariantsWithLocations = geneVariants.map(d => {
                    const geneLocation = geneLocations[d.geneId];
                    return {
                        ...d,
                        canonicalTranscript: geneLocation.canonicalTranscript
                    }
                })
                return {
                    genes: genesWithLocations,
                    variants,
                    leadVariants,
                    diseases,
                    geneVariants: geneVariantsWithLocations,
                    variantLeadVariants,
                    leadVariantDiseases,
                    maxGwasPValue
                }
            })
        },
        locusTable: (_, { chromosome, start, end, g2VMustHaves, g2VScore, r2, gwasPValue, selectedId, selectedType, offset, limit }) => {
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
                    const [geneId, variantId] = selectedId.split('-');
                    selectedSql = `AND (gene_id = "${geneId}" AND ld_snp_rsID = "${variantId}")`;
                    break;
                case 'variantLeadVariant':
                    const [varId, leadVariantId] = selectedId.split('-');
                    selectedSql = `AND (ld_snp_rsID = "${varId}" AND gwas_snp = "${leadVariantId}")`;
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
                ld_snp_rsID as variantId,
                GRCH38_chrom as variantChromosome,
                GRCh38_pos as variantPosition,
                gwas_snp as leadVariantId,
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
        },
        diseaseTable: (_, { efoId, offset, limit }) => {
            const params = {$efoId: efoId, $offset: offset, $limit: limit};

            // rows
            const rowsSql = `
            SELECT
                "index",
                gene_id as geneId,
                gene_symbol as geneSymbol,
                GRCH38_gene_chrom as geneChromosome,
                GRCh38_gene_pos as geneTss,
                ld_snp_rsID as variantId,
                GRCH38_chrom as variantChromosome,
                GRCh38_pos as variantPosition,
                gwas_snp as leadVariantId,
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
        }
    },
};

// merge schema with resolvers
const schema = makeExecutableSchema({ typeDefs, resolvers });

// create express app
const app = express();
app.use(cors())
app.use('/graphql', bodyParser.json(), graphqlExpress({ schema }))
app.use('/graphiql', graphiqlExpress({ endpointURL: '/graphql' }));

// start
const server = app.listen(4000, () => {
    const host = server.address().address;
    const port = server.address().port;
    console.log('Listening at http://%s:%s', host, port);
});
