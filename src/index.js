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
        locus: (_, { chromosome, start, end }) => {
            const params = {$chromosome: chromosome, $start: start, $end: end};
            const oldWhere = `
            WHERE
                (GRCh38_gene_chrom=$chromosome AND GRCh38_gene_pos>=$start AND GRCh38_gene_pos<=$end)
                OR (GRCh38_chrom=$chromosome AND GRCh38_pos>=$start AND GRCh38_pos<=$end)
            `;
            const commonWhere = `
            WHERE
                (GRCh38_gene_start IS NOT NULL)
                AND (GRCh38_gene_end IS NOT NULL)
                AND (GRCh38_gwas_snp_pos IS NOT NULL)
                AND (
                    (GRCh38_gene_chrom=$chromosome AND GRCh38_gene_start>=$start AND GRCh38_gene_start<=$end)
                    OR (GRCh38_gene_chrom=$chromosome AND GRCh38_gene_end>=$start AND GRCh38_gene_end<=$end)
                    OR (GRCh38_chrom=$chromosome AND GRCh38_pos>=$start AND GRCh38_pos<=$end)
                    OR (GRCh38_gwas_snp_chrom=$chromosome AND GRCh38_gwas_snp_pos>=$start AND GRCh38_gwas_snp_pos<=$end))
            `;

            // genes
            const genesSql = `
            SELECT DISTINCT
                gene_id as id,
                gene_symbol as symbol,
                GRCH38_gene_chrom as chromosome,
                GRCh38_gene_pos as tss,
                GRCh38_gene_start as start,
                GRCh38_gene_end as end
            FROM processed
            ${commonWhere}
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
            SELECT DISTINCT
                ld_snp_rsID as id,
                GRCH38_chrom as chromosome,
                GRCh38_pos as position
            FROM processed
            ${commonWhere}
            `;
            const variantsQuery = db.all(variantsSql, params);

            // lead variants
            const leadVariantsSql = `
            SELECT DISTINCT
                gwas_snp as id,
                GRCh38_gwas_snp_chrom as chromosome,
                GRCh38_gwas_snp_pos as position
            FROM processed
            ${commonWhere}
            `;
            const leadVariantsQuery = db.all(leadVariantsSql, params);

            // diseases
            const diseasesSql = `
            SELECT DISTINCT
                disease_efo_id as id,
                disease_name as name
            FROM processed
            ${commonWhere}
            `;
            const diseasesQuery = db.all(diseasesSql, params);

            // geneVariants
            const geneVariantsSql = `
            SELECT DISTINCT
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
            ${commonWhere}
            `;
            const geneVariantsQuery = db.all(geneVariantsSql, params)

            // variantLeadVariants
            const variantLeadVariantsSql = `
            SELECT DISTINCT
                (ld_snp_rsID || "-" || gwas_snp) AS id,
                ld_snp_rsID as variantId,
                GRCH38_chrom as variantChromosome,
                GRCh38_pos as variantPosition,
                gwas_snp as leadVariantId,
                GRCh38_gwas_snp_chrom as leadVariantChromosome,
                GRCh38_gwas_snp_pos as leadVariantPosition,
                r2
            FROM processed
            ${commonWhere}
            `;
            const variantLeadVariantsQuery = db.all(variantLeadVariantsSql, params)

            // leadVariantDiseases
            const leadVariantDiseasesSql = `
            SELECT DISTINCT
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
            ${commonWhere}
            `;
            const leadVariantDiseasesQuery = db.all(leadVariantDiseasesSql, params)


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
            ]).then(([genes, variants, leadVariants, diseases, geneVariants, variantLeadVariants, leadVariantDiseases, geneLocations]) => {
                console.log(geneVariants)
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
                    leadVariantDiseases
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
