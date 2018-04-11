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

import ensemblClient from './ensemblClient';

// open a connection to the database
let db = new sqlite3.Database('postgap.db', sqlite3.OPEN_READONLY, err => {
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

            // genes
            const genesSql = `
            SELECT DISTINCT
                gene_id as id,
                gene_symbol as symbol,
                GRCH38_gene_chrom as chromosome,
                GRCh38_gene_pos as tss
            FROM raw
            WHERE
                (GRCh38_gene_chrom=$chromosome AND GRCh38_gene_pos>=$start AND GRCh38_gene_pos<=$end)
                OR (GRCh38_chrom=$chromosome AND GRCh38_pos>=$start AND GRCh38_pos<=$end)
            `;
            const genesQuery = db.all(genesSql, params);

            // gene locations (from Ensembl)
            const geneLocationsQuery = genesQuery.then(genes => {
                // retrieve canonical transcript for each gene from Ensembl API
                const geneIds = genes.map(d => d.id);

                // look up in cache object
                const cachedGeneIds = geneIds.filter(d => geneCache[d]);
                const notCachedGeneIds = geneIds.filter(d => !geneCache[d]);
                console.log(`geneCache: ${cachedGeneIds.length} cached, ${notCachedGeneIds.length} uncached`)
                const cachedGenes = {}
                cachedGeneIds.forEach(d => {
                    cachedGenes[d] = geneCache[d];
                });

                if (notCachedGeneIds.length === 0) {
                    // all in the cache
                    return cachedGenes;
                }

                return ensemblClient.fetchGenes(notCachedGeneIds)
                .then(notCachedGenes => {
                    // update cache
                    geneCache = {
                        ...geneCache,
                        ...notCachedGenes
                    }

                    // return merged
                    return {
                        ...cachedGenes,
                        ...notCachedGenes
                    }
                });
            });

            // variants
            const variantsSql = `
            SELECT DISTINCT
                ld_snp_rsID as id,
                GRCH38_chrom as chromosome,
                GRCh38_pos as position
            FROM raw
            WHERE
                (GRCh38_gene_chrom=$chromosome AND GRCh38_gene_pos>=$start AND GRCh38_gene_pos<=$end)
                OR (GRCh38_chrom=$chromosome AND GRCh38_pos>=$start AND GRCh38_pos<=$end)
            `;
            const variantsQuery = db.all(variantsSql, params);

            // lead variants
            const leadVariantsSql = `
            SELECT DISTINCT
                gwas_snp as id
            FROM raw
            WHERE
                (GRCh38_gene_chrom=$chromosome AND GRCh38_gene_pos>=$start AND GRCh38_gene_pos<=$end)
                OR (GRCh38_chrom=$chromosome AND GRCh38_pos>=$start AND GRCh38_pos<=$end)
            `;
            const leadVariantsQuery = db.all(leadVariantsSql, params);

            // lead variant locations (from Ensembl)
            const leadVariantsLocationsQuery = leadVariantsQuery.then(leadVariants => {
                // retrieve chromosome/position for each lead variant from Ensembl API
                const leadVariantIds = leadVariants.map(d => d.id);

                // look up in cache object
                const cachedLeadVariantIds = leadVariantIds.filter(d => leadVariantCache[d]);
                const notCachedLeadVariantIds = leadVariantIds.filter(d => !leadVariantCache[d]);
                console.log(`leadVariantCache: ${cachedLeadVariantIds.length} cached, ${notCachedLeadVariantIds.length} uncached`)
                const cachedLeadVariants = {}
                cachedLeadVariantIds.forEach(d => {
                    cachedLeadVariants[d] = leadVariantCache[d];
                });

                if (notCachedLeadVariantIds.length === 0) {
                    // all in the cache
                    return cachedLeadVariants;
                }

                return ensemblClient.fetchVariants(notCachedLeadVariantIds)
                .then(notCachedLeadVariants => {
                    // update cache
                    leadVariantCache = {
                        ...leadVariantCache,
                        ...notCachedLeadVariants
                    }

                    // return merged
                    return {
                        ...cachedLeadVariants,
                        ...notCachedLeadVariants
                    }
                });
            });

            // diseases
            const diseasesSql = `
            SELECT DISTINCT
                disease_efo_id as id,
                disease_name as name
            FROM raw
            WHERE
                (GRCh38_gene_chrom=$chromosome AND GRCh38_gene_pos>=$start AND GRCh38_gene_pos<=$end)
                OR (GRCh38_chrom=$chromosome AND GRCh38_pos>=$start AND GRCh38_pos<=$end)
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
                VEP as vep,
                GTEx as gtex,
                PCHiC as pchic,
                Fantom5 as fantom5,
                DHS as dhs,
                Nearest as nearest
            FROM raw
            WHERE
                (GRCh38_gene_chrom=$chromosome AND GRCh38_gene_pos>=$start AND GRCh38_gene_pos<=$end)
                OR (GRCh38_chrom=$chromosome AND GRCh38_pos>=$start AND GRCh38_pos<=$end)
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
                r2
            FROM raw
            WHERE
                (GRCh38_gene_chrom=$chromosome AND GRCh38_gene_pos>=$start AND GRCh38_gene_pos<=$end)
                OR (GRCh38_chrom=$chromosome AND GRCh38_pos>=$start AND GRCh38_pos<=$end)
            `;
            const variantLeadVariantsQuery = db.all(variantLeadVariantsSql, params)

            // leadVariantDiseases
            const leadVariantDiseasesSql = `
            SELECT DISTINCT
                (gwas_snp || "-" || disease_efo_id) AS id,
                gwas_snp as leadVariantId,
                disease_efo_id as efoId,
                disease_name as efoName,
                gwas_pvalue as gwasPValue,
                gwas_odds_ratio as gwasOddsRatio,
                gwas_beta as gwasBeta,
                gwas_study as gwasStudy,
                gwas_pmid as gwasPMId,
                gwas_size as gwasSize
            FROM raw
            WHERE
                (GRCh38_gene_chrom=$chromosome AND GRCh38_gene_pos>=$start AND GRCh38_gene_pos<=$end)
                OR (GRCh38_chrom=$chromosome AND GRCh38_pos>=$start AND GRCh38_pos<=$end)
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
                leadVariantsLocationsQuery
            ]).then(([genes, variants, leadVariants, diseases, geneVariants, variantLeadVariants, leadVariantDiseases, geneLocations, leadVariantsLocations]) => {
                const genesWithLocations = genes.map(d => {
                    const geneLocation = geneLocations[d.id];
                    return {
                        ...d,
                        start: geneLocation.start,
                        end: geneLocation.end,
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
                const leadVariantsWithLocations = leadVariants.map(d => {
                    const leadVariantLocation = leadVariantsLocations[d.id];
                    return {
                        ...d,
                        chromosome: leadVariantLocation.chromosome,
                        position: leadVariantLocation.position
                    }
                })
                const variantLeadVariantsWithLocations = variantLeadVariants.map(d => {
                    const leadVariantLocation = leadVariantsLocations[d.leadVariantId];
                    return {
                        ...d,
                        leadVariantChromosome: leadVariantLocation.chromosome,
                        leadVariantPosition: leadVariantLocation.position
                    }
                })
                const leadVariantDiseasesWithLocations = leadVariantDiseases.map(d => {
                    const leadVariantLocation = leadVariantsLocations[d.leadVariantId];
                    return {
                        ...d,
                        leadVariantChromosome: leadVariantLocation.chromosome,
                        leadVariantPosition: leadVariantLocation.position
                    }
                })
                
                return {
                    genes: genesWithLocations,
                    variants,
                    leadVariants: leadVariantsWithLocations,
                    diseases,
                    geneVariants: geneVariantsWithLocations,
                    variantLeadVariants: variantLeadVariantsWithLocations,
                    leadVariantDiseases: leadVariantDiseasesWithLocations
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
            FROM raw
            WHERE disease_efo_id=$efoId
            LIMIT $limit
            OFFSET $offset
            `;
            const rowsQuery = db.all(rowsSql, params);

            // total rows
            const totalRowsSql = `
            SELECT COUNT(*) as total
            FROM raw
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
