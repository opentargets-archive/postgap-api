import fs from 'fs';
import path from 'path';
import express from 'express';
import graphqlHTTP from 'express-graphql';
import { buildSchema } from 'graphql';
import { makeExecutableSchema } from 'graphql-tools';
import sqlite3 from 'sqlite3';
import {promisify} from 'bluebird';

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

            // wait for all queries and return composite object
            return Promise.all([
                genesQuery,
                variantsQuery,
                diseasesQuery
            ]).then(([genes, variants, diseases]) => {
                return {
                    genes,
                    variants,
                    diseases
                }
            })
        }
    },
};

// merge schema with resolvers
const schema = makeExecutableSchema({ typeDefs, resolvers });

const app = express();
app.use('/', graphqlHTTP(req => {
    const startTime = Date.now();
    return {
        schema,
        graphiql: true,
        extensions: ({ document, variables, operationName, result }) => ({
            timeTaken: Date.now() - startTime,
        })
    };
}));
app.listen(4000);
