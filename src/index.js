import fs from 'fs';
import path from 'path';
import express from 'express';
import cors from 'cors';
import responseTime from 'response-time';
import bodyParser from 'body-parser';
import { graphqlExpress, graphiqlExpress } from 'apollo-server-express';
import { buildSchema } from 'graphql';
import { makeExecutableSchema } from 'graphql-tools';
import sqlite3 from 'sqlite3';
import { promisify } from 'bluebird';
import { addMiddleware } from 'graphql-add-middleware';
import locus, { resolveLocus } from './resolvers/locus';
import locusTable from './resolvers/locusTable';
import diseaseTable from './resolvers/diseaseTable';

// open a connection to the database
let db = new sqlite3.Database('postgap.20180324.db', sqlite3.OPEN_READONLY, err => {
    if (err) {
        console.error(err.message)
    } else {
        console.log('Connected to POSTGAP database.')
    }
})

// sqlite3 methods do not support promises out of the box
db.all = promisify(db.all);
db.get = promisify(db.get);

// sqlite3 emits a profile event whenever sql is executed
// use this to log query execution times (where the sql is commented
// with a query type)
db.on('profile', (sql, ms) => {
    const queryTypeRegex = /SQL_QUERY_TYPE=(\w*)/
    const matches = sql.match(queryTypeRegex);
    if (matches) {
        const queryType = matches[1];
        console.log('sql-execution-time', queryType, ms)
    }
});

// total number of genes is ~20,000,
// so cache ensembl calls on an object
let geneLocationsCache = {};
const geneLocationsSql = `
-- SQL_QUERY_TYPE=loadGeneLocations
SELECT DISTINCT
    gene_id as geneId,
    description,
    strand,
    canonical_transcript as canonicalTranscript
FROM gene
`;
db.all(geneLocationsSql).then(genes => {
    // create a lookup cache
    genes.forEach(d => {
        geneLocationsCache[d.geneId] = {
            geneId: d.geneId,
            description: d.description,
            forwardStrand: (d.strand === 1),
            canonicalTranscript: JSON.parse(d.canonicalTranscript)
        }
    });
});

// load the schema
const schemaFile = path.join(__dirname, 'schema.gql');
const typeDefs = fs.readFileSync(schemaFile, 'utf8');

// specify the resolution methods for allowed query
const resolvers = {
    Query: {
        locus: resolveLocus,
        locusTable,
        diseaseTable
    },
    DrawableLocus: locus
};

// merge schema with resolvers
const schema = makeExecutableSchema({ typeDefs, resolvers });


// setup cors options
const corsOptions = {
    preflightContinue: false
};

// context for resolvers
const context = {
    db,
    geneLocationsCache
};

// create express app
const app = express();
app.use('/graphql', responseTime(), cors(corsOptions), bodyParser.json(), graphqlExpress({ schema, context }))
app.use('/graphiql', graphiqlExpress({ endpointURL: '/graphql' }));

// start
const server = app.listen(4000, () => {
    const host = server.address().address;
    const port = server.address().port;
    console.log('Listening at http://%s:%s', host, port);
});
