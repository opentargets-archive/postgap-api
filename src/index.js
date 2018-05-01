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
import locus from './resolvers/locus';
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
        locus,
        locusTable,
        diseaseTable
    },
};

// merge schema with resolvers
const schema = makeExecutableSchema({ typeDefs, resolvers });

// timing helper (per resolver)
const timeResolver = async (rootValue, args, context, info, next) => {
    const start = Date.now();
    const result = await next();
    const end = Date.now();
    console.log(`Query.locus took: ${end - start}ms`);
    return result;
};

// time specific resolver calculations
addMiddleware(schema, 'Query.locus', timeResolver);
addMiddleware(schema, 'Query.locus.genes', timeResolver);

// setup cors options
const corsOptions = {
    preflightContinue: false
};

// context for resolvers
const context = {
    db
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
