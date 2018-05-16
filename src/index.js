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
import morgan from 'morgan';

import locus, { resolveLocus } from './resolvers/locus';
import locusTable from './resolvers/locusTable';
import diseaseTable from './resolvers/diseaseTable';
import geneLocation from './resolvers/geneLocation';
import variantLocation from './resolvers/variantLocation';
import search from './resolvers/search';

const DB_FILENAME = 'postgap.20180514.db';

// open a connection to an in-memory database
let db = new sqlite3.Database(':memory:', sqlite3.OPEN_READWRITE, err => {
    if (err) {
        console.error(err.message)
    } else {
        // sqlite3 methods do not support promises out of the box
        db.all = promisify(db.all);
        db.get = promisify(db.get);
        db.exec = promisify(db.exec);
        db.run = promisify(db.run);

        // temporarily attach the file db and copy into memory
        db.run(`ATTACH DATABASE "${DB_FILENAME}" AS filedb`)
        .then(() => console.log('Attached file db.'))
        .then(() => {
            return db.all('SELECT type,name,sql FROM filedb.sqlite_master;').then(data => {
                const buildSql = data.map(d => d.sql).join(';');
                const tables = data.filter(d => d.type === 'table')
                    // .filter(d => ['gene', 'lead_variant', 'chr_7', 'processed'].indexOf(d.name) >= 0)
                    .map(d => d.name);

                const copyTables = () => Promise.all(
                    tables.map(table => {
                        return db.exec(`INSERT INTO ${table} SELECT * FROM filedb.${table};`)
                            .then(() => console.log(`Copied all rows of ${table} into memory.`))
                    })
                );
                return db.exec(buildSql).then(() => console.log('Created tables in memory.'))
                .then(copyTables);
            })
        })
        .then(() => {
            return db.run(`DETACH DATABASE filedb`);
        })
        .then(() => console.log('Detached file db.'))
        ;
    }
});

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

// cache gene locations on object
let geneLocationsCache = JSON.parse(fs.readFileSync('genes.json','utf8'));


// load the schema
const schemaFile = path.join(__dirname, 'schema.gql');
const typeDefs = fs.readFileSync(schemaFile, 'utf8');

// specify the resolution methods for allowed query
const resolvers = {
    Query: {
        locus: resolveLocus,
        locusTable,
        diseaseTable,
        geneLocation,
        variantLocation,
        search
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
app.use('/graphql', morgan('tiny'), responseTime(), cors(corsOptions), bodyParser.json(), graphqlExpress({ schema, context }))
app.use('/graphiql', graphiqlExpress({ endpointURL: '/graphql' }));
app.get('/', (req, res) => res.send('All is well at the root!'))


// start
const server = app.listen(4000, () => {
    const host = server.address().address;
    const port = server.address().port;
    console.log('Listening at http://%s:%s', host, port);
});
