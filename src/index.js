import express from 'express';
import graphqlHTTP from 'express-graphql';
import { buildSchema } from 'graphql';
import sqlite3 from 'sqlite3';
import {promisify} from 'bluebird';

let db = new sqlite3.Database('postgap.db', sqlite3.OPEN_READONLY, err => {
    if (err) {
        console.error(err.message)
    } else {
        console.log('Connected to POSTGAP database.')
    }
})

db.all = promisify(db.all);

const books = [
    {
        title: 'Harry Potter',
        author: 'JK Rowling'
    },
    {
        title: 'Jurassic Park',
        author: 'Michael Crichton'
    }
];

const schema = buildSchema(`
type Query {
    hello: String
    books: [Book]
    locus(chromosome: String, start: Int, end: Int): [String]
}
type Book {
    title: String
    author: String
}
`);

const rootValue = {
    hello: () => 'Hello World!',
    books: () => books,
    locus: ({ chromosome, start, end }) => {
        const sql = `SELECT * FROM raw WHERE (GRCh38_gene_chrom=? AND GRCh38_gene_pos>=? AND GRCh38_gene_pos<=?) OR (GRCh38_chrom=? AND GRCh38_pos>=? AND GRCh38_pos<=?)`;
        const params = [chromosome, start, end, chromosome, start, end];
        return db.all(sql,params).then(rows => {
            // process rows here
            return [`${rows.length}`];
        });
    }
};

const app = express();
app.use('/', graphqlHTTP({
    schema,
    rootValue,
    graphiql: true
}));
app.listen(4000);
