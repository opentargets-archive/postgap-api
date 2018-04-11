# POSTGAP API
This is currently an experimental application to build a GraphQL API to serve POSTGAP data specifically for the POSTGAP web application. The main motivations are to improve performance and to enable more flexible queries.

## Usage
### Build db step
Assumption: You have downloaded the latest POSTGAP output file (eg `postgap.20180324.txt.gz`) to the root of a clone of this repo.

TODO: Setup requirements.txt and add instructions here

Transform the data to meet the Open Targets requirements (outputs `<input_file>.transformed`):
```
python3 ./data-prep/transform.py postgap.20180324.txt.gz
```

Build the `sqlite3` database (outputs `postgap.db`):
```
python3 ./data-prep/build.py postgap.20180324.txt.gz.transformed
```

### GraphQL server
Assumption: You have a recent `node` environment set up and `yarn` installed.

To run the GraphQL server in development mode:
```
yarn install
yarn run start
```

### Docker
You can build the GraphQL server as a docker image.
```
docker build -t <name>:<version> .
```

To run, you might then do the following.
```
docker run -p 4000:4000 -d --name=<container-name> <name>:<version>
```


## Inputs
This application should build on a POSTGAP flat file. This primary file should then be decorated with the following additional columns of data.
Essential:
* OT G2V score, which is calculated by a known algorithm.
* GR38 chromosome and position of lead variant
* GR38 start, end and strand of each gene's canonical transcript
* GR38 exons (start and end) of each gene's canonical transcript (expressed as a JSON string)
Stretch goals:
* GR38 PCHiC, Fantom5, DHS raw data blocks (transformed from BED files)
* GTEx tissue expression values per gene

## Outputs
This application should be run as a standalone docker container.

## Queries
See the GraphQL schema for full details, but the main queries are:
* Get rows specific to a disease (indexed by `efoId`). Allow pagination and sorting by any field.
* Get rows specific to a locus (indexed by `chromosome`, `start`, `end`). Allow pagination and sorting by any field. Allow filtering by connection values (`otG2VScore`, `VEP`, `r2`, `gwasPValue`, ...)
* Get drawable entities specific to a locus (indexed by `chromosome`, `start`, `end`). Allow filtering by connection values (`otG2VScore`, `VEP`, `r2`, `gwasPValue`, ...). The drawable entities are:
  * Unique genes, variants, lead variants and diseases within the locus (or linked to from another entity within the locus).
  * Unique gene-variants, variant-lead variants, lead variant-diseases within the locus (or linked to from another entity within the locus).
