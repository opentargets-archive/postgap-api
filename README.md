# POSTGAP API
This is currently an experimental application to build a GraphQL API to serve POSTGAP data specifically for the POSTGAP web application. The main motivations are to improve performance and to enamble more flexible queries.

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
