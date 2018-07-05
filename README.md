# POSTGAP API
This is currently an experimental application to build a GraphQL API to serve POSTGAP data specifically for the POSTGAP web application. The main motivations are to improve performance and to enable more flexible queries.

## Development

### Build container with google cloudbuilder

To use the google cloud builder...
```sh
gcloud auth login
gcloud config set project open-targets-eu-dev
gcloud container builds submit --config cloudbuild.yaml .
```


### Build db step

Setup python requirements using [pipenv](https://docs.pipenv.org/):
```sh
pip3 install --user pipenv #or `brew install pipenv` or `apt install pipenv`
pipenv install
```

Transform the data to meet the Open Targets requirements (outputs `<input_file>.transformed.gz`):
```
pipenv run python data-prep/transform.py --sample
```
Note: the command above will download a 1GB file and might take a while to execute.
If you want to run it on a local file, use instead:
```sh
pipenv run python data-prep/transform.py -f <postgapfile.<version>.txt.gz>
```

Build the `sqlite3` database (outputs `postgap.<dateversion>.db`). You may need to `brew install mysql` first.
```
pipenv run python data-prep/build.py /data-prep/postgap.<version>.txt.transformed.gz
```

If this step is succesful you should obtain a `postgap.<version>.db` file in
your root directory. Now just symlink it to `postgap.db` for everything to work:
```sh
ln -s postgap.<version>.db postgap.db
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

## Deployment
To run an image built locally in the cloud, first push to google cloud:
```
# run once
gcloud auth configure-docker
# then
docker tag <LOCAL_IMAGE>:<LOCAL_TAG> eu.gcr.io/open-targets-eu-dev/postgap-api:<TAG_NAME>
gcloud config set project opent-targets-eu-dev
docker push eu.gcr.io/open-targets-eu-dev/postgap-api:<TAG_NAME>
```

Update the tag in `postgap-api.yml` to `TAG_NAME`.

There is a kubernetes cluster set up on google cloud under the project `open-targets-eu-dev`. Replace the existing development cluster. To deploy the new image, run the following:
```
kubectl delete deploy/postgap-api
kubectl create -f postgap-api.yml
```

**NOTE**: This deployment process is subject to change in the near future. It should eventually be triggered by Circle CI.

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
