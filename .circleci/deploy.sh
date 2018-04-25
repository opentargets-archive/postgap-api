gcloud compute ssl-certificates create wild-opentargets-io --certificate '*.opentargets.io.chained.crt' --private-key '*.opentargets.io.key' --project open-targets-eu-dev

kubectl run postgapapi --image=eu.gcr.io/open-targets-eu-dev/postgap-api:24052018 --port=4000

kubectl expose deployment postgapapi --target-port=4000 --type=NodePort

kubectl create -f ingress.yml