#!/usr/bin/env bash
kubectl delete deploy/postgapapi
kubectl run postgapapi --image=eu.gcr.io/open-targets-eu-dev/postgap-api:$1 --port=4000

