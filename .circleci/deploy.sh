#!/usr/bin/env bash
kubectl delete deploy/postgapapi
kubectl run postgapapi --image=eu.gcr.io/open-targets-eu-dev/postgap-api:20180426sample --port=4000

