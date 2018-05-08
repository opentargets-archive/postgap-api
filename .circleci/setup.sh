#!/usr/bin/env bash
kubectl expose deployment postgapapi --target-port=4000 --type=NodePort
kubectl create -f ingress.yml