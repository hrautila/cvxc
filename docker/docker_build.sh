#!/bin/bash

version=$(git describe | sed 's/^[a-z]*//')
base=ubuntu
name=cvxc

if [ "$1" != "" ]; then
    base=$1
fi

DOCKERFILE="Dockerfile.${base}"
if [ ! -f $DOCKERFILE ]; then
    echo "Error: $DOCKERFILE not found" >&2
    exit 1
fi

TAG="${name}:${version}"

docker build -f $DOCKERFILE -t $TAG .
docker tag $TAG ${name}:latest

