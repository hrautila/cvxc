#!/bin/bash

version=0.1
base=alpine
name=cvxc

if [ "$1" != "" ]; then
    base=$1
fi

DOCKERFILE="Dockerfile.${base}"
if [ ! -f $DOCKERFILE ]; then
    echo "Error: $DOCKERFILE not found" >&2
    exit 1
fi

TAG="${name}:${version}-${base}"

docker build -f $DOCKERFILE -t $TAG .

