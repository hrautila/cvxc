#!/bin/bash

images=$(docker images | awk '$1 ~ /cvxc/ {printf "%s:%s\n", $1, $2}')

if [ "$images" = "" ]; then
    echo "No cvxc images. Pull one!" >&2
    exit 1
fi

docker_image=$(echo "$images" | tr ' ' '\n' | awk -F: '$2 ~ /latest/ {print $0}')
if [ "$docker_image" = "" ]; then
    docker_image=$(echo "$images" | tr ' ' '\n' | head -1)
fi

istty="-t"
if [ ! -t 0 ]; then
    istty=
fi

docker run --rm -i $istty --volume $PWD:/var/run $docker_image $*

