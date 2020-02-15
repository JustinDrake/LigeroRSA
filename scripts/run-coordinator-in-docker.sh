#!/bin/bash

docker run --network="host" -it ligero:build ./coordinator_full_protocol "$@"
