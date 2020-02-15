#!/bin/bash

docker run --network="host" -it ligero:build ./party_full_protocol "$@"
