#!/bin/bash

docker run -v `pwd`:/opt/src --network="host" ligero/ceremony/cmake ./scripts/rebuild-and-run-tests.sh
