#!/bin/bash

#docker build -f Dockerfile -t ligero/ceremony/cmake .
docker build -f Dockerfile.build -t ligero:build .
