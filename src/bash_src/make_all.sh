#!/bin/bash

cd ..
module load gcc/8.2.0 openmpi
make distributed distributed2 hybrid
make sequential multithreaded