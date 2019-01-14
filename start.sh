#!/bin/bash

mpicxx main.cpp
mpirun -np 4 ./a.out

