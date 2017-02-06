#!/bin/sh

echo 'Running Navier-stokes Solver'
make -C source/

source/./result

rm source/result