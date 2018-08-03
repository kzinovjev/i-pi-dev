#!/bin/bash

rm pids
rm -rf /tmp/ipi*
rm -rf results
mkdir results
cd results
mpiexec -n 8 ../walker.sh &
echo $! >> ../pids
