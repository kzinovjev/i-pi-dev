#!/bin/bash
IPI= #path to i-pi
DRIVER= #path to sander_driver

rm -rf results
mkdir results
cd results

$IPI ../ecDHFR.xml &
sleep 5
for i in `seq 1 8`;do
    $DRIVER ../ecDHFR.in ../ecDHFR.prmtop ../ecDHFR.rst7 localhost 21142 0 &
done
