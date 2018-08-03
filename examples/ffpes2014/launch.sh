#!/bin/bash
IPI=/ssd/i-pi-ghts/bin/i-pi

rm -rf results
mkdir results
cd results
cp ../metoh.xml xml

$IPI xml
