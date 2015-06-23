#!/bin/bash

# Input should be the STAR *.Log.final.out file
input=$1

grep '|' $input | \
sed 's/^\s\+//g' | \
sed 's/\s\+|\s\+/\t/g'
