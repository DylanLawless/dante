#!/bin/bash

cd ../data

find . -type f -name '*tsv' -exec wc -l {} + | sed 's/[[:blank:]]\+/ /g' | sort -k2 > count_data.md

