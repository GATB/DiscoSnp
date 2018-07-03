#! /bin/bash

grep ^# $1 > clustersup$2_$1
grep -v ^# $1 | awk '{print $1"_"$0}' | awk -F '_' '{if ($4 <= '"$2"') print $0}' | cut -d "_" -f 9- >> clustersup$2_$1


