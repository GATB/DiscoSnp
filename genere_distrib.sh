#!/bin/sh
#************

# 3 arguments for the version number such that X1.X2.X3

if [ "$#" -ne 3 ]; then
echo "requires 3 arguments for the version number such that X1.X2.X3"
exit 0
fi

rm -rf build
mkdir build
cd build/
cmake -DMAJOR=$1 -DMINOR=$2 -DPATCH=$3 ..


make delivery -j 2

## Pour tester ce qu'il y a dans les archives avant d'envoyer sur la forge
#make package -j 2
#make package_source -j 2
