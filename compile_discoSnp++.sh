#*****************************************************************************
#   DiscoMore: discovering polymorphism from raw unassembled NGS reads
#   A tool from the GATB (Genome Assembly Tool Box)
#   Copyright (C) 2014  INRIA
#   Authors: P.Peterlongo, E.Drezen
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as
#  published by the Free Software Foundation, either version 3 of the
#  License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#*****************************************************************************

#!/bin/bash


MY_PATH="`( cd \"$MY_PATH\" && pwd )`"  # absolutized and normalized
if [ -z "$MY_PATH" ] ; then
# error; for some reason, the path is not accessible
# to the script (e.g. permissions re-evaled after suid)
exit 1  # fail
fi
echo -e "\tCompiling discoSnp, from directory "$MY_PATH


if [ -f $MY_PATH/build/Makefile ]; then
echo -e File $MY_PATH/build/Makefile exists, we use it
cd build
else
mkdir build
cd build
cmake ..
fi

make -j8
if [ $? -eq 0 ]
then
echo "Compilation done, have fun! "

exit
fi

cd ..
