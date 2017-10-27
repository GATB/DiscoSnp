#!/bin/bash
#*****************************************************************************
#   discoSnpRad: discovering polymorphism from raw unassembled RADSEQ NGS reads
#   A tool from the GATB (Genome Assembly Tool Box)
#   Copyright (C) 2014  INRIA
#   Authors: J. Gauthier, C. Mouden, C. Lemaitre, P. Peterlongo
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


EDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
function myrealpath { echo $(cd $(dirname $1); pwd)/$(basename $1); }
cmd="$EDIR/run_discoSnp++.sh "$@" -x -t -e"
echo "I run discoSnp with following command line: " ${cmd}
${cmd}

#TODO: call automatically the clustering +radseq filters. 




