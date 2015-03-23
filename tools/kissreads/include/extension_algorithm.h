/*****************************************************************************
 *   discoSnp++: discovering polymorphism from raw unassembled NGS reads
 *   A tool from the GATB (Genome Assembly Tool Box)
 *   Copyright (C) 2014  INRIA
 *   Authors: P.Peterlongo, E.Drezen
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

/*
 * extension_algorithm.h
 *
 *  Created on: 15 sept. 2010
 *      Author: ppeterlo
 */

#ifndef EXTENSION_ALGORITHM_H_
#define EXTENSION_ALGORITHM_H_

#include<fragment_info.h>
#include<stdio.h>
#include <zlib.h>



float read_mapping (char ** read_file_names,
                      int read_set_id,
                      const int k,
                      const int min_coverage,
                      p_fragment_info * all_starters,
                      int qual,
                      FILE * sam_out,
                      int subst_allowed,
                      const int minimal_read_overlap
                      );


void set_read_coherency(p_fragment_info * all_starters, const int nb_events_per_set, const int nb_fragment_per_event, const char paired, const int min_coverage, const int read_set_id);
#endif /* EXTENSION_ALGORITHM_H_ */
