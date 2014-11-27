/*****************************************************************************
 *   DiscoMore: discovering polymorphism from raw unassembled NGS reads
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
 * outputs.h
 *
 *  Created on: 27 oct. 2010
 *      Author: ppeterlo
 */

#ifndef OUTPUTS_H_
#define OUTPUTS_H_
void print_results_2_paths_per_event(FILE * coherent_out, FILE * uncoherent_out,   const p_fragment_info * results_against_set, int number_of_read_sets, int nb_events_per_set, int qual);
void print_results_invs(FILE * coherent_out, FILE * uncoherent_out,   const p_fragment_info * results_against_set, int number_of_read_sets, int nb_events_per_set, int qual);
void print_generic_results(FILE * coherent_out, FILE * uncoherent_out,  const p_fragment_info * results_against_set, const int number_of_read_sets, int nb_events_per_set, int qual);
#endif /* OUTPUTS_H_ */
