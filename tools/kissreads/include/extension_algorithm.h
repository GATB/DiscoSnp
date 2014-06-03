/**
 * Copyright INRIA , contributors Peterlongo
 * pierre.peterlongo@inria.fr
 *
 *
 * This software is a computer program whose purpose is to detect the
 * presence of a sequence in a set of NGS reads, and to compute its average quality and coverage
 *
 * This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software.  You can  use,
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability.
 
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or
 * data to be ensured and,  more generally, to use and operate it in the
 * same conditions as regards security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 */

/*
 * extension_algorithm.h
 *
 *  Created on: 15 sept. 2010
 *      Author: ppeterlo
 */

#ifndef EXTENSION_ALGORITHM_H_
#define EXTENSION_ALGORITHM_H_

#include <extending_fragment.h>
#include<fragment_info.h>
#include<stdio.h>



float read_coherence (FILE * reads_file,
                      char * reads_file_name,
                      const int k,
                      const int min_coverage,
                      p_fragment_info * all_starters,
                      const int number_starters,
                      int read_file,
                      int qual,
                      int nb_events_per_set,
                      int nb_fragment_per_event,
                      FILE * sam_out,
                      int subst_allowed,
                      const char no_subsutitution_on_central_position);

#endif /* EXTENSION_ALGORITHM_H_ */
