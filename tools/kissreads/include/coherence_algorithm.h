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
 * coherence_algorithm.h
 *
 *  Created on: 15 sept. 2010
 *      Author: peterlongo chikhi
 */

#ifndef COHERENCE_ALGORITHM_H_
#define COHERENCE_ALGORITHM_H_


char read_coherent_generic (const int pwi, const char * fragment, const char * read, const int subst_allowed);
char read_coherent_SNP(const int pwi, const char * fragment, const char * read, const int subst_allowed, const char * SNP_position);

#endif /* COHERENCE_ALGORITHM_H_ */
