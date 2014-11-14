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
 * couple.h
 *
 *  Created on: 17 sept. 2010
 *      Author: ppeterlo
 */
//#include<tree.h>
#ifndef COUPLE_H_
#define COUPLE_H_

typedef struct {
	int a;
	int b;
} couple, * p_couple;


typedef struct {
	char* a;
	int b;
} charint_couple, * p_charint_couple;

typedef struct {
	void *a;
	void *b;
} pointers_couple, * p_pointers_couple;


//typedef struct {
//	p_node t;
//	int b;
//} node_couple, * p_node_couple; // used to store the info of a seed maping one position (b) of a fragment
//

p_couple create_couple(int a,int b);
void free_couple (const void * v_c);

p_charint_couple create_charint_couple(char* a,int b);
void free_charint_couple (const void * v_c);

p_pointers_couple create_pointers_couple(void* a, void * b);
void free_pointers_couple (const void * v_c);
//
////p_node_couple create_node_couple(p_node t, int pos);//, int fragment_id);
//void free_node_couple( const void * v_sc);

#endif /* COUPLE_H_ */
