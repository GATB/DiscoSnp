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
 * couple.c
 *
 *  Created on: 25 oct. 2010
 *      Author: ppeterlo
 */
#include <couple.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <commons.h>
p_couple create_couple(int a,int b){
	p_couple cpl = malloc(sizeof(couple));
	test_alloc(cpl);
	cpl->a=a;
	cpl->b=b;
	return cpl;
}


void free_couple (const void * v_c){
	couple * c = (couple *) v_c;
	if(c != NULL){
		free(c);
	}
}


p_charint_couple create_charint_couple(char* a,int b){
	p_charint_couple cpl = malloc(sizeof(charint_couple));
	test_alloc(cpl);
	cpl->a=a;
	cpl->b=b;
	return cpl;
}


void free_charint_couple (const void * v_c){
	charint_couple * c = (charint_couple *) v_c;
	if(c != NULL){
		free(c);
	}
}


p_pointers_couple create_pointers_couple(void* a, void * b){
	p_pointers_couple cpl = malloc(sizeof(pointers_couple));
	test_alloc(cpl);
	cpl->a=a;
	cpl->b=b;
	return cpl;
}

void free_pointers_couple (const void * v_c){
	pointers_couple * c = (pointers_couple *) v_c;
	if(c != NULL){
		free(c);
	}
}
//
//p_node_couple create_node_couple(p_node t, int pos){//, int fragment_id){
//	p_node_couple pnc = (p_node_couple) malloc(sizeof(node_couple)); test_alloc(pnc);
//	pnc->t=t;
//	pnc->b=pos;
////	pnc->fragment_id=fragment_id;
//	return pnc;
//}
//
//void free_node_couple( const void * v_nc){
//	node_couple * nc = (node_couple *) v_nc;
//	if(nc != NULL) {
//		free(nc);
//	}
//}
