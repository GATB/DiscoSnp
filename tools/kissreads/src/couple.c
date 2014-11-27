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
