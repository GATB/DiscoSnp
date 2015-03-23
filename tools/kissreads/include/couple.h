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
