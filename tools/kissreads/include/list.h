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
 * list.h
 *
 *  Created on: 16 sept. 2010
 *      Author: ppeterlo
 */

#ifndef LIST_H_
#define LIST_H_
#include<stdio.h>
#include<stdlib.h>
#include<couple.h>
#include<commons.h>
#include <string.h>


typedef struct CELL {
   void * val;
   struct CELL *prox;
} cell, * p_cell;

typedef struct {
	int size;
	p_cell first;
} list;

list *list_create(void);
void list_add(list *l, void* val);



void list_free(const void *list, void (*specific_free)(void *));
void list_of_generic_free(const void *v_list);
void list_of_generic_empty(const void * v_list);


int numberInList (list *l);

/////////////////////////////////////////////////////////////////////////////////
///////// For speed purpose, here is a list dedicated to int only ///////////////
/////////////////////////////////////////////////////////////////////////////////

typedef struct CELLINT {
    int val;
    struct CELLINT *prox;
} cellint, * p_cellint;

typedef struct {
	int size;
	p_cellint first;
} listint;

listint *listint_create(void);
void listint_add(listint *l, int val);


char listint_contains(listint *l, int val);

void listint_free(const void *v_list);
void listint_empty(const void * v_list);


int numberInListint (listint *l);


#endif /* LIST_H_ */
