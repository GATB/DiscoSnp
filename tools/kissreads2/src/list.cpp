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
 * list.c
 *
 *  Created on: 16 sept. 2010
 *      Author: ppeterlo
 */
#include "list.h"

int numberInList (mylist *l)
{
	return l->size;
}

mylist* list_create(void){
	mylist * l = (mylist *)malloc(sizeof(mylist)); test_alloc(l);
	l->size=0;
	l->first=NULL;
	return l;
}

void list_add(mylist *l, void * val)
{
	cell * newcell = (cell *)malloc(sizeof(cell)); test_alloc(newcell);
	newcell->val = val;

	l->size++;
	newcell->prox = l->first;
	l->first=newcell;
}

/**
 * free a list and its generic content: free all cells, and specifically free their val using a function given as argument.
 */
void list_free(const void *v_list, void (*specific_free)(void *))
{
	mylist * p_list = (mylist *) v_list;
	cell *aux;
	cell *l = p_list->first;

	while (l != NULL)
	{
		aux = l;
		l = l->prox;
		specific_free(aux->val);
		free(aux);
	}
    free(p_list);
}

/**
 * empty a list and its generic content: free all cells, and  free their val using free().
 */
void list_of_generic_empty(const void * v_list)
{
	mylist * p_list = (mylist *) v_list;
	cell *aux;
	if(!p_list) return;
	cell * l = p_list->first;
    
	while (l != NULL)
	{
//        printf("freeing value %d\n", *(int*)l->val); //DEB
		if(l->val != NULL) free(l->val);
		aux = l;
		l = l->prox;
		free(aux);
	}
    p_list->size=0;
	p_list->first=NULL;
}

/**
 * free a list and its generic content: free all cells, and  free their val using free().
 */
void list_of_generic_free(const void * v_list)
{
	list_of_generic_empty(v_list);
	free((mylist *) v_list);
}

/*
 * for debugging purposes..
 */
void list_print(mylist *l)
{
	cell *c=l->first;
	while (c != NULL)
	{
		c=c->prox;
	}
}

char list_of_int_contains(mylist *l, int val){
    cell *c=l->first;
	while (c != NULL)
	{
		if ((*(int*)c->val)==val) return 1;
        c=c->prox;
	}
    return 0;
}



/////////////////////////////////////////////////////////////////////////////////
///////// For speed purpose, here is a list dedicated to int only ///////////////
/////////////////////////////////////////////////////////////////////////////////

listint *listint_create(void){
	listint * l = (listint *)malloc(sizeof(listint)); test_alloc(l);
	l->size=0;
	l->first=NULL;
	return l;
}

void listint_add(listint *l, int val)
{
	cellint * newcell = (cellint *)malloc(sizeof(cellint)); test_alloc(newcell);
	newcell->val = val;
    
	l->size++;
	newcell->prox = l->first;
	l->first=newcell;
}



bool listint_contains(listint *l, int val){
    
    
    
    cellint *c=l->first;
    while (c != NULL)
    {
        if ((c->val)==val) return true;
        c=c->prox;
    }
      return false;
}


void listint_empty(const void *v_list){
    
	listint * p_list = (listint *) v_list;
	cellint *aux;
	if(!p_list) return;
	cellint * l = p_list->first;
    
	while (l != NULL)
	{
		aux = l;
		l = l->prox;
		free(aux);
	}
    p_list->size=0;
	p_list->first=NULL;
}
void listint_free(const void * v_list)
{
	listint_empty(v_list);
	free((listint *) v_list);
}



int numberInListint (listint *l){
	return l->size;
}

