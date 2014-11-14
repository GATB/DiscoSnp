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
 * list.c
 *
 *  Created on: 16 sept. 2010
 *      Author: ppeterlo
 */
#include "list.h"

int numberInList (list *l)
{
	return l->size;
}

list* list_create(void){
	list * l = (list *)malloc(sizeof(list)); test_alloc(l);
	l->size=0;
	l->first=NULL;
	return l;
}

void list_add(list *l, void * val)
{
	cell * new = (cell *)malloc(sizeof(cell)); test_alloc(new);
	new->val = val;

	l->size++;
	new->prox = l->first;
	l->first=new;
}

/**
 * free a list and its generic content: free all cells, and specifically free their val using a function given as argument.
 */
void list_free(const void *v_list, void (*specific_free)(void *))
{
	list * p_list = (list *) v_list;
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
	list * p_list = (list *) v_list;
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
	free((list *) v_list);
}

/*
 * for debugging purposes..
 */
void list_print(list *l)
{
	cell *c=l->first;
	while (c != NULL)
	{
		c=c->prox;
	}
}

char list_of_int_contains(list *l, int val){
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
	cellint * new = (cellint *)malloc(sizeof(cellint)); test_alloc(new);
	new->val = val;
    
	l->size++;
	new->prox = l->first;
	l->first=new;
}



char listint_contains(listint *l, int val){
    cellint *c=l->first;
    while (c != NULL)
    {
        if ((c->val)==val) return 1;
        c=c->prox;
    }
    return 0;
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

