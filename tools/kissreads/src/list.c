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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "list.h"
#include<commons.h>

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
		/*
		 * (rayan) FIXME: what is this line doing here? i'm disabling it for now
		 */
		//if(aux->val != NULL) free(aux->val);
		free(aux);
	}
}

/**
 * free a list and its generic content: free all cells, and  free their val using free().
 */
void list_of_generic_free(const void * v_list)
{
	list * p_list = (list *) v_list;
	cell *aux;
	if(!p_list) return;
	cell * l = p_list->first;

	while (l != NULL)
	{
		if(l->val != NULL) free(l->val);
		l->val=NULL;
		aux = l;
		l = l->prox;
		free(aux);
	}
	free(p_list);
}

void list_del(list *l, cell *c)
{
	cell *cur = l->first;
	if (cur == c)
	{
		cell *next=cur->prox;
		free(c);
		l->first=next;
		return;
	}
	while (cur != NULL)
	{
		if (cur->prox == c)
		{
			cur->prox=c->prox;
			free(c);
			return;
		}
		if (cur->prox == NULL)
			return;
		cur=cur->prox;
	}
}

/*
 * for debugging purposes..
 */
void list_print(list *l)
{
	cell *c=l->first;
	while (c != NULL)
	{
		//printf("cell %d\n",((couple*)c->val)->b);
		//printf("cell %X\n",(c->val));
		//		printf("cell %X\n",(c));
		c=c->prox;
	}
}
