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
 * tree.c
 *
 *  Created on: 2 déc. 2010
 *      Author: ppeterlo
 */


#include <tree.h>
#include <commons.h>
#include <stdlib.h>
#include <string.h>
#include <hash.h>
#include <fragment_info.h>
#include <fragment_index.h>
#include<limits.h>
#include<assert.h>
#include <string.h>


int global_max_deep;

#define INIT_TREE_VALUES(fragment)	\
tree t = (tree) malloc (sizeof(node)); test_alloc(t); \
t->to_be_indexed=1; \
t->children=NULL; \
t->Mapped=NULL; \
t->fragment=fragment; \
t->out_degree=0; \
t->in_degree=0; \
t->simplified=0; \
t->merged=0; \
t->pruned=0; \
t->kmers_indexation_done=0; \
t->global_coverage=-1; \
t->had_been_extended=0;


tree tree_create_tree (char * fragment){
	INIT_TREE_VALUES(fragment)
	t->fragments_already_used=hash_create(1); test_alloc(t->fragments_already_used);
	return t;
}
tree tree_create_child (char * fragment, hash_t reads_already_used){
	INIT_TREE_VALUES(fragment)
	t->fragments_already_used=reads_already_used;
	return t;
}

void tree_add_child (tree t, char * child_fragment){
	p_node child_node = tree_create_child(child_fragment, t->fragments_already_used);
	child_node->in_degree++;
//	if(t->children==NULL){
//		t->children = (p_node *)malloc(sizeof(p_node)); test_alloc(t->children);
//		t->children[0]=child_node;
//		t->out_degree=1;
//		return;
//	}
//	else{
		t->children = (p_node *)realloc(t->children, sizeof(p_node)*(t->out_degree+1)); test_alloc(t->children);
		t->children[t->out_degree++]=child_node;
//	}
}

void tree_add_child_node (tree t, p_node child_node){
	child_node->in_degree++;
//	if(t->children==NULL){
//		t->children = (p_node *)malloc(sizeof(p_node)); test_alloc(t->children);
//		t->children[0]=child_node;
//		t->out_degree=1;
//		return;
//	}
//	else{
		t->children = (p_node *)realloc(t->children, sizeof(p_node)*(t->out_degree+1)); test_alloc(t->children);
		t->children[t->out_degree++]=child_node;
//	}
}


void tree_free_tree (tree t){
	if (t==NULL) return;
	int i;
	for(i=0;i<t->out_degree;i++) tree_free_tree(t->children[i]);
	free(t->fragment);
	free(t->children);
	free(t->Mapped);
	free(t);
}

//
//int main_test_tree(int argc, char **argv) {
//
//	tree t=tree_create_tree(strdup(""));
//	tree_add_child(t, strdup("premier fils"));
//	tree_add_child(t, strdup("second fils"));
//	tree_add_child(t->children[0], strdup("premier fils du premier fils"));
//	tree_add_child(t->children[0], strdup("second fils du premier fils"));
//	tree_print_tree(t);
//
//	printf("recursive all paths: \n");
//	tree_print_all_paths(t, "prefix_", 0, 0, stdout, 0, 12);
//
//	tree_free_tree(t);
//	return 1;
//}


