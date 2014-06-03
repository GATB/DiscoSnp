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
 * tree.h
 *
 *  Created on: 2 déc. 2010
 *      Author: ppeterlo
 */
#include<stdio.h>
#include <fragment_info.h>

#ifndef TREE_H_
#define TREE_H_

typedef struct _node {
	char * fragment; // fragment of text extending something
	int to_be_indexed; // used to know if a node is reused or if it's just created. Alsoo used as id for nodes during outputs
	int out_degree;
	int in_degree;
	char merged;     // boolean to know if the node was already merged
	char simplified; // boolean to know if the node was already simplified.
	char pruned; // boolean to know if the fragment of the node was already pruned
	char had_been_extended; // boolean to know if a node had been extended (with or without success).
	char kmers_indexation_done; // boolean to know if we already indexed the kmer present in the fragment for computing coverage
	int global_coverage; // support of each kmer present in the fragment divided by |fragment|-k+1. Initially fixed to -1 to know it has not been computed.
	int * locus_coverage; // gives the coverage, position by position. At each position, store the kmer that started at this position. Note that position < k exist as the coverage is computed while k-1 overlap existed with previous node.
	struct _node ** children; // array of children
	hash_t fragments_already_used; // reads already used for this inital consensus fragment
	hash_t Mapped;             // couples read -> position used during extension algorithms: reads mapping on the extending fragment
} *p_node, node, *tree;



tree tree_create_tree (char * fragment);
//tree tree_create_tree (char * fragment);
tree tree_create_child (char * fragment, hash_t reads_already_used);
void tree_add_child (tree t, char * child_fragment);
void tree_add_child_node (tree t, p_node child_node);
#endif /* TREE_H_ */
