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
 * commons.h
 *
 *  Created on: 17 sept. 2010
 *      Author: ppeterlo
 */

#ifndef COMMONS_H_
#define COMMONS_H_


//#define GET_ONLY_UPPER_CHARS // can be used for analysing outputs of kissnp where the extension is in upper case while the 2k+1 snp is in upper case. On wants only to analyse the 2k+1 snp
#include<stdio.h>
#include <stdint.h>
#include<zlib.h>
//int artificial_overlap;
char comp ['t'+1];
char nuc [4];
char standard_fasta;
char verbose;
char silent;
char whole_coverage;
char quality;
int kmer_size;
int size_seeds;

#ifdef INPUT_FROM_KISSPLICE
int min_overlap;
int countingOption;
#endif
int average_size_reads;
int nb_event_sets;
int number_of_read_sets;
int size_before_reads_starting; // see fragment_info.h
char * anykmer;
char **  sort_strings (char ** strings, int number);
void * mymalloc(const int size);
void * mycalloc(const int size, const int size_2);
char * mystrdup (const char *s1);
void print_rev_comp(char s[], gzFile out);
void revcomp(char s[], int len);
void rev(char s[], int len);
void init_static_variables(const int k);
char * to_upper (char  * word);
char * to_lower (char  * word);
char * format_comment(char * raw_comment);
int get_next_sequence_and_comments (gzFile file, char * sequence, char * comment, char * line);
int get_next_sequence_and_comments_for_starters (char * sequence, char * comment, const char input_only_upper, char * line);
int get_next_fasta_sequence (gzFile file, char * read, char * line);
int get_next_sequence_and_comments_for_fastq (gzFile file, char * sequence, char * comment, char * quality, char * line);
int get_next_sequence_for_fastq (gzFile file, char * read, char * quality, char * line);
int number_of_sequences_in_file(gzFile file, char * line);
gzFile file;

uint64_t sum_memory;
uint64_t sum_memory_strdup;

uint64_t  mask_nbseed ;
uint64_t  mask_offset_seed; 
unsigned int nbits_nbseeds;

// macro to test if a variable is null (i.e., a malloc failed)
#define test_alloc( variable) {	if(variable == NULL){		fprintf(stderr,"cannot allocate memory for variable %s, exit\n",#variable);		exit(1);	}}

typedef  uint64_t kmer_type;
kmer_type  codeSeed(const char *seq);

//saturated increment of unsigned char
#define Sinc8(a)  ((a == 0xFF ) ? 0xFF : a++)
#define Sinc24(a)  ((a == 0xFFFFFF ) ? 0xFFFFFF : a++)

#endif /* COMMONS_H_ */
