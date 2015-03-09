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
 * commons.h
 *
 *  Created on: 17 sept. 2010
 *      Author: ppeterlo
 */

#ifndef COMMONS_H_
#define COMMONS_H_


//#define GET_ONLY_UPPER_CHARS // can be used for analysing outputs of kissnp where the extension is in lower case while the 2k+1 snp is in upper case. On wants only to analyse the 2k+1 snp
#include<stdio.h>
#include <stdint.h>
#include<zlib.h>
typedef  uint64_t kmer_type;
//int artificial_overlap;
char comp ['t'+1];
char nuc [4];
char standard_fasta;
char silent;
char quality;
char only_print;
int kmer_size;
int size_seeds;
int minimal_read_overlap;
kmer_type mask_code_seed;

int valid_character(const char c);
int average_size_reads;
int nb_event_sets;
int number_of_read_sets;
int size_before_reads_starting; // see fragment_info.h
char * anykmer;
char **  sort_strings (char ** strings, int number);
void * mymalloc(const int size);
void * mycalloc(const int size, const int size_2);
char * mystrdup (const char *s1);
void print_rev_comp(char s[], FILE* out);
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

kmer_type  codeSeed(const char *seq);
kmer_type  updateCodeSeed(const char *seq, kmer_type *x);

//saturated increment of unsigned char
#define Sinc8(a)  ((a == 0xFF ) ? 0xFF : a++)
#define Sinc24(a)  ((a == 0xFFFFFF ) ? 0xFFFFFF : a++)

#endif /* COMMONS_H_ */
