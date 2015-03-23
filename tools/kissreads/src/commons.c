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
 * commons.c
 *
 *  Created on: 17 sept. 2010
 *      Author: ppeterlo
 */

#include<commons.h>
#include <ctype.h>
#include<stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h> // Added by Pierre Peterlongo on 10/09/2012.
#include <stdint.h>
#include<math.h>


#include<inttypes.h> // DEBUG

//#define GET_ONLY_UPPER_CHARS // can be used for analysing outputs of kissnp where the extension is in upper case while the 2k+1 snp is in upper case. On wants only to analyse the 2k+1 snp

#define MAX_SIZE_LINE 1048576
//char line[MAX_SIZE_LINE];

extern kmer_type mask_code_seed;

static int
cmpstringp(const void *p1, const void *p2)
{
	/* Les arguments de cette fonction sont des "pointeurs de
              pointeurs sur des caractères", mais les arguments de
              strcmp(3) sont des "pointeurs sur des caractères", d’où
              le forçage de type et l’utilisation de l’astérisque */

	return strcmp(* (char * const *) p1, * (char * const *) p2);
}



char ** sort_strings (char ** strings, int number)
{
  qsort(strings,number, sizeof(char *), cmpstringp);
  return strings;
}


char * mystrdup (const char *s1)
{
    sum_memory_strdup+=strlen(s1);
    //printf("Strdup %i \n", strlen(s1));

    return strdup(s1);
}
void * mymalloc(const int size){
  sum_memory+=size;
 // printf("Allocate %d, sum=%lli  = %lli Mo\n", size, sum_memory, sum_memory /1024LL/1024LL);
  return malloc(size);
}

void * mycalloc(const int size, const int size_2){
  sum_memory+=size*size_2;
  //printf("Allocate %d, sum=%lli  = %lli Mo\n", size*size_2, sum_memory, sum_memory /1024LL/1024LL);
  return calloc(size,size_2);
}


void init_static_variables(const int k){
	int i;
	for (i=0;i<'T'+1;i++) comp[i]=i; // for other iupac alphabet letters
	comp['A']='T';
	comp['T']='A';
	comp['C']='G';
	comp['G']='C';
	comp['a']='t';
	comp['t']='a';
	comp['c']='g';
	comp['g']='c';

	nuc[0]='A';
	nuc[1]='C';
	nuc[2]='G';
	nuc[3]='T';

	anykmer = (char *) malloc(k+1); test_alloc(anykmer);
    
    
    mask_code_seed=1; // don't know why but 1<<(2*k)  does not work with k>32. This is why I made this stupid loop/
    int z;
    for (z=0;z<(2*k);z++){
        mask_code_seed = mask_code_seed<<1;
    }
    
    mask_code_seed = mask_code_seed-1;
    
}

int valid_character(const char c){
	if(c=='A' || c=='C' || c=='G' || c=='T') return 1;
	return 0;
}

void print_rev_comp(char s[], FILE* out){
	int i;
	const int len=strlen(s);
	for(i=len-1;i>-1;i--)fprintf(out, "%c",comp[(int)s[i]]);
}

void revcomp(char s[], int len)
{
	int i;
	char t;
	for (i=0;i<len/2;i++)
	{
		t=s[i];
		s[i] = comp [(int)(s[len-i-1])];
		s[len-i-1] = comp [(int)(t)];
	}
	if (len%2==1)
		s[len/2]=comp[(int)(s[len/2])];

}

void rev(char s[], int len)
{
	int i;
	char t;
	for (i=0;i<len/2;i++)
	{
		t=s[i];
		s[i] = s[len-i-1];
		s[len-i-1] = t;
	}
	if (len%2==1)
		s[len/2]=s[len/2];
}

char * to_upper (char  * word){
	int i=0;
	while(word[i]!='\0') {word[i]=toupper(word[i]); i++;}
	return word;
}

char * to_lower (char  * word){
	int i=0;
	while(word[i]!='\0') {word[i]=tolower(word[i]); i++;}
	return word;
}

//char line[1048576];

/**
 * removes eventual first '>' character
 * remove eventual last '\n' character
 * duplicate the raw_comment
 */
char * format_comment(char * raw_comment){
	if(raw_comment[0]=='>') raw_comment++;
	int l = strlen(raw_comment);
	if(raw_comment[l-1]=='\n') raw_comment[l-1]='\0';
	return strdup(raw_comment);
}

int get_next_sequence_and_comments_for_starters_fasta (char * sequence, char * comment, const char input_only_upper, char * line){
	char *rv;
	char *p;
	int nextchar=0;
	rv=gzgets(file,(char *)comment,MAX_SIZE_LINE);// read comment ('>read00xxxx...\n')
	

	if(rv == NULL) return 0;
	do{
	  rv=gzgets(file,(char *)sequence,MAX_SIZE_LINE); //
	}while(sequence[0]=='>');


	p = (char *)strchr((char*)sequence, '\n');
	if (p) *p = '\0';
	p = (char *)strchr((char*)sequence, '\r');
	if (p) *p = '\0';

	nextchar=gzgetc(file); // cheat, reads the next '>' character in order to induce EOF

	while (nextchar!='>' && !gzeof(file))
	{
		gzseek(file, -1, SEEK_CUR);
		rv=gzgets(file,(char *)line,MAX_SIZE_LINE);// read comment ('>read00xxxx...\n')
		rv = strchr(line, '\n'); // find the last \n char
		if(rv) *rv = '\0';       // change it into \0
		rv = strchr(line, '\r'); // find the last \r char
		if(rv) *rv = '\0';       // change it into \0

		strcat(sequence, line); // concat the restult in the sequence

		nextchar=gzgetc(file); // cheat, reads the next '>' character in order to induce EOF
	}
	gzseek(file, -1, SEEK_CUR); // Go back to previous read character
    
//#ifndef GET_ONLY_UPPER_CHARS
    if(!input_only_upper)
        to_upper(sequence);
//#endif
	//printf("return sequence %s\n",sequence);
//	free(line);
	return strlen(sequence); // readlen
}


int get_next_sequence_and_comments_for_starters_fastq (char * sequence, char * comment, const char input_only_upper, char * line){

	char *rv, *qv;
	char *p;
	char quality[96000];
	//does not work if the sequence is written on several lines
//	char * line = malloc(sizeof(char)*1048576);

	rv=gzgets(file,(char *)comment,MAX_SIZE_LINE);// read comment1 ('@read00xxxx...\n')

	if(rv == NULL) return 0;
	do{
	  rv=gzgets(file, (char *)sequence,MAX_SIZE_LINE); //
	}while(sequence[0]=='@');

	qv=gzgets(file, (char *)line,MAX_SIZE_LINE);// read comment2 ('+read00xxxx...\n')
	if(qv == NULL) return 0;
	qv=gzgets(file, (char *)quality,MAX_SIZE_LINE); //

	p = (char *)strchr((char*)sequence, '\n');
	if (p) *p = '\0';
	p = (char *)strchr((char*)sequence, '\r');
	if (p) *p = '\0';

	//nextchar=gzgetc(file); // cheat, reads the next '>' character in order to induce EOF
//#ifndef GET_ONLY_UPPER_CHARS
    if(!input_only_upper)
        to_upper(sequence);
//#endif
	//printf("return sequence %s\n",sequence);
//	free(line);
	return strlen(sequence); // readlen
}

int get_next_sequence_and_comments_for_starters (char * sequence, char * comment, const char input_only_upper, char * line){
  char nextchar=gzgetc(file);
  gzseek(file, -1, SEEK_CUR); // Go back to previous read character
  if(nextchar=='@') return  get_next_sequence_and_comments_for_starters_fastq(sequence, comment, input_only_upper,line);
  if(nextchar=='>') return  get_next_sequence_and_comments_for_starters_fasta(sequence, comment, input_only_upper,line);
  gzgets(file, sequence,MAX_SIZE_LINE);
  fprintf(stderr,"could not determine if the file is fasta or fastq in line %s, exit\n", sequence);
  exit(1);
}

int get_next_sequence_and_comments_for_fastq (gzFile file, char * sequence, char * comment, char * quality, char * line){
	char *rv, *qv;
	char *p;
	//does not work if the sequence is written on several lines
//char * line = malloc(sizeof(char)*1048576);

	rv=gzgets(file,(char *)comment,MAX_SIZE_LINE);// read comment1 ('@read00xxxx...\n')

	if(rv == NULL) return 0;
	do{
	  rv=gzgets(file, (char *)sequence,MAX_SIZE_LINE); //
	}while(sequence[0]=='@');

	qv=gzgets(file, (char *)line,MAX_SIZE_LINE);// read comment2 ('+read00xxxx...\n')
	if(qv == NULL) return 0;
	qv=gzgets(file, (char *)quality,MAX_SIZE_LINE); //

	p = (char *)strchr((char*)sequence, '\n');
	if (p) *p = '\0';
	p = (char *)strchr((char*)sequence, '\r');
	if (p) *p = '\0';

	//	nextchar=gzgetc(file); // cheat, reads the next '>' character in order to induce EOF

	to_upper(sequence);
	//printf("return sequence %s\n",sequence);
//free(line);
	return strlen(sequence); // readlen
}

int get_next_fasta_sequence (gzFile file, char * sequence , char * line){
	char *rv;
	char *p;
	int nextchar=0;
//char * line = malloc(sizeof(char)*1048576);
	
 rv=gzgets(file, (char *)line,MAX_SIZE_LINE);// read comment ('>read00xxxx...\n')
	if(rv == NULL) return -1;
	do{
	  rv=gzgets(file,(char *)sequence,MAX_SIZE_LINE); //
	}while(sequence[0]=='>');

	
	p = (char *)strchr((char*)sequence, '\n');
	if (p) *p = '\0';
	p = (char *)strchr((char*)sequence, '\r');
	if (p) *p = '\0';

	nextchar=gzgetc(file); // cheat, reads the next '>' character in order to induce EOF

	while (nextchar!='>' && !gzeof(file))
	{
	  
	  gzseek(file, -1, SEEK_CUR);
	  rv=gzgets(file, (char *)line,MAX_SIZE_LINE);// read a line
	  rv = strchr(line, '\n'); // find the last \n char
	  if(rv) *rv = '\0';       // change it into \0
	  rv = strchr(line, '\r'); // find the last \r char
	  if(rv) *rv = '\0';       // change it into \0
	  
	  strcat(sequence, line); // concat the result in the sequence
	  
	  nextchar=gzgetc(file); // cheat, reads the next '>' character in order to induce EOF
	}
    // gzseek(file, -1, SEEK_CUR); // Go back to previous read character
	to_upper(sequence);
	return strlen(sequence); // readlen
}
int get_next_sequence_for_fastq (gzFile file, char * sequence, char * quality, char * line){
	char *rv, *qv;
	char *p;
	//does not work if the sequence is written on several lines
//char * line = malloc(sizeof(char)*MAX_SIZE_LINE);

	rv=gzgets(file,(char *)line,MAX_SIZE_LINE);// read comment ('@read00xxxx...\n')
	if(rv == NULL) return -1;
	do{
		rv=gzgets(file,(char *)sequence,MAX_SIZE_LINE); //
	}while(sequence[0]=='@');

	qv=gzgets(file,(char *)line,MAX_SIZE_LINE);// read comment ('+read00xxxx...\n')
	if(qv == NULL) return -1;
	qv=gzgets(file,(char *)quality,MAX_SIZE_LINE); //

	p = (char *)strchr((char*)sequence, '\n');
	if (p) *p = '\0';
	p = (char *)strchr((char*)sequence, '\r');
	if (p) *p = '\0';

	//	nextchar=gzgetc(file); // cheat, reads the next '@' character in order to induce EOF

	to_upper(sequence);
	//printf("return sequence %s\n",sequence);
//free(line);
	return strlen(sequence); // readlen
}

int number_of_sequences_in_fasta_file(gzFile file, char * line){
	int sequences_number=0;
	int previous_is_a_comment=0;
	gzrewind(file);
//	char line[1048576];
//char * line = malloc(sizeof(char)*1048576);
	do{
		if(gzgets(file,line,MAX_SIZE_LINE) == NULL) break;
		if(line[0]=='>'){
			if(!previous_is_a_comment){
				previous_is_a_comment=1;
				sequences_number++;
			}
		}
		else previous_is_a_comment=0;
	}
	while(1);
	gzrewind(file);
//free(line);
	return sequences_number;
}

int number_of_sequences_in_fastq_file(gzFile file, char * line){
	int sequences_number=0;
	int previous_is_a_comment=0;
	int line_number=0;
	gzrewind(file);
//	char line[1048576];
	do{
		if(gzgets(file,line,MAX_SIZE_LINE) == NULL) break;
		line_number++;
		if(line[0]=='@' && line_number%4 == 1){
			//this does not work if the sequence is given on several lines, is this possible for fastq???
			if(!previous_is_a_comment){
				previous_is_a_comment=1;
				sequences_number++;
			}
		}
		else previous_is_a_comment=0;
	}
	while(1);
	gzrewind(file);
	return sequences_number;
}

int number_of_sequences_in_file(gzFile file, char * line){

	gzrewind(file);
	if(gzgets(file,line,MAX_SIZE_LINE) == NULL) return 0;
	gzrewind(file);
	
	if(line[0]=='@') return number_of_sequences_in_fastq_file(file,line);
	if(line[0]=='>') return number_of_sequences_in_fasta_file(file,line);
	
	fprintf(stderr,"Error, read file starts with line %s, witch is not fastq nor fasta format\n", line);
	return 0;
}


// binary code of any character coded on 2 bits.
// Among them: N or n or G or G=11 --- A or a = 00 --- C or c = 01 --- T or t = 10
inline int NT2int(const char nt)
{
    return (nt>>1)&3;
    
}



// update a code of a seed with a new character O(1)
kmer_type  updateCodeSeed(const char *seq, kmer_type *x) // update of a seed (shift and adding a new character)
{
//    printf("seed2 %" PRIu64 " mask %" PRIu64 "\n", *x, mask_code_seed);
    *x = (*x)*4 + NT2int(seq[size_seeds-1]); // add the code of the new nucleotid
    
    *x = *x & mask_code_seed; // remove the leftmost couple of bits
//    printf("seed3 %" PRIu64 " mask %" PRIu64 "\n", *x, mask_code_seed);
    return *x;
}


// transform a character seed into a code seed O(size seed)
kmer_type codeSeed(const char *seq) // initialisation of a seed
{
    int i;
    kmer_type x=0;
    for (i=0; i<size_seeds; ++i)
    {
        x = x<<2;
        x+=NT2int(seq[i]);
        
    }

    return x;
}




