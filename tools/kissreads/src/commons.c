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

//#define GET_ONLY_UPPER_CHARS // can be used for analysing outputs of kissnp where the extension is in upper case while the 2k+1 snp is in upper case. On wants only to analyse the 2k+1 snp

#define MAX_SIZE_LINE 1048576
//char line[MAX_SIZE_LINE];

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
}

void print_rev_comp(char s[], gzFile out){
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
	  rv=gzgets(file, (char *)line,MAX_SIZE_LINE);// read comment ('>read00xxxx...\n')
	  rv = strchr(line, '\n'); // find the last \n char
	  if(rv) *rv = '\0';       // change it into \0
	  rv = strchr(line, '\r'); // find the last \r char
	  if(rv) *rv = '\0';       // change it into \0
	  
	  strcat(sequence, line); // concat the restult in the sequence
	  
	  nextchar=gzgetc(file); // cheat, reads the next '>' character in order to induce EOF
	}
	//	gzseek(file, -1, SEEK_CUR); // Go back to previous read character
	to_upper(sequence);
//free(line);
	return strlen(sequence); // readlen
}
int get_next_sequence_for_fastq (gzFile file, char * sequence, char * quality, char * line){
	char *rv, *qv;
	char *p;
	//does not work if the sequence is written on several lines
//char * line = malloc(sizeof(char)*MAX_SIZE_LINE);

	rv=gzgets(file,(char *)line,MAX_SIZE_LINE);// read comment ('@read00xxxx...\n')
	if(rv == NULL) return 0;
	do{
		rv=gzgets(file,(char *)sequence,MAX_SIZE_LINE); //
	}while(sequence[0]=='@');

	qv=gzgets(file,(char *)line,MAX_SIZE_LINE);// read comment ('+read00xxxx...\n')
	if(qv == NULL) return 0;
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
//	char line[1048576];
	if(gzgets(file,line,MAX_SIZE_LINE) == NULL) return 0;
	gzrewind(file);
	
	if(line[0]=='@') return number_of_sequences_in_fastq_file(file,line);
	if(line[0]=='>') return number_of_sequences_in_fasta_file(file,line);
	
	fprintf(stderr,"Error, read file starts with line %s, witch is not fastq nor fasta format\n", line);
	return 0;
}

int NT2int(char nt)
{
    int i;
    i = nt;
    i = (i>>1)&3; // that's quite clever, guillaume.
    return i;
}


kmer_type  codeSeed(const char *seq)
{
    int i;
    kmer_type x;
    
    x=0;
    for (i=0; i<size_seeds; ++i)
    {
        x = x*4 + NT2int(seq[i]);
    }
    return x;
}




