#ifdef ED

/**
 * Copyright INRIA and ENS, contributors Peterlongo and Chikhi
 * pierre.peterlongo@inria.fr
 * rayan.chikhi@irisa.fr
 *
 * This software is a computer program whose purpose is to detect the
 * presence of a starter in a set of NGS reads, and to provide the neighbor
 * in case of success..
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

#include "commons.h"
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern size_t	strlen (const char*);
char comp['t'+1];
char nuc [4];
char keyi2nucleotide_array ['T'+1];
int size_seeds=25; // USED ONLY FOR MAPPING, might be different from the kmer_size used in minia. Size_seeds is no more than 32.
char * prefix_trashable=(char *)"trashmeplease";

unsigned int nbits_nbseeds;
uint64_t  mask_nbseed ;
uint64_t  mask_offset_seed; 


/**
 * redefined this fonction for mac users, absent from string.h
 */
char * strndup (const char *s, size_t n)
{
  char *result;
  size_t len = strlen (s);

  if (n < len)
    len = n;

  result = (char *) malloc (len + 1);
  if (!result)
    return 0;

  result[len] = '\0';
  return (char *) memcpy (result, s, len);
}


static int cmpstringp(const void *p1, const void *p2)
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


void init_static_variables(){
  
	int i;
	for (i=0;i<'T'+1;i++) {
		comp[i]=i; // for other iupac alphabet letters
	}
	printf("Init static variables\n");
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

	keyi2nucleotide_array['A']=0;
        keyi2nucleotide_array['C']=1;
        keyi2nucleotide_array['G']=2;
        keyi2nucleotide_array['T']=3;

}

void print_rev_comp(char s[], FILE * out){
	int i;
	const int len=strlen(s);
	for(i=len-1;i>-1;i--)fprintf(out, "%c",comp[(int)s[i]]);
}

void revcomp(char *s, int len)
{
	int i;
	char t;
	for (i=0;i<len/2;i++)
	{

	  t=s[i];
	  s[i] = comp[(int)(s[len-i-1])];
	  s[len-i-1] = comp[(int)(t)];
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
}

char * to_upper (char  * word){
	int i=0;
//	printf("before upper=%s %c\n", word, word[0]);
	while(word[i]!='\0')
		{
//			printf("%d %c %d", i, word[i], word[i]);
			word[i]=toupper(word[i]);
//			printf(" -> %c |", word[i]);
			i++;
		}
//	printf("after  upper=%s\n", word);
	return word;
}

char * to_lower (char  * word){ 
	int i=0;
	while(word[i]!='\0') {word[i]=tolower(word[i]); i++;}
	return word;
}

char line[1048576];

// TODO: REMOVE THIS AND CHECK ALL POTENTIAL PROJECTS USING THIS FUNCTION
int last_index_of(char * string, char character){
	int i, pos=-1;
	for(i=0;i<strlen(string);i++) if(string[i]==character) pos=i;
	return pos;
}

uint64_t codeSeed_uint64_t(const char *seq)
{
    int i;
    uint64_t x;
    
    x=0;
    for (i=0; i<size_seeds; ++i)
    {
        x = x*4 + NT2int(seq[i]);
    }
    return x;
}

int valid_character(const char c){
	if(c=='A' || c=='C' || c=='G' || c=='T') return 1;
	return 0;
}


#endif
