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
 * commons.h
 *
 *  Created on: 17 sept. 2010
 *      Author: ppeterlo
 */

#ifndef COMMONS_H_
#define COMMONS_H_


#define NBITS_OFFSET_SEED 40
// macro to test if a variable is null (i.e., a malloc failed)
#define test_alloc( variable) {	if(variable == NULL){		fprintf(stderr,"cannot allocate memory for variable %s, exit\n",#variable);		exit(1);	}}
#define Sinc8(a)  ((a == 0xFF ) ? 0xFF : a++)
#define Sinc24(a)  ((a == 0xFFFFFF ) ? 0xFFFFFF : a++)


//#define GET_ONLY_UPPER_CHARS // can be used for analysing outputs of kissnp where the extension is in lower case while the 2k+1 snp is in upper case. On wants only to analyse the 2k+1 snp
#include<stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <gatb/gatb_core.hpp>


#define MAX_SIZE_LINE 16777215
// KMERS ARE LIMITED TO 32 Nucleotides
typedef  uint64_t kmer_type;

class GlobalValues {
public:
    
    unsigned char *comp; //['t'+1];
    unsigned char *nuc; //[4];
    // OPTIONS
    bool standard_fasta;
    // BINARY ALPHABET
    int size_seeds;
    int index_stride;
    int minimal_read_overlap;
    bool compute_genotypes;
    kmer_type mask_code_seed;
    int number_of_read_sets;
    int subst_allowed;
    vector<int> min_coverage; //minimal coverage per read set
    
    uint64_t  mask_nbseed;
    uint64_t  mask_offset_seed;
    unsigned int nbits_nbseeds;//saturated increment of unsigned char
    bool radseq_option;
    bool phasing;
    
    GlobalValues(){
        
        
        comp=(unsigned char *)malloc(sizeof(unsigned char)*('t'+1)); //test_alloc(comp)
        
        nuc=(unsigned char *)malloc(sizeof(unsigned char)*4); //test_alloc(nuc);
        
        
        for (int i=0;i<'T'+1;i++) comp[i]=i; // for other iupac alphabet letters
        comp[(unsigned char)'A']='T';
        comp[(unsigned char)'T']='A';
        comp[(unsigned char)'C']='G';
        comp[(unsigned char)'G']='C';
        comp[(unsigned char)'a']='t';
        comp[(unsigned char)'t']='a';
        comp[(unsigned char)'c']='g';
        comp[(unsigned char)'g']='c';
        
        nuc[0]='A';
        nuc[1]='C';
        nuc[2]='G';
        nuc[3]='T';
        
        mask_offset_seed = (1ULL << (NBITS_OFFSET_SEED)) -1 ;
        nbits_nbseeds = 8*sizeof(uint64_t)- NBITS_OFFSET_SEED ;
        mask_nbseed  = ( 1ULL <<  (uint64_t) nbits_nbseeds  ) -1 ;
        
        
    };
    
    void set_mask_code_seed(){
        mask_code_seed=1; // don't know why but 1<<(2*k)  does not work with k>32. This is why I made this stupid loop/
        for (int z=0;z<(2*size_seeds);z++){
            mask_code_seed = mask_code_seed<<1;
        }
        mask_code_seed = mask_code_seed-1;
    }
    
    void revcomp(char s[]);
    void rev(char s[]);
    int valid_character(const char c);
    kmer_type  codeSeed(const char *seq);
    kmer_type  updateCodeSeed(const char *seq, kmer_type *x);
    char * format_comment(const char * raw_comment);
} ;





static int NT2int(char nt);



#endif /* COMMONS_H_ */
