/*****************************************************************************
 *   discoSnp++: discovering polymorphism from raw unassembled NGS reads
 *   A tool from the GATB (Genome Assembly Tool Box)
 *   Copyright (C) 2020  INRIA
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




void init_static_variables(const int k){
    
}

void GlobalValues::revcomp(char s[])
{
	int i;
    const int len = strlen(s);
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

void GlobalValues::rev(char s[])
{
	int i;
    const int len = strlen(s);
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



/** Transform a nucleotide in ASCII form into an integer form as:
 *     - A=0
 *     - C=1
 *     - T=2
 *     - G=3
 * \param[in] nt : the nucleotide in ASCII
 * \return the translated nucleotide */
static int NT2int(char nt)  {  return (nt>>1)&3;  }

// update a code of a seed with a new character O(1)
kmer_type  GlobalValues::updateCodeSeed(const char *seq, kmer_type *x) // update of a seed (shift and adding a new character)
{
    *x = (*x)*4 + NT2int(seq[size_seeds-1]); // add the code of the new nucleotid
    
    *x = *x & mask_code_seed; // remove the leftmost couple of bits
    return *x;
}


// transform a character seed into a code seed O(size seed)
kmer_type GlobalValues::codeSeed(const char *seq) // initialisation of a seed
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




