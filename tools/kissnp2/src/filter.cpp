//Copyright inria / irisa (2013)
//
//
//raluca.uricaru@gmail.com
//pierre.peterlongo@inria.fr
//
//This software is a computer program whose purpose is to call SNPs from NGS reads.
//
//This software is governed by the CeCILL license under French law and
//abiding by the rules of distribution of free software.  You can  use,
//modify and/ or redistribute the software under the terms of the CeCILL
//license as circulated by CEA, CNRS and INRIA at the following URL
//"http://www.cecill.info".
//
//As a counterpart to the access to the source code and  rights to copy,
//modify and redistribute granted by the license, users are provided only
//with a limited warranty  and the software's author,  the holder of the
//economic rights,  and the successive licensors  have only  limited
//liability.
//
//In this respect, the user's attention is drawn to the risks associated
//with loading,  using,  modifying and/or developing or reproducing the
//software by the user in light of its specific status of free software,
//that may mean  that it is complicated to manipulate,  and  that  also
//therefore means  that it is reserved for developers  and  experienced
//professionals having in-depth computer knowledge. Users are therefore
//encouraged to load and test the software's suitability as regards their
//requirements in conditions enabling the security of their systems and/or
//data to be ensured and,  more generally, to use and operate it in the
//same conditions as regards security.
//
//The fact that you are presently reading this means that you have had
//knowledge of the CeCILL license and that you accept its terms.


#include "filter.h"

#include "../minia/Kmer.h"
#include "Kmer_for_kissnp2.h"


/**
 * Methode permettant d'appliquer le filtre Low Complexity
 * Les zones ininteressantes des sequences sont mises en lettres minuscules
 * \param la sequence
 * \param lenseq la longueur de la sequence
 * \return 1 si l'operation s'est deroulee correctement
 */
int filterLowComplexity (char* data, int lenseq, int threshold)
{
	int DUSTSCORE[64];
	int i, j, k, s, m;

	for (i=0; i<64; i++) DUSTSCORE[i]=0;

	for (j=2; j<lenseq; ++j)
	{
		k = data[j-2]*16 + data[j-1]*4 + data[j];
		++DUSTSCORE[k];
	}

	s=0;
	for (i=0; i<64; ++i)
	{
		m = DUSTSCORE[i];
		s = s + m*(m-1)/2;
	}

    /*if (s>=threshold) return s;
    else return -1;*/
    return s;
}

/**
 * Methode permettant d'appliquer le filtre Low Complexity sur les deux chemins d'un SNP
 **/
int filterLowComplexity2Paths (char* seq1, char *seq2, int lenseq, int threshold)
{
	int DUSTSCORE1[64], DUSTSCORE2[64];
	int i, j, k, s1, m, s2;
    
	for (i=0; i<64; i++)
    {
        DUSTSCORE1[i]=0;
        DUSTSCORE2[i]=0;
        
    }
    
	for (j=2; j<lenseq; ++j)
	{
		k = seq1[j-2]*16 + seq1[j-1]*4 + seq1[j];
		++DUSTSCORE1[k];
        k = seq2[j-2]*16 + seq2[j-1]*4 + seq2[j];
		++DUSTSCORE2[k];
	}
    
	s1=0;s2=0;
	for (i=0; i<64; ++i)
	{
		m = DUSTSCORE1[i];
		s1 = s1 + (m*(m-1))/2;
        m = DUSTSCORE2[i];
		s2 = s2 + (m*(m-1))/2;
	}
    
    /*if (s1+s2>=threshold) return s1+s2;
     else return -1;*/
    return s1+s2;
}
