/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
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

#include <Filter.hpp>

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
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

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
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
