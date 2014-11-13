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
static int NT2int(char nt)  {  return (nt>>1)&3;  }

int filterLowComplexity2Paths (const std::string& seq1, const std::string& seq2)
{
    int DUSTSCORE1[64], DUSTSCORE2[64];
    int i, j, k, s1, m, s2;

    for (i=0; i<64; i++)
    {
        DUSTSCORE1[i]=0;
        DUSTSCORE2[i]=0;

    }

    size_t lenseq = seq1.size();
    for (j=2; j<lenseq; ++j)
    {
        k = NT2int(seq1[j-2])*16 + NT2int(seq1[j-1])*4 + NT2int(seq1[j]);
        ++DUSTSCORE1[k];
        k = NT2int(seq2[j-2])*16 + NT2int(seq2[j-1])*4 + NT2int(seq2[j]);
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

    return s1+s2;
}
