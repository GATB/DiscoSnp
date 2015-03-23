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

#include <Filter.hpp>
#include <gatb/gatb_core.hpp>
using namespace std;

static int NT2int(char nt)  {  return (nt>>1)&3;  }


/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  : True if the sequence is of high complexity else return false
 ** REMARKS : TODO: comment and validate this function.
 *********************************************************************/
bool filterLowComplexityPath(const std::string& seq){
    int DUSTSCORE[64]; // all tri-nucleotides
    for (int i=0; i<64; i++) DUSTSCORE[i]=0;
    size_t lenseq = seq.size();
    for (int j=2; j<lenseq; ++j) ++DUSTSCORE[NT2int(seq[j-2])*16 + NT2int(seq[j-1])*4 + NT2int(seq[j])];
    int m,s=0;
    
    for (int i=0; i<64; ++i)
    {
        m = DUSTSCORE[i];
        s  += (m*(m-1))/2;
    }
    
    return s<((lenseq-2)/4 * (lenseq-6)/4)/2;
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  : True if the two sequences are of high complexity else return false
 ** REMARKS : The smaller the score, the higher the complexity
 *********************************************************************/
bool filterLowComplexity2Paths (const std::string& seq1, const std::string& seq2)
{
    return filterLowComplexityPath(seq1) && filterLowComplexityPath(seq2);
}
