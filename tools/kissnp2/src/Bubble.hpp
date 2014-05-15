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

#ifndef _TOOL_BUBBLE_HPP_
#define _TOOL_BUBBLE_HPP_

/********************************************************************************/
#include <gatb/gatb_core.hpp>
#include <vector>
/********************************************************************************/

struct Bubble
{
    Bubble (const Graph& graph, const Node& starter)
        : graph(graph), sizeKmer(graph.getKmerSize()),
          start1(starter), start2(starter),
          path1(2*graph.getKmerSize()-1), path2(2*graph.getKmerSize()-1)
    {
        /** We copy the kmer in path1 and path2*/
        for (size_t i=0; i<sizeKmer; i++)  {  path1[i] = path2[i] = graph.getNT(starter,i);  }
    }

    /** */
    bool mutate  (char nt)
    {
        start2 = graph.mutate (start1, sizeKmer-1, (Nucleotide)nt);
        path2[sizeKmer-1] = nt;

        return graph.contains(start2);
    }

    /** */
    void set (size_t idx, Nucleotide nt)
    {
        path1[sizeKmer-1+idx] = nt;
        path2[sizeKmer-1+idx] = nt;
    }

    void finish (const Node& stop1, const Node& stop2)  { this->stop1=stop1; this->stop2=stop2; }

    const Graph& graph;

    Node  start1;
    Node  start2;

    Node  stop1;
    Node  stop2;

    size_t sizeKmer;

    std::vector<char> path1;
    std::vector<char> path2;
};

/********************************************************************************/

#endif /* _TOOL_BUBBLE_HPP_ */

