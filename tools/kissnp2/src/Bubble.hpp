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
    Bubble (const Graph& graph)
        : graph(graph), sizeKmer(graph.getKmerSize()),
          path1(2*graph.getKmerSize()-1), path2(2*graph.getKmerSize()-1), index(0), score(0)
    {
    }

    /** */
    void init (const Node& starter)
    {
        pathsSet = false;
        start1 = starter;
    }

    /** */
    bool mutate  (char nt)
    {
        /** We build the mutated node from the given nucleotide. */
        start2 = graph.mutate (start1, sizeKmer-1, (Nucleotide)nt);

        /** We check whether the mutation is in the graph or not. If not, direct return. */
        if (graph.contains(start2) == false) { return false; }

        /** We may have to configure both paths if not set. Note the optimization with the boolean. */
        if (pathsSet==false)
        {
            graph.copy (start1,path1.data());
            graph.copy (start2,path2.data());
            pathsSet=true;
        }

        /** We update the last mutated nucleotide of the mutation branch. */
        path2[sizeKmer-1] = nt;

        /** The mutation exists in the graph, so we return true. */
        return true;
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
    bool pathsSet;

    int nucleotideRight;
    int nucleotideLeft;

    Path extensionRight;
    Path extensionLeft;

    int where_to_extend;

    int score;

    size_t index;

    Sequence seq1;
    Sequence seq2;

    size_t divergenceLeft;
    size_t divergenceRight;
};

/********************************************************************************/

#endif /* _TOOL_BUBBLE_HPP_ */

