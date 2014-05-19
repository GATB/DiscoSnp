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
    Bubble (const Graph& graph) : graph(graph)  {}

    const Graph& graph;

    Node begin[2];
    Node end  [2];

    int score;


    size_t index;
    int where_to_extend;

    int closureRight;
    int closureLeft;

    Path extensionRight;
    Path extensionLeft;

    size_t divergenceLeft;
    size_t divergenceRight;

    Sequence seq1;
    Sequence seq2;
};

/********************************************************************************/

#endif /* _TOOL_BUBBLE_HPP_ */

