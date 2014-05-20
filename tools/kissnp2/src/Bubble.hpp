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

/** Forward declaraions. */
class Kissnp2;

/********************************************************************************/

class BubbleFinder
{
public:
    BubbleFinder (Kissnp2& tool, const Graph& graph);

    BubbleFinder (const BubbleFinder& bf);

    ~BubbleFinder ();

    void operator() (const Node& node);

private:

    Kissnp2&     tool;
    const Graph& graph;
    size_t       sizeKmer;

    /** Current Bubble instance built by this BubbleFinder instance. */
    Bubble bubble;

    /** */
    Terminator* _terminator;
    void setTerminator (Terminator* terminator)  { SP_SETATTR(terminator); }

    /** */
    Traversal* _traversal;
    void setTraversal (Traversal* traversal) { SP_SETATTR(traversal); }

    /** Start a bubble detection from a given node. This node is mutated (its last nucleotide)
     * in a second node, and so this couple of nodes is set as the starting branch of a potential
     * bubble.
     * \param[in] node : the starting node. */
    void start (Bubble& bubble, const Node& node);

    /** Extension of a bubble by testing extensions from both branches of the bubble.
     * \param[in] pos : position of the nucleotide to be added to both branches.
     */
    void expand (
        int pos,
        Bubble& bubble,
        const Node& node1,
        const Node& node2,
        const Node& previousNode1,
        const Node& previousNode2
    );

    /** Extend the bubble to the left/right with a small assembly part of the de Bruijn graph.
     * \return true if the bubble has been extended, false otherwise. */
    bool extend (Bubble& bubble);

    /** Finish the bubble, ie output the pair of sequences in the output bank.
     * \param[in] bubble: bubble to be dumped in the output bank
     */
    void finish (Bubble& bubble);

    /** Check whether new node is similar to two previous nodes.
     * \param[in] previous : previous node
     * \param[in] current : current node
     * \param[in] next : next node
     * \return true if next node is different to current and previous nodes.*/
    bool checkNodesDiff (const Node& previous, const Node& current, const Node& next) const;

    /** Check whether the first kmer of the first path is smaller than the first kmer
     * of the revcomp(first path), this should avoid repeated SNPs
     * \param[in] path : branch of a bubble.
     * \return true if first path is lower than last reverse path. */
    bool checkPath (Bubble& bubble) const;

    /** Check bubble according to user choice for branching.
     * \param[in] node 1 : bubble branch last node
     * \param[in] node 2 : bubble branch last node
     * \return true if bubble is ok */
    bool checkBranching (const Node& node1, const Node& node2) const;

    /** Check complexity for a bubble. Also returns a complexity score.
     * \param[in] path1 : branch of the bubble
     * \param[in] path2 : branch of the bubble
     * \param[out] score : complexity score
     * \return true if the complexity is ok
     */
    bool checkLowComplexity (Bubble& bubble) const;

    /** Fill a Sequence object for a given branch of a bubble.
     * \param[in] path : branch of the bubble.
     * \param[in] type : used to set the comment part of the sequence (likely 'higher' or 'lower')
     * \param[in] score : score for the bubble.
     * \param[in] where_to_extend : got from 'extend' method.
     * \param[in] seqIndex : index of the sequence (more exactly index for the pair of sequences)
     * \param[out] seq : sequence to be filled
     */
    void buildSequence (Bubble& bubble, size_t pathIdx, const char* type, Sequence& seq);

    /** */
    bool two_possible_extensions_on_one_path (const Node& node) const;
    bool two_possible_extensions (Node node1, Node node2) const;
};

/********************************************************************************/

#endif /* _TOOL_BUBBLE_HPP_ */

