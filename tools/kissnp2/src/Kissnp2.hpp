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

#ifndef _TOOL_KISSNP2_HPP_
#define _TOOL_KISSNP2_HPP_

/********************************************************************************/
#include <gatb/gatb_core.hpp>
#include <Bubble.hpp>
/********************************************************************************/

/********************************************************************************/
/** \brief Kissnp2 that looks for SNP
 *
 * The BubbleFinder class aims to find SNP in a provided de Bruijn graph.
 * The output is a bank with pairs of sequences defining a bubble.
 *
 * This class may be inherited by refining some methods; for instance, the pairs
 * of sequences in the output bank may have some left/right extensions from the
 * found bubbles.
 */
class Kissnp2 : public Tool
{
public:

    /** Constructor. */
    Kissnp2 ();

    /** Destructor. */
    ~Kissnp2 ();

    /** Implementation of Algorithm::execute method. */
    void execute ();

    /** Start a bubble detection from a given node. This node is mutated (its last nucleotide)
     * in a second node, and so this couple of nodes is set as the starting branch of a potential
     * bubble.
     * \param[in] node : the starting node. */
    void start (Bubble& bubble, const Node& node);

    /** Getter for the graph.
     * \return the de Bruign graph */
    const Graph& getGraph() const  { return graph; }

protected:

    /** We define a typedef for the nucleotide integer values: A=0, C=1, T=2, G=3 */
    typedef char NT;
    typedef NT* PATH;

    /** */
    void configure ();

    /** Extend the bubble to the left/right with a small assembly part of the de Bruijn graph.
     * \return -1 not closed, 0 no unique extension, 1 only left, 2 only right, 3 both */
    bool extend (Bubble& bubble);

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

    /** */
    void finish (Bubble& bubble);

    /** De Bruijn graph, likely built by the 'dbgh5' binary from GATB-CORE. */
    Graph graph;

    /** Shortcut attribute for the kmer size of the de Bruijn graph. */
    size_t sizeKmer;

    /** Threshold (computed from the kmer size). */
    int threshold;

    /** Output bank of the bubbles (as a pair of sequences). Note here: we use the IBank
     * interface here, and not a specific implementation (like BankFasta), so we could
     * deal with different kinds of banks. */
    IBank* _outputBank;
    void setOutputBank (IBank* outputBank)  { SP_SETATTR(outputBank); }

    bool low;

    /* authorised_branching =
    *   0: branching forbidden in any path
    *   1: same branching on both path forbidden (i.e. 2 distinct nucleotides may be used in both paths for extension)
    *   2: no restriction on branching */
    int authorised_branching;

    bool print_extensions;
    int  min_size_extension;

    /** We need a synchronizer for dumping the sequences into the output bank. */
    ISynchronizer* _synchronizer;
    void setSynchronizer (ISynchronizer* synchronizer)  { SP_SETATTR(synchronizer); }

    /** Statistics abouth the bubbles lookup. */
    size_t nb_bubbles;
    size_t nb_bubbles_high;
    size_t nb_bubbles_low;

    size_t nb_where_to_extend[4];

    /** */
    enum TraversalKind { NONE=0, UNITIG=1, CONTIG=2 };
    TraversalKind traversalKind;

    /** */
    gatb::core::debruijn::impl::Terminator* _terminator;
    void setTerminator (gatb::core::debruijn::impl::Terminator* terminator)  { SP_SETATTR(terminator); }

    /** */
    gatb::core::debruijn::impl::Traversal* _traversal;
    void setTraversal (gatb::core::debruijn::impl::Traversal* traversal) { SP_SETATTR(traversal); }

    /** Fill a Sequence object for a given branch of a bubble.
     * \param[in] path : branch of the bubble.
     * \param[in] type : used to set the comment part of the sequence (likely 'higher' or 'lower')
     * \param[in] score : score for the bubble.
     * \param[in] where_to_extend : got from 'extend' method.
     * \param[in] seqIndex : index of the sequence (more exactly index for the pair of sequences)
     * \param[out] seq : sequence to be filled
     */
    void buildSequence (Bubble& bubble, size_t pathIdx, const char* type, Sequence& seq);

    bool two_possible_extensions_on_one_path (const Node& node) const;
    bool two_possible_extensions (Node node1, Node node2) const;

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
};

/********************************************************************************/

#endif /* _TOOL_KISSNP2_HPP_ */

