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
    void start (const Node& node);

    /** Getter for the graph.
     * \return the de Bruign graph */
    const Graph& getGraph() const  { return graph; }

protected:

    /** */
    void configure ();

    /** Extend the bubble to the left/right with a small assembly part of the de Bruijn graph.
     * This method is virtual and may be refined in sub classes to have specific bubble extensions.
     * \return -1 not closed, 0 no unique extension, 1 only left, 2 only right, 3 both */
    virtual int extend (const char * path1, const char * path2, char * path1_c, char * path2_c);

    /** */
    virtual void addExtraCommentInformation (
        std::stringstream& ss,
        std::pair<char*,int> left_extension,
        std::pair<char*,int> right_extension,
        int where_to_extend
    );

    /** Extension of a bubble by testing extensions from both branches of the bubble.
     * \param[in] pos : position of the nucleotide to be added to both branches.
     */
    void expand (
        int pos,
        char* path1,
        char* path2,
        const Node& node1,
        const Node& node2,
        const Node& previousNode1,
        const Node& previousNode2
    );

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

    enum ExtensionMode {  NONE=0, UNITIG=1, CONTIG=2 };
    ExtensionMode extensionMode;

    bool print_extensions;
    int  min_size_extension;

    /** We need a synchronizer for dumping the sequences into the output bank. */
    ISynchronizer* _synchronizer;
    void setSynchronizer (ISynchronizer* synchronizer)  { SP_SETATTR(synchronizer); }

    /** Statistics abouth the bubbles lookup. */
    size_t nb_bubbles;
    size_t nb_bubbles_high;
    size_t nb_bubbles_low;

    /** Fill a Sequence object for a given branch of a bubble.
     * \param[in] path : branch of the bubble.
     * \param[in] type : used to set the comment part of the sequence (likely 'higher' or 'lower')
     * \param[in] score : score for the bubble.
     * \param[in] where_to_extend : got from 'extend' method.
     * \param[in] seqIndex : index of the sequence (more exactly index for the pair of sequences)
     * \param[out] seq : sequence to be filled
     */
    void retrieveSequence (char* path, const char* type, int score, int where_to_extend, size_t seqIndex, Sequence& seq);

    bool two_possible_extensions_on_one_path (const Node& node) const;
    bool two_possible_extensions (Node node1, Node node2) const;

    /** Check whether new node is similar to two previous nodes.
     * \param[in] previous : previous node
     * \param[in] current : current node
     * \param[in] next : next node
     * \return true if next node is different to current and previous nodes.*/
    bool checkKmersDiff (const Node& previous, const Node& current, const Node& next) const;

    /** Check whether the first kmer of the first path is smaller than the first kmer
     * of the revcomp(first path), this should avoid repeated SNPs
     * \param[in] path : branch of a bubble.
     * \return true if first path is lower than last reverse path. */
    bool checkPath (char* path) const;

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
    bool checkLowComplexity (char* path1, char* path2, int& score) const;
};

/********************************************************************************/

#endif /* _TOOL_KISSNP2_HPP_ */

