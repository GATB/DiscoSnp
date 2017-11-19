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

#ifndef _TOOL_BUBBLE_HPP_
#define _TOOL_BUBBLE_HPP_

/********************************************************************************/
#include <gatb/gatb_core.hpp>
#include <string>
/********************************************************************************/

/** We define string constants for command line options. */
#define STR_MAX_AMBIGOUS_INDELS             "-max_ambigous_indel"
#define STR_DISCOSNP_LOW_COMPLEXITY         "-l"
#define STR_DISCOSNP_AUTHORISED_BRANCHING   "-b"
#define STR_DISCOSNP_TRAVERSAL_UNITIG       "-t"
#define STR_DISCOSNP_TRAVERSAL_CONTIG       "-T"
#define STR_KISSNP2_COVERAGE_FILE_NAME      "-coverage_file"
#define STR_KISSNP2_DONT_OUTPUT_FIRST_COV   "-dont_output_first_coverage"

#define STR_MAX_INDEL_SIZE                  "-D"
#define STR_MAX_POLYMORPHISM                "-P"
#define STR_MAX_SYMMETRICAL_CROSSROADS      "-max_symmetrical_crossroads"
#define STR_RADSEQ                          "-x"


/********************************************************************************/

/** We define a structure holding all the information about a bubble. */
struct Bubble
{
    // Here is a description of the attributes defining a bubble.
    //
    //    begin      end
    //   ------X   X------    <-- branch1   (X is the common nucleotide between the 2 nodes)
    //   ------Y   Y------    <-- branch2   (Y is the common nucleotide between the 2 nodes)
    //
    //  => the bubble is defined by two branches
    //
    //  => a branch is defined by two nodes [begin,end] that overlap on one nucleotide
    //     so the path length defined by these two nodes is (2*kmerSize - 1)
    Node begin[2];
    Node end  [2];
    
    
    // A Bubble may represent a SNP or a deletion (currently).
    // SNP case:
    // ---A
    //    A---
    //   and
    // ---T
    //    T---
    // We have 2 upper nodes (begin[0] and end[0]) overlaping by 1 nucleotides
    // We have 2 lower nodes (begin[1] and end[1]) overlaping by 1 nucleotides
    // Thus in the SNP case, size_overlap[0] == size_overlap[1] == 1
    
    // In the DELETION case:
    //   123456
    // 123IiiiY456 (I is the insertion here of length 5)
    // The output is : one path:
    // 1234 (k=4)
    //   3456
    // other path:
    // 123I and Y456
    // In this case we store
    //      begin[0]=1234 end[0]=3456 size_overlap[0]=2
    //  and begin[1]=123I end[1]=Y456 size_overlap[1]=0 and we store the insertion removing hte first and last character (iii in this example).
    
    
    // strings storing central information coming after the first node of a path.
    std::string extended_string[2];
    
    std::string polymorphism_type;
    
    // is this bubble of high complexity.
    bool high_complexity;
    
    /** indicates if the bubble passes the complexity filter (high complexity or don't care about complexity) */
    bool acceptable_complexity;
    
    
    
    /** indicates when a predicted bubble is finished,  if it is canonical */
    bool isCanonical;
    
    /** indicates that a bubble had been closed. Maybe it was not dumped if the bubble was not canonical */
    bool closed_bubble;
    
    int type; // 0 = isolated SNP, 1 = isolated insertion, ... other to come
    
    
    int final_nb_polymorphism; // number of SNPs in a bubble, could be one (isolated SNP or an insertion) or more (an indel+n SNPs (n+1)) or n SNPs (n)
    
    // Index of the bubble
    size_t index;
    
    // Information about the bubble extension:  0=nothing, 1=left only, 2=right only, 3=both
    int where_to_extend;
    
    // Closing right nucleotide of the bubble
    Nucleotide closureRight;
    
    // Closing left nucleotide of the bubble
    Nucleotide closureLeft;
    
    // Unitig/Contig path on the right of the bubble.
    Path extensionRight;
    
    // Unitig/Contig path on the left of the bubble.
    Path extensionLeft;
    
    // Length of the right and left unitigs
    //  - if the traversal is a simple traversal, then this length is equal to the length of the extension
    //  - else (if the traversal is a monument traversal), then this length is equal to the starting position of the first bubble (if exist))
    size_t divergenceLeft;
    size_t divergenceRight;
    
    // Sequence instances to be configured and then dumped in the output bank.
    Sequence seq1;
    Sequence seq2;
};


/********************************************************************************/

/** \brief class that tries to build a bubble from a starting node
 *
 * This class does all the job for expanding a bubble from a starting node if possible and
 * potentially extend it with right and left unitigs/contigs.
 *
 * This class is intended to be instantiated N times, one per thread.
 *
 * The starting node is provided to the operator() (so, one can see this class as a functor).
 */
class BubbleFinder
{
public:
    
    /** We define a structure gathering information during bubble detection. */
    struct Stats
    {
        Stats ()  : nb_bubbles(0), nb_bubbles_snp(0), nb_bubbles_snp_high(0), nb_bubbles_snp_low(0), nb_bubbles_del(0), nb_bubbles_del_high(0), nb_bubbles_del_low(0)  { memset (nb_where_to_extend_snp, 0, sizeof(nb_where_to_extend_snp)); memset (nb_where_to_extend_del, 0, sizeof(nb_where_to_extend_del)); }
        
        size_t nb_bubbles;
        
        size_t nb_bubbles_snp;
        size_t nb_bubbles_snp_high;
        size_t nb_bubbles_snp_low;
        size_t nb_where_to_extend_snp[4];
        
        
        
        size_t nb_bubbles_del;
        size_t nb_bubbles_del_high;
        size_t nb_bubbles_del_low;
        size_t nb_where_to_extend_del[4];
    };
    
    /** Constructor. */
    BubbleFinder (IProperties* props, const Graph& graph, Stats& stats);
    
    /** Copy constructor. It will be used when cloning N times the instance by the dispatcher
     * \param[in] bf : instance to be copied.*/
    BubbleFinder (const BubbleFinder& bf);
    
    /** Destructor. */
    ~BubbleFinder ();
    
    /** Starting method that gets a node as argument and tries to build a bubble from it.
     * NOTE: defined as a template, because we can use either Node or BranchingNode instances
     * as starting nodes. See also 'start' method, which is template too, with too possible
     * specializations: one for Node, one for BranchingNode
     * \param[in] node : the starting node. */
    template<typename T>
    void operator() (const T& node)
    {
        /** We start the SNP in both directions (forward and reverse). */
        start (bubble, node);
        start (bubble, graph.reverse(node));
    }
    
    /** Get a properties object with the configuration of the finder. */
    IProperties* getConfig () const;
    
    /** Constants */
    static const char* STR_BFS_MAX_DEPTH;
    static const char* STR_BFS_MAX_BREADTH;
protected:
    
    /** */
    const Graph& graph;
    
    /** Statistics about the bubbles lookup. */
    Stats& stats;
    
    /** Current Bubble instance built by this BubbleFinder instance. */
    Bubble bubble;
    
    /** Shortcut attribute for the kmer size of the de Bruijn graph. */
    size_t sizeKmer;
    
    /** Maximal number of polymorphism per bubble. Isolated = zero **/
    int max_polymorphism;
    
    /** Max deletion size **/
    int max_indel_size;
    
    /** Max indel size **/
    int max_indel_ambiguity;
    
    bool accept_low; // Option set: do we accept low complexity bubbles
    bool accept_truncated_bubbles; //CHARLOTTE Option set: we accept truncated bubbles (for radseq sequencing)
    
    std::queue<std::pair<Node,std::string> > breadth_first_queue;
    
    
    int max_recursion_depth;
    int current_recursion_depth;
    
    int max_depth;   // for unitigs/contigs extensions
    int max_breadth; // for unitigs/contigs extensions
    
    /* authorised_branching =
     *   0: branching forbidden in any path
     *   1: same branching on both path forbidden (i.e. 2 distinct nucleotides may be used in both paths for extension)
     *   2: no restriction on branching */
    int authorised_branching;
    
    
    /** In b 2: maximaml number of symetrically branches traversed while walking the bubble**/
    int max_sym_branches;
    
    
    
    /** Gives the kind of traversal to be done at left/right of the bubble. */
    TraversalKind traversalKind;
    
    /** Output bank of the bubbles (as a pair of sequences). Note here: we use the IBank
     * interface here, and not a specific implementation (like BankFasta), so we could
     * deal with different kinds of banks. */
    IBank* _outputBank;
    void setOutputBank (IBank* outputBank)  { SP_SETATTR(outputBank); }
    
    /** We need a synchronizer for dumping the sequences into the output bank. */
    ISynchronizer* _synchronizer;
    void setSynchronizer (ISynchronizer* synchronizer)  { SP_SETATTR(synchronizer); }
    
    /** Terminator for marking branching nodes (used by the Traversal instance) */
    BranchingTerminator* _terminator;
    void setTerminator (BranchingTerminator* terminator)  { SP_SETATTR(terminator); }
    
    /** Used for computing unitigs or contigs (according to traversal kind choice) at the left
     * and right of the bubble. */
    Traversal* _traversal;
    void setTraversal (Traversal* traversal) { SP_SETATTR(traversal); }
    
    /** Start a bubble detection from a given node. This node is mutated (its last nucleotide) in a
     * second node, and so this couple of nodes is set as the starting branch of a potential bubble.
     * NOTE: defined as template, with 2 specializations: one for Node, one for BranchingNode (see cpp file)
     * \param[in] node : the starting node. */
    template<typename T>
    void start (Bubble& bubble, const T& node);
    
    /** Extension of a bubble given two nextNodes to be tested.
     *
     */
    bool expand_heart(
                      const int nb_polymorphism,
                      Node& nextNode1,
                      Node& nextNode2,
                      Node& node1,
                      Node& node2,
                      Node& previousNode1,
                      Node& previousNode2,
                      std::string local_extended_string1,
                      std::string local_extended_string2,
                      int sym_branches,
                      int stack_size,
                      const int dissyetrical_end_size
                      );
    
    /**  Extension of a single single node. Extension is non branching and stops after max_depth path. Returns true if the created path if smaller or equal to max_depth.
     *
     */
    bool expand_one_simple_path (
                                 Node& node,                       // Node currently tested
                                 string& local_extended_string,    // add nucleotides to this string
                                 const int max_depth,              // maximal size of the path
                                 int & size_extension
    );
    
    /** Extension of a bubble by testing extensions from both branches of the bubble.
     *
     */
    bool expand (const int nb_polymorphism,
                 Node node1, // In case of indels, this node is the real extended one, but we keep it at depth 1
                 Node node2, // In case of indels, this node is not extended (depth 1)
                 Node previousNode1,
                 Node previousNode2,
                 std::string local_extended_string1,
                 std::string local_extended_string2,
                 int sym_branches,
                 int stack_size);
    
    /** Extend the bubble to the left/right with a small assembly part of the de Bruijn graph.
     * \return true if the bubble has been extended, false otherwise. */
    void extend ();
    
    /** Finish the bubble, ie output the pair of sequences in the output bank.
     * \param[in] bubble: bubble to be dumped in the output bank
     */
    void finish (const int dissyetrical_end_size);
    
    /** Check whether new node is similar to two previous nodes.
     * \param[in] previous : previous node
     * \param[in] current : current node
     * \param[in] next : next node
     * \return true if next node is different to current and previous nodes.*/
    bool checkNodesDiff (Node& previous, Node& current, Node& next) const;
    
    /** Check whether the first kmer of the first path is smaller than the first kmer
     * of the revcomp(first path), this avoids repeated SNPs
     * \param[in] path : branch of a bubble.
     * \return set isCanonical to true if first path is lower than last reverse path. */
    void checkPath ();
    
    /** Check bubble according to user choice for branching.
     * \param[in] node 1 : bubble branch last node
     * \param[in] node 2 : bubble branch last node
     * \return true if bubble is ok */
    bool checkBranching (Node& node1, Node& node2,  int & sym_branches) const;
    
    /** Check that indel bubbles respect the maximal size of the position ambiguity */
    bool checkRepeatSize (string &extension1, string &extension2) const;
    
    /** Check complexity for a bubble.
     * \param[in] path1 : branch of the bubble
     * \param[in] path2 : branch of the bubble
     * set acceptable_complexity to  true if the complexity is ok or if we accept low complexity bubbles.
     */
    void checkLowComplexity ();
    
    /** Fill a Sequence object for a given branch of a bubble.
     * \param[in] path : branch of the bubble.
     * \param[in] type : used to set the comment part of the sequence (likely 'higher' or 'lower')
     * \param[in] score : score for the bubble.
     * \param[in] where_to_extend : got from 'extend' method.
     * \param[in] seqIndex : index of the sequence (more exactly index for the pair of sequences)
     * \param[out] seq : sequence to be filled
     */
    void buildSequence (size_t pathIdx, const char* type, Sequence& seq, std::string polymorphism_comments);
    
    /** */
    bool two_possible_extensions_on_one_path (Node& node) const;
    bool two_possible_extensions (Node node1, Node node2) const;
private:
    bool recursive_indel_prediction(
                                    int extended_path_id,
                                    std::string tried_extension,
                                    Node current,
                                    size_t insert_size,
                                    const char end_insertion);
    void start_snp_prediction();
    void start_indel_prediction();
};

#endif /* _TOOL_BUBBLE_HPP_ */

