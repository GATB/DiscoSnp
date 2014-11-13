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

#include <Bubble.hpp>
#include <Filter.hpp>

using namespace std;

#define DEBUG(a) //  a

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
BubbleFinder::BubbleFinder (IProperties* props, const Graph& graph, Stats& stats)
    : graph(graph), stats(stats), _terminator(0), _traversal(0), _outputBank(0), _synchronizer(0)
{
    assert (props != 0);

    /** We retrieve the kmer size. */
    sizeKmer = graph.getKmerSize();

    /** We set the complexity threshold. */
    threshold  = (sizeKmer/2-2)*(sizeKmer/2-3);

    /** We set attributes according to user choice. */
    low                  = props->get    (STR_DISCOSNP_LOW_COMPLEXITY) != 0;
    authorised_branching = props->getInt (STR_DISCOSNP_AUTHORISED_BRANCHING);

    /** We set the traversal kind. */
    traversalKind = Traversal::NONE;
    if (props->get(STR_DISCOSNP_TRAVERSAL_UNITIG) != 0)  { traversalKind = Traversal::UNITIG; }
    if (props->get(STR_DISCOSNP_TRAVERSAL_CONTIG) != 0)  { traversalKind = Traversal::CONTIG; }

    /** We set the name of the output file. */
    stringstream ss;
    ss << props->getStr(STR_URI_OUTPUT)  << "_k_" << sizeKmer  << "_c_" << graph.getInfo().getInt("abundance");
    ss << ".fa";

    /** We set the output file. So far, we force FASTA output usage, but we could make it configurable. */
    setOutputBank (new BankFasta (ss.str()));

    /** We need a synchronizer for dumping high/low sequences into the output bank in an atomic way
     * (ie avoids potential issues with interleaved high/low sequences in multithread execution). */
    setSynchronizer (System::thread().newSynchronizer());
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
BubbleFinder::BubbleFinder (const BubbleFinder& bf)
    :  graph(bf.graph), stats(bf.stats), _terminator(0), _traversal(0), _outputBank(0), _synchronizer(0)
{
    sizeKmer             = bf.sizeKmer;
    threshold            = bf.threshold;
    low                  = bf.low;
    authorised_branching = bf.authorised_branching;
    traversalKind        = bf.traversalKind;

    /** Copy by reference (not by value). */
    setOutputBank   (bf._outputBank);
    setSynchronizer (bf._synchronizer);

    /** NOT A TRUE COPY: each instance created by this constructor will have its own
     *  Traversal/Terminator instances. */
    setTerminator (new BranchingTerminator(graph));
    setTraversal  (Traversal::create (traversalKind, graph, *_terminator));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
BubbleFinder::~BubbleFinder ()
{
    setOutputBank   (0);
    setSynchronizer (0);
    setTerminator   (0);
    setTraversal    (0);
}


/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<>
void BubbleFinder::start (Bubble& bubble, const BranchingNode& node)
{
    DEBUG ((cout << "[BubbleFinder::start] BRANCHING NODE " << graph.toString(node) << endl));

    /** We compute the successors of the node. */
    Graph::Vector<Node> successors = graph.successors<Node> (node);

    for (size_t i=0; i<successors.size(); i++)
    {
        bubble.begin[0] = successors[i];

        // In case two or more branching nodes lead to the same extensions : (b1 -> s1 and b1 -> s2 and b2 -> s1 and b2 -> s2), then
        // we need to construct the bubble s1... s2... from only one of the two branching nodes b1 or b2.
        // We chose the smallest from all the possible starting nodes to start a bubble.
        Graph::Vector<Node> predecessors = graph.predecessors<Node> (successors[i]);
        if (predecessors.size()>1)
        {
            for (size_t k=0; k<predecessors.size(); k++) { if (predecessors[k].kmer < node.kmer) { return; } }
        }

        for (size_t j=i+1; j<successors.size(); j++)
        {
            bubble.begin[1] = successors[j];
            expand (1, bubble, bubble.begin[0], bubble.begin[1], Node(~0), Node(~0));
        }
    }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void BubbleFinder::expand (
    int pos,
    Bubble& bubble,
    const Node& node1,
    const Node& node2,
    const Node& previousNode1,
    const Node& previousNode2
)
{
    /** A little check won't hurt. */
    assert (pos <= sizeKmer-1);

    /** We may have to stop the extension according to the branching mode. */
    if (checkBranching(node1,node2) == false)  { return; }

    /** We get the common successors of node1 and node2. */
    /** Returns the successors of two nodes, ie with the same transition nucleotide from both nodes. */
    Graph::Vector < pair<Node,Node> > successors = graph.successors<Node> (node1, node2);
    
    /** A second little check won't hurt. */
    /** We should not have several extensions possible unless authorised_branching==2 */
    assert(authorised_branching==2 || successors.size==1);

    /** We loop over the successors of the two nodes. */
    for (size_t i=0; i<successors.size(); i++)
    {
        /** Shortcuts. */
        Node& nextNode1 = successors[i].first;
        Node& nextNode2 = successors[i].second;

        /** We check whether the new nodes are different from previous ones. */
        bool checkPrevious =
            checkNodesDiff (previousNode1, node1, nextNode1) &&
            checkNodesDiff (previousNode2, node2, nextNode2);

        if (!checkPrevious)  { continue; }

        /************************************************************/
        /**                   RECURSION CONTINUES                  **/
        /************************************************************/
        if (pos < sizeKmer-1)
        {
            DEBUG((cout<<pos<<" nextNode1.value "<<graph.toString(nextNode1)<<" nextNode2.value "<<graph.toString(nextNode2)<<endl));
            /** We call recursively the method (recursion on 'pos'). */
            expand (pos+1, bubble,  nextNode1, nextNode2,  node1, node2);

            /** There's only one branch to expand if we keep non branching SNPs only, therefore we can safely stop the for loop */
            if ( authorised_branching==0 || authorised_branching==1 )   {  break; }
        }

        /************************************************************/
        /**                   RECURSION FINISHED                   **/
        /************************************************************/
        else
        {
            DEBUG((cout<<"last "<<pos<<" nextNode1.value "<<graph.toString(nextNode1)<<" nextNode2.value "<<graph.toString(nextNode2)<<endl));
            /** We check the branching properties of the next kmers. */
            if (checkBranching(nextNode1, nextNode2)==false)  { return; }
            
            /** We check that the two last kmers can be right extended, leading to at least a common successor */
            /** As we start a bubble from a branching node, we also have to symetrically check than the bubble is right closed */
            if (graph.successors<Node> (nextNode1, nextNode2).size()<1) { return; }
            /** We finish the bubble with both last nodes. */
            bubble.end[0] = nextNode1;
            bubble.end[1] = nextNode2;

            /** We check several conditions (the first path vs. its revcomp and low complexity). */
            if (checkPath(bubble)==true && checkLowComplexity(bubble)==true)
            {
                /** We extend the bubble on the left and right (unitigs or contigs). */
                if (extend (bubble) == true)
                {
                    /** We got all the information about the bubble, we finish it. */
                    finish (bubble);
                }
            }
        }

    } /* end of for (size_t i=0; i<successors.size(); i++) */
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
bool BubbleFinder::extend (Bubble& bubble)
{
    Nucleotide closureLeft  = NUCL_UNKNOWN;
    Nucleotide closureRight = NUCL_UNKNOWN;

    /** We may have to extend the bubble according to the user choice. */
    if (traversalKind != Traversal::NONE)
    {
        /** We ask for the predecessors of the first node and successors of the last node. */
        Graph::Vector<Node> successors   = graph.successors<Node>   (bubble.end[0]);
        Graph::Vector<Node> predecessors = graph.predecessors<Node> (bubble.begin[0]);

        /** We need to reset branching nodes between extensions in case of overlapping extensions. */
        _terminator->reset ();

        /** If unique, we keep the left/right extensions. */
        if (successors.size()==1)
        {
            /** We compute right extension of the node. */
            closureRight = graph.getNT (successors  [0], sizeKmer-1);
            _traversal->traverse (successors[0], DIR_OUTCOMING, bubble.extensionRight);
            bubble.divergenceRight = _traversal->getBubbles().empty() ? bubble.extensionRight.size() : _traversal->getBubbles()[0].first;
        }

        if (predecessors.size()==1)
        {
            /** We compute left extension of the node. */
            closureLeft  = graph.getNT (predecessors[0], 0);
            _traversal->traverse (graph.reverse(predecessors[0]), DIR_OUTCOMING, bubble.extensionLeft);
            bubble.divergenceLeft = _traversal->getBubbles().empty() ? bubble.extensionLeft.size() : _traversal->getBubbles()[0].first;
        }
    }

    /** We return a code value according to left/right extensions status. */
         if (closureLeft==NUCL_UNKNOWN && closureRight==NUCL_UNKNOWN)  { bubble.where_to_extend = 0; }
    else if (closureLeft!=NUCL_UNKNOWN && closureRight==NUCL_UNKNOWN)  { bubble.where_to_extend = 1; }
    else if (closureLeft==NUCL_UNKNOWN && closureRight!=NUCL_UNKNOWN)  { bubble.where_to_extend = 2; }
    else if (closureLeft!=NUCL_UNKNOWN && closureRight!=NUCL_UNKNOWN)  { bubble.where_to_extend = 3; }

    bubble.closureLeft  = closureLeft;
    bubble.closureRight = closureRight;

    return true;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void BubbleFinder::finish (Bubble& bubble)
{
    /** We set the bubble index. NOTE: We have to make sure it is protected against concurrent
     * accesses since we may be called here from different threads. */
    bubble.index = __sync_add_and_fetch (&(stats.nb_bubbles), 1);

    /** We build two Sequence objects from the information of the bubble. */
    buildSequence (bubble, 0, "higher", bubble.seq1);
    buildSequence (bubble, 1, "lower",  bubble.seq2);

    /** We have to protect the sequences dump wrt concurrent accesses. We use a {} block with
     * a LocalSynchronizer instance with the shared ISynchonizer of the Kissnp2 class. */
    {
        LocalSynchronizer sync (_synchronizer);

        /** We insert the two sequences into the output bank. */
        _outputBank->insert (bubble.seq1);
        _outputBank->insert (bubble.seq2);

        /** Stats update (in concurrent access protection block). */
        stats.nb_where_to_extend[bubble.where_to_extend] ++;

        if (bubble.score < threshold)  { stats.nb_bubbles_high++; }
        else                           { stats.nb_bubbles_low++;  }
    }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
bool BubbleFinder::two_possible_extensions_on_one_path (const Node& node) const
{
    return graph.indegree(node)>=2 || graph.outdegree(node)>=2;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
bool BubbleFinder::two_possible_extensions (Node node1, Node node2) const
{
    return
        graph.successors<Edge> (node1, node2).size() >= 2  ||
        graph.successors<Edge> (graph.reverse (node1),graph.reverse (node2)).size() >= 2;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void BubbleFinder::buildSequence (Bubble& bubble, size_t pathIdx, const char* type, Sequence& seq)
{
    stringstream commentStream;

    /** We build the comment for the sequence. */
    commentStream << "SNP_" << type << "_path_" << bubble.index << "|" << (bubble.score >= threshold ? "low" : "high");

    /** We may have extra information for the comment. */
    if (traversalKind == Traversal::UNITIG)
    {
        commentStream << "|left_unitig_length_";
        commentStream << (bubble.where_to_extend%2==1 ? (bubble.extensionLeft.size() +1) : 0);
        commentStream << "|right_unitig_length_";
        commentStream << (bubble.where_to_extend>1    ? (bubble.extensionRight.size()+1) : 0);
    }
    if (traversalKind == Traversal::CONTIG)
    {
        commentStream << "|left_unitig_length_";
        commentStream << (bubble.where_to_extend%2==1 ? (bubble.divergenceLeft +1) : 0);
        commentStream << "|right_unitig_length_";
        commentStream << (bubble.where_to_extend>1    ? (bubble.divergenceRight+1) : 0);

        commentStream << "|left_contig_length_";
        commentStream << (bubble.where_to_extend%2==1 ? (bubble.extensionLeft.size() +1) : 0);
        commentStream << "|right_contig_length_";
        commentStream << (bubble.where_to_extend>1    ? (bubble.extensionRight.size()+1) : 0);
    }

    /** We assign the comment of the sequence. */
    seq.setComment (commentStream.str());

    size_t lenLeft  = bubble.extensionLeft.size ();
    size_t lenRight = bubble.extensionRight.size ();
    size_t len      = (2*sizeKmer-1);

    if (bubble.closureLeft  != NUCL_UNKNOWN)  { len += 1 + lenLeft;  }
    if (bubble.closureRight != NUCL_UNKNOWN)  { len += 1 + lenRight; }

    /** We resize the sequence data if needed. Note: +1 for ending '\0'
     * NOTE: we use resize if we need more space, setSize otherwise. */
    if (seq.getData().size() < len+1)  {  seq.getData().resize  (len+1); }
    else                               {  seq.getData().setSize (len+1); }

    char* output = seq.getDataBuffer();

    /** We add the left extension if any. Note that we use lower case for extensions. */
    if (bubble.closureLeft != NUCL_UNKNOWN)
    {
        for (size_t i=0; i<lenLeft; i++)  {  *(output++) = tolower(ascii (reverse(bubble.extensionLeft [lenLeft-i-1])));  }
        *(output++) = tolower(ascii(bubble.closureLeft));
    }

    /** We add the bubble path. */
    string begin = graph.toString (bubble.begin[pathIdx]);
    string end   = graph.toString (bubble.end[pathIdx]);

    for (size_t i=0; i<sizeKmer-1; i++)  {  *(output++) = begin[i];  }
    for (size_t i=0; i<sizeKmer;   i++)  {  *(output++) = end  [i];  }

    /** We add the right extension if any. Note that we use lower case for extensions. */
    if (bubble.closureRight != NUCL_UNKNOWN)
    {
        *(output++) =  tolower(ascii(bubble.closureRight));
        for (size_t i=0; i<lenRight; i++)  {  *(output++) = tolower(ascii (bubble.extensionRight[i]));  }
    }

    /** We add a null terminator for the strings. */
    *(output++) = '\0';
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
bool BubbleFinder::checkNodesDiff (const Node& previous, const Node& current, const Node& next) const
{
    return (next.kmer != current.kmer) && (next.kmer != previous.kmer);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
bool BubbleFinder::checkPath (Bubble& bubble) const
{
    /** We test whether the first kmer of the first path is smaller than
     * the first kmer of the revcomp(first path), this should avoid repeated SNPs */
    return graph.toString (bubble.begin[0])  <  graph.toString (graph.reverse(bubble.end[0]));
}

/*********************************************************************
** METHOD  : BubbleFinder::checkBranching
** PURPOSE : Checks the branching properties of a couple of nodes. If authorised_branching==0: no branching is authorized in any of the two nodes.
**           If authorised_branching==1: no symetrical branching authorized from the two nodes (ie. N1-->A and N1-->B and N2-->A and N2-->B not authorized)
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
bool BubbleFinder::checkBranching (const Node& node1, const Node& node2) const
{
    // stop the extension if authorised_branching==0 (not branching in any path) and any of the two paths is branching
    if (authorised_branching==0 && (two_possible_extensions_on_one_path(node1) || two_possible_extensions_on_one_path(node2)))
    {
        return false;
    }

    // stop the extension if authorised_branching==1 (not branching in both path) and both the two paths are branching
    if (authorised_branching==1 && two_possible_extensions (node1, node2))
    {
        return false;
    }

    return true;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
bool BubbleFinder::checkLowComplexity (Bubble& bubble) const
{
    string path1 = graph.toString (bubble.begin[0]).substr(0, sizeKmer-1) + graph.toString (bubble.end[0]);
    string path2 = graph.toString (bubble.begin[1]).substr(0, sizeKmer-1) + graph.toString (bubble.end[1]);

    /** We compute the low complexity score of the two paths. */
    bubble.score = filterLowComplexity2Paths (path1, path2);

    return (bubble.score < threshold || (bubble.score>=threshold && low));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
IProperties* BubbleFinder::getConfig () const
{
    IProperties* props = new Properties();

    /** We aggregate information for user. */
    props->add (0, "config",   "");
    props->add (1, "kmer_size",   "%d", sizeKmer);
    props->add (1, "threshold",   "%d", threshold);
    props->add (1, "auth_branch", "%d", authorised_branching);
    props->add (1, "low",         "%d", low);
    props->add (1, "traversal",   "%s", Traversal::getName(traversalKind));

    return props;
}
