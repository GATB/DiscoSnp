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

#include <Kissnp2.hpp>
#include <Filter.hpp>

using namespace std;

/********************************************************************************
 *
 * A QUICK OVERVIEW OF THE CLASS...
 *
 * This class implements the detection of SNP in a provided de Bruijn graph.
 *
 * It is implemented as a subclass of Tool, so we get all the facilities of the
 * Tool class. We can therefore see two main parts here:
 *
 * 1) constructor: we define all the command line parameters available
 *
 * 2) 'execute' method: 'main' method of the class where the main loop (ie. iteration
 *    over nodes of the graph) is done. Each node (and its revcomp) is used as the first
 *    branch of a bubble, the second branch being a mutated node of the first node.
 *    This bubble [node1,node2] is then expanded kmerSize if possible. If still ok, a few
 *    checks are done to avoid redundant computations.
 *
 * A 'configure' method is called at the beginning of 'execute' in order to configure
 * all the attributes (with command line parameters information)
 *
 * There are several 'checkXXX' methods used to avoid unwanted (or duplicate) bubbles.
 *
********************************************************************************/

/** We define string constants for command line options. */
static const char* STR_DISCOSNP_LOW_COMPLEXITY       = "-l";
static const char* STR_DISCOSNP_AUTHORISED_BRANCHING = "-b";
static const char* STR_DISCOSNP_TRAVERSAL_UNITIG     = "-t";
static const char* STR_DISCOSNP_TRAVERSAL_CONTIG     = "-T";
static const char* STR_DISCOSNP_EXTENSION_SIZE       = "-e";

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Kissnp2::Kissnp2 ()
    : Tool ("Kissnp2"), _outputBank(0), _synchronizer(0), nb_bubbles(0), nb_bubbles_high(0), nb_bubbles_low(0),
      low(false), authorised_branching(0),  print_extensions(true),
      sizeKmer (0), min_size_extension(-1), threshold(0),
      traversalKind(NONE), _terminator(0), _traversal(0)
{
    /** We add options known by kissnp2. */
    getParser()->push_front (new OptionNoParam  (STR_DISCOSNP_LOW_COMPLEXITY,       "conserve low complexity SNPs",     false));
    getParser()->push_front (new OptionOneParam (STR_DISCOSNP_AUTHORISED_BRANCHING, "branching mode\n"
            "\t0: forbid SNPs for wich any of the two paths is branching (high precision, low recall)\n"
            "\t1: forbid SNPs for wich the two paths are branching (e.g. the two paths can be created either with a 'A' or a 'C' at the same position (default value)\n"
            "\t2: No limitation on branching (low precision, high recall)",  false, "1"));
    getParser()->push_front (new OptionOneParam (STR_DISCOSNP_EXTENSION_SIZE,       "extend found SNPs  and conserve only those whose min(left and right extension) is bigger or equal to length",  false, "-1"));
    getParser()->push_front (new OptionNoParam  (STR_DISCOSNP_TRAVERSAL_UNITIG,     "extend found and stop at first polymorphism (strict extension=unitigs) SNPs. Uncompatible with -T",  false));
    getParser()->push_front (new OptionNoParam  (STR_DISCOSNP_TRAVERSAL_CONTIG,     "extend found and stop at large polymorphism (extension=contigs) SNPs. Uncompatible with -t",  false));
    getParser()->push_front (new OptionOneParam (STR_URI_OUTPUT,                    "output name",                      true));
    getParser()->push_front (new OptionOneParam (STR_URI_INPUT,                     "input file (likely a hdf5 file)",  true));

    /** Attribute initialization. */
    memset (nb_where_to_extend, 0, sizeof(nb_where_to_extend));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Kissnp2::~Kissnp2 ()
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
void Kissnp2::configure ()
{
#if 1
BankFasta::setDataLineSize (100000);
#endif

    /** We load the graph from the provided uri. */
    graph = Graph::load (getInput()->getStr(STR_URI_INPUT));

    /** We retrieve the kmer size. */
    sizeKmer = getGraph().getKmerSize();

    /** We set the threshold. */
    threshold  = (sizeKmer/2-2)*(sizeKmer/2-3);

    /** We set attributes according to user choice. */
    low                  = getInput()->get    (STR_DISCOSNP_LOW_COMPLEXITY) != 0;
    authorised_branching = getInput()->getInt (STR_DISCOSNP_AUTHORISED_BRANCHING);
    min_size_extension   = getInput()->getInt (STR_DISCOSNP_EXTENSION_SIZE);

    /** We set the traversal kind. */
    if (getInput()->get(STR_DISCOSNP_TRAVERSAL_UNITIG) != 0)  { traversalKind = UNITIG; }
    if (getInput()->get(STR_DISCOSNP_TRAVERSAL_CONTIG) != 0)  { traversalKind = CONTIG; }

    if (traversalKind != NONE)
    {
        setTerminator (new BranchingTerminator(graph));
        setTraversal  (Traversal::create (traversalKind==UNITIG ? "unitig" : "monument", graph, *_terminator));
    }

    /** We set the name of the output file. */
    stringstream ss;
    ss << getInput()->getStr(STR_URI_OUTPUT)  << "_k_" << sizeKmer  << "_c_" << getGraph().getInfo().getInt("abundance");
    if (min_size_extension > -1)  { ss << "_e_" << min_size_extension; }
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
/** We define a functor that only calls Kissnp2::start in the context of thread call
 * made through a Dispatcher object.
 * NOTE: it would be much simpler to use a lambda expression but we must be sure that this
 * code will compile with old compilers, so explicit functor struct is mandatory :(
 */
struct FunctorExecute
{
    Kissnp2& finder;
    Bubble   bubble;

    FunctorExecute (Kissnp2& finder) : finder(finder), bubble(finder.getGraph())  {}

    void operator() (const Node& node)
    {
        /** We start the SNP in both directions (forward and reverse). */
        finder.start (bubble, node);
        finder.start (bubble, finder.getGraph().reverse(node));
    }
};

/** */
void Kissnp2::execute ()
{
    /** We configure the object with command line arguments. */
    configure ();

    /** We get an iterator over the nodes of the graph. */
    ProgressGraphIterator<Node,ProgressTimer> it (getGraph().iterator<Node>(), "nodes");

    /** THIS IS THE MAIN LOOP... We launch the iteration over all the nodes of the graph. */
    IDispatcher::Status status = getDispatcher()->iterate (it, FunctorExecute(*this));

    /** We aggregate information for user. */
    getInfo()->add (1, "config",   "");
    getInfo()->add (2, "kmer_size",      "%d", sizeKmer);
    getInfo()->add (2, "threshold",      "%d", threshold);
    getInfo()->add (2, "auth_branch",    "%d", authorised_branching);
    getInfo()->add (2, "low",            "%d", low);
    getInfo()->add (2, "traversal_kind", "%d", traversalKind);

    getInfo()->add (1, "bubbles",   "");
    getInfo()->add (2, "nb",      "%lu", nb_bubbles);
    getInfo()->add (2, "nb_high", "%lu", nb_bubbles_high);
    getInfo()->add (2, "nb_low",  "%lu", nb_bubbles_low);
    getInfo()->add (2, "extensions",  "");
    getInfo()->add (3, "none",       "%d", nb_where_to_extend[0]);
    getInfo()->add (3, "left",       "%d", nb_where_to_extend[1]);
    getInfo()->add (3, "right",      "%d", nb_where_to_extend[2]);
    getInfo()->add (3, "left|right", "%d", nb_where_to_extend[3]);

    getInfo()->add (1, "time", "");
    getInfo()->add (2, "find", "%d", status.time);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Kissnp2::start (Bubble& bubble, const Node& node)
{
    /** We get the mutations of the given node at position sizeKmer-1.
     * IMPORTANT: the third argument (set to 1) tells that the allowed nucleotide
     * variants must be greater than the nucleotide at position (sizeKmer-1) of the given node.
     * => We try all the possible extensions that were not previously tested (clever :-)) */
    Graph::Vector<Node> mutations = graph.mutate (node, sizeKmer-1, 1);

    /** We initialize the first path of the bubble. */
    bubble.begin[0] = node;

    /** We loop over all mutations. */
    for (size_t i=0; i<mutations.size(); i++)
    {
        /** We initialize the second path of the bubble. */
        bubble.begin[1] = mutations[i];

        /** We try to expand this new putative bubble. */
        expand (1, bubble, node, mutations[i], Node(~0), Node(~0));
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
void Kissnp2::expand (
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
    Graph::Vector < pair<Node,Node> > successors = getGraph().successors<Node> (node1, node2);

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
            /** We call recursively the method (recursion on 'pos'). */
            expand (pos+1, bubble,  nextNode1, nextNode2,  node1, node2);

            /** There's only one branch to expand if we keep non branching SNPs only, therefore we can safely stop the for loop */
            if ( authorised_branching==0 || authorised_branching==1 )   {  break;  }
        }

        /************************************************************/
        /**                   RECURSION FINISHED                   **/
        /************************************************************/
        else
        {
            /** We check the branching properties of the next kmers. */
            if (checkBranching(nextNode1, nextNode2)==false)  { return; }

            /** We finish the bubble with both last nodes. */
            bubble.end[0] = nextNode1;
            bubble.end[1] = nextNode2;

            /** We check several conditions (the first path vs. its revcomp and low complexity). */
            if (checkPath(bubble)==true && checkLowComplexity(bubble)==true)
            {
                /** We extend the bubble on the left and right. */
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
void Kissnp2::finish (Bubble& bubble)
{
    /** We set the bubble index. NOTE: We have to make sure it is protected against concurrent
     * accesses since we may be called here from different threads. */
    bubble.index = __sync_add_and_fetch (&nb_bubbles, 1);

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
        nb_where_to_extend[bubble.where_to_extend] ++;

        if (bubble.score < threshold)  { nb_bubbles_high++; }
        else                           { nb_bubbles_low++;  }
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
bool Kissnp2::extend (Bubble& bubble)
{
    static const int BAD = -1;

    int closureLeft  = BAD;
    int closureRight = BAD;

    /** We may have to extend the bubble according to the user choice. */
    if (traversalKind != NONE)
    {
        /** We ask for the predecessors of the first node and successors of the last node. */
        Graph::Vector<Node> predecessors = graph.predecessors<Node> (bubble.begin[0]);
        Graph::Vector<Node> successors   = graph.successors<Node>   (bubble.end[0]);

        /** If unique, we keep the left/right extensions. */
        if (predecessors.size()==1)  { closureLeft  = graph.getNT (predecessors[0], 0);           }
        if (successors.size()  ==1)  { closureRight = graph.getNT (successors  [0], sizeKmer-1);  }

        _terminator->reset ();

        /** We compute left extension of the node. */
        _traversal->traverse (graph.reverse(predecessors[0]), DIR_OUTCOMING, bubble.extensionLeft);
        bubble.divergenceLeft = _traversal->getBubbles().empty() ? 0 : _traversal->getBubbles()[0].first;

        /** We compute right extension of the node. */
        _traversal->traverse (successors[0], DIR_OUTCOMING, bubble.extensionRight);
        bubble.divergenceRight = _traversal->getBubbles().empty() ? 0 : _traversal->getBubbles()[0].first;

        /** We return a code value according to left/right extensions status. */
             if (closureLeft==BAD && closureRight==BAD)  { bubble.where_to_extend = 0; }
        else if (closureLeft!=BAD && closureRight==BAD)  { bubble.where_to_extend = 1; }
        else if (closureLeft==BAD && closureRight!=BAD)  { bubble.where_to_extend = 2; }
        else if (closureLeft!=BAD && closureRight!=BAD)  { bubble.where_to_extend = 3; }
    }

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
bool Kissnp2::two_possible_extensions_on_one_path (const Node& node) const
{
    return getGraph().indegree(node)>=2 || getGraph().outdegree(node)>=2;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
bool Kissnp2::two_possible_extensions (Node node1, Node node2) const
{
    return
        getGraph().successors<Edge> (node1, node2).size() >= 2  ||
        getGraph().successors<Edge> (getGraph().reverse (node1),getGraph().reverse (node2)).size() >= 2;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Kissnp2::buildSequence (Bubble& bubble, size_t pathIdx, const char* type, Sequence& seq)
{
    stringstream commentStream;

    /** We build the comment for the sequence. */
    commentStream << "SNP_" << type << "_path_" << bubble.index << "|" << (bubble.score >= threshold ? "low" : "high");

    /** We may have extra information for the comment. */
    if (traversalKind == UNITIG)
    {
        commentStream << "|left_unitig_length_";
        commentStream << (bubble.where_to_extend%2==1 ? (bubble.extensionLeft.size() +1) : 0);
        commentStream << "|right_unitig_length_";
        commentStream << (bubble.where_to_extend>1    ? (bubble.extensionRight.size()+1) : 0);
    }
    if (traversalKind == CONTIG)
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

    static const int BAD = -1;

    size_t len = (2*sizeKmer-1) + bubble.extensionLeft.size() + bubble.extensionRight.size();

    if (bubble.closureLeft  != BAD)  { len++; }
    if (bubble.closureRight != BAD)  { len++; }

    /** We resize the sequence data if needed. Note: +1 for ending '\0'*/
    if (seq.getData().size() < len+1)  {  seq.getData().resize (len+1); }

    char* output = seq.getDataBuffer();

    size_t lenLeft  = bubble.extensionLeft.size ();
    for (size_t i=0; i<lenLeft; i++)  {  *(output++) = tolower(ascii (reverse(bubble.extensionLeft [lenLeft-i-1])));  }

    /** We add the left extension if any. Note that we use lower case for extensions. */
    if (bubble.closureLeft != BAD)   {  *(output++) = tolower(bin2NT[bubble.closureLeft]);  }

    /** We add the bubble path. */
    string begin = graph.toString (bubble.begin[pathIdx]);
    string end   = graph.toString (bubble.end[pathIdx]);

    for (size_t i=0; i<sizeKmer-1; i++)  {  *(output++) = begin[i];  }
    for (size_t i=0; i<sizeKmer; i++)    {  *(output++) = end  [i];  }

    /** We add the right extension if any. Note that we use lower case for extensions. */
    if (bubble.closureRight != BAD)  {  *(output++) =  tolower(bin2NT[bubble.closureRight]);  }

    size_t lenRight  = bubble.extensionRight.size ();
    for (size_t i=0; i<lenRight; i++)  {  *(output++) = tolower(ascii (bubble.extensionRight[i]));  }

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
bool Kissnp2::checkNodesDiff (const Node& previous, const Node& current, const Node& next) const
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
bool Kissnp2::checkPath (Bubble& bubble) const
{
    /** We test whether the first kmer of the first path is smaller than
     * the first kmer of the revcomp(first path), this should avoid repeated SNPs */
    return graph.toString (bubble.begin[0])  <  graph.toString (graph.reverse(bubble.end[0]));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
bool Kissnp2::checkBranching (const Node& node1, const Node& node2) const
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
bool Kissnp2::checkLowComplexity (Bubble& bubble) const
{
    string path1 = graph.toString (bubble.begin[0]).substr(0, sizeKmer-1) + graph.toString (bubble.end[0]);
    string path2 = graph.toString (bubble.begin[1]).substr(0, sizeKmer-1) + graph.toString (bubble.end[1]);

    /** We compute the low complexity score of the two paths. */
    bubble.score = filterLowComplexity2Paths (path1, path2);

    return (bubble.score < threshold || (bubble.score>=threshold && low));
}
