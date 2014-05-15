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
 * A QUICK OVERVIEW...
 *
 * This class implements the detection of SNP in a provided de Bruijn graph.
 *
 * It is implemented as a subclass of Tool, so we get all the facilities of the
 * Tool class. We can therefore see two main parts here:
 *
 * 1) constructor: we define all the command line parameters available
 *
 * 2) 'execute' method: 'main' method of the class where the main loop (ie. iteration
 *    over nodes of the graph) is done.
 *
 * A 'configure' method is called at the beginning of 'execute' in order to configure
 * all the attributes (with command line parameters information)
 *
********************************************************************************/

/** We define string constants for command line options. */
static const char* STR_DISCOSNP_LOW_COMPLEXITY       = "-l";
static const char* STR_DISCOSNP_AUTHORISED_BRANCHING = "-b";
static const char* STR_DISCOSNP_EXTENSION_MODE       = "-t";
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
      low(false), authorised_branching(0), extensionMode (NONE), print_extensions(true),
      sizeKmer (0), min_size_extension(-1), threshold(0)
{
    /** We add options known by kissnp2. */
    getParser()->push_front (new OptionNoParam  (STR_DISCOSNP_LOW_COMPLEXITY,       "conserve low complexity SNPs",     false));
    getParser()->push_front (new OptionOneParam (STR_DISCOSNP_AUTHORISED_BRANCHING, "branching mode\n"
            "\t0: forbid SNPs for wich any of the two paths is branching (high precision, low recall)\n"
            "\t1: forbid SNPs for wich the two paths are branching (e.g. the two paths can be created either with a 'A' or a 'C' at the same position (default value)\n"
            "\t2: No limitation on branching (low precision, high recall)",  false, "1"));
    getParser()->push_front (new OptionOneParam (STR_DISCOSNP_EXTENSION_SIZE,       "extend found SNPs  and conserve only those whose min(left and right extension) is bigger or equal to length",  false, "-1"));
    getParser()->push_front (new OptionOneParam (STR_DISCOSNP_EXTENSION_MODE,       "extension mode: 0=none, 1=unitig, 2=contig",  false, "0"));
    getParser()->push_front (new OptionOneParam (STR_URI_OUTPUT,                    "output name",                      true));
    getParser()->push_front (new OptionOneParam (STR_URI_INPUT,                     "input file (likely a hdf5 file)",  true));
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

    /** We set the extension mode. */
    int mode = getInput()->getInt (STR_DISCOSNP_EXTENSION_MODE);
    if ( !(mode>=0 && mode<=2) )  { throw Exception ("Bad parameter for extension mode (%d)", mode); }
    extensionMode = (ExtensionMode) mode;

    /** We set the name of the output file. */
    stringstream ss;
    ss << getInput()->getStr(STR_URI_OUTPUT)  << "_k_" << sizeKmer  << "_c_" << getGraph().getInfo().getInt("abundance");
    if (min_size_extension > -1)  { ss << "_e_" << min_size_extension; }
    ss << ".fa";

    /** We set the output file. So far, we force FASTA output usage, but we could make it configurable. */
    setOutputBank (new BankFasta (ss.str()));

    /** We need a synchronizer for dumping high/low sequences into the output file in an atomic way
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
 * NOTE: it would much simpler to use a lambda expression but we must be sure that this
 * code will compile with old compilers, so explicit functor struct is mandatory :(
 */
struct FunctorExecute
{
    Kissnp2& finder;

    FunctorExecute (Kissnp2& finder) : finder(finder)  {}

    void operator() (const Node& node)
    {
        /** We start the SNP in both directions (forward and reverse). */
        finder.start (node);
        finder.start (finder.getGraph().reverse(node));
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

    /** We aggregate information. */
    getInfo()->add (1, "config",   "");
    getInfo()->add (2, "kmer_size",   "%d", sizeKmer);
    getInfo()->add (2, "threshold",   "%d", threshold);
    getInfo()->add (2, "auth_branch", "%d", authorised_branching);
    getInfo()->add (2, "low",         "%d", low);
    getInfo()->add (2, "extend_mode", "%d", extensionMode);

    getInfo()->add (1, "bubbles",   "");
    getInfo()->add (2, "nb",      "%lu", nb_bubbles);
    getInfo()->add (2, "nb_high", "%lu", nb_bubbles_high);
    getInfo()->add (2, "nb_low",  "%lu", nb_bubbles_low);

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
void Kissnp2::start (const Node& starter)
{
    /** We start a new bubble from the given starter node. */
    Bubble bubble (graph, starter);

    /** We get the last nucleotide of the starting node. */
    int nt = graph.getNT (starter, sizeKmer-1);

    /** We try all the possible extensions that were not previously tested (clever :-)) */
    for (int i=nt+1; i<4; i++)
    {
        /** We check whether the bubble can be mutated. */
        if (bubble.mutate (i) == true)
        {
            /** We open a new putative bubble. */
            expand (1, bubble, bubble.start1, bubble.start2, Node(~0), Node(~0));
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
void Kissnp2::expand (
    int pos,
    Bubble& bubble,
    const Node& node1,
    const Node& node2,
    const Node& previousNode1,
    const Node& previousNode2
)
{
    Sequence seq1, seq2;

    /** A little check won't hurt. */
    if (pos > sizeKmer-1)  { throw Exception("Bad recursion in Kissnp2..."); }

    /** We may have to stop the extension according to the branching mode. */
    if (checkBranching(node1,node2) == false)  { return; }

    /** We get the common successors of node1 and node2. */
    Graph::Vector < pair<Edge,Edge> > successors = getGraph().successors<Edge> (node1, node2);

    /** We loop over the successors of the two nodes. */
    for (size_t i=0; i<successors.size(); i++)
    {
        /** Shortcuts. */
        pair<Edge,Edge>& p = successors[i];
        Node nextNode1 = p.first.to;
        Node nextNode2 = p.second.to;
        Nucleotide nt  = p.first.nt;

        /** We check whether the new nodes are different from previous ones. */
        bool checkPrevious =
            checkKmersDiff (previousNode1, node1, nextNode1) &&
            checkKmersDiff (previousNode2, node2, nextNode2);

        if (!checkPrevious)  { continue; }

        /** We add the current nucleotide to the bubble paths. */
        bubble.set (pos, nt);

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
            /** We finish the bubble with both last nodes. */
            bubble.finish (nextNode1, nextNode2);

            /** We check the first path vs. its revcomp. */
            if (checkPath(bubble)==true)
            {
                int score=0;

                /** We check the branching properties of the next kmers. */
                if (checkBranching(nextNode1, nextNode2)==false)  { return; }

                /** We check the low complexity of paths (we also get the score). */
                if (checkLowComplexity (bubble, score) == true)
                {
                    char path1_c[2*sizeKmer+2], path2_c[2*sizeKmer+2]; // +2 stands for the \0 character

                    /** We update statistics. NOTE: We have to make sure it is protected against concurrent
                     * accesses since we may be called here from different threads. */
                    size_t currentNbBubbles = __sync_add_and_fetch (&(nb_bubbles), 1 );

                    if (score < threshold)  { __sync_add_and_fetch (&(nb_bubbles_high), 1 );  }
                    else                    { __sync_add_and_fetch (&(nb_bubbles_low),  1 );  }

                    /** We may have to close/extend the bubble. */
                    int where_to_extend = extend (bubble, path1_c, path2_c);

                    /** We retrieve the two sequences. */
                    retrieveSequence (path1_c, "higher", score, where_to_extend, currentNbBubbles, seq1);
                    retrieveSequence (path2_c, "lower",  score, where_to_extend, currentNbBubbles, seq2);

                    /** We have to protect the sequences dump wrt concurrent accesses. We use a {} block with
                     * a LocalSynchronizer instance with the shared ISynchonizer of the Kissnp2 class. */
                    {
                        LocalSynchronizer sync (_synchronizer);

                        /** We insert the two sequences into the output bank. */
                        _outputBank->insert (seq1);
                        _outputBank->insert (seq2);
                    }
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
int Kissnp2::extend (Bubble& bubble, char* path1_c, char* path2_c)
{
    static const int BAD = -1;

    int leftExtension  = BAD;
    int rightExtension = BAD;

    /** We may have to extend the bubble according to the user choice. */
    if (extensionMode != NONE)
    {
        Graph::Vector<Node> predecessors = graph.predecessors<Node> (bubble.start1);
        Graph::Vector<Node> successors   = graph.successors<Node>   (bubble.stop1);

        if (predecessors.size()==1)  { leftExtension  = graph.getNT (predecessors[0], 0);           }
        if (successors.size()  ==1)  { rightExtension = graph.getNT (successors  [0], sizeKmer-1);  }
    }

    /** We add the left extension if any. Note that we use lower case for extensions. */
    if (leftExtension != BAD)   {  *(path1_c++) = *(path2_c++) = tolower(bin2NT[leftExtension]);  }

    /** We add the bubble paths. */
    for (size_t i=0; i<2*sizeKmer-1; i++)
    {
        *(path1_c++) = bin2NT[bubble.path1[i]];
        *(path2_c++) = bin2NT[bubble.path2[i]];
    }

    /** We add the right extension if any. Note that we use lower case for extensions. */
    if (rightExtension != BAD)  {  *(path1_c++) = *(path2_c++) = tolower(bin2NT[rightExtension]);  }

    /** We add a null terminator for the strings. */
    *(path1_c++) = *(path2_c++) = '\0';

    /** We return a code value according to left/right extensions status. */
         if (leftExtension==BAD && rightExtension==BAD)  { return 0; }
    else if (leftExtension!=BAD && rightExtension==BAD)  { return 1; }
    else if (leftExtension==BAD && rightExtension!=BAD)  { return 2; }
    else if (leftExtension!=BAD && rightExtension!=BAD)  { return 3; }

    throw Exception ("Bad extension status in kissnp2");
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
void Kissnp2::retrieveSequence (PATH path, const char* type, int score, int where_to_extend, size_t seqIndex, Sequence& seq)
{
    //Here, path have 2k, 2k-1 or 2k+1 size
    int size_path = strlen(path);

    pair<char*, int> right_contig;
    pair<char*, int> left_contig;

    stringstream commentStream;

    /** We build the comment for the sequence. */
    commentStream << "SNP_" << type << "_path_" << seqIndex << "|" << (score >= threshold ? "low" : "high");

    /** We may have extra information for the comment. */
    addExtraCommentInformation (commentStream, left_contig, right_contig, where_to_extend);

    /** We assign the comment of the sequence. */
    seq.setComment (commentStream.str());

    int start = 0;
    int stop  = strlen(path);

    /** We resize the sequence data. */
    seq.getData().resize (stop-start);

    /** We fill the data buffer of the sequence. */
    for (int i=start; i<stop;i++)  { seq.getDataBuffer()[i] = path[i]; }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
// prints the headers depending of the extension choice.
void Kissnp2::addExtraCommentInformation (
    std::stringstream& ss,
    std::pair<char*,int> left_extension,
    std::pair<char*,int> right_extension,
    int where_to_extend)
{
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
bool Kissnp2::checkKmersDiff (const Node& previous, const Node& current, const Node& next) const
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
    PATH p1 = bubble.path1.data();
    PATH p2 = bubble.path1.data() + 2*sizeKmer - 2;

    static char tableDirect[]  = { 'A', 'C', 'T', 'G' };
    static char tableReverse[] = { 'T', 'G', 'A', 'C' };

    /** Here, we do the comparison between the two kmers 'on the fly': we don't get
     * a string representation of both kmers and then use strcmp for instance.
     * Instead we test each nucleotide (first translated back in ASCII) */
    for (size_t i=0; i<sizeKmer; i++)
    {
        char n1 = tableDirect  [p1[ i]];
        char n2 = tableReverse [p2[-i]];

             if (n1 < n2)  { return true;  }
        else if (n1 > n2)  { return false; }
    }

    return false;
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
bool Kissnp2::checkLowComplexity (Bubble& bubble, int& score) const
{
    /** We compute the low complexity score of the two paths. */
    score = filterLowComplexity2Paths (bubble.path1.data(), bubble.path2.data(), 2*sizeKmer-1, threshold);

    return (score < threshold || (score>=threshold && low));
}
