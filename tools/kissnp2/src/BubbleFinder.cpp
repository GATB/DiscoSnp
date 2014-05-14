/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
 *   Authors: R.Chikhi, G.Rizk, E.Drezen
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

/********************************************************************************/
// We include required definitions
/********************************************************************************/

#include <BubbleFinder.hpp>
#include <Filter.hpp>
using namespace std;

/********************************************************************************/

#define DEBUG(a)   printf a

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
BubbleFinder<span>::BubbleFinder (const Graph& graph, IProperties* input)
    : Algorithm ("bubble", -1, input),
      graph(graph), sizeKmer(graph.getKmerSize()),
      _outputBank(0), _synchronizer(0),
      nb_bubbles(0), nb_bubbles_high(0), nb_bubbles_low(0)
{
    /** We set the threshold. */
    threshold  = (sizeKmer/2-2)*(sizeKmer/2-3);

    /** We set attributes according to user choice. */
    low                  = getInput()->getInt (STR_DISCOSNP_LOW_COMPLEXITY);
    authorised_branching = getInput()->getInt (STR_DISCOSNP_AUTHORISED_BRANCHING);
    min_size_extension   = getInput()->getInt (STR_DISCOSNP_WITH_EXTENSION);

    /** We set the name of the output file. */
    stringstream ss;
    ss << getInput()->getStr(STR_URI_OUTPUT)  << "_k_" << sizeKmer  << "_c_" << graph.getInfo().getInt("abundance");
    if (min_size_extension > -1)  { ss << "_e_" << min_size_extension; }
    ss << ".fa";

    /** We set the output file. */
    setOutputBank (new BankFasta (ss.str()));

    extend_snps        = false;
    print_extensions   = true;
    strict_extension   = true;

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
template<size_t span>
BubbleFinder<span>::~BubbleFinder ()
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
/** We define a functor that only calls BubbleFinder::start in the context of thread call
 * made through a Dispatcher object.
 * NOTE: it would much simpler to use a lambda expression but we must be sure that this
 * code will compile with old compilers, so explicit functor struct is mandatory :(
 */
template<size_t span>
struct FunctorExecute
{
    BubbleFinder<span>& finder;

    FunctorExecute (BubbleFinder<span>& finder) : finder(finder)  {}

    void operator() (const Node& node)
    {
        /** We start the SNP in both directions (forward and reverse). */
        finder.start (node);
        finder.start (finder.getGraph().reverse(node));
    }
};

/** */
template<size_t span>
void BubbleFinder<span>::execute ()
{
    // We get an iterator for all nodes of the graph.
    ProgressGraphIterator<Node,ProgressTimer> it (graph.iterator<Node>(), "nodes");

    /** We launch the iteration over all the nodes of the graph. */
    IDispatcher::Status status = getDispatcher()->iterate (it, FunctorExecute<span>(*this));

    // We aggregate information
    getInfo()->add (0, "config",   "");
    getInfo()->add (1, "kmer_size",   "%d", sizeKmer);
    getInfo()->add (1, "threshold",   "%d", threshold);
    getInfo()->add (1, "auth_branch", "%d", authorised_branching);
    getInfo()->add (1, "low",         "%d", low);

    getInfo()->add (0, "bubbles",   "");
    getInfo()->add (1, "nb",      "%lu", nb_bubbles);
    getInfo()->add (1, "nb_high", "%lu", nb_bubbles_high);
    getInfo()->add (1, "nb_low",  "%lu", nb_bubbles_low);

    getInfo()->add (0, "time", "");
    getInfo()->add (1, "find", "%d", status.time);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
void BubbleFinder<span>::start (const Node& node1)
{
    char path1[2*sizeKmer-1], path2[2*sizeKmer-1];

    /** We copy the kmer in path1 and path2*/
#if 0
    for (size_t i=0; i<sizeKmer; i++)  {  path1[i] = path2[i] = graph.getNT(node1,i);  }
#else
    kmer_type kmer = node1.kmer.get<kmer_type>();
    if (node1.strand == STRAND_FORWARD) {  for (size_t i=0; i<sizeKmer; i++)  {  path1[i] = path2[i] = kmer[sizeKmer-1-i];  } }
    else                                {  for (size_t i=0; i<sizeKmer; i++)  {  path1[i] = path2[i] = reverse((Nucleotide)kmer[i]);    } }
#endif

    /** We try all the possible extensions that were not previously tested (clever :-)) */
    for (int i=path1[sizeKmer-1]+1; i<4; i++)
    {
        /** We mutate the last nucleotide (index 'sizeKmer-1') of the current node. */
        Node node2 = graph.mutate (node1, sizeKmer-1, (Nucleotide)i);

        /** We check whether the mutated kmer belongs to the graph. */
        if (graph.contains(node2) == true)
        {
            /** We update the path2 with the current mutation nucleotide. */
            path2[sizeKmer-1] = i;

            Node previousNode1 (~0);
            Node previousNode2 (~0);

            /** We open a new putative bubble. */
            expand (1, path1, path2, node1, node2, previousNode1, previousNode2);
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
template<size_t span>
void BubbleFinder<span>::expand (
    int pos,
    char* path1,
    char* path2,
    const Node& node1,
    const Node& node2,
    const Node& previousNode1,
    const Node& previousNode2
)
{
    Sequence seq1, seq2;

    /** A little check won't hurt. */
    if (pos > sizeKmer-1)  { throw Exception("Bad recursion in BubbleFinder..."); }

    /** We may have to stop the extension according to the branching mode. */
    if (checkBranching(node1,node2) == false)  { return; }

    /** We get the common successors of node1 and node2. */
    Graph::Vector < pair<Edge,Edge> > successors = graph.successors<Edge> (node1, node2);

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
        path1[sizeKmer-1+pos] = nt;
        path2[sizeKmer-1+pos] = nt;

        /************************************************************/
        /**                   RECURSION CONTINUES                  **/
        /************************************************************/
        if (pos < sizeKmer-1)
        {
            /** We call recursively the method (recursion on 'pos'). */
            expand (pos+1, path1, path2,  nextNode1, nextNode2,  node1, node2);

            /** There's only one branch to expand if we keep non branching SNPs only, therefore we can safely stop the for loop */
            if ( authorised_branching==0 || authorised_branching==1 )   {  break;  }
        }

        /************************************************************/
        /**                   RECURSION FINISHED                   **/
        /************************************************************/
        else
        {
            /** We check the first path vs. its revcomp. */
            if (checkPath(path1)==true)
            {
                int score=0;

                /** We check the branching properties of the next kmers. */
                if (checkBranching(nextNode1, nextNode2)==false)  { return; }

                /** We check the low complexity of paths (we also get the score). */
                if (checkLowComplexity (path1, path2, score) == true)
                {
                    char path1_c[2*sizeKmer+2], path2_c[2*sizeKmer+2]; // +2 stands for the \0 character

                    /** We update statistics. NOTE: We have to make sure it is protected against concurrent
                     * accesses since we may be called here from different threads. */
                    size_t currentNbBubbles = __sync_add_and_fetch (&(nb_bubbles), 1 );

                    if (score < threshold)  { __sync_add_and_fetch (&(nb_bubbles_high), 1 );  }
                    else                    { __sync_add_and_fetch (&(nb_bubbles_low),  1 );  }

                    /** We may have to close the SNP. */
                    int where_to_extend = extend (path1, path2, path1_c, path2_c);

                    /** We retrieve the two sequences. */
                    retrieveSequence (path1_c, "higher", score, where_to_extend, currentNbBubbles, seq1);
                    retrieveSequence (path2_c, "lower",  score, where_to_extend, currentNbBubbles, seq2);

                    /** We have to protect the sequences dump wrt concurrent accesses. We use a {} block with
                     * a LocalSynchronizer instance with the shared ISynchonizer of the BubbleFinder class. */
                    {
                        LocalSynchronizer sync (_synchronizer);

                        /** We insert the two sequences into the output bank. */
                        _outputBank->insert (seq1);
                        _outputBank->insert (seq2);
                    }
                }
            }
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
template<size_t span>
int BubbleFinder<span>::extend (const char* path1, const char* path2, char* path1_c, char* path2_c)
{
    /** This implementation doesn't close the snp => only reports paths. */
    size_t i=0;
    for(i=0; i<2*sizeKmer-1; i++)
    {
        path1_c[i]=bin2NT[path1[i]];
        path2_c[i]=bin2NT[path2[i]];
    }
    path1_c[i]='\0';
    path2_c[i]='\0';

    return 0;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
bool BubbleFinder<span>::two_possible_extensions_on_one_path (const Node& node) const
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
template<size_t span>
bool BubbleFinder<span>::two_possible_extensions (Node node1, Node node2) const
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
template<size_t span>
void BubbleFinder<span>::retrieveSequence (char* path, const char* type, int score, int where_to_extend, size_t seqIndex, Sequence& seq)
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

    // left: first nucleotide is an extension
    if (where_to_extend%2==1) {  start++;  }

    // right: last nucleotide is an extension
    if (where_to_extend>1)    {  stop--;   }

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
template<size_t span>
void BubbleFinder<span>::addExtraCommentInformation (
    std::stringstream& ss,
    std::pair<char*,int> left_extension,
    std::pair<char*,int> right_extension,
    int where_to_extend)
{
#if 0

    if (extend_snps && strict_extension)
    {
        if(where_to_extend == 0) { ss << "|left_unitig_length_0|right_unitig_length_0"; }

        if(where_to_extend == 1) { ss <<  "|left_unitig_length_%d|right_unitig_length_0", (int)strlen(left_extension.first)-sizeKmer+1); //+1 because of close_snp

        if(where_to_extend == 2)
            fprintf(file,"|left_unitig_length_0|right_unitig_length_%d", (int)strlen(right_extension.first)-sizeKmer+1);//+1 because of close_snp

        if(where_to_extend == 3)
            fprintf(file,"|left_unitig_length_%d|right_unitig_length_%d", (int)strlen(left_extension.first)-sizeKmer+1,(int)strlen(right_extension.first)-sizeKmer+1);//+1 because of close_snp
    }

    if (extend_snps && !strict_extension)
    {
        if(where_to_extend == 0)
        {
            fprintf(file,"|left_unitig_length_0|right_unitig_length_0");
            fprintf(file,"|left_contig_length_0|right_contig_length_0");
        }
        if(where_to_extend == 1)
        {
            fprintf(file,"|left_unitig_length_%d|right_unitig_length_0", left_extension.second-sizeKmer+1); //+1 because of close_snp
            fprintf(file,"|left_contig_length_%d|right_contig_length_0", (int)strlen(left_extension.first)-sizeKmer+1); //+1 because of close_snp
        }
        if(where_to_extend == 2)
        {
            fprintf(file,"|left_unitig_length_0|right_unitig_length_%d", right_extension.second-sizeKmer+1);//+1 because of close_snp
            fprintf(file,"|left_contig_length_0|right_contig_length_%d", (int)strlen(right_extension.first)-sizeKmer+1);//+1 because of close_snp
        }
        if(where_to_extend == 3)
        {
            fprintf(file,"|left_unitig_length_%d|right_unitig_length_%d", left_extension.second-sizeKmer+1,right_extension.second-sizeKmer+1);//+1 because of close_snp
            fprintf(file,"|left_contig_length_%d|right_contig_length_%d", (int)strlen(left_extension.first)-sizeKmer+1,(int)strlen(right_extension.first)-sizeKmer+1);//+1 because of close_snp
        }
    }
#endif
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
bool BubbleFinder<span>::checkKmersDiff (const Node& previous, const Node& current, const Node& next) const
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
template<size_t span>
bool BubbleFinder<span>::checkPath (char* path) const
{
    /** We test whether the first kmer of the first path is smaller than
     * the first kmer of the revcomp(first path), this should avoid repeated SNPs */
    char* p1 = path;
    char* p2 = path + 2*sizeKmer - 2;

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
template<size_t span>
bool BubbleFinder<span>::checkBranching (const Node& node1, const Node& node2) const
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
template<size_t span>
bool BubbleFinder<span>::checkLowComplexity (char* path1, char* path2, int& score) const
{
    /** We compute the low complexity score of the two paths. */
    score = filterLowComplexity2Paths (path1, path2, 2*sizeKmer-1, threshold);

    return (score < threshold || (score>=threshold && low));
}

/*********************************************************************
#######  #     #  #######  #######  #     #   #####   ###  #######  #     #
#         #   #      #     #        ##    #  #     #   #   #     #  ##    #
#          # #       #     #        # #   #  #         #   #     #  # #   #
#####       #        #     #####    #  #  #   #####    #   #     #  #  #  #
#          # #       #     #        #   # #        #   #   #     #  #   # #
#         #   #      #     #        #    ##  #     #   #   #     #  #    ##
#######  #     #     #     #######  #     #   #####   ###  #######  #     #
*********************************************************************/

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
BubbleFinderWithExtension<span>::BubbleFinderWithExtension (const Graph& graph, IProperties* input)
    : BubbleFinder<span> (graph,input)
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
template<size_t span>
BubbleFinderWithExtension<span>::~BubbleFinderWithExtension ()
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
template<size_t span>
int BubbleFinderWithExtension<span>::extend (const char* path1, const char* path2, char* path1_c, char* path2_c)
{
}

/********************************************************************************/
template class BubbleFinder <KSIZE_1>;
template class BubbleFinder <KSIZE_2>;
template class BubbleFinder <KSIZE_3>;
template class BubbleFinder <KSIZE_4>;
