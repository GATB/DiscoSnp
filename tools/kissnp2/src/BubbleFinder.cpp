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
      modelMin(graph.getKmerSize(), KMER_MINIMUM), modelDirect(graph.getKmerSize(), KMER_DIRECT),
      _outputBank(0), _synchronizer(0),
      nb_bubbles(0), nb_bubbles_high(0), nb_bubbles_low(0)
{
    threshold  = (sizeKmer/2-2)*(sizeKmer/2-3);
    setOutputBank (new BankFasta ((getInput()->getStr(STR_URI_OUTPUT)+string(".fa")).c_str()));

    low                  = getInput()->getInt (STR_LOW_COMPLEXITY);
    authorised_branching = getInput()->getInt (STR_AUTHORISED_BRANCHING);

    extend_snps        = false;
    print_extensions   = true;
    min_size_extension = -1;
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
 * code will compile with old compilers :(
 */
template<size_t span>
struct FunctorExecute
{
    BubbleFinder<span>& finder;
    FunctorExecute (BubbleFinder<span>& finder) : finder(finder) {}
    void operator() (const Node& node)  {  finder.start (node);  }
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
void BubbleFinder<span>::start (const Node& node)
{
    // Shortcuts
    kmer_type kmer1 = node.kmer.get<kmer_type>();

    char path1[2*sizeKmer-1], path2[2*sizeKmer-1];

    for (int direction=0; direction<=1; direction++)
    {
        /** We start bubble for kmer in reverse form. */
        if (direction == 1)  {  kmer1 = modelMin.reverse(kmer1);  }

        /** We copy the kmer in path1 and path2*/
        for (size_t i=0; i<sizeKmer; i++)  {  path1[i] = path2[i] = kmer1[sizeKmer-1-i];  }

        /** We try all the possible extensions that were not previously tested (clever :-)) */
        for (int i=kmer1[0]+1; i<4; i++)
        {
            /** We mutate the last nucleotide of the current node. */
            kmer_type kmer2 = modelDirect.mutate (kmer1, 0, i);

            // BETTER NOT USE IT : avoid some double SNPs: do only A-C, A-G couples instead of G-T and C-T.
            // The problem comes from A-T and C-G that are still double.

            /** We check whether the mutated kmer belongs to the graph. */
            Node::Value val (min (kmer2, modelMin.reverse(kmer2)));
            if (graph.contains(val) == true) // the tried kmer is indexed.
            {
                /** We update the path2 with the current mutation nucleotide. */
                path2[sizeKmer-1] = i;

                kmer_type previous_kmer1 (~0);
                kmer_type previous_kmer2 (~0);

                /** We open a new putative bubble. */
                expand (1, path1, path2, kmer1, kmer2, previous_kmer1, previous_kmer2);
            }
        }

    }  // for (direction...)
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
    char* path1, char* path2,
    const kmer_type& kmer1, const kmer_type& kmer2,
    const kmer_type& previous_kmer1,
    const kmer_type& previous_kmer2
)
{
    kmer_type next_kmer1, next_kmer2;
    Sequence seq1, seq2;

    if (pos > sizeKmer-1)  { return; }

    /** We may have to stop the extension according to the branching mode. */
    if (checkBranching(kmer1,kmer2) == false)  { return; }

    for (int nt=0; nt<4; nt++)
    {
        next_kmer1 = modelMin.codeSeedRight (kmer1, nt, Data::INTEGER);
        next_kmer2 = modelMin.codeSeedRight (kmer2, nt, Data::INTEGER);

        Node::Value val1 (next_kmer1);
        Node::Value val2 (next_kmer2);

        if (graph.contains(val1) && graph.contains(val2))
        {
            /** We check whether the new kmers are different from previous ones. */
            bool checkPrevious =
                checkKmersDiff (previous_kmer1, kmer1, next_kmer1) &&
                checkKmersDiff (previous_kmer2, kmer2, next_kmer2);

            if (!checkPrevious)  { continue; }

            /** We add the current nucleotide to the bubble paths. */
            path1[sizeKmer-1+pos] = nt;
            path2[sizeKmer-1+pos] = nt;

            /************************************************************/
            /**                   RECURSION CONTINUES                  **/
            /************************************************************/
            if (pos < sizeKmer-1)
            {
                next_kmer1 = modelDirect.codeSeedRight (kmer1, nt, Data::INTEGER);
                next_kmer2 = modelDirect.codeSeedRight (kmer2, nt, Data::INTEGER);

                /** We call recursively the method (recursion on 'pos'). */
                expand (pos+1, path1, path2,  next_kmer1, next_kmer2,  kmer1, kmer2);

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
                    next_kmer1 = modelDirect.codeSeedRight (kmer1, nt, Data::INTEGER);
                    next_kmer2 = modelDirect.codeSeedRight (kmer2, nt, Data::INTEGER);
                    if (checkBranching(next_kmer1, next_kmer2)==false)  { return; }

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
                        int where_to_extend = close_snp (path1, path2, path1_c, path2_c);

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

        } /* end of if (graph.contains(val1) && graph.contains(val2)... */

    } /* end of for(nt=0; nt<4; nt++) */
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
int BubbleFinder<span>::close_snp (const char* path1, const char* path2, char* path1_c, char* path2_c)
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
bool BubbleFinder<span>::two_possible_extensions_on_one_path (kmer_type kmer) const
{
    Node::Value val(kmer);
    Node node (val);
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
bool BubbleFinder<span>::two_possible_extensions (kmer_type kmer1, kmer_type kmer2) const
{
    kmer_type next_kmer1, next_kmer2;
    bool already_an_extension = false;

    for (int d=0; d<=1; d++)
    {
        if (d == 1)
        {
            kmer1 = modelMin.reverse (kmer1);
            kmer2 = modelMin.reverse (kmer2);
        }

        for(int nt=0; nt<4; nt++)
        {
            next_kmer1 = modelMin.codeSeedRight (kmer1, nt, Data::INTEGER);
            next_kmer2 = modelMin.codeSeedRight (kmer2, nt, Data::INTEGER);

            Node::Value val1 (next_kmer1);
            Node::Value val2 (next_kmer2);

            if (graph.contains (val1) && graph.contains(val2))
            {
                if (!already_an_extension)  {  already_an_extension = true;  }
                else                        {  return true;                  }
            }
        }
        already_an_extension = false;
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
bool BubbleFinder<span>::checkKmersDiff (const kmer_type& previous, const kmer_type& current, const kmer_type& next) const
{
    /** We have to get the min/revcomp of current and previous because they may not be the minimal,
     * which is mandatory to compare to the 'next' kmer. */
    kmer_type previousMin = min (previous, modelMin.reverse (previous));
    kmer_type currentMin  = min (current,  modelMin.reverse (current));

    return (next != currentMin) && (next != previousMin);
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

    kmer_type kmerFirst = modelDirect.codeSeed (path + 0,          Data::INTEGER);
    kmer_type kmerLast  = modelDirect.codeSeed (path + sizeKmer-1, Data::INTEGER);

    kmerLast = modelMin.reverse (kmerLast);

    return modelMin.toString(kmerFirst).compare (modelMin.toString (kmerLast)) < 0;
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
bool BubbleFinder<span>::checkBranching (kmer_type kmer1, kmer_type kmer2) const
{
    // stop the extension if authorised_branching==0 (not branching in any path) and any of the two paths is branching
    if (authorised_branching==0 && (two_possible_extensions_on_one_path(kmer1) || two_possible_extensions_on_one_path(kmer2)))
    {
        return false;
    }

    // stop the extension if authorised_branching==1 (not branching in both path) and both the two paths are branching
    if (authorised_branching==1 && two_possible_extensions (kmer1, kmer2))
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

/********************************************************************************/
template class BubbleFinder <KSIZE_1>;
template class BubbleFinder <KSIZE_2>;
template class BubbleFinder <KSIZE_3>;
template class BubbleFinder <KSIZE_4>;
