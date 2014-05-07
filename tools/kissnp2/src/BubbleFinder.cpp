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
BubbleFinder<span>::BubbleFinder (
    const Graph& graph,
    const char* SNP_file_name, int low, int authorised_branching,
    bool extend_snps, int min_size_extension,bool print_extensions, bool strict_extension
)
    : Algorithm ("bubble"),

      graph(graph), modelMin(graph.getKmerSize(), KMER_MINIMUM), modelDirect(graph.getKmerSize(), KMER_DIRECT),

      low(low), authorised_branching(authorised_branching),

      nb_bubbles(0), nb_bubbles_high(0), nb_bubbles_low(0),

      sizeKmer(graph.getKmerSize()),
nbKmer2(0), SKIP1(0), SKIP2(0), SKIP3(0), checksumPrint1(0),
      extend_snps(extend_snps), min_size_extension(min_size_extension),
      print_extensions (print_extensions), strict_extension (strict_extension)
{
    threshold  = (sizeKmer/2-2)*(sizeKmer/2-3);

    SNP_file = fopen (SNP_file_name,"w");
    if (SNP_file == NULL)  {  throw Exception ("Cannot open file %s, exit\n",SNP_file_name);  }
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
    if (SNP_file)  { fclose(SNP_file); }
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
void BubbleFinder<span>::execute ()
{
    // We get an iterator for all nodes of the graph.
    ProgressGraphIterator<Node,ProgressTimer> it (graph.iterator<Node>(), "nodes");

    // NOTE: we use a {} block to get the execution time of the loop
    {
        TIME_INFO (getTimeInfo(), "find");

        // We loop each node of the graph
        for (it.first(); !it.isDone(); it.next())  {  start (it.item());  }
    }

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

    getInfo()->add (0, getTimeInfo().getProperties("time"));

    stringstream ss;
    getInfo()->add (0, "checksums",   "");
    ss.str(""); ss << checksumStart1;   getInfo()->add (1, "checksumStart1",  "%s", ss.str().c_str());
    ss.str(""); ss << checksumStart2;   getInfo()->add (1, "checksumStart2",  "%s", ss.str().c_str());
    ss.str(""); ss << checksumStart3;   getInfo()->add (1, "checksumStart3",  "%s", ss.str().c_str());
    ss.str(""); ss << checksumExpand1;  getInfo()->add (1, "checksumExpand1", "%s", ss.str().c_str());
    ss.str(""); ss << checksumExpand2;  getInfo()->add (1, "checksumExpand2", "%s", ss.str().c_str());
    ss.str(""); ss << checksumGraph;    getInfo()->add (1, "checksumGraph",   "%s", ss.str().c_str());
    ss.str(""); ss << checksumPrint1;   getInfo()->add (1, "checksumPrint1",  "%s", ss.str().c_str());
    getInfo()->add (1, "nbKmer2", "%d", nbKmer2);
    getInfo()->add (1, "SKIP1",   "%d", SKIP1);
    getInfo()->add (1, "SKIP2",   "%d", SKIP2);
    getInfo()->add (1, "SKIP3",   "%d", SKIP3);
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
    kmer_type kmer1   = node.kmer.get<kmer_type>();

    int direction;
    char path1[2*sizeKmer-1], path2[2*sizeKmer-1];
    char p1[sizeKmer+1], p2[sizeKmer+1], p3[sizeKmer+1], p4[sizeKmer+1];

    for (direction=0; direction<=1; direction++)
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

checksumStart1 += kmer1;
checksumStart2 += kmer2;
checksumStart3 += min (kmer2, modelMin.reverse(kmer2));


            // BETTER NOT USE IT : avoid some double SNPs: do only A-C, A-G couples instead of G-T and C-T.
            // The problem comes from A-T and C-G that are still double.

            /** We check whether the mutated kmer belongs to the graph. */
            Node::Value val (min (kmer2, modelMin.reverse(kmer2)));
            if (graph.contains(val) == true) // the tried kmer is indexed.
            {
//cout << "------------->  GRAPH kmer2=" << kmer2 << "  " << model.toString(kmer2) << endl;
checksumGraph += kmer2;
nbKmer2++;
                /** We update the path2 with the current mutation. */
                path2[sizeKmer-1] = i;

                strcpy(p1,"\0");
                strcpy(p3,"\0");

                /** We setup strings for further process. */
                modelMin.toString (kmer1, p2);
                modelMin.toString (kmer2, p4);

                /** We open a new putative bubble. */
                expand (0, path1, path2, kmer1, kmer2, 1, p1, p2, p3, p4);
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
    int direction,
    char* path1, char* path2,
    kmer_type kmer1, kmer_type kmer2,
    int pos,
    char* p1, char* p2, char* p3, char* p4
)
{
    kmer_type next_kmer1, next_kmer2;
    int nt,nt2,i,score=0;
    char p_aux1[sizeKmer+1], p_aux2[sizeKmer+1];

    if (direction == 1)
    {
        kmer1 = modelMin.reverse (kmer1);
        kmer2 = modelMin.reverse (kmer2);
    }

    if (pos > sizeKmer-1)  { return; }

    /** We may have to stop the extension according to the branching mode. */
    if (checkBranching(kmer1,kmer2) == false)  { return; }

    for(nt=0; nt<4; nt++)
    {
        next_kmer1 = modelMin.codeSeedRight (kmer1, nt, Data::INTEGER);
        next_kmer2 = modelMin.codeSeedRight (kmer2, nt, Data::INTEGER);

        modelMin.toString (next_kmer1, p_aux1);
        modelMin.toString (next_kmer2, p_aux2);

        Node::Value val1 (next_kmer1);
        Node::Value val2 (next_kmer2);

        if (graph.contains(val1) && graph.contains(val2)
            && strcmp(p_aux1,p1) && strcmp(p_aux1,p2) && strcmp(p_aux2,p3) && strcmp(p_aux2,p4)
        )
        {
            strcpy(p1, p2);
            strcpy(p3, p4);
            strcpy(p2, p_aux1);
            strcpy(p4, p_aux2);

            if ( direction == 0 )
            {
                //if everything is ok then we only use direction = 0
                path1[sizeKmer-1+pos]=nt;
                path2[sizeKmer-1+pos]=nt;
            }
            else
            {
                //shouldn't go in here
                path1[sizeKmer-1+pos]=nt+4;
                path2[sizeKmer-1+pos]=nt+4;
            }

            //if we finished the bubble and need to close it
            if (pos == sizeKmer-1)
            {
                //TEST whether the first kmer of the first path is smaller than the first kmer of the revcomp(first path), this should avoid repeated SNPs
#if 0
                char first_kmer[sizeKmer+1], first_kmer_rev[sizeKmer+1];
                for ( i=0; i<sizeKmer; i++ )
                    first_kmer[i] = bin2NT[path1[i]];
                first_kmer[sizeKmer]='\0';

                for ( i=sizeKmer-1; i<2*sizeKmer-1; i++ )
                    first_kmer_rev[i-sizeKmer+1] = bin2NT[path1[i]];
                first_kmer_rev[sizeKmer]='\0';

                revcomp(first_kmer_rev, sizeKmer);

                if (strcmp(first_kmer, first_kmer_rev)<0)
#endif

                {
                    /** We check the branching properties of the next kmers. */
                    next_kmer1 = modelDirect.codeSeedRight (kmer1, nt, Data::INTEGER);
                    next_kmer2 = modelDirect.codeSeedRight (kmer2, nt, Data::INTEGER);
                    if (checkBranching(next_kmer1, next_kmer2)==false)  { return; }

                    nb_bubbles++;

                    // TODO: avoid this test if low is not required.
                    score = filterLowComplexity2Paths (path1, path2, 2*sizeKmer-1, threshold);
                    if ( score < threshold || (score>=threshold && low))
                    {
                        char path1_c[2*sizeKmer+2], path2_c[2*sizeKmer+2]; // +2 stands for the \0 character
                        if (score < threshold)
                            nb_bubbles_high++;
                        else
                            nb_bubbles_low++;

                        //do not close snps and output only 2k-1 paths if extend_snps is not true
                        int where_to_extend=0;
                        if (!extend_snps)
                        {
                            for(i=0;i<2*sizeKmer-1;i++)
                            {
                                path1_c[i]=bin2NT[path1[i]];
                                path2_c[i]=bin2NT[path2[i]];
                            }
                            path1_c[i]='\0';
                            path2_c[i]='\0';
                        }
#if 0
                        else
                            where_to_extend = close_snp(path1, path2, path1_c, path2_c);
#endif

                        print_sequence_and_eventually_contigs (path1_c, path2_c, score, where_to_extend); //2k-1, 2k or 2k+1 paths
                    }
                }

            }
            else //bubble not yet finished
            {
                next_kmer1 = modelDirect.codeSeedRight (kmer1, nt, Data::INTEGER);
                next_kmer2 = modelDirect.codeSeedRight (kmer2, nt, Data::INTEGER);

checksumExpand1 += next_kmer1;
checksumExpand2 += next_kmer2;

                modelMin.toString (next_kmer1, p_aux1);
                modelMin.toString (next_kmer2, p_aux2);

                /** We call recursively the method. */
                expand (direction, path1, path2, next_kmer1, next_kmer2, pos+1, p1, p2, p3, p4);

                //there's only one branch to expand if we keep non branching SNPs only, therefore we can safely stop the for loop
                if ( authorised_branching==0 || authorised_branching==1 )
                {
                    SKIP3++;
                    break;
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
#if 1
/**
 * Extends if necessary left and/or right parts of the bubble.
 * Prints results
 * path1 and path2: character pathes (no binary format). Of length 2k-1, 2k, or 2k+1, depending ont the results of "close_snps"
 * score=complexity score
 * where_to_extend : 0=nothing, 1=left only, 2=right only, 3=both
 */
template<size_t span>
void BubbleFinder<span>::print_sequence_and_eventually_contigs (
    char* path1,
    char* path2,
    const int score,
    int where_to_extend
)
{
checksumPrint1 += 1;


    int i;
#if 0
    terminator->reset(); // need to reset branching kmers between extensions in case of overlapping extensions.
#endif
    FILE * file = SNP_file;

    //Here, path1 and path2 may have 2k, 2k-1 or 2k+1 size
    int size_path = strlen(path1);

    pair<char*, int> right_contig;
    pair<char*, int> left_contig;

#if 0
    if(extend_snps)
    {
        if(strict_extension)    IterativeExtensions::traversal_type = IterativeExtensions::SimplePaths; // strict
        else                    IterativeExtensions::traversal_type = IterativeExtensions::Monument;    // contigs

        if (where_to_extend == 2 || where_to_extend == 3)
        {
            right_contig=IterativeExtensions::extend_a_snp(path1,10000); // path1 or path2 provides the same result (the extended kmer is not involved)
        }

        if (where_to_extend == 1 || where_to_extend == 3)
        {
            revcomp(path1, size_path); // get the reverse complement of the sequence
            left_contig=IterativeExtensions::extend_a_snp(path1,10000); // get the left extension of this sequence
            revcomp(left_contig.first, strlen(left_contig.first)); // put back the left contig in the right order
            revcomp(path1, size_path); // put back the sequence in the right order.
        }

        if(min_size_extension >-1)
        {
            if(where_to_extend<3 || (int)strlen(left_contig.first)-sizeKmer<min_size_extension || (int)strlen(right_contig.first)-sizeKmer<min_size_extension) {
                //            printf("Don't output this, too short\n");
                free(left_contig.first);
                free(right_contig.first);
                // FIXME: do I have to free the pairs themselves ?
                return;
            }
        }
    }
#endif


    fprintf(file, ">SNP_higher_path_%lu|", nb_bubbles);
    if ( score >= threshold )   { fprintf(file,"low");  }
    else                        { fprintf(file,"high"); }

    prints_contig_informations (file, left_contig, right_contig, extend_snps, strict_extension, where_to_extend);
    fprintf(file,"\n");

#if 0
    if (print_extensions && extend_snps)
    {
        // prints the left extension
        if(where_to_extend%2==1) for(i=0;i<strlen(left_contig.first)-sizeKmer;i++) fprintf(file,"%c", tolower(left_contig.first[i]));

        // change the case of first and/or last character of the central node
        if (where_to_extend%2==1) // left: first nucleotide is an extension
            path1[0]=tolower(path1[0]);

        if (where_to_extend>1) // right: last nucleotide is an extension
            path1[size_path-1]=tolower(path1[size_path-1]);

        fprintf(file,"%s",path1);
        if(extend_snps && print_extensions && where_to_extend>1) for(i=sizeKmer;i<strlen(right_contig.first);i++) fprintf(file,"%c", tolower(right_contig.first[i]));
    }
    else
#endif
    {
        int start = 0;
        int stop = strlen(path1);
        if(!print_extensions)
        {
            if(where_to_extend%2==1) // left: first nucleotide is an extension
                start++;
            if (where_to_extend>1) // right: last nucleotide is an extension
                stop--;
        }
        for(i=start;i<stop;i++) fprintf(file,"%c", path1[i]);
    }
    fprintf(file,"\n");


    fprintf(file, ">SNP_lower_path_%lu|", nb_bubbles);
    if ( score >= threshold )
        fprintf(file,"low");
    else
        fprintf(file,"high");

    prints_contig_informations(file, left_contig, right_contig, extend_snps, strict_extension, where_to_extend);

    fprintf(file,"\n");

#if 0
    if (print_extensions && extend_snps)
    {
        if(print_extensions && where_to_extend%2==1) for(i=0;i<strlen(left_contig.first)-sizeKmer;i++) fprintf(file,"%c", tolower(left_contig.first[i]));
        // change the case of first and/or last character of the central node
        if (where_to_extend%2==1) // left: first nucleotide is an extension
            path2[0]=tolower(path2[0]);
        if (where_to_extend>1) // right: last nucleotide is an extension
            path2[size_path-1]=tolower(path2[size_path-1]);
        fprintf(file,"%s",path2);
        if(extend_snps && print_extensions && where_to_extend>1) for(i=sizeKmer;i<strlen(right_contig.first);i++) fprintf(file,"%c", tolower(right_contig.first[i]));
    }
    else
#endif
    {
        int start = 0;
        int stop = strlen(path2);
        if(!print_extensions)
        {
            if(where_to_extend%2==1) // left: first nucleotide is an extension
                start++;
            if (where_to_extend>1) // right: last nucleotide is an extension
                stop--;
        }
        for(i=start;i<stop;i++) fprintf(file,"%c", path2[i]);
    }
    fprintf(file,"\n");

#if 0
    if(extend_snps && where_to_extend%2==1)
    {
        free(left_contig.first);
    }
    if(extend_snps && where_to_extend>1)
    {
        free(right_contig.first);
    }
#endif
}
#endif

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
void BubbleFinder<span>::prints_contig_informations (
    FILE* file,
    std::pair<char*,int> left_extension,
    std::pair<char*,int> right_extension,
    bool extend_snps,
    bool strict_extension,
    int where_to_extend)
{
    if (extend_snps && strict_extension)
    {
        if(where_to_extend == 0)
            fprintf(file,"|left_unitig_length_0|right_unitig_length_0");
        if(where_to_extend == 1)
            fprintf(file,"|left_unitig_length_%d|right_unitig_length_0", (int)strlen(left_extension.first)-sizeKmer+1); //+1 because of close_snp
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

/********************************************************************************/
template class BubbleFinder <KSIZE_1>;
template class BubbleFinder <KSIZE_2>;
template class BubbleFinder <KSIZE_3>;
template class BubbleFinder <KSIZE_4>;
