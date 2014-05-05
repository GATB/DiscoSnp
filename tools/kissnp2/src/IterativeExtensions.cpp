//Copyright inria / irisa (2013)
//
//
//raluca.uricaru@gmail.com
//pierre.peterlongo@inria.fr
//
//This software is a computer program whose purpose is to call SNPs from NGS reads.
//
//This software is governed by the CeCILL license under French law and
//abiding by the rules of distribution of free software.  You can  use,
//modify and/ or redistribute the software under the terms of the CeCILL
//license as circulated by CEA, CNRS and INRIA at the following URL
//"http://www.cecill.info".
//
//As a counterpart to the access to the source code and  rights to copy,
//modify and redistribute granted by the license, users are provided only
//with a limited warranty  and the software's author,  the holder of the
//economic rights,  and the successive licensors  have only  limited
//liability.
//
//In this respect, the user's attention is drawn to the risks associated
//with loading,  using,  modifying and/or developing or reproducing the
//software by the user in light of its specific status of free software,
//that may mean  that it is complicated to manipulate,  and  that  also
//therefore means  that it is reserved for developers  and  experienced
//professionals having in-depth computer knowledge. Users are therefore
//encouraged to load and test the software's suitability as regards their
//requirements in conditions enabling the security of their systems and/or
//data to be ensured and,  more generally, to use and operate it in the
//same conditions as regards security.
//
//The fact that you are presently reading this means that you have had
//knowledge of the CeCILL license and that you accept its terms.

#include "IterativeExtensions.h"

//#define DONTMARK

AssocPairedSet *pairedBranchingKmers;

// this type bears similarity with Traversal.h:kmer_strand_nt but I wanted to rename "nt" to "depth" for clarity
struct kmer_strand_depth {
    kmer_type kmer;
    int strand;
    int depth;
    kmer_strand_depth(kmer_type kmer, int strand, int depth) : kmer(kmer), strand(strand), depth(depth) {}
    bool operator<(const kmer_strand_depth &other) const { // needed for comparisons inside a list
        if (kmer != other.kmer)
            return (kmer < other.kmer);
        if (depth != other.depth)
            return depth < other.depth;
        return (strand < other.strand);
    }
};

/*
 * our assembly graph is connected by (k-1)-overlaps,
 * so this function is used to make sure we see each (k-1)-overlap in at most one right extremity 
 */
bool compare_and_mark_last_k_minus_one_mer(string node, set<kmer_type> &kmers_set)
{
    kmer_type kmer_fw, kmer_rc;

    sizeKmer--;
    kmer_type kmer = extractKmerFromRead( (char *)node.c_str(), node.size() - sizeKmer, &kmer_fw, &kmer_rc, false);
    int strand = (kmer == kmer_rc);
//    char kmer_seq[100];
//    code2seq(kmer_fw,kmer_seq); // convert starting kmer to nucleotide seq
//    //printf("checking marking kmer: %s\n",kmer_seq);
    sizeKmer++;

    if (kmers_set.find(kmer) != kmers_set.end())
         return true;

    kmers_set.insert(kmer); 
    return false;
}

    
// default extension modes (do not change)
IterativeExtensions::Traversal_type IterativeExtensions::traversal_type = Monument;
IterativeExtensions::When_to_stop_extending IterativeExtensions::when_to_stop_extending = Until_max_depth;
bool IterativeExtensions::dont_output_first_nucleotide = false; 


void IterativeExtensions::construct_linear_seqs(string L, int max_depth, string output_file)
{
    construct_linear_seqs(L, string(),max_depth,output_file,2,1000000000);
}

/*
 * return the contig which starts with L, as well as all ths contigs that follow him, up to max_depth. 
 * results go to a fasta file (output_file)
 *
 * requires that:
 *  global variables terminator, bloo1 and false_positives are already constructed by minia 
 *
 * outputs:
 *  a set of contigs in the output file */
void IterativeExtensions::construct_linear_seqs(string L,string R, int max_depth, string output_file , int verb, int max_nodes,  parcours_t search_mode, bool swf )

{
    bool debug = verb>=2 ? 1 : 0;
    debug = false; // force undebug mode
   // debug = true;
        kmer_type kmer;
    char kmer_seq[sizeKmer+1];
    FILE *linear_seqs_file;
    Traversal *traversal;
    if (traversal_type == IterativeExtensions::SimplePaths)
        traversal = new SimplePathsTraversal(bloo1,false_positives,terminator);
    else if (traversal_type == IterativeExtensions::Monument)
        traversal = new MonumentTraversal(bloo1,false_positives,terminator);
    traversal->set_maxlen(max_depth);
    traversal->set_max_depth(500);
    traversal->set_max_breadth(20);

    long long nbNodes = 0;
    long long totalnt=0;
    long long contig_len =0;
    long long maxlen=1000000;
    char *right_traversal = (char *) malloc(maxlen*sizeof(char));
    char *node          = (char *) malloc(maxlen*sizeof(char));
    node [0] = '\0';
    
    if(!output_file.empty())
    {
        linear_seqs_file = fopen((char * )output_file.c_str(),"w");
    }
    else {
        linear_seqs_file = NULL;
        IterativeExtensions::when_to_stop_extending = IterativeExtensions::After_first_contig; // sequence
    }
    
    STARTWALL(nodes);

    // heuristics: start from the last kmer of L
    
    vector < kmer_strand_depth > kmers_to_traverse;
    kmer_type kmer_fw, kmer_rc;

    kmer = extractKmerFromRead( (char *)L.c_str(), L.size() - sizeKmer, &kmer_fw, &kmer_rc, false);
    int strand = (kmer == kmer_rc);
    kmer_strand_depth first_kmer (kmer,strand,0);
    kmers_to_traverse.push_back(first_kmer);

#ifndef DONTMARK
    set<kmer_type> already_extended_from;
 //   compare_and_mark_last_k_minus_one_mer(L, already_extended_from); // mark first kmer to never extend from it again, // L will be marked at first iteration below
#endif

    while (kmers_to_traverse.size() > 0) // min_depth is max_gap_length here
    {
        
        kmer_strand_depth ksd (0,0,0);
        if (search_mode == PROFONDEUR)
        {
            
             ksd = kmers_to_traverse.back();
             kmers_to_traverse.pop_back();
        }
        else if (search_mode == LARGEUR)
        {
            ksd = kmers_to_traverse.front();
            kmers_to_traverse.erase (kmers_to_traverse.begin());

        }


        
        kmer = ksd.kmer;
        int strand = ksd.strand;
        int depth = ksd.depth;

        code2seq(kmer,kmer_seq); // convert starting kmer to nucleotide seq
        if (strand == 1)
            revcomp_sequence(kmer_seq,sizeKmer);

        if (debug)
            printf(" --- iteration: kmer %s%s (of lenght %d) depth %d nbNodes explored %lli ---\n",kmer_seq,strand?" (internally rc)":"", (int)strlen(kmer_seq), depth,nbNodes);

        // right extension
        int len_right = traversal->traverse(kmer,right_traversal,strand);

        if (debug)
            printf("right traversal %i = %s\n", len_right,right_traversal);
        // save the node
        strcpy(node,kmer_seq);//               + starting_kmer if we ask for the de bruijn graph or for a text output
        strcat(node,right_traversal);//           + right_traversal
	
        int node_len=len_right+sizeKmer;
	//TODO: watch the reason of wrong length with depth = 0
	
	fprintf(linear_seqs_file,">%lli__len__%i__depth__%i\n",nbNodes,node_len,depth);

	  if (depth == 0 && dont_output_first_nucleotide )
        /* we need this in mapsembler: 
          // the first used kmer should be a k-1 mer. Indeed, if a first kmer is extracted from a sequence :
          /// -------------****** (k=6), then this node is the one linked to a new one starting with ******, thus with an overlap of k and not k-1. */
	   fprintf(linear_seqs_file,"%s\n",node+1);
       else
	   fprintf(linear_seqs_file,"%s\n",node);
	
            
        nbNodes++;
        totalnt+=node_len;
        
        if (debug)
        {
          //  printf ("[L=%s] assembled a %d bp node at depth %d\n",L.c_str(),node_len, depth);
            printf("(node seq: %s)\n",node);
        }
	   
        // if we only want 1 extension, stop now
        if (when_to_stop_extending == IterativeExtensions::After_first_contig)
                      break;

        if(swf)
        {
            char * found = NULL;
            found =strstr(node, R.c_str());
            if(found!=NULL && depth > sizeKmer)
            {
               // printf("swf STOP \n");

                break;
            }
        }

        if(totalnt > max_nodes) //GR stop when too complex  huum when to stop ?
            break;
        

        // if max depth reached, don't extend that one
        if (depth + totalnt > max_depth)
            continue;


#ifndef DONTMARK
        // make sure this is the only time we see this (k-1)-overlap
        bool already_seen = compare_and_mark_last_k_minus_one_mer(node, already_extended_from);
        if (already_seen)
            continue;
#endif


        kmer = extractKmerFromRead( node, node_len - sizeKmer, &kmer_fw, &kmer_rc, false);
        strand = (kmer == kmer_rc);

        if (debug)
        {
            code2seq(kmer,kmer_seq);
            if (strand == 1)
                revcomp_sequence(kmer_seq,sizeKmer);
            printf("... restarting from kmer %s %s\n",kmer_seq,strand?"(internally rc)":"");
        }

        // continue extending from immediately overlapping kmers
        // there may be just one 1 possibility (there was in-branching)
        
        for(int test_nt=0; test_nt<4; test_nt++) 
        {
            int current_strand = strand;
            kmer_type current_kmer = next_kmer(kmer,test_nt, &current_strand);

            if(bloo1->contains(current_kmer) && !false_positives->contains(current_kmer)){
                kmers_to_traverse.push_back(kmer_strand_depth (current_kmer,current_strand,depth + len_right +1 )); // ou plutot depth + len_right +1 (+1 = la nt ajoutee ici) (et pas node_len)  ?

                if (debug)
                {
                    code2seq(current_kmer,kmer_seq); // convert starting kmer to nucleotide seq
                    if (current_strand == 1)
                        revcomp_sequence(kmer_seq,sizeKmer);
                    printf("... --> found kmer immediately after: %s %s \n",kmer_seq,current_strand?"(internally rc)":"");
                }
            }
        } 

    }

    //delete terminator;
    //delete traversal;
    free(right_traversal);
    //SolidKmers->close();
    free(node);
    fclose(linear_seqs_file);
 }


void IterativeExtensions::construct_linear_seqs_paired(string L, int max_depth, string output_file)
{
    construct_linear_seqs_paired(L, string(),max_depth,output_file,2,1000000000);
}

/*
 * return the contig which starts with L, as well as all ths contigs that follow him, up to max_depth.
 * results go to a fasta file (output_file)
 *
 * requires that:
 *  global variables terminator, bloo1 and false_positives are already constructed by minia
 *
 * outputs:
 *  a set of contigs in the output file */
void IterativeExtensions::construct_linear_seqs_paired(string L,string R, int max_depth, string output_file , int verb, int max_nodes)
{
    
    ///code for paired kmer
    int number_of_kmers_memorized = 1000; //should be greater than insert length, but maybe too much can be detrimental ?
    kmer_type  *  kmer_history = (kmer_type *) malloc (sizeof(kmer_type)*number_of_kmers_memorized);
    memset(kmer_history,0,sizeof(kmer_type)*number_of_kmers_memorized);
    int store_index =0;
    ///
    
    bool debug = verb>=2 ? 1 : 0;
    char bin2NT[4]    = {'A','C','T','G'};
    char bin2NTrev[4] = {'T','G','A','C'};
    char binrev[4]    = {2,3,0,1};
    
    
    kmer_type kmer;
    char kmer_seq[sizeKmer+1];
    FILE *linear_seqs_file;
    Traversal *traversal;
    if (traversal_type == IterativeExtensions::SimplePaths)
        traversal = new SimplePathsTraversal(bloo1,false_positives,terminator);
    else if (traversal_type == IterativeExtensions::Monument)
        traversal = new MonumentTraversal(bloo1,false_positives,terminator);
    traversal->set_maxlen(max_depth);
    traversal->set_max_depth(500);
    traversal->set_max_breadth(20);
    
    long long nbNodes = 0;
    long long totalnt=0;
    long long contig_len =0;
    long long maxlen=1000000;
    char *right_traversal = (char *) malloc(maxlen*sizeof(char));
    char *node          = (char *) malloc(maxlen*sizeof(char));
    node [0] = '\0';
    
    if(!output_file.empty())
    {
        linear_seqs_file = fopen((char * )output_file.c_str(),"w");
    }
    else {
        linear_seqs_file = NULL;
        IterativeExtensions::when_to_stop_extending = IterativeExtensions::After_first_contig; // sequence
    }
    
    STARTWALL(nodes);
    
    // heuristics: start from the last kmer of L
    
    vector < kmer_strand_depth > kmers_to_traverse;
    kmer_type kmer_fw, kmer_rc;
    
    kmer = extractKmerFromRead( (char *)L.c_str(), L.size() - sizeKmer, &kmer_fw, &kmer_rc, false);
    int strand = (kmer == kmer_rc);
    kmer_strand_depth first_kmer (kmer,strand,0);
    kmers_to_traverse.push_back(first_kmer);
    
#ifndef DONTMARK
    set<kmer_type> already_extended_from;
    compare_and_mark_last_k_minus_one_mer(L, already_extended_from); // mark first kmer to never extend from it again
#endif
    
    while (kmers_to_traverse.size() > 0) // min_depth is max_gap_length here
    {
        kmer_strand_depth ksd = kmers_to_traverse.back();
        kmers_to_traverse.pop_back();
        kmer = ksd.kmer;
        int strand = ksd.strand;
        int depth = ksd.depth;
        
        code2seq(kmer,kmer_seq); // convert starting kmer to nucleotide seq
        if (strand == 1)
            revcomp_sequence(kmer_seq,sizeKmer);
        
        if (debug)
            printf("iteration: kmer %s%s (of lenght %zu) depth %d nbNodes explored %lli\n",kmer_seq,strand?" (internally rc)":"", strlen(kmer_seq), depth,nbNodes);
        
        // right extension
        int len_right = traversal->traverse(kmer,right_traversal,strand);
        if (debug)
            printf("right traversal = %s\n", right_traversal);
        // save the node
        strcpy(node,kmer_seq);//               + starting_kmer if we ask for the de bruijn graph or for a text output
        strcat(node,right_traversal);//           + right_traversal
        int node_len=len_right+sizeKmer;
        fprintf(linear_seqs_file,">%lli__len__%i__depth__%i\n",nbNodes,node_len,depth);
        
        if (depth == 0 && dont_output_first_nucleotide )
        /* we need this in mapsembler:
         // the first used kmer should be a k-1 mer. Indeed, if a first kmer is extracted from a sequence :
         /// -------------****** (k=6), then this node is the one linked to a new one starting with ******, thus with an overlap of k and not k-1. */
	   fprintf(linear_seqs_file,"%s\n",node+1);
       else
	   fprintf(linear_seqs_file,"%s\n",node);
        
        
        
        ////memorize previous kmers  in the kmer_history, it is a circular buffer
        kmer_type kmer, graine, graine_revcomp ;
        for (int ii=0; ii<node_len-sizeKmer+1; ii++)
        {
            kmer = extractKmerFromRead(right_traversal,ii,&graine,&graine_revcomp);
            kmer_history[store_index]=kmer;
            store_index = (store_index + 1) % number_of_kmers_memorized ;
        }
        
        nbNodes++;
        totalnt+=node_len;
        
        if (debug)
        {
            printf ("[L=%s] assembled a %d bp node at depth %d\n",L.c_str(),node_len, depth);
            printf("(node seq: %s)\n",node);
        }
        
        // if we only want 1 extension, stop now
        if (when_to_stop_extending == IterativeExtensions::After_first_contig)
            break;
        
#ifdef STOPWHENFOUND
        char * found = NULL;
        found =strstr(node, R.c_str());
        if(found!=NULL)
            break;
#endif
        
        if(nbNodes > max_nodes) //GR stop when too complex  huum when to stop ?
            break;
        
        // if max depth reached, don't extend that one and clear kmer history, a new branch is explored
        if (depth + totalnt > max_depth)
        {
            memset(kmer_history,0,sizeof(kmer_type)*number_of_kmers_memorized);
            store_index = 0;
            continue;
        }
        
#ifndef DONTMARK
        // make sure this is the only time we see this (k-1)-overlap
        bool already_seen = compare_and_mark_last_k_minus_one_mer(node, already_extended_from);
        if (already_seen)
            continue;
#endif
        
        //this is the last kmer of the contig
        kmer = extractKmerFromRead( node, node_len - sizeKmer, &kmer_fw, &kmer_rc, false);
        strand = (kmer == kmer_rc);
        
        if (debug)
        {
            code2seq(kmer,kmer_seq);
            if (strand == 1)
                revcomp_sequence(kmer_seq,sizeKmer);
            printf("restarting from kmer %s %s\n",kmer_seq,strand?"(internally rc)":"");
        }
        
        // continue extending from immediately overlapping kmers
        // there may be just one 1 possibility (there was in-branching)
        
        //we are at the end of a contig, we push the starting point of all possible branches in a stack of kmers 
        //todo : push only kmers that are consistent with  pairing information
        //use kmer_history and the pairedBranchingKmers object for that
        
        //Raluca
        //first, compute the degree of the kmer named "kmer"
        int degree = 0;
        int current_strand;
        kmer_type current_kmer, aux_kmer;
        int nt1=-1, nt2=-1;
        bool ok_nt1=false, ok_nt2=false;
        
        for(int test_nt=0; test_nt<4; test_nt++)
        {
            current_strand = strand;
            current_kmer = next_kmer(kmer,test_nt, &current_strand);
                
            if(bloo1->contains(current_kmer) && !false_positives->contains(current_kmer))
            {
                degree++;
                nt1==-1 ? nt1 = test_nt : nt2 = test_nt;
            }
        }

        if(debug)
        {
            printf("degre %i \n",degree);
        }
        //second, depending on the degree and on whether we find the current_kmer in the history, we decide to push it or not
        if ( degree == 0 ) // this path ends, a new one will be explored, clear kmer history
        {
            memset(kmer_history,0,sizeof(kmer_type)*number_of_kmers_memorized);
            store_index = 0;
        }
        else
        if ( degree == 1 ) // not branching, is this possible knowing that we're at the end of a contig?
        {
            //if it is possible, then simply push the kmer obtained with nt1
            ok_nt1 = true;
            ok_nt2 = false;
        }
        else
        if (degree == 2 )
        {
            //we have the 2 extension, nt1 and nt2
            //first, get the entries in pairedBranchingKmers for kmer named "kmer"
            pair_nt_kmer_t val;
            int is_present = pairedBranchingKmers->get(kmer, &val);
            
            if (!is_present) // strange
            {
                //simply push the two possible branches
                ok_nt1 = ok_nt2 = true;
            }
            else
            {
                if(debug)
                {
                    printf("is_present,   %c  %c \n",bin2NT[nt1],bin2NT[nt2]);
                }
                
                int aux_nt;
                //search each of the two nucleotides in val of kmer
                for (int i=1; i<=2; i++)
                {
                    aux_kmer = -1;
                    i==1 ? aux_nt = bin2NT[nt1] : aux_nt = bin2NT[nt2];
                    
                    if(debug)
                    {
//                        char seq[100];
//                        code2seq(aux_kmer, seq);
                        printf("loop %i   aux_nt %c \n",i,aux_nt);
                    }
                    
                    if ( val.nk1.nt == aux_nt )
                        aux_kmer = val.nk1.prev_kmer;
                    else
                    if ( val.nk2.nt == aux_nt )
                        aux_kmer = val.nk2.prev_kmer;
                    else
                    {
                        i==1 ? ok_nt1 = false : ok_nt2 = false;
                        if ( i==2 && !ok_nt1 && !ok_nt2 )
                            break;
                    }
                    
                    if(debug)
                    {
                        printf("loop %i   ok_nt1 %i ok_nt2 %i  \n",i,ok_nt1,ok_nt2);
                    }
                    
                    //if one of the nt1 and nt2 corresponds, then search the corresponding prev_kmer (aux_kmer) in the history, if the prev_kmer is not present then ok_nt = false
                    if ( aux_kmer != -1 )
                    {
                        for ( int j=0; j<number_of_kmers_memorized; j++)
                            if ( kmer_history[j] == aux_kmer)
                            {
                                i==1 ? ok_nt1 = true : ok_nt2 = true;
                               
                                if(debug)
                                {
                                    char seq[100];
                                    code2seq(aux_kmer, seq);
                                    printf("found kmer   ok_nt1 %i ok_nt2 %i     %s \n",ok_nt1,ok_nt2,seq);
                                }
                                break;
                            }
                    
                        // if both ok_nt false, then push the two of them
                        if ( i==2 && !ok_nt1 && !ok_nt2 )
                            ok_nt1 = ok_nt2 = true;
                    }
                }
                    
            }
        }
        else //degree== 3 or 4
        {
            
            
            continue;
        }
        
        
        
        // if both ok_nt false, then push the two of them
        if ( !ok_nt1 && !ok_nt2 )
            ok_nt1 = ok_nt2 = true;
        
        
        
        
        //else ???
            // we cannot decide among the 3 or 4 knowing that our paired history keeps only 2 possible branchings per kmer TRUE? so we push them all?
        
        if ( ok_nt1 )
        {
            current_strand = strand;
            current_kmer = next_kmer(kmer, nt1, &current_strand);
                
            kmers_to_traverse.push_back(kmer_strand_depth (current_kmer,current_strand,depth + node_len));
            if (debug)
            {
                code2seq(current_kmer,kmer_seq); // convert starting kmer to nucleotide seq
                if (current_strand == 1)
                    revcomp_sequence(kmer_seq,sizeKmer);
                printf("push kmer immediately after: %s %s \n",kmer_seq,current_strand?"(internally rc)":"");
            }
        }
        
        if ( ok_nt2 )
        {
            current_strand = strand;
            current_kmer = next_kmer(kmer, nt2, &current_strand);
                
            kmers_to_traverse.push_back(kmer_strand_depth (current_kmer,current_strand,depth + node_len));
            if (debug)
            {
                code2seq(current_kmer,kmer_seq); // convert starting kmer to nucleotide seq
                if (current_strand == 1)
                    revcomp_sequence(kmer_seq,sizeKmer);
                printf("push kmer immediately after: %s %s \n",kmer_seq,current_strand?"(internally rc)":"");
            }
        }
        
    }
    
    free(kmer_history);
    //delete terminator;
    //delete traversal;
    free(right_traversal);
    //SolidKmers->close();
    free(node);
    fclose(linear_seqs_file);
}


/* requires:
 *  bloo1 and false_positives are already constructed by minia
 *
 * outputs:
 *  a pair containing: 
 *  1/ the extension
 *  2/ the length of the first unitig composing this sequence.
 *    - if the traversal is a simple traversal, then this length is equal to the length of the extension
 *    - else (if the traversal is a monument traversal), then this length is equal to the starting position of the first bubble (if exist)) 
 */
pair<char *, int> IterativeExtensions::extend_a_snp(string L, int max_depth)
{
//    printf("\nExtended sequ = %s\n", L.c_str()); // DEB
//    bool debug = true;
    bool debug = false;
    kmer_type kmer;
    char kmer_seq[sizeKmer+1];
    Traversal *traversal;
    if (traversal_type == IterativeExtensions::SimplePaths)
        traversal = new SimplePathsTraversal(bloo1,false_positives,terminator);
    else if (traversal_type == IterativeExtensions::Monument)
        traversal = new MonumentTraversal(bloo1,false_positives,terminator);
    traversal->set_maxlen(max_depth);
    traversal->set_max_depth(500);
    traversal->set_max_breadth(20);
    
    long long nbNodes = 0;
    long long totalnt=0;
    long long contig_len =0;
    long long maxlen=1000000;
    char *right_traversal = (char *) malloc(maxlen*sizeof(char));
    char *node          = (char *) malloc(maxlen*sizeof(char));
    node [0] = '\0';
    IterativeExtensions::when_to_stop_extending = IterativeExtensions::After_first_contig; // sequence
    
    
    STARTWALL(nodes);
    
    // heuristics: start from the last kmer of L
    
    kmer_type kmer_fw, kmer_rc;
    
    kmer = extractKmerFromRead( (char *)L.c_str(), L.size() - sizeKmer, &kmer_fw, &kmer_rc, false);
    int strand = (kmer == kmer_rc);
    
    code2seq(kmer,kmer_seq); // convert starting kmer to nucleotide seq
//    printf("Extended kmer = %s strand = %d\n", kmer_seq, strand); // DEB
    if (strand == 1)
        revcomp_sequence(kmer_seq,sizeKmer);
    
//    if (debug)
//        printf("iteration: kmer %s%s depth %d\n",kmer_seq,strand?" (internally rc)":"", 0);
    
    // right extension
    int len_right = traversal->traverse(kmer,right_traversal,strand);
    
    
//    if(debug) printf("right traversal = %s\n", right_traversal);
    // save the node
    strcpy(node,kmer_seq);//               + starting_kmer if we ask for the de bruijn graph or for a text output
    strcat(node,right_traversal);//           + right_traversal
    free(right_traversal);
    int first_divergence=strlen(node);
    if (traversal_type == IterativeExtensions::Monument && traversal->bubbles_positions.size()>0) {
        first_divergence=traversal->bubbles_positions[0].first+sizeKmer;
    }
    return std::make_pair(node,first_divergence);
}

