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


//
//  SNP.cpp
//

#ifndef ASSERTS
#define NDEBUG // disable asserts, they're computationnally intensive
#endif

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <algorithm> // for max

#include "SNP.h"
#include "../minia/Terminator.h"
#include "../minia/Traversal.h" // for extensions()
#include "filter.h"
#include "commons.h"
#include "IterativeExtensions.h"

extern int sizeKmer;
extern int threshold;
char bin2NT[8] = {'A','C','T','G','A','C','T','G'};
//#define DEBUG
FILE * SNP_file;
//bool already_an_extension;

using namespace::std;

unsigned long int nb_bubbles=0, nb_bubbles_high=0, nb_bubbles_low=0;

unsigned char Bubble::branching_structure(kmer_type graine)
{
    kmer_type next_kmer;
    unsigned char result = 0;
    int nt;
    int strand;
    
    for(nt=0; nt<4; nt++)
    {
        // forward right extensions
        strand=0;
        next_kmer = next_kmer_new(graine,nt,strand);
        if(bloom_solid_kmers->contains(next_kmer) && !debloom->contains(next_kmer)){
            result|=1<<nt;
        }
        
        // reverse right extensions
        strand=1;
        next_kmer = next_kmer_new(graine,nt,strand);
        if(bloom_solid_kmers->contains(next_kmer) && !debloom->contains(next_kmer)){
            result|=1<<(nt+4);
        }
    }
    return result;
}


//extend two 2k-1 paths with a new nucleotide left and a new nucleotide right.
//returns  0: no unique extension neither left nor right, 1 unique only left, 2 unique only right, 3: both
int Bubble::close_snp(const char * path1, const char * path2, char *path1_c, char *path2_c){
    kmer_type kmer_fw, kmer_rc;
    int strand;
    int i;
    kmer_type kmer;
    char right_extension=-1, left_extension=-1;
    
    
    char * seq1 = (char *) malloc((2*sizeKmer)*sizeof(char));
    for(i=0;i<2*sizeKmer-1;i++) seq1[i]=bin2NT[path1[i]];
    seq1[2*sizeKmer-1]='\0';
    
    //RIGHT EXTENSION
    kmer = extractKmerFromRead(seq1, strlen(seq1) - sizeKmer, &kmer_fw, &kmer_rc, false);
    strand = (kmer == kmer_rc);
    
    bool left_extensible;
    bool right_extensible;
    
    int nb_right_extensions = 0;
    for(int test_nt=0; test_nt<4; test_nt++)
    {
        int current_strand = strand;
        kmer_type current_kmer = next_kmer(kmer,test_nt, &current_strand);
        if(bloo1->contains(current_kmer) && !false_positives->contains(current_kmer)){
            right_extension=test_nt;
            nb_right_extensions++;
            //break;
        }
    }
    
    if (nb_right_extensions == 1) // one right extension found, extend to the right
        right_extensible=true;
    else
        right_extensible=false;
    
    
    // LEFT EXTENSION
    int nb_left_extensions = 0;
    revcomp(seq1, strlen(seq1));
    kmer = extractKmerFromRead(seq1, strlen(seq1) - sizeKmer, &kmer_fw, &kmer_rc, false);
    strand = (kmer == kmer_rc);
    for(int test_nt=0; test_nt<4; test_nt++)
    {
        int current_strand = strand;
        kmer_type current_kmer = next_kmer(kmer,test_nt, &current_strand);
        if(bloo1->contains(current_kmer) && !false_positives->contains(current_kmer)){
            left_extension=test_nt;
            nb_left_extensions++;
            //break;
        }
    }
    revcomp(seq1, strlen(seq1));
    
    if (nb_left_extensions == 1 ) // one left extension found, extend to the left
        left_extensible = true;
    else
        left_extensible=false;
    
    free(seq1);
    //returns 0 no extension neither left nor right, 1 unique only left, 2 unique only right, unique 3 both
    
    // Add the new results. (revcomp of left extension and rigth extension).
    // put the found bubble in the middle of the 2k+1 length array
    // DO NOT MODIFY path1 AND path2 as we still need them outside this method
    if (!left_extensible && !right_extensible){//no extension, 2k-1 paths
        for(i=0;i<2*sizeKmer-1;i++) {
            path1_c[i]=bin2NT[path1[i]];
            path2_c[i]=bin2NT[path2[i]];
        }
        path1_c[i]='\0';
        path2_c[i]='\0';
        return 0;
    }
    if (left_extensible && !right_extensible){//2k path, only left
        path1_c[0]=bin2NT[revcomp_int(left_extension)];
        path2_c[0]=bin2NT[revcomp_int(left_extension)];
        for(i=1;i<=2*sizeKmer-1;i++) {
            path1_c[i]=bin2NT[path1[i-1]];
            path2_c[i]=bin2NT[path2[i-1]];
        }
        path1_c[i]='\0';
        path2_c[i]='\0';
        
        return 1;
        
    }
    if (!left_extensible && right_extensible){//2k path, only right
        for(i=0;i<2*sizeKmer-1;i++) {
            path1_c[i]=bin2NT[path1[i]];
            path2_c[i]=bin2NT[path2[i]];
        }
        path1_c[i]=bin2NT[right_extension];
        path2_c[i]=bin2NT[right_extension];
        i++;
        path1_c[i]='\0';
        path2_c[i]='\0';
        
        return 2;
        
    }
    
    if (left_extensible && right_extensible) {
        path1_c[0]=bin2NT[revcomp_int(left_extension)];
        
        path2_c[0]=path1_c[0];
        for(i=1;i<2*sizeKmer;i++) {
            path1_c[i]=bin2NT[path1[i-1]];
            path2_c[i]=bin2NT[path2[i-1]];
        }
        path1_c[i]=bin2NT[right_extension];
        path2_c[i]=bin2NT[right_extension];
        i++;
        path1_c[i]='\0';
        path2_c[i]='\0';
        
        return 3;
    }
    
    assert(1==0); // should not come here
}


// prints the headers depending of the extension choice.
inline void prints_contig_informations(FILE * file, std::pair<char*,int> left_extension, std::pair<char*,int> right_extension, bool extend_snps, bool strict_extension, int where_to_extend){
    if (extend_snps && strict_extension){
        if(where_to_extend == 0)
            fprintf(file,"|left_unitig_length_0|right_unitig_length_0");
        if(where_to_extend == 1)
            fprintf(file,"|left_unitig_length_%d|right_unitig_length_0", (int)strlen(left_extension.first)-sizeKmer+1); //+1 because of close_snp
        if(where_to_extend == 2)
            fprintf(file,"|left_unitig_length_0|right_unitig_length_%d", (int)strlen(right_extension.first)-sizeKmer+1);//+1 because of close_snp
        if(where_to_extend == 3)
            fprintf(file,"|left_unitig_length_%d|right_unitig_length_%d", (int)strlen(left_extension.first)-sizeKmer+1,(int)strlen(right_extension.first)-sizeKmer+1);//+1 because of close_snp
    }
    
    if(extend_snps && !strict_extension){
        
        if(where_to_extend == 0){
            fprintf(file,"|left_unitig_length_0|right_unitig_length_0");
            fprintf(file,"|left_contig_length_0|right_contig_length_0");
        }
        if(where_to_extend == 1){
            fprintf(file,"|left_unitig_length_%d|right_unitig_length_0", left_extension.second-sizeKmer+1); //+1 because of close_snp
            fprintf(file,"|left_contig_length_%d|right_contig_length_0", (int)strlen(left_extension.first)-sizeKmer+1); //+1 because of close_snp
        }
        if(where_to_extend == 2){
            fprintf(file,"|left_unitig_length_0|right_unitig_length_%d", right_extension.second-sizeKmer+1);//+1 because of close_snp
            fprintf(file,"|left_contig_length_0|right_contig_length_%d", (int)strlen(right_extension.first)-sizeKmer+1);//+1 because of close_snp
        }
        if(where_to_extend == 3){
            fprintf(file,"|left_unitig_length_%d|right_unitig_length_%d", left_extension.second-sizeKmer+1,right_extension.second-sizeKmer+1);//+1 because of close_snp
            fprintf(file,"|left_contig_length_%d|right_contig_length_%d", (int)strlen(left_extension.first)-sizeKmer+1,(int)strlen(right_extension.first)-sizeKmer+1);//+1 because of close_snp
        }
    }
}


/**
 * Extends if necessary left and/or right parts of the bubble.
 * Prints results
 * path1 and path2: character pathes (no binary format). Of length 2k-1, 2k, or 2k+1, depending ont the results of "close_snps"
 * score=complexity score
 * where_to_extend : 0=nothing, 1=left only, 2=right only, 3=both
 */
void Bubble::print_sequence_and_eventually_contigs(char * path1, char * path2, const int score, int where_to_extend){
    int i;
    terminator->reset(); // need to reset branching kmers between extensions in case of overlapping extensions.
    FILE * file;
    //Here, path1 and path2 may have 2k, 2k-1 or 2k+1 size
    int size_path = strlen(path1);
    
    
    pair<char*, int> right_contig;
    pair<char*, int> left_contig;
    if(extend_snps){// && where_to_extend > 0){
        if(strict_extension)
            IterativeExtensions::traversal_type = IterativeExtensions::SimplePaths; // strict
        else
            IterativeExtensions::traversal_type = IterativeExtensions::Monument; // contigs
        if (where_to_extend == 2 || where_to_extend == 3)
            right_contig=IterativeExtensions::extend_a_snp(path1,10000); // path1 or path2 provides the same result (the extended kmer is not involved)
        if (where_to_extend == 1 || where_to_extend == 3)
        {
            revcomp(path1, size_path); // get the reverse complement of the sequence
            left_contig=IterativeExtensions::extend_a_snp(path1,10000); // get the left extension of this sequence
            revcomp(left_contig.first, strlen(left_contig.first)); // put back the left contig in the right order
            revcomp(path1, size_path); // put back the sequence in the right order.
        }
        
        if(min_size_extension >-1){
            if(where_to_extend<3 || (int)strlen(left_contig.first)-sizeKmer<min_size_extension || (int)strlen(right_contig.first)-sizeKmer<min_size_extension) {
                //            printf("Don't output this, too short\n");
                free(left_contig.first);
                free(right_contig.first);
                // FIXME: do I have to free the pairs themselves ?
                return;
            }
        }
    }
    
    file = SNP_file;
    
    //fprintf(file, ">SNP_higher_path_%lu|score_%d|", nb_bubbles, score);
    fprintf(file, ">SNP_higher_path_%lu|", nb_bubbles);
    if ( score >= threshold )
        fprintf(file,"low");
    else
    {
        fprintf(file,"high");
    }
    
    
    prints_contig_informations(file, left_contig, right_contig, extend_snps, strict_extension, where_to_extend);
    fprintf(file,"\n");
    
    if (print_extensions && extend_snps){
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
    else{
        int start = 0;
        int stop = strlen(path1);
        if(!print_extensions){
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
    
    if (print_extensions && extend_snps){
        if(print_extensions && where_to_extend%2==1) for(i=0;i<strlen(left_contig.first)-sizeKmer;i++) fprintf(file,"%c", tolower(left_contig.first[i]));
        // change the case of first and/or last character of the central node
        if (where_to_extend%2==1) // left: first nucleotide is an extension
            path2[0]=tolower(path2[0]);
        if (where_to_extend>1) // right: last nucleotide is an extension
            path2[size_path-1]=tolower(path2[size_path-1]);
        fprintf(file,"%s",path2);
        if(extend_snps && print_extensions && where_to_extend>1) for(i=sizeKmer;i<strlen(right_contig.first);i++) fprintf(file,"%c", tolower(right_contig.first[i]));
    }
    else{
        int start = 0;
        int stop = strlen(path2);
        if(!print_extensions){
            if(where_to_extend%2==1) // left: first nucleotide is an extension
                start++;
            if (where_to_extend>1) // right: last nucleotide is an extension
                stop--;
        }
        for(i=start;i<stop;i++) fprintf(file,"%c", path2[i]);
    }
    fprintf(file,"\n");
    
    if(extend_snps && where_to_extend%2==1){
        free(left_contig.first);
    }
    if(extend_snps && where_to_extend>1){
        free(right_contig.first);
    }
}


// determines if a kmer is branching or not, function copied (I know, it's horrible) from Terminator.cpp in minia
//bool Bubble::is_branching(kmer_type kmer)
//{
//    bool reversed = false;
//    if (kmer > revcomp(kmer)) {
//        kmer = revcomp(kmer);
//        reversed = true;
//    }
//
//    // cannot really be optimized, because most kmers will be non-branching, hence computing branching_structure() takes optimal time
//    int nb_forward_links = 0, nb_reverse_links = 0;
//    int i;
//    unsigned char branching = branching_structure(kmer);
//
//    for (i=0;i<4;i++)
//        nb_forward_links += (branching>>i)&1;
//
//    for (i=4;i<8;i++)
//        nb_reverse_links += (branching>>i)&1;
//
//    if(reversed) kmer = revcomp(kmer); // go back to the good value.
//    return !(nb_forward_links == 1 && nb_reverse_links == 1);
//
//}


//void Bubble::expand_BubblePierre(int direction, char* path1, char* path2, kmer_type kmer1, kmer_type kmer2, int pos, char *p1, char *p2, char *p3, char *p4, int low, int authorised_branching){
//    kmer_type next_kmer1, next_kmer2;
//    int nt,nt2,i,score;
//    char p_aux1[sizeKmer+1], p_aux2[sizeKmer+1];
//
//
//
//    if (pos <= sizeKmer-1)
//    {
//        for(nt=0; nt<4; nt++)
//        {
//            next_kmer1 = next_kmer_new(kmer1,nt,direction);
//            next_kmer2 = next_kmer_new(kmer2,nt,direction);
//
//            if(authorised_branching &&
//               (is_branching(next_kmer1) || is_branching(next_kmer2))
//               ) continue; // one of the two new kmers is branching, we continue.
//
//            code2seq(next_kmer1, p_aux1);
//            code2seq(next_kmer2, p_aux2);
//#ifdef DEBUG
//            if ( !(strcmp(p_aux1,p1) && strcmp(p_aux1,p2) && strcmp(p_aux2,p3) && strcmp(p_aux2,p4)) )
//            {
//                if ( !strcmp(p_aux1,p1) )
//                    printf("\n1st branch: 1st and 3rd equal kmers %s\n", p_aux1);
//                if ( !strcmp(p_aux1,p2) )
//                    printf("\n1st branch: 2nd and 3rd equal kmers %s\n", p_aux1);
//                if ( !strcmp(p_aux2,p3) )
//                    printf("\n2nd branch: 1st and 3rd equal kmers %s\n", p_aux2);
//                if ( !strcmp(p_aux2,p4) )
//                    printf("\n2nd branch: 2nd and 3rd equal kmers %s\n", p_aux2);
//            }
//#endif
//            if(bloom_solid_kmers->contains(next_kmer1) && !debloom->contains(next_kmer1) && bloom_solid_kmers->contains(next_kmer2) && !debloom->contains(next_kmer2)
//               && strcmp(p_aux1,p1) && strcmp(p_aux1,p2) && strcmp(p_aux2,p3) && strcmp(p_aux2,p4))
//            {
//
//                strcpy(p1, p2);
//                strcpy(p3, p4);
//                strcpy(p2, p_aux1);
//                strcpy(p4, p_aux2);
//
//                if ( direction == 0 ){
//                    //if everything is ok then we only use direction = 0
//                    path1[sizeKmer-1+pos]=nt;
//                    path2[sizeKmer-1+pos]=nt;
//                }
//                else
//                {
//                    //shouldn't go in here
//                    path1[sizeKmer-1+pos]=nt+4;
//                    path2[sizeKmer-1+pos]=nt+4;
//                }
//
//                //if we finished the bubble and need to close it
//                if ( pos == sizeKmer-1 )
//                {
//                    //TEST whether the first kmer of the first path is smaller than the first kmer of the revcomp(first path), this should avoid repeated SNPs
//                    char first_kmer[sizeKmer+1], first_kmer_rev[sizeKmer+1];
//                    for ( i=0; i<sizeKmer; i++ )
//                        first_kmer[i] = bin2NT[path1[i]];
//                    first_kmer[sizeKmer]='\0';
//
//                    for ( i=sizeKmer-1; i<2*sizeKmer-1; i++ )
//                        first_kmer_rev[i-sizeKmer+1] = bin2NT[path1[i]];
//                    first_kmer_rev[sizeKmer]='\0';
//                    revcomp(first_kmer_rev, sizeKmer);
//
//                    if (strcmp(first_kmer, first_kmer_rev)<0)
//                    {
//                        nb_bubbles++;
//
//                        score = filterLowComplexity2Paths(path1, path2, 2*sizeKmer-1, threshold);
//
//                        if ( score < threshold || (score>=threshold && low))
//                        {
//                            char path1_c[2*sizeKmer+1], path2_c[2*sizeKmer+1];
//                            if(close_snp(path1, path2, path1_c, path2_c)) // PRINT A SNP Only if we can close it right and left
//                            {
//                                if (score < threshold)
//                                    nb_bubbles_high++;
//                                else
//                                    nb_bubbles_low++;
//
//                                print_sequence_and_eventually_contigs(path1_c, path2_c, score);
//                            }
//                        }
//                    }
//                }
//                else //bubble not yet finished
//                {
//
//                    next_kmer1 = next_kmer_new_norevcomp(kmer1,nt,direction);
//                    next_kmer2 = next_kmer_new_norevcomp(kmer2,nt,direction);
//                    expand_Bubble(direction, path1, path2, next_kmer1, next_kmer2, pos+1, p1, p2, p3, p4, low, authorised_branching);
//                }
//            }
//        }
//    }
//}
bool Bubble::two_possible_extensions_on_one_path(kmer_type kmer){
    kmer_type next_kmer;
    bool already_an_extension;
    for (int d=0; d<=1; d++)
    {
        already_an_extension = false;
        for(int nt=0; nt<4; nt++)
        {
            next_kmer = next_kmer_new(kmer,nt,d);
            if(bloom_solid_kmers->contains(next_kmer) && !debloom->contains(next_kmer)){
                if (already_an_extension)
                    return true;
                already_an_extension = true;
            }
        }
    }
    return false;
}

bool Bubble::two_possible_extensions(kmer_type kmer1, kmer_type kmer2){
    kmer_type next_kmer1, next_kmer2;
    char p_aux1[sizeKmer+1], p_aux2[sizeKmer+1];
    bool already_an_extension = false;
    
    for (int d=0; d<=1; d++)
    {
        
        for(int nt=0; nt<4; nt++)
        {
            next_kmer1 = next_kmer_new(kmer1,nt,d);
            next_kmer2 = next_kmer_new(kmer2,nt,d);
            
            code2seq(next_kmer1, p_aux1);
            code2seq(next_kmer2, p_aux2);
            
            if(bloom_solid_kmers->contains(next_kmer1) && !debloom->contains(next_kmer1) && bloom_solid_kmers->contains(next_kmer2) && !debloom->contains(next_kmer2)
               )//&& strcmp(p_aux1,p1) && strcmp(p_aux1,p2) && strcmp(p_aux2,p3) && strcmp(p_aux2,p4))
            {
                if (!already_an_extension)
                {
#ifdef DEBUG
                    printf("\nFirst branching in direction %d with %c at pos %d\n", d, bin2NT[nt], pos);
#endif
                    already_an_extension = true;
                }
                else
                {
#ifdef DEBUG
                    code2seq(kmer1, p_aux1);
                    code2seq(kmer2, p_aux2);
                    printf("\nDon't expand, branching in direction %d with %c for %s and %s at pos %d\n", d, bin2NT[nt], p_aux1, p_aux2, pos);
#endif
                    return true;
                }
            }
        }
        already_an_extension = false;
    }
    return false;
}


/****************************************************************************************************************
 //the 2 functions that follow build bubbles starting from the first node in the bubble (after the switching node)
 *****************************************************************************************************************/
/**
 * ...
 * authorised_branching =
 *  - 0: branching forbiden in any path
 *  - 1: same branching on both path forbiden (i.e. 2 disctinct nucelotides may be used in both paths for extension)
 *  - 2: no restriction on branching
 */
void Bubble::expand_Bubble(int direction, char* path1, char* path2, kmer_type kmer1, kmer_type kmer2, int pos, char *p1, char *p2, char *p3, char *p4, int low, int authorised_branching){
    kmer_type next_kmer1, next_kmer2;
    int nt,nt2,i,score;
    char p_aux1[sizeKmer+1], p_aux2[sizeKmer+1];
    
    
    //  ---x--- (2k-1)
    // 0123456
    //    *    (size kmer -1)
    
    
    if (pos <= sizeKmer-1)
    {
        
        // stop the extension if authorised_branching==0 (not branching in any path) and any of the two paths is branching
        if (authorised_branching==0 && (two_possible_extensions_on_one_path(kmer1) || two_possible_extensions_on_one_path(kmer2))) return;
            
        // stop the extension if authorised_branching==1 (not branching in both path) and both the two paths are branching
        if (authorised_branching==1 && two_possible_extensions
            (kmer1, kmer2)) return;
        
        
        
        for(nt=0; nt<4; nt++)
        {
            next_kmer1 = next_kmer_new(kmer1,nt,direction);
            next_kmer2 = next_kmer_new(kmer2,nt,direction);
            code2seq(next_kmer1, p_aux1);
            code2seq(next_kmer2, p_aux2);
            
#ifdef DEBUG
            printf("\nTry %d at pos %d\n for kmers %s %s\n", nt, pos, p_aux1, p_aux2);
            if ( !(strcmp(p_aux1,p1) && strcmp(p_aux1,p2) && strcmp(p_aux2,p3) && strcmp(p_aux2,p4)) )
            {
                if ( !strcmp(p_aux1,p1) )
                    printf("\n1st branch: 1st and 3rd equal kmers %s\n", p_aux1);
                if ( !strcmp(p_aux1,p2) )
                    printf("\n1st branch: 2nd and 3rd equal kmers %s\n", p_aux1);
                if ( !strcmp(p_aux2,p3) )
                    printf("\n2nd branch: 1st and 3rd equal kmers %s\n", p_aux2);
                if ( !strcmp(p_aux2,p4) )
                    printf("\n2nd branch: 2nd and 3rd equal kmers %s\n", p_aux2);
            }
#endif
            if(bloom_solid_kmers->contains(next_kmer1) && !debloom->contains(next_kmer1) && bloom_solid_kmers->contains(next_kmer2) && !debloom->contains(next_kmer2)
               && strcmp(p_aux1,p1) && strcmp(p_aux1,p2) && strcmp(p_aux2,p3) && strcmp(p_aux2,p4))
            {
                strcpy(p1, p2);
                strcpy(p3, p4);
                strcpy(p2, p_aux1);
                strcpy(p4, p_aux2);
                
                if ( direction == 0 ){
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
                if ( pos == sizeKmer-1 )
                {
                    //TEST whether the first kmer of the first path is smaller than the first kmer of the revcomp(first path), this should avoid repeated SNPs
                    
                    char first_kmer[sizeKmer+1], first_kmer_rev[sizeKmer+1];
                    for ( i=0; i<sizeKmer; i++ )
                        first_kmer[i] = bin2NT[path1[i]];
                    first_kmer[sizeKmer]='\0';
                    
                    for ( i=sizeKmer-1; i<2*sizeKmer-1; i++ )
                        first_kmer_rev[i-sizeKmer+1] = bin2NT[path1[i]];
                    first_kmer_rev[sizeKmer]='\0';
                    revcomp(first_kmer_rev, sizeKmer);
                    
                    if (strcmp(first_kmer, first_kmer_rev)<0)
                    {
                        
                        
                        
                        // stop the extension if authorised_branching==0 (not branching in any path) and any of the two paths is branching on last kmers
                        if (authorised_branching==0 && (two_possible_extensions_on_one_path(next_kmer_new_norevcomp(kmer1,nt,direction)) || two_possible_extensions_on_one_path(next_kmer_new_norevcomp(kmer2,nt,direction)))) return;
                        // stop the extension if authorised_branching==1 (not branching in both path) and both the two paths are branching on last kmers
                        if(authorised_branching==1 && two_possible_extensions(next_kmer_new_norevcomp(kmer1,nt,direction), next_kmer_new_norevcomp(kmer2,nt,direction))) return;
                        
                        nb_bubbles++;
                        
                        // TODO: avoid this test if low is not required.
                        score=filterLowComplexity2Paths(path1, path2, 2*sizeKmer-1, threshold);
                        
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
                                for(i=0;i<2*sizeKmer-1;i++) {
                                    path1_c[i]=bin2NT[path1[i]];
                                    path2_c[i]=bin2NT[path2[i]];
                                }
                                path1_c[i]='\0';
                                path2_c[i]='\0';
                                
                            }
                            else
                                where_to_extend = close_snp(path1, path2, path1_c, path2_c);
                            
                            print_sequence_and_eventually_contigs(path1_c, path2_c, score, where_to_extend); //2k-1, 2k or 2k+1 paths
                            
                        }
                    }
                    
                }
                else //bubble not yet finished
                {
                    
                    next_kmer1 = next_kmer_new_norevcomp(kmer1,nt,direction);
                    next_kmer2 = next_kmer_new_norevcomp(kmer2,nt,direction);
                    code2seq(next_kmer1, p_aux1);
                    code2seq(next_kmer2, p_aux2);
#ifdef DEBUG
                    printf("\nExpanded with %d at pos %d\n for kmers %s %s\n", nt, pos, p_aux1, p_aux2);
#endif
                    expand_Bubble(direction, path1, path2, next_kmer1, next_kmer2, pos+1, p1, p2, p3, p4, low, authorised_branching);
#ifdef DEBUG
                    code2seq(next_kmer1, p_aux1);
                    code2seq(next_kmer2, p_aux2);
                    printf("\nReturn to expansion with %d at pos %d\n for kmers %s %s\n", nt, pos, p_aux1, p_aux2);
#endif
                    //there's only one branch to expand if we keep non branching SNPs only, therefore we can safely stop the for loop
                    if ( authorised_branching==0 || authorised_branching==1 ) break;
                }
            }
        }
    }
    
}


//tries to start a bubble from the current node, first node containing a SNP
void Bubble::start_Bubble(kmer_type kmer1, int low, int authorised_branching)
{
    char seq1[256], seq2[256];//, seq_aux[26];
    int i,j, direction;
    char path1[2*sizeKmer-1], path2[2*sizeKmer-1];
    kmer_type kmer2, kmer_aux;
    char p1[sizeKmer+1], p2[sizeKmer+1], p3[sizeKmer+1], p4[sizeKmer+1];
    
    
#ifdef DEBUG
    {
        printf("\nInitial kmer ");
        code2seq(kmer1, seq1);
        for (j=0;j<=sizeKmer-1;j++)
            printf("%c",seq1[j]);
    }
#endif
    for (direction=0; direction<=1; direction++)
    {
        if ( direction == 1 )
        {
            // start bubble for kmer in reverse form
            kmer1 = revcomp(kmer1);
        }
        code2seq(kmer1, seq1);
        code2seq(kmer1, seq2);
        
        for (i=0; i<sizeKmer-1; i++)
        {
            path1[i] = NT2int(seq1[i]);
            path2[i] = NT2int(seq1[i]);
        }
        path1[sizeKmer-1] = NT2int(seq1[sizeKmer-1]);
        
        for (i=NT2int(seq1[sizeKmer-1])+1;i<4;i++) // try all the possible extensions that were not previously tested (clever :-))
        {
            seq2[sizeKmer-1]=bin2NT[i];
            
            //BETTER NOT USE IT
            //avoid some double SNPs: do only A-C, A-G couples instead of G-T and C-T. The problem comes from A-T and C-G that are still double.
            /*if ((seq1[sizeKmer-1]=='G' || seq1[sizeKmer-1]=='C') && bin2NT[i]=='T')
             continue;*/
            
            kmer2 = codeSeed(seq2);
            
#ifdef DEBUG
            {
                printf("\nTry %c ", bin2NT[i]);
                
                for (j=0;j<=sizeKmer-1;j++)
                    printf("%c",seq2[j]);
                printf(" ");
                for (j=0;j<=sizeKmer-1;j++)
                    printf("%c",seq1[j]);
            }
#endif
            if (bloom_solid_kmers->contains(min(kmer2,revcomp(kmer2))) && !debloom->contains(min(kmer2,revcomp(kmer2)))) // the tried kmer is indexed.
            {
                path2[sizeKmer-1] = NT2int(seq2[sizeKmer-1]);
                
                strcpy(p1,"\0");
                strcpy(p3,"\0");
                code2seq(kmer1, p2);
                code2seq(kmer2, p4);
#ifdef DEBUG
                printf("\nCouple with %c \n", bin2NT[i]);
#endif
                //                    already_an_extension = false;
                expand_Bubble(0, path1, path2, kmer1, kmer2, 1, p1, p2, p3, p4, low, authorised_branching); // we open a new putative bubble
                
#ifdef DEBUG
                printf("\n Starting Bubble\n");
                if ( direction == 1)
                    printf("reverse\n");
                for (j=1;j<=sizeKmer;j++)
                    printf("%c",bin2NT[path1[j]]);
                printf("\n");
                for (j=1;j<=sizeKmer;j++)
                    printf("%c",bin2NT[path2[j]]);
#endif
                
            }
            else
            {
#ifdef DEBUG
                printf("\ndoes not contain it");
#endif
            }
        }
    }
    
}


void Bubble::find_bubbles(const char * SNP_file_name, int low, int authorised_branching)
{
    kmer_type kmer;
    printf("\nSearch for bubbles with threshold %d\n", threshold);
    
    SNP_file=fopen(SNP_file_name,"w");
    if(SNP_file == NULL){
        fprintf(stderr,"cannot open file %s, exit\n",SNP_file_name);
        exit(1);
    }
    
    printf("Find Bubbles:\n");
    
    /*char path1[51], path2[51], seq_aux[52];
     int i;
     strcpy(seq_aux, "TGTGTGCACACAAGAGTGGCCAGGGTGGAGCCTCGGGATCGGTTGGATAGT");
     for (i=0; i<2*sizeKmer+1; i++)
     {
     path1[i] = NT2int(seq_aux[i]);
     }
     printf("\n%d = score of %s\n",filterLowComplexity(path1, 51, 90), seq_aux );
     strcpy(seq_aux, "TGTGTGCACACAAGAGTGGCCAGGGGGGAGCCTCGGGATCGGTTGGATAGT");
     for (i=0; i<2*sizeKmer+1; i++)
     {
     path2[i] = NT2int(seq_aux[i]);
     }
     printf("\n%d = score of %s\n",filterLowComplexity(path2, 51, 90), seq_aux );
     printf("\n%d = score total\n",filterLowComplexity2Paths(path1, path2, 51, 90));
     
     
     strcpy(seq_aux, "ACTATCCAACCGATCCCGAGGCTCCACCCTGGCCACTCTTGTGTGCACACA");
     for (i=0; i<2*sizeKmer+1; i++)
     {
     path1[i] = NT2int(seq_aux[i]);
     }
     printf("\n%d = score of %s\n",filterLowComplexity(path1, 51, 90),seq_aux );
     strcpy(seq_aux, "ACTATCCAACCGATCCCGAGGCTCCCCCCTGGCCACTCTTGTGTGCACACA");
     for (i=0; i<2*sizeKmer+1; i++)
     {
     path2[i] = NT2int(seq_aux[i]);
     }
     printf("\n%d = score of %s\n",filterLowComplexity(path2, 51, 90),seq_aux );
     printf("\n%d = score total\n",filterLowComplexity2Paths(path1, path2, 51, 90));
     
     exit(0);
     
     char start_kmer [] = "AACAACAGGTGCTGGAGAGGATGTGGAGAAA";
     kmer = codeSeed(start_kmer);
     start_Bubble(kmer, low, authorised_branching);*/
    
    unsigned long int nb=0;
    off_t nbel = SolidKmers->nb_elements();
    while (SolidKmers->read_element(&kmer))
    {
        if(nb%1000==0)
        {printf("%c \t %llu %%",13,(nb*100)/nbel);
        }
        nb++;
        //start_SwitchingNode(kmer);
        start_Bubble(kmer, low, authorised_branching);
    }
    
    printf("%c \t 100 %%\n", 13);
    fclose(SNP_file);
    
    printf("\nFound %lu bubbles (not branching if demanded and with no restriction on closing). Among these, we select the closing bubbles with non-branching kmers (if demanded), from which %lu are high complexity bubbles and %lu are low complexity bubbles (by default, this is 0 if low bubbles do not have to be printed)\n", nb_bubbles, nb_bubbles_high, nb_bubbles_low);
}


/*void Bubble::read_bubble_file(char *bubble_filename)
 {
 FILE *bubble_file;
 char line[512], line2[512];
 char SNP[2*sizeKmer-1];
 int i, score, hc_bubbles=0, lc_bubbles=0;
 printf("kissnp is opeing file %s\n", bubble_filename);
 bubble_file = fopen(bubble_filename,"r");
 if (bubble_file == NULL)
 {
 printf("kissnp error opening file: %s\n",bubble_filename);
 exit(1);
 }
 
 while(fgets(line, 512, bubble_file) != NULL){
 if ( strstr(line, "Bubble") != NULL)
 {
 strcpy(line2, line);
 if ( fgets(line, 512, bubble_file) != NULL )
 {
 for (i=0; i<2*sizeKmer-1; i++)
 SNP[i]=NT2int(line[i+1]);
 
 score = filterLowComplexity(SNP, 2*sizeKmer-1,threshold);
 if ( score >= 2*threshold )
 {
 printf("Low complexity bubble with score %d\n", score);
 lc_bubbles++;
 }
 else
 {
 printf("%s", line2);
 printf("High complexity bubble with score %d\n", score);
 hc_bubbles++;
 printf("%s\n", line);
 }
 
 }
 }
 }
 printf("%d high complexity bubbles and %d low complexity bubbles, with threshold %d", hc_bubbles, lc_bubbles, threshold);
 fclose(bubble_file);
 }*/
