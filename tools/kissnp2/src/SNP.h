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
//  SNP.h

#ifndef SNP_H_INCLUDED
#define SNP_H_INCLUDED


#include <stdlib.h>
#include <inttypes.h>
#include <stdint.h>
#include <cmath> // for log2f


#include "../minia/Bloom.h"
#include "../minia/Kmer.h"
#include "Kmer_for_kissnp2.h"
#include "../minia/Set.h"
#include "../minia/Bank.h"



class Bubble{

protected:

    BinaryBank *SolidKmers;
    Bloom *bloom_solid_kmers;
    Set *debloom;
    bool extend_snps;
    bool print_extensions;
    bool strict_extension; // true: strict extension, fales: contig extension
    int min_size_extension;
    
    /**
     * Extends if necessary left and/or right parts of the bubble.
     * Prints results
     * path1 and path2: character pathes (no binary format). Of length 2k-1, 2k, or 2k+1, depending ont the results of "close_snps"
     * score=complexity score
     * where_to_extend : 0=nothing, 1=left only, 2=right only, 3=both
     */
    void print_sequence_and_eventually_contigs(char * path1, char * path2, const int score, int where_to_extend); //extension==0 no extension, =1 only left, =2 only right, =3 both
    unsigned char branching_structure(kmer_type graine);
    bool is_branching(kmer_type kmer);
    bool two_possible_extensions_on_one_path(kmer_type);
    int close_snp(const char * path1, const char * path2, char * path1_c, char * path2_c); //returns -1 not closed, 0 no unique extension, 1 only left, 2 only right, 3 both
public:
    void find_bubbles(const char * SNP_file_name, int low, int authorised_branching);
    void start_SwitchingNode(kmer_type graine);
    void start_Bubble(kmer_type graine, int low, int authorised_branching);
    bool two_possible_extensions(kmer_type kmer1, kmer_type kmer2);
//    void expand_Bubble_till_SN(int direction, char* path1, char* path2, kmer_type kmer1, kmer_type kmer2, int pos);
    void expand_Bubble(int direction, char* path1, char* path2, kmer_type kmer1, kmer_type kmer2, int pos, char* p1, char* p2, char* p3, char* p4, int low, int no_branching);
    bool next(kmer_type *kmer);
    void read_bubble_file(char *bubble_filename);

    Bubble(BinaryBank *given_SolidKmers, Bloom *given_bloom, Set *given_debloom, bool extend_snps, int min_size_extension,bool print_extensions, bool strict_extension) : SolidKmers(given_SolidKmers), bloom_solid_kmers(given_bloom), debloom(given_debloom), extend_snps(extend_snps), min_size_extension(min_size_extension), print_extensions (print_extensions), strict_extension (strict_extension) { }
    Bubble() { }
};


#endif // SNP_H_INCLUDED
