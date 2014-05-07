/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
 *   Authors: P.Peterlongo, R.Uricary, R.Chikhi, G.Rizk, E.Drezen
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

#ifndef _BUBBLE_HPP_
#define _BUBBLE_HPP_

/********************************************************************************/

#include <cstdlib>
#include <cmath>
#include <gatb/gatb_core.hpp>

/********************************************************************************/
template<size_t span=KSIZE_1>
class BubbleFinder : public Algorithm
{
public:

    /** Shortcuts. */
    typedef typename Kmer<span>::Model Model;
    typedef typename Kmer<span>::Type  kmer_type;

    /** Constructor. */
    BubbleFinder (
        const Graph& graph,
        const char* SNP_file_name, int low, int authorised_branching,
        bool extend_snps, int min_size_extension, bool print_extensions, bool strict_extension
    );

    /** Destructor. */
    ~BubbleFinder ();

    /** */
    void execute ();


protected:

    /** */
    void start (const Node& node);

    /** */
    void expand(int direction, char* path1, char* path2, kmer_type kmer1, kmer_type kmer2, int pos, char* p1, char* p2, char* p3, char* p4);

    /** */
    bool next(kmer_type *kmer);

    /** */
    void read_bubble_file(char *bubble_filename);

    const Graph& graph;
    Model        modelMin;
    Model        modelDirect;

    size_t sizeKmer;
    int threshold;

    int low;
    int authorised_branching;

//    bool authorised_branching;

    size_t nb_bubbles;
    size_t nb_bubbles_high;
    size_t nb_bubbles_low;

    kmer_type checksumStart1;
    kmer_type checksumStart2;
    kmer_type checksumStart3;

    kmer_type checksumExpand1;
    kmer_type checksumExpand2;

    kmer_type checksumGraph;
    size_t    nbKmer2;

    size_t checksumPrint1;

    size_t SKIP1;
    size_t SKIP2;
    size_t SKIP3;

    bool extend_snps;
    bool print_extensions;
    bool strict_extension; // true: strict extension, fales: contig extension
    int  min_size_extension;
    
    FILE* SNP_file;

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

    bool two_possible_extensions_on_one_path(kmer_type) const;
    bool two_possible_extensions(kmer_type kmer1, kmer_type kmer2) const;


    int close_snp(const char * path1, const char * path2, char * path1_c, char * path2_c); //returns -1 not closed, 0 no unique extension, 1 only left, 2 only right, 3 both

    void prints_contig_informations (
        FILE* file,
        std::pair<char*,int> left_extension,
        std::pair<char*,int> right_extension,
        bool extend_snps,
        bool strict_extension,
        int where_to_extend);

    bool checkBranching (kmer_type kmer1, kmer_type kmer2) const;
};

/********************************************************************************/

#endif // _BUBBLE_HPP_
