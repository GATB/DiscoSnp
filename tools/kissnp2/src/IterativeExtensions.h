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

#ifndef ITERATIVE_EXTENSION_H
#define ITERATIVE_EXTENSION_H
#include <assert.h>
#include <utility>      // std::pair, std::get
#ifdef MINIA_IS_IN_PARENT_FOLDER

#include "../minia/Utils.h"
#include "../minia/Terminator.h"
#include "../minia/Traversal.h"
#include "../minia/Kmer.h"
#include "../minia/Bank.h"
#include "../minia/Debloom.h"

#else

#include "minia/Utils.h"
#include "minia/Terminator.h"
#include "minia/Traversal.h"
#include "minia/Kmer.h"
#include "minia/Bank.h"
#include "minia/Debloom.h"

#endif

extern int64_t genome_size;
extern Bloom * bloo1;
extern BinaryBank * SolidKmers;
extern BranchingTerminator * terminator;
//Raluca
extern AssocPairedSet *pairedBranchingKmers ;

typedef enum {
    LARGEUR,
    PROFONDEUR,
    NB_PARCOURS
} parcours_t;


class IterativeExtensions {
    public:
        
    IterativeExtensions();
    static void construct_linear_seqs(string L,string R, int max_depth, string output_file, int verb, int  max_nodes, parcours_t search_mode = PROFONDEUR,  bool swf = false );
    static void construct_linear_seqs(string L, int max_depth, string output_file);
    static void construct_linear_seqs_paired(string L,string R, int max_depth, string output_file , int verb, int max_nodes);
    static void construct_linear_seqs_paired(string L, int max_depth, string output_file);



    // controls behavior of construct_linear_seqs
    enum Traversal_type { SimplePaths, Monument };
    enum When_to_stop_extending { After_first_contig, Until_max_depth };
    static Traversal_type traversal_type; // initialized in .cpp
    static When_to_stop_extending when_to_stop_extending;    
    static bool dont_output_first_nucleotide;
    
    // extra function added by pierre
    static pair<char *, int>  extend_a_snp(string L, int max_depth);
};

#endif
