//
//  Kissreads2.h
//  discoSnp_GATB
//
//  Created by Pierre Peterlongo on 03/07/15.
//  Copyright (c) 2015 Pierre Peterlongo. All rights reserved.
//

#ifndef __discoSnp_GATB__Kissreads2__
#define __discoSnp_GATB__Kissreads2__

#include <iostream>
// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <commons.h>
#include <hash.h>
#include <outputs.h>
#include <fragment_info.h>
#include <fragment_index.h>
#include <extension_algorithm.h>


/** We define string constants for command line options. */
#define STR_KISSREADS_MAX_HAMMING          "-hamming"
#define STR_KISSREADS_GENOTYPE             "-genotype"
#define STR_KISSREADS_OUTPUT_FASTA         "-output_fasta"
#define STR_KISSREADS_SIZE_SEEDS           "-size_seeds"
#define STR_KISSREADS_INDEX_STRIDE         "-index_stride"
#define STR_KISSREADS_SIZE_K               "-k"
#define STR_KISSREADS_COVERAGE_FILE_NAME   "-coverage_file"
//#define STR_KISSREADS_MIN_COVERAGE         "-abundance-min"
#define STR_URI_OUTPUT_COHERENT            "-co"
#define STR_URI_OUTPUT_UNCOHERENT          "-unco"
#define STR_URI_READS_INPUT                "-reads"
#define STR_URI_PREDICTION_INPUT           "-predictions"
#define STR_RADSEQ                         "-x"
#define STR_PHASING                        "-phasing"

/** \brief Tool class that looks for SNP
 *
 * The Kissnp2 is the front class for SNP detection in a provided de Bruijn graph.
 * The output is a bank with pairs of sequences defining a bubble.
 */
class Kissreads2 : public Tool
{
public:
    
    /** Constructor. */
    Kissreads2 ();
    
    /** Implementation of Tool::execute method. */
    void execute ();
    
};

#endif /* defined(__discoSnp_GATB__Kissreads2__) */
