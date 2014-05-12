#ifdef ED

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

#ifndef ASSERTS
#define NDEBUG // disable asserts, they're computationnally intensive
#endif

#include <stdio.h>
#include <assert.h>
#include <algorithm> // for min
#include "../minia/Kmer.h"
#include "Kmer_for_kissnp2.h"
#include "../minia/lut.h"

using namespace std;


kmer_type next_kmer_new(kmer_type graine, int added_nt, int strand)
{
    assert(added_nt<4);
    assert((strand == NULL) || (strand<2));

    kmer_type new_graine;
    kmer_type temp_graine;

    if (strand == 1)// the kmer we're extending is actually a revcomp sequence in the bidirected debruijn graph node
        temp_graine = revcomp(graine);
    else
        temp_graine = graine;

    new_graine = (((temp_graine) * 4 )  + added_nt) & kmerMask;
    //new_graine = (((graine) >> 2 )  + ( ((kmer_type)added_nt) << ((sizeKmer-1)*2)) ) & kmerMask; // previous kmer
    kmer_type revcomp_new_graine = revcomp(new_graine);

    return min(new_graine,revcomp_new_graine);
}

kmer_type next_kmer_new_norevcomp(kmer_type graine, int added_nt, int strand)
{
    assert(added_nt<4);
    assert((strand == NULL) || (strand<2));

    kmer_type new_graine;
    kmer_type temp_graine;

    if (strand == 1)// the kmer we're extending is actually a revcomp sequence in the bidirected debruijn graph node
        temp_graine = revcomp(graine);
    else
        temp_graine = graine;

    new_graine = (((temp_graine) * 4 )  + added_nt) & kmerMask;
    //new_graine = (((graine) >> 2 )  + ( ((kmer_type)added_nt) << ((sizeKmer-1)*2)) ) & kmerMask; // previous kmer
    kmer_type revcomp_new_graine = revcomp(new_graine);

    return new_graine;
}

#endif
