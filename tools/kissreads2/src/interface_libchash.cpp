/*****************************************************************************
 *   discoSnp++: discovering polymorphism from raw unassembled NGS reads
 *   A tool from the GATB (Genome Assembly Tool Box)
 *   Copyright (C) 2020  INRIA
 *   Authors: P.Peterlongo, E.Drezen
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

/*
 * interface between libchash and the generic hash functions used in kissreads2
 */

//#include <assert.h>
#include <libchash.h>
#include <hash.h>

#define debug_libchash_access

/*
 * prepare the key for the hash table:
 * we only use uppercase keys
 */




hash_t hash_create_binarykey(){
    /**
     *  arguments for AllocateHashTable():
     *
     *  cchKey: if it's a positive number, then each key is a
     *                fixed-length record of that length.  If it's 0,
     *                the key is assumed to be a \0-terminated string.
     *        fSaveKey: normally, you are responsible for allocating
     *                  space for the key.  If this is 1, we make a
     *                  copy of the key for you.
     */
    return (hash_t)AllocateHashTable(sizeof(kmer_type),1);
}



 


inline void get_offset_and_nb_from_sinfo(hash_val  sinfo, uint64_t * offset_seed,uint64_t * nb_seeds, GlobalValues &gv )
{
    *nb_seeds = sinfo & gv.mask_nbseed;
    *offset_seed = (sinfo >> gv.nbits_nbseeds) & gv.mask_offset_seed ;
}

inline void set_offset_and_nb_into_sinfo(hash_val * sinfo, uint64_t  offset_seed,uint64_t  nb_seeds, GlobalValues &gv  )
{
    hash_val new_val = 0 ;
    new_val =  (offset_seed << gv.nbits_nbseeds ) + nb_seeds ;
    *sinfo = new_val;
}

void iterate_and_fill_offsets(hash_t map, GlobalValues &gv ){
    
    uint64_t offset_courant=0;
   // seedinfo sinfo;
    hash_val sinfo;
    uint64_t offset_seed;
    uint64_t nb_seeds;
    //iterate over hash table
    
    uint64_t nbdiffseeds= 0;
    HTItem *bck;
    
    for ( bck = HashFirstBucket((struct HashTable*)map); bck; bck = HashNextBucket((struct HashTable*)map) )
    {
        sinfo = (hash_val) (bck->data);
        get_offset_and_nb_from_sinfo(sinfo, &offset_seed, &nb_seeds, gv);

        set_offset_and_nb_into_sinfo(&sinfo, offset_courant, 0, gv);
        offset_courant += nb_seeds ;

        bck->data = sinfo;
        
        nbdiffseeds ++;
    }
    
}

void hash_fill_kmer_index(hash_t map, const kmer_type * key, std::pair <uint64_t, int > * seed_table, const uint64_t fragment_id, const int position_on_fragment, GlobalValues &gv ){
    HTItem *res;

    hash_val sinfo;
    uint64_t offset_seed;
    uint64_t nb_seeds;
    uint64_t indexseed;
    
    std::pair <uint64_t, int > new_couple;
    new_couple.first = fragment_id;
    new_couple.second = position_on_fragment;
    
    
    res=HashFind((struct HashTable*)map, *key);

    assert(res!=NULL); //key should be present
    
    sinfo = (hash_val) (res->data);
    get_offset_and_nb_from_sinfo(sinfo, &offset_seed, &nb_seeds, gv);
    indexseed = offset_seed +  nb_seeds ;
    


    seed_table[indexseed]=new_couple;
    
    Sinc24(nb_seeds) ;
    set_offset_and_nb_into_sinfo(&sinfo, offset_seed, nb_seeds,gv);
    
    res->data = sinfo;
    
    
}

int get_seed_info(hash_t map, const kmer_type * key, uint64_t * offset_seed, uint64_t * nb_seeds, GlobalValues &gv ){
    
    const HTItem *res;
//    hash_val sinfo;
    
    
    res=HashFind((struct HashTable*)map, *key);
    
    if(res!=NULL)
    {
//        sinfo = (hash_val) (res->data);
        get_offset_and_nb_from_sinfo((hash_val) (res->data), offset_seed, nb_seeds, gv);
        
        return 1;
    }
    return 0;
    
}



void hash_incr_kmer_count(hash_t map, const kmer_type * key, GlobalValues& gv){
    HTItem *res;
    uint64_t offset_seed;
    uint64_t nb_seeds;
    hash_val sinfo;
    

    res=HashFindOrInsert( (struct HashTable*)map, *key,(ulong) 0); // find or insert 0
    
    sinfo = (hash_val) (res->data);
    get_offset_and_nb_from_sinfo(sinfo, &offset_seed, &nb_seeds,gv);


    Sinc24(nb_seeds) ; // saturated inc to 24 bit
    
    set_offset_and_nb_into_sinfo(&sinfo, offset_seed, nb_seeds,gv);
    res->data = sinfo;

    
}


