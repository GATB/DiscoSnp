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
//#include <libchash.h>
#include <xhash.h>
#include <commons.h>

#define debug_libchash_access

/*
 * prepare the key for the hash table:
 * we only use uppercase keys
 */



xhash xhash_create_seed_index(){
    xhash x;
    xh_init_uint64_t(&x, sizeof(uint64_t));
    return x;
}

//hash_t hash_create_binarykey(){
//    /**
//     *  arguments for AllocateHashTable():
//     *
//     *  cchKey: if it's a positive number, then each key is a
//     *                fixed-length record of that length.  If it's 0,
//     *                the key is assumed to be a \0-terminated string.
//     *        fSaveKey: normally, you are responsible for allocating
//     *                  space for the key.  If this is 1, we make a
//     *                  copy of the key for you.
//     */
//    return (hash_t)AllocateHashTable(sizeof(kmer_type),1);
//}



 


inline void get_offset_and_nb_from_sinfo(uint64_t  sinfo, uint64_t * offset_seed,uint64_t * nb_seeds, GlobalValues &gv )
{
    *nb_seeds = sinfo & gv.mask_nbseed;
    *offset_seed = (sinfo >> gv.nbits_nbseeds) & gv.mask_offset_seed ;
}

inline void set_offset_and_nb_into_sinfo(uint64_t * sinfo, uint64_t  offset_seed,uint64_t  nb_seeds, GlobalValues &gv  )
{
    uint64_t new_val = 0 ;
    new_val =  (offset_seed << gv.nbits_nbseeds ) + nb_seeds ;
    *sinfo = new_val;
}

void iterate_and_fill_offsets(xhash * map, GlobalValues &gv ){
    
    uint64_t offset_courant=0;
   // seedinfo sinfo;
    uint64_t sinfo;
    uint64_t offset_seed;
    uint64_t nb_seeds;
    //iterate over hash table
    
    uint64_t nbdiffseeds= 0;
//    HTItem *bck;
    
    xh_entry *it;
    for (it = xh_begin(map); it != NULL; it = xh_next(it)) {
        sinfo = xh_val(it, uint64_t);
        get_offset_and_nb_from_sinfo(sinfo, &offset_seed, &nb_seeds, gv);

        set_offset_and_nb_into_sinfo(&sinfo, offset_courant, 0, gv);
        offset_courant += nb_seeds ;

//        bck->data = sinfo;
        // replace old value for the key:
        xh_put_uint64_t(map, xh_key(it, uint64_t), &sinfo);
        
        nbdiffseeds ++;
    }
    
}

void hash_fill_kmer_index(xhash * map, const kmer_type * key, std::pair <uint64_t, int > * seed_table, const uint64_t fragment_id, const int position_on_fragment, GlobalValues &gv ){
//    HTItem *res;

    uint64_t sinfo;
    uint64_t offset_seed;
    uint64_t nb_seeds;
    uint64_t indexseed;
    
    std::pair <uint64_t, int > new_couple;
    new_couple.first = fragment_id;
    new_couple.second = position_on_fragment;
    
    xh_entry *e = xh_get_uint64_t(map, *key); // todo change from * key to & key in the declaration
    assert(e!=NULL); //key should be present
    sinfo = xh_val(e, uint64_t);

    get_offset_and_nb_from_sinfo(sinfo, &offset_seed, &nb_seeds, gv);
    indexseed = offset_seed +  nb_seeds ;
    


    seed_table[indexseed]=new_couple;
    
    Sinc24(nb_seeds) ;
    set_offset_and_nb_into_sinfo(&sinfo, offset_seed, nb_seeds,gv);
    
    xh_put_uint64_t(map, *key, &sinfo);
//    res->data = sinfo;
    
    
}

int get_seed_info(xhash * map, const kmer_type * key, uint64_t * offset_seed, uint64_t * nb_seeds, GlobalValues &gv ){
    
//    const HTItem *res;
    
    xh_entry *e;
    e = xh_get_uint64_t(map, *key);
//    if (e == NULL)
    
    
//    res=HashFind((struct HashTable*)map, *key);
    
    if(e!=NULL)
    {
        get_offset_and_nb_from_sinfo(xh_val(e, uint64_t), offset_seed, nb_seeds, gv);
        
        return 1;
    }
    return 0;
    
}



void hash_incr_kmer_count(xhash * map, const kmer_type * key, GlobalValues& gv){
//    HTItem *res;
    uint64_t offset_seed;
    uint64_t nb_seeds;
    uint64_t sinfo;
    

//    res=HashFindOrInsert( (struct HashTable*)map, *key,(ulong) 0); // find or insert 0
    xh_entry *e;
    e = xh_get_uint64_t(map, *key);
    if (e == NULL){
        int zero = 0;
        xh_put_uint64_t(map, *key, &zero);
    }
    e = xh_get_uint64_t(map, *key);
    
    
    sinfo = xh_val(e, uint64_t);
  
    get_offset_and_nb_from_sinfo(sinfo, &offset_seed, &nb_seeds,gv);


    Sinc24(nb_seeds) ; // saturated inc to 24 bit
    
    set_offset_and_nb_into_sinfo(&sinfo, offset_seed, nb_seeds,gv);
    xh_put_uint64_t(map, *key, &sinfo);
//    res->data = sinfo;

    
}


