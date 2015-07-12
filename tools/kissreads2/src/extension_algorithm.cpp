/*****************************************************************************
 *   discoSnp++: discovering polymorphism from raw unassembled NGS reads
 *   A tool from the GATB (Genome Assembly Tool Box)
 *   Copyright (C) 2014  INRIA
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
 * extension_algorithm.c
 *
 *  Created on: 16 sept. 2010
 *      Author: ppeterlo
 */

#include <extension_algorithm.h>


//#define DEBUG_MAPPING
//#define DEBUG_QUALITY
#define min(a, b) ((a) < (b) ? (a) : (b))



void feed_coherent_positions(vector<FragmentInfo*> starters, const int starter_id, const int start, const int length_read, string quality, int start_on_read, int read_set_id, GlobalValues& gv){
    
    
    int i, start_on_starter, stop_on_starter;
	if(start<0) start_on_starter=0;
    else start_on_starter=start;
    
    FragmentInfo* the_starter=starters[starter_id];
    FragmentInfo* the_reference_starter = starters[2*(starter_id/2)]; // In case of snps, only the upper path starter contains informations such as the positions of the SNPs. This is the reference?
	
    if(start+length_read<the_starter->upperCaseSequence.size()) stop_on_starter=start+length_read;
    else stop_on_starter=the_starter->upperCaseSequence.size();
    
    
	the_starter->number_mapped_reads[read_set_id]++;

    if ( quality.length()>0 ){
        if (the_reference_starter->nbOfSnps>0) {
            //            printf("there are %d SNPs\n", the_reference_starter->nbOfSnps);
            int snp_id;
            for(snp_id=0;snp_id<the_reference_starter->nbOfSnps;snp_id++){ // we only add the qualities of the mapped SNPs
                //                printf("the_reference_starter->SNP_positions[%d]=%d \n", snp_id, the_reference_starter->SNP_positions[snp_id]); //DEB
                i=the_reference_starter->SNP_positions[snp_id];
                //                printf(" %d \n", start_on_read + i - start_on_starter); //DEB
                the_starter->sum_qualities[read_set_id] += (unsigned int) quality[start_on_read + i - start_on_starter];
                the_starter->nb_mapped_qualities[read_set_id] += 1;
            }
        }
        else{ // we sum all qualities and divide by the number of positions
            int sum_temp=0;
            int denom=0;
            for(i=start_on_starter;i<stop_on_starter;i++) { // to avoid to increase too much the the_starter->sum_qualities array, we add the average quality of the whole read.
                denom+=1;
                sum_temp+=(unsigned int) quality[start_on_read + i - start_on_starter];
            }
            if(denom>0){
                the_starter->sum_qualities[read_set_id] += (unsigned int)sum_temp/denom;
                the_starter->nb_mapped_qualities[read_set_id] += 1;
            }
        }
    }
    
    
    
#ifdef KMER_SPANNING
    // the position i is contained into a kmer fully contained into only 1 mapped read, return 1
    // for doing this we stored on each position of the fragment the number of k-mers starting at this position that fully belong to a read that was mapped
	
    //  -------------------------------------------------------------- starter
    //            ************************
    //                       <-----k----->
    //  00000000001111111111110000000000000000000000000000000000000000 the_starter->local_coverage[read_file_id]
	if(start+length_read-gv.minimal_read_overlap<the_starter->upperCaseSequence.size()) stop_on_starter=start+length_read-gv.minimal_read_overlap;
    else stop_on_starter=the_starter->upperCaseSequence.size();
#endif
    
	for(i=start_on_starter;i<stop_on_starter;i++) Sinc8(the_starter->local_coverage[read_set_id][i]);
	
    
    
}




/**
 *  | pwi (may be negative)
 *     --------------  fragment
 *  *************      read
 *
 *  Tests if the overlapping part between read and fragment do not have more than <code>subst_allowed</code> substitions
 * In case of SNPs, we need to avoid any substitution on the central fragment position (the one containing the SNP)
 * Thus in this function, we return 0 if any substitution occurs on this central position, whatever the number of substitution_seen
 *  returns 1 if true between read and fragment, 0 else
 */
char constrained_read_coherent(const int pwi, const char * fragment, const char * read, const int subst_allowed, const char * SNP_positions){
	int substitution_seen=0; // number of seen substitutions for now
    
	int pos_on_read, pos_on_fragment; // where to start
    
    //    print_mapping(pwi,fragment,read); //DEB
    
	/*
	 *  | pwi (negative)
	 *     --------------  fragment
	 *  *************      read
	 *     | we start here
	 */
	if(pwi<0) {
		pos_on_fragment=0;
		pos_on_read=-pwi;
	}
    
	/*
	 *        | pwi (positive)
	 *     --------------  fragment
	 *        *************      read
	 *        | we start here
	 */
	else{
		pos_on_fragment=pwi;
		pos_on_read=0;
	}
    
    //const int snp_pos = strlen(fragment)/2;
    char snp_pos = SNP_positions[0];
    
    int id_array_SNP_position=0;
    
	// walk the read and the fragment together, detecting substitutions.
	// stop if the number of substitution is too high
	while(fragment[pos_on_fragment]!='\0' && read[pos_on_read]!='\0'){
        
        if (pos_on_fragment>snp_pos)
        {
            id_array_SNP_position++;
            snp_pos = SNP_positions[id_array_SNP_position];
        }
		if (fragment[pos_on_fragment]!=read[pos_on_read] &&
            fragment[pos_on_fragment]!='*' &&
            fragment[pos_on_fragment]!='?' &&
            fragment[pos_on_fragment]!='N'){ // one subsitution
			substitution_seen++;
			if(substitution_seen>subst_allowed) return 0; // too much subsitutions
            if(pos_on_fragment==snp_pos) {
                return 0; // substition should not be on the snp
            }
		}
		pos_on_fragment++;
		pos_on_read++;
	}
	return 1;
}


// For each position in the starter:
// find the kmer supported by the smallest number of reads, e.g.:
//                                i
//  ------------------------------X------------------------------- starter
// a         ************************
// b         **************************
// c                    *******************
// d                       *************************
//                     <-----k-----> : 2 reads (a and b)
//                      <-----k-----> : 3 reads (a,b and b)
//                       <-----k-----> : 2 reads (b and c)
//                        <-----k-----> : 1 reads (c)
//                         <-----k-----> : 2 reads (c and d)
//                            ...
// the position i is contained into a kmer fully contained into only 1 mapped read, return 1
// for doing this we stored on each position of the fragment the number of k-mers starting at this position that fully belong to a read that was mapped.
int minimal_kmer_coverage(FragmentInfo the_starter, int read_file_id, GlobalValues& gv){
    
    int i, val_min=INT_MAX;
    const  int stopi=the_starter.upperCaseSequence.size();
    for(i=0;i<stopi;i++){ // for each position on the read
        val_min=min(val_min, the_starter.local_coverage[read_file_id][i]);
    }
    return val_min;
}


/**
 * Performs the first extension of the algorithm:
 * For each starter:
 *  - verify which reads maps on it with at most subst_allowed substitutions
 *  - for fragments that are fully covered by reads (each position is covered at least once):
 *    1) add their id in a list (that is returned by the function)
 *    2) detect the extending reads (at least "min_extension_coverage_depth" reads) that enable to extend right the starter
 *    3) store this extending reads in the structure of the fragments (fragment_info)
 
 * returns the average read length
 */
void ReadMapper::map_all_reads_from_a_file (
                                GlobalValues & gv,
                                FragmentIndex& index
                                ){

    
    
    //////////////////////////////////////////////////////////////////////////
	/////////////// read all reads - storing those coherent with reads ///////
	//////////////////////////////////////////////////////////////////////////
	char * starter;
    uint64_t offset_seed;
    uint64_t nb_occurrences;
    
	// working variables
	int read_len, i, ii,  pwi, stop, read_coherence;
    long int  read_number=0;
    
	// map of starter -> position (for each read and direction, stores the starter and position already tested.)
	// book space for the read that going to be read.
    
	   
    listint * tested_starters = listint_create(); // given a read, stores all starter ids on which a seed was seen. Enables to free quickly list of tested pwis of each tested starters.
	
    
	// read all the reads
    // We create an iterator over this bank.
    Iterator<Sequence>* it = inputBank->iterator();
    LOCAL (it);
    // We loop over sequences.
    for (it->first(); !it->isDone(); it->next())
    {
        read_number++;
        // Shortcut
        Sequence& seq = it->item();
        char *read = strdup(seq.toString().c_str());
        read_len = strlen(read);
        char * quality = strdup(seq.getQuality().c_str());
        
        const int minimal_pwi = gv.minimal_read_overlap - seq.getDataSize();
        
        // The read must overlap the fragment with at least minimal_read_overlap positions.
        // here is the first position on which the read may map :
        //        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  starter
        //     **********       read (10)
        //        <-----> minimal_read_overlap (7)
        //     <-> -pwi (-3)
        // -pwi+minimal_read_overlap <= |read|
        // -pwi <= |read|-minimal_read_overlap
        // pwi >= minimal_read_overlap-|read|
        // pwi >= 7-10 = -3
        // minimal_pwi = minimal_read_overlap-|read|
        
        
		stop = read_len-gv.size_seeds+1;
		// read all seeds present on the read:
        int direction;
        kmer_type coded_seed;
        char toinit=1;
        //        char validSeed;
        for(direction=0;direction<2;direction++){ // try the two possible directions of the read
            toinit=1; // we have to init a new seed
            for (i=0;i<stop;i++){ // for all possible seed on the read
                if(toinit) {
                    coded_seed=gv.codeSeed(read+i); // init the seed
                    toinit=0;
                }
                else { // previous seed was correct, we extend it.
                    coded_seed=gv.updateCodeSeed(read+i,&coded_seed); // utpdate the previous seed with a
                }
                if(get_seed_info(index.seeds_count,&coded_seed,&offset_seed,&nb_occurrences,gv)){
                    // for each occurrence of this seed on the starter:
                    for (ii=offset_seed; ii<offset_seed+nb_occurrences; ii++) {
                        couple * value = &(index.seed_table[ii]);
                        const char * starter = index.all_predictions[value->a]->upperCaseSequence.c_str();
                        
#ifdef DEBUG_MAPPING
                        cout<<"seed = "<<read+i<<"in "<<starter<<" pos "<<value->b<<starter+value->b<<endl;//DEB
#endif
                        

                        if(index.all_predictions[value->a]->mapped_with_current_read[read_set_id]) {
                            continue; // a mach was already found with this starter, we consider only one match per couple (starter/read)
                        }
                        pwi = value->b-i; // starting position of the read on the starter.
                        if (listint_contains(index.all_predictions[value->a]->tested_pwis_with_current_read[read_set_id],pwi)) {
                            continue; // this reads was already (unsuccessfuly) tested with this starter at this position. No need to try it again.
                        }
                        if (numberInListint (index.all_predictions[value->a]->tested_pwis_with_current_read[read_set_id])==0){ // this is the first time we meet this starter with this read, we add it in the list of met starters (in order to free the pwis lists of each encounted starters for this read)
                            listint_add(tested_starters,value->a);
                        }
                        listint_add(index.all_predictions[value->a]->tested_pwis_with_current_read[read_set_id],pwi); // We store the fact that this read was already tested at this position on this starter.
                        
                        
                        
                        // overview general situation:
                        
                        //        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  starter
                        //        <---------> b
                        //                   [--------]                     seed
                        //             ******************************       read
                        //             <----> i
                        //        <---> pwi
                        
                        const int maximal_pwi = strlen(starter)-gv.minimal_read_overlap;
                        
//                        cout<<minimal_pwi<<"<"<<pwi<<"<"<<maximal_pwi<<"???"<< endl; //DEB
                        
                        if (pwi<minimal_pwi) {
                            continue; // this read to not overlap enough with the starter.
                        }
                        if (pwi > maximal_pwi) {
                            continue; // this read to not overlap enough with the starter.
                        }
                        
                        //        ;;;;;;;;;;;  starter (11)
                        //             ******************************       read
                        //             <----> minimal_read_overlap (6)
                        //        <---> pwi (5)
                        // |starter| <= pwi+minimal_read_overlap
                        
                        read_coherence=0;
                        read_coherence = constrained_read_coherent(pwi, starter, read, gv.subst_allowed, index.all_predictions[value->a-value->a%2]->SNP_positions);

                        if(read_coherence == 1){ // tuple read starter position is read coherent
                            
                            index.all_predictions[value->a]->mapped_with_current_read[read_set_id]=true;
                            
                            
#ifdef DEBUG_MAPPING
                            printf("SUCCESS %d %d \n", pwi, value->a);
                            cout<<pwi<<" "<<index.all_predictions[value->a]->upperCaseSequence<<" "<<read<<endl; //DEB
#endif
                            feed_coherent_positions(index.all_predictions, value->a , pwi, (int)strlen(read), quality, i, read_set_id, gv);
                            
                        } // end tuple read starter position is read coherent
                    }
                } // end all infos for the current seed
            } // end all seeds of the read
            gv.revcomp(read);
            
            gv.rev (quality);
            
            // free the lists of tested positions:
            cellint * tested_starter=tested_starters->first;
            while (tested_starter != NULL)
            {
                listint_empty(index.all_predictions[tested_starter->val]->tested_pwis_with_current_read[read_set_id]);
                index.all_predictions[tested_starter->val]->mapped_with_current_read[read_set_id] = false;
                //                printf("freeing starter %d\n", *(int *)tested_starter->val); //DEB
                tested_starter=tested_starter->prox;
            }
            listint_empty(tested_starters);
            
            
        } // end both directions
        
        free(read);
        free(quality);
	}// end all reads of the file
    
    if(!gv.silent) printf("\r %ld reads treated\n", read_number);  // print progress
    
  
    listint_free(tested_starters);
    
    

}





void ReadMapper::set_read_coherency(GlobalValues& gv, FragmentIndex index){
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////// for each starter: check those fully coherent and store left and right reads covering them ///////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    int starter_id;
	for (starter_id=0;starter_id < index.all_predictions.size();starter_id++){
        index.all_predictions[starter_id]->set_read_coherent(read_set_id,gv);
	} // end all fragments
	
	
    
    
}



    





