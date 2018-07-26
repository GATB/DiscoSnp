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


// Replaces SNP_higher_path_3780|P_1:30_A/G|high|nb_pol_1 by 3780h
string parse_variant_id(string in){
    string res="";
    // push the values while they are in [0-9]
    for (char c : in){
        if (c>=int('0') && c<=int('9'))
            res=res+c;
        if (c=='|')
            break;
    }
    
    res=res+in[4]; // h or l.
    return res;
}


//feed_coherent_positions(index.all_predictions, value->a , pwi, (int)strlen(read), quality, seed_position, read_set_id, gv);

void feed_coherent_positions(vector<FragmentInfo*> & predictions, const int prediction_id, const int pwi, const int length_read, string quality, int read_set_id, GlobalValues& gv){
    int start_on_prediction, stop_on_prediction;
    int start_on_read;
    /*
     *  | pwi (negative)
     *     --------------  fragment
     *  *************      read
     *     | we start here
     */
    if(pwi<0) {
        start_on_prediction=0;
        start_on_read=-pwi;
    }
    
    /*
     *        | pwi (positive)
     *     --------------  fragment
     *        *************      read
     *        | we start here
     */
    else{
        start_on_prediction=pwi;
        start_on_read=0;
    }
    
    int i;
    
    
    FragmentInfo* the_prediction=predictions[prediction_id];
    FragmentInfo* the_reference_prediction = predictions[2*(prediction_id/2)]; // In case of snps, only the upper path prediction contains informations such as the positions of the SNPs. This is the reference?
    
    if(pwi+length_read<the_prediction->upperCaseSequence.size()) stop_on_prediction=pwi+length_read;
    else stop_on_prediction=the_prediction->upperCaseSequence.size();
    
    
    __sync_fetch_and_add ( & the_prediction->number_mapped_reads[read_set_id],1);
    
    if ( quality.length()>0 ){
        if (the_reference_prediction->nbOfSnps>0) {
            int snp_id;
            for(snp_id=0;snp_id<the_reference_prediction->nbOfSnps;snp_id++){ // we only add the qualities of the mapped SNPs
                i=the_reference_prediction->SNP_positions[snp_id];
                if (start_on_read + i - start_on_prediction>=0 && start_on_read + i - start_on_prediction < length_read){
                    the_prediction->sum_qualities[read_set_id] += (unsigned int) quality[start_on_read + i - start_on_prediction];
                    the_prediction->nb_mapped_qualities[read_set_id] += 1;
                }
            }
        }
        else{ // we sum all qualities and divide by the number of positions
            int sum_temp=0;
            int denom=0;
            for(i=start_on_prediction;i<stop_on_prediction;i++) { // to avoid to increase too much the the_prediction->sum_qualities array, we add the average quality of the whole read.
                denom+=1;
                sum_temp+=(unsigned int) quality[start_on_read + i - start_on_prediction];
            }
            if(denom>0){
                the_prediction->sum_qualities[read_set_id] += (unsigned int)sum_temp/denom;
                the_prediction->nb_mapped_qualities[read_set_id] += 1;
            }
        }
    }
    
    
    
#ifdef KMER_SPANNING
    // the position i is contained into a kmer fully contained into only 1 mapped read, return 1
    // for doing this we stored on each position of the fragment the number of k-mers starting at this position that fully belong to a read that was mapped
    
    //  -------------------------------------------------------------- prediction
    //            ************************
    //                       <-----k----->
    //  00000000001111111111110000000000000000000000000000000000000000 the_prediction->local_coverage
    if(pwi+length_read-gv.minimal_read_overlap<the_prediction->upperCaseSequence.size()) stop_on_prediction=pwi+length_read-gv.minimal_read_overlap;
    else stop_on_prediction=the_prediction->upperCaseSequence.size();
#endif
    
    for(i=start_on_prediction;i<stop_on_prediction;i++) Sinc8(the_prediction->local_coverage[i]);
    
    
    
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
bool constrained_read_mappable(const int pwi, const char * fragment, const char * read, const int subst_allowed, const char * SNP_positions, const int seed_position_on_read, const int size_seed){
    int substitution_seen=0; // number of seen substitutions for now
    int pos_on_read, pos_on_fragment; // where to start
    
    //       print_mapping(pwi,fragment,read); //DEB
    
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
    
    int snp_pos = SNP_positions[0];
    
    int id_array_SNP_position=0;
    
    
    // walk the read and the fragment together, detecting substitutions.
    // stop if the number of substitution is too high
    while(fragment[pos_on_fragment]!='\0' && read[pos_on_read]!='\0'){
        // we know that the seed has a perfect match. We can skip this positions.
        // TODO: test this latter, i've found a valgrind error (28 oct 2015)
        //        if (pos_on_read==seed_position_on_read){
        //            pos_on_fragment+=size_seed;
        //            pos_on_read+=size_seed;
        //        }
        if (pos_on_fragment>snp_pos)
        {
            id_array_SNP_position++;
            snp_pos = SNP_positions[id_array_SNP_position];
        }
        if (fragment[pos_on_fragment]!=toupper(read[pos_on_read]) &&
            fragment[pos_on_fragment]!='*' &&
            fragment[pos_on_fragment]!='?' &&
            fragment[pos_on_fragment]!='N'){ // one subsitution
            substitution_seen++;
            if(substitution_seen>subst_allowed) return false; // too much subsitutions
            if(pos_on_fragment==snp_pos) {
                return false; // substition should not be on the snp
            }
        }
        pos_on_fragment++;
        pos_on_read++;
    }
    return true;
}


// For each position in the prediction:
// find the kmer supported by the smallest number of reads, e.g.:
//                                i
//  ------------------------------X------------------------------- prediction
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
//int minimal_kmer_coverage(FragmentInfo the_prediction, int read_file_id, GlobalValues& gv){

//    int i, val_min=INT_MAX;
//    const  int stopi=the_prediction.upperCaseSequence.size();
//    for(i=0;i<stopi;i++){ // for each position on the read
//        val_min=min(val_min, the_prediction.local_coverage[read_file_id][i]);
//    }
//    return val_min;
//}


// We define a functor that will be cloned by the dispatcher
struct Functor
{
    //    ISynchronizer* synchro;    fstream& file;
    
    map<u_int64_t, set<u_int64_t> >  tested_prediction_and_pwis;          // stores for this read, the pwi positions tested for each prediction.
    set<u_int64_t> mapped_prediction_as_set;                              // stores for this read, the succesfully mapped predictions - enables a quick existance testing of an element
    //    list<u_int64_t> mapped_prediction_as_list;                            // stores for this read, the succesfully mapped predictions - conserve the predictions order.
    
    
    
    GlobalValues & gv;
    FragmentIndex& index;
    const int read_set_id;
    u_int64_t * number_of_mapped_reads;
    map<string,int> & phased_variants;
    
    
    
    
    Functor (GlobalValues & gv, FragmentIndex& index, const int read_set_id, u_int64_t * number_of_mapped_reads, map<string,int> & phased_variants) : gv(gv), index(index), read_set_id(read_set_id), number_of_mapped_reads(number_of_mapped_reads), phased_variants(phased_variants){}
    
    
    map<int,int64_t> core_mapping(char *read, char * quality){
        map<int,int64_t> pwi_and_mapped_predictions;                        // stores for this reads the succesfully mapped predictions together with their pwi.
        
        const uint64_t read_len = strlen(read);
        
        const int minimal_pwi = gv.minimal_read_overlap - read_len;//seq.getDataSize();
        uint64_t offset_seed;
        uint64_t nb_occurrences;
        
        
        
        // The read must overlap the fragment with at least minimal_read_overlap positions.
        // here is the first position on which the read may map :
        //        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  prediction
        //     **********       read (10)
        //        <-----> minimal_read_overlap (7)
        //     <-> -pwi (-3)
        // -pwi+minimal_read_overlap <= |read|
        // -pwi <= |read|-minimal_read_overlap
        // pwi >= minimal_read_overlap-|read|
        // pwi >= 7-10 = -3
        // minimal_pwi = minimal_read_overlap-|read|
        const int stop = read_len-gv.size_seeds+1;
        kmer_type coded_seed;
        
        // for both dirrections of the read
        for(int direction=0;direction<2;direction++) // try the two possible directions of the read
        {
            // read all seeds present on the read:
            for (int seed_position=0;seed_position<stop;seed_position++){ // for all possible seed on the read
                if(seed_position==0) {
                    coded_seed=gv.codeSeed(read+seed_position); // init the seed
                }
                else { // previous seed was correct, we extend it.
                    coded_seed=gv.updateCodeSeed(read+seed_position,&coded_seed); // utpdate the previous seed
                }
                
                if(get_seed_info(index.seeds_count,&coded_seed,&offset_seed,&nb_occurrences,gv)){
                    // for each occurrence of this seed on the prediction:
                    for (int occurrence_id=offset_seed; occurrence_id<offset_seed+nb_occurrences; occurrence_id++) {
                        couple * value = &(index.seed_table[occurrence_id]);
                        if (mapped_prediction_as_set.count(value->a)!=0) {
                            continue; // This prediction was already mapped with this read.
                        }
                        
                        
                        
                        
                        // shortcut
                        set<u_int64_t> & tested_positions = tested_prediction_and_pwis[value->a];
                        
                        // get the corresponding prediction sequence
                        const char * prediction = index.all_predictions[value->a]->upperCaseSequence.c_str();
                        
#ifdef DEBUG_MAPPING
                        cout<<"seed = "<<read+seed_position<<"in "<<prediction<<" pos "<<value->b<<prediction+value->b<<endl;//DEB
#endif
                        
                        
                        
                        const int pwi = value->b-seed_position; // starting position of the read on the prediction.
                        if (tested_positions.count(pwi) != 0) continue; // this reads was already tested with this prediction at this position. No need to try it again.
                        tested_positions.insert(pwi); // We store the fact that this read was already tested at this position on this prediction.
                        
                        
                        
                        
                        
                        
                        // overview general situation:
                        
                        //        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  prediction
                        //        <---------> b
                        //                   [--------]                     seed
                        //             ******************************       read
                        //             <----> i
                        //        <---> pwi
                        
                        const int maximal_pwi = strlen(prediction)-gv.minimal_read_overlap;
                        
                        
                        if (pwi<minimal_pwi) {
                            continue; // this read to not overlap enough with the prediction.
                        }
                        if (pwi > maximal_pwi) {
                            continue; // this read to not overlap enough with the prediction.
                        }
                        //        ;;;;;;;;;;;  prediction (11)
                        //             ******************************       read
                        //             <----> minimal_read_overlap (6)
                        //        <---> pwi (5)
                        // |prediction| <= pwi+minimal_read_overlap
                        
                        
                        
                        const bool is_read_mapped = constrained_read_mappable(pwi, prediction, read, gv.subst_allowed, index.all_predictions[value->a-value->a%2]->SNP_positions, seed_position, gv.size_seeds);
                        
#ifdef DEBUG_MAPPING
                        if (is_read_mapped) cout<<endl<<read<<" mapped on "<<prediction<<" "<<value->a<<" pos "<<pwi<<endl;
#endif
                        
                        if(is_read_mapped){ // tuple read prediction position is read coherent
                            __sync_fetch_and_add (number_of_mapped_reads, 1);
                            
//    #ifdef PHASING
                            if(gv.phasing){
                            mapped_prediction_as_set.insert     (value->a);     // This prediction whould not be mapped again with the same read
                            // currently the phasing works better with SNPs, as boths paths of  an indel may be mapped by a same read
                            if (index.all_predictions[value->a]->nbOfSnps !=0){      // If this is not an indel (todo phase also indels)
                                //                            mapped_prediction_as_list.push_back (value->a);     // Store the prediction in the list (remaining the order)
                                int sign=direction==0?1:-1;
                                
                                if (direction == 0){
                                    if (pwi_and_mapped_predictions.find(pwi) == pwi_and_mapped_predictions.end())  pwi_and_mapped_predictions[pwi] = sign*value->a;
                                    // TODO what if this read maps already a variant at the same position ?
                                }
                                else{
                                    
                                    //        ;;;;;;;;;;;  prediction (11)
                                    //             ******************************       read (30)
                                    //             <----> minimal_read_overlap (6)
                                    //        <---> pwi (5)
                                    //                   <----------------------> rc_pwi (-24)
                                    // we have : |read_overlap| = |prediction|-pwi
                                    // we have : rc_pwi = |read_overlap|-|read|
                                    // we have : rc_pwi = |prediction|-pwi-read
                                    // rc_pwi = 11-5-30 = -24.
                                    
                                    // Validation with pwi<0:
                                    /*
                                     *  <-> pwi (-3)
                                     *     --------------  prediction (14)
                                     *  ***********        read (11)
                                     *             <----> rc_pwi=6
                                     * |prediction|-pwi-read = 14-(-3)-11 = 6 (CQFD :))
                                     */
                                    const int rc_pwi = strlen(prediction) - pwi - read_len;
                                    if (pwi_and_mapped_predictions.find(rc_pwi) == pwi_and_mapped_predictions.end())  pwi_and_mapped_predictions[rc_pwi] = sign*value->a;
                                    // TODO what if this read maps already a variant at the same position ?
                                    ///
                                }
                            }
                            }
                            ////// END PHASING
//#endif //PHASING
                            
#ifdef DEBUG_MAPPING
                            printf("SUCCESS %d %d \n", pwi, value->a);
                            cout<<pwi<<" "<<index.all_predictions[value->a]->upperCaseSequence<<" "<<read<<endl; //DEB
#endif
                            feed_coherent_positions(index.all_predictions, value->a , pwi, (int)strlen(read), quality, read_set_id, gv);
                            
                        } // end tuple read prediction position is read coherent
                    }
                } // end all infos for the current seed
            } // end all seeds of the read
            
            // Clean the temp mapping positions and tested pwi //
            for (std::map<u_int64_t, set<u_int64_t> > ::iterator it=tested_prediction_and_pwis.begin(); it!=tested_prediction_and_pwis.end(); ++it){
                it->second.clear();
            }
            tested_prediction_and_pwis.clear();
            mapped_prediction_as_set.clear();
            
            // Return the read if necessary //
            if (direction == 0){
                gv.revcomp(read);
                gv.rev (quality);
            }
            
        } // end both directions
        return pwi_and_mapped_predictions;
    }
    
    void operator() (std::pair<Sequence,Sequence>& pair){
        // Shortcut
        char *read1 = strdup(pair.first.toString().c_str());
        char * quality1 = strdup(pair.first.getQuality().c_str());
        char *read2 = strdup(pair.second.toString().c_str());
        char * quality2 = strdup(pair.second.getQuality().c_str());
        
        map<int,int64_t> pwi_and_mapped_predictions1 = core_mapping(read1, quality1);
        map<int,int64_t> pwi_and_mapped_predictions2 = core_mapping(read2, quality2);
        
        // clear (if one still have to check the reverse complement of the read) or free (else) the list of int for each prediction_id on which we tried to map the current read
        
        /////// PHASING
//        #ifdef PHASING
        if (gv.phasing){
        if ((pwi_and_mapped_predictions1.size() + pwi_and_mapped_predictions2.size())>1){                                            // If two or more variants mapped by the same read
            string phased_variant_ids ="";                                          // Create a string containing the (lexicographically) ordered set of variant ids.
            
            for (map<int,int64_t>::iterator it=pwi_and_mapped_predictions1.begin(); it!=pwi_and_mapped_predictions1.end(); ++it){
                //            for (set<pair<int,u_int64_t>> ::iterator it=pwi_and_mapped_predictions.begin(); it!=pwi_and_mapped_predictions.end(); ++it){
                // TODO: optimize this
                //                const int pwi = it->first;
                int64_t var_id =it->second;
                string sign = "";
                if (var_id<0){
                    sign ="-";
                    var_id=-var_id;
                }
                string phased_variant_id = sign+parse_variant_id(index.all_predictions[var_id]->sequence.getComment());
                phased_variant_ids = phased_variant_ids+phased_variant_id+';';
            }
            phased_variant_ids += ' ';
            
            for (map<int,int64_t>::iterator it=pwi_and_mapped_predictions2.begin(); it!=pwi_and_mapped_predictions2.end(); ++it){
                //            for (set<pair<int,u_int64_t>> ::iterator it=pwi_and_mapped_predictions.begin(); it!=pwi_and_mapped_predictions.end(); ++it){
                // TODO: optimize this
                //                const int pwi = it->first;
                int64_t var_id =it->second;
                string sign = "";
                if (var_id<0){
                    sign ="-";
                    var_id=-var_id;
                }
                string phased_variant_id = sign+parse_variant_id(index.all_predictions[var_id]->sequence.getComment());
                phased_variant_ids = phased_variant_ids+phased_variant_id+';';
            }
            
            // Associate this string to the number of times it is seen when mapping this read set
            if (phased_variants.find(phased_variant_ids) == phased_variants.end())  phased_variants[phased_variant_ids] = 1;
            else                                                                    phased_variants[phased_variant_ids] = phased_variants[phased_variant_ids]+1;
        }
        
        pwi_and_mapped_predictions1.clear();
        pwi_and_mapped_predictions2.clear();
        }
//#endif // Phasing
        /////// END PHASING
        
        
        
        free(read1);
        free(quality1);
        free(read2);
        free(quality2);
        
    }
    
    
    
    
    void operator() (Sequence& seq)
    {
        // Shortcut
        char *read = strdup(seq.toString().c_str());
        char * quality = strdup(seq.getQuality().c_str());
        
        map<int,int64_t> pwi_and_mapped_predictions = core_mapping(read, quality);
        
        // clear (if one still have to check the reverse complement of the read) or free (else) the list of int for each prediction_id on which we tried to map the current read
        
        /////// PHASING
        if (gv.phasing){
        if (pwi_and_mapped_predictions.size()>1){                                   // If two or more variants mapped by the same read
            string phased_variant_ids ="";                                          // Create a string containing the (lexicographically) ordered set of variant ids.
            
            for (map<int,int64_t>::iterator it=pwi_and_mapped_predictions.begin(); it!=pwi_and_mapped_predictions.end(); ++it){
                int64_t var_id =it->second;
                string sign = "";
                if (var_id<0){
                    sign ="-";
                    var_id=-var_id;
                }
                string phased_variant_id = sign+parse_variant_id(index.all_predictions[var_id]->sequence.getComment());
                
                //DEBUG
//                cout<<"phased_variant_id        "<<phased_variant_id<<endl;
//                cout<<"from sequence:           "<<index.all_predictions[var_id]->sequence.getComment()<<endl;
//                cout<<"parsed from sequence:    "<<parse_variant_id(index.all_predictions[var_id]->sequence.getComment())<<endl;
//                cout<<it->first<<" "<<it->second<<endl;
//                    phased_variant_id+="_"+index.all_predictions[var_id]->upperCaseSequence; //DEBUG
                //ENDDEBUG
                phased_variant_ids = phased_variant_ids+phased_variant_id+';';
                
            }
            // Associate this string to the number of times it is seen when mapping this read set
            if (phased_variants.find(phased_variant_ids) == phased_variants.end())  phased_variants[phased_variant_ids] = 1;
            else                                                                    phased_variants[phased_variant_ids] = phased_variants[phased_variant_ids]+1;
        }
        
        pwi_and_mapped_predictions.clear();
        }
        /////// END PHASING

        
        free(read);
        free(quality);
        
    }
};

// We a define a functor that will be called during bank parsing
struct ProgressFunctor : public IteratorListener  {  void inc (u_int64_t current)   {  std::cout << ".";  } };

/**
 * Performs the first extension of the algorithm:
 * For each read:
 *  - map it on prediction(s)
 *  - add information to mapped predictions
 * returns the number of mapped reads
 */
u_int64_t ReadMapper::map_all_reads_from_a_file (
                                                 GlobalValues & gv,
                                                 FragmentIndex& index,
                                                 const int read_set_id
                                                 ){
    //////////////////////////////////////////////////////////////////////////
    /////////////// read all reads - storing those coherent with reads ///////
    //////////////////////////////////////////////////////////////////////////
    
    
    
    
    
    u_int64_t number_of_mapped_reads = 0;
    map<string,int> phased_variants;
    
    // Few tests for finding pair of banks.
    //    cout <<inputBank->getId()<<" "<<inputBank->getCompositionNb()<<endl;
    //    for (int subBankId=0; subBankId<inputBank->getCompositionNb(); subBankId++){
    //        cout<<"sub " <<inputBank->getIdNb(subBankId)<<endl;
    //        IBank* subbank = inputBank->getBanks()[subBankId];
    //        cout<<"subank id "<<subbank->getId()<<endl;
    //    }
    
    const std::vector<IBank*>& subbanks = inputBank->getBanks();
    // Test if a bank is composed of two read files.
    if (inputBank->getCompositionNb()==2 && subbanks[0]->getCompositionNb()==1 && subbanks[1]->getCompositionNb()==1){ // PAIRED END
        IBank* bank1 =subbanks[0]; LOCAL(bank1);
        IBank* bank2 =subbanks[1]; LOCAL(bank2);
        PairedIterator<Sequence> *  itPair  = new  PairedIterator<Sequence> (bank1->iterator(), bank2->iterator());
        LOCAL(itPair);
        
        ProgressIterator< std::pair <Sequence, Sequence>> prog_iter (itPair, Stringify::format ("Mapping pairend read set %d", read_set_id).c_str(), bank1->estimateNbItems());
        Dispatcher(nbCores,2047).iterate (prog_iter, Functor(gv, index, read_set_id, &number_of_mapped_reads, phased_variants));
    }
    
    else{ // SINGLE END
        // We create a sequence iterator for the bank with progress information
        ProgressIterator<Sequence> iter (*inputBank, Stringify::format ("Mapping read set %d", read_set_id).c_str());
        Dispatcher(nbCores,2047).iterate (iter, Functor(gv, index, read_set_id, &number_of_mapped_reads, phased_variants));
    }
    
    // PHASING:
    if (gv.phasing){
    stringstream phasingFileName;
    phasingFileName<<"phased_alleles_read_set_id_"<<(read_set_id+1)<<".txt";
    cout<<"print in phasing information in "<<phasingFileName.str()<<endl;
    ofstream phasingFile (phasingFileName.str());
    phasingFile <<"#"<<inputBank->getId()<<endl;
    for (map<string,int>::iterator it=phased_variants.begin(); it!=phased_variants.end(); ++it)
        phasingFile << it->first << " => " << it->second << '\n';
    phasingFile.close();
    }
    // ENDPHASING
    
    
    return number_of_mapped_reads;
}


void ReadMapper::set_read_coherency(GlobalValues& gv, FragmentIndex index){
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////// for each prediction: check those fully coherent and store left and right reads covering them ///////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    int prediction_id;
    for (prediction_id=0;prediction_id < index.all_predictions.size();prediction_id++){
        index.all_predictions[prediction_id]->set_read_coherent(read_set_id,gv);
    } // end all fragments
}








