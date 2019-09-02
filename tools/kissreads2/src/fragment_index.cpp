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
 * fragment_index.c
 *
 *  Created on: 16 sept. 2010
 *      Author: ppeterlo
 */

#include<fragment_index.h>


//#define DEBUG_INDEXING

int line_num(FILE * f)
{
	rewind(f);
	char c;
	int lines = 0;
	while((c = fgetc(f)) != EOF) if(c == '\n') lines++;
	if(c != '\n') lines++;
	rewind(f);
	return lines;
}



char *  strdup_upper_case(char * in){
    // count number of upper case letters in "in"
    int count =0;
    unsigned long i;
    for(i=0;i<strlen(in);i++) if(in[i]>='A' && in[i]<='Z') count++;
    char * temp = (char *) malloc(sizeof(char)*(count+1)); test_alloc(temp);
    unsigned long  j=0;
    for(i=0;i<strlen(in);i++) if(in[i]>='A' && in[i]<='Z') temp[j++]=in[i];
    temp[j]='\0';
    return temp;
}


char * strdup_first_lower(char * in){
    // count number of first lower case letters in "in"
    int count =0;
    unsigned long  i;
    for(i=0;i<strlen(in);i++)
        if(in[i]>='a' && in[i]<='z') count++;
        else break;
    char * temp = (char *) malloc(sizeof(char)*(count+1)); test_alloc(temp);
    unsigned long  j=0;
    for(i=0;i<strlen(in);i++)
        if(in[i]>='a' && in[i]<='z') temp[j++]=in[i];
        else break;
    temp[j]='\0';
    return temp;
}

char * strdup_last_lower(char * in){
    // count number of first lower case letters in "in"
    int count =0;
    int i;
    for(i=strlen(in)-1;i>=0;i--)
        if(in[i]>='a' && in[i]<='z') {
            count++;
        }
        else break;
    
    char * temp = (char *) malloc(sizeof(char)*(count+1)); test_alloc(temp);
    int j=count-1;
    for(i=strlen(in)-1;i>=0;i--)
        if(in[i]>='a' && in[i]<='z'){
            temp[j--]=in[i];
        }
        else break;
    temp[count]='\0';
    return temp;
}


void FragmentIndex::empty_coverage(){
    unsigned long prediction_id;
    for (prediction_id=0;prediction_id < all_predictions.size();prediction_id++){
        for(unsigned long z=0;z<all_predictions[prediction_id]->upperCaseSequence.size(); z++)
            all_predictions[prediction_id]->local_coverage[z]=(unsigned char)0;
    } // end all fragments
}

// read and store all fragments presents in the pointed file.
// index by seeds of length k all these fragments.
// each fragment is stored twice: one direct, one reverse complement.
void FragmentIndex::index_predictions (BankFasta inputBank, GlobalValues& gv){
	kmer_type coded_seed;
	int i,stop;
    uint64_t total_seeds = 0 ;
    
        
    


#ifdef DEBUG_INDEXING
    printf("indexing predictions, allocating memory for storing info about %d read sets\n",  gv.number_of_read_sets);
#endif
	// each fragment has an id. each seed pobints to couples (id, position). and the fragment is then found thanks to all_fragment[id];
    
    // We create an iterator over this bank.
    Iterator<Sequence>* it = inputBank.iterator();
    LOCAL (it);
    
   
    
    
    // First loop over the sequences: count seeds occurrences
    for (it->first(); !it->isDone(); it->next())
    {
        FragmentInfo * currentFragment = new FragmentInfo(it->item(), gv.number_of_read_sets);
		// read all the seeds present on the fragment
        const char * w = currentFragment->upperCaseSequence.c_str();
        stop=strlen(w)-gv.size_seeds+1;
		for (i=0;i<stop;i+= gv.index_stride){
                coded_seed=gv.codeSeed(w+i); // init the seed (as seeds are not consecutives)
                hash_incr_kmer_count(seeds_count,&coded_seed, gv);
                total_seeds++;
		}
        all_predictions.push_back(currentFragment);
	}
    
#ifdef DEBUG_INDEXING
    printf("%d seeds\n", total_seeds);
#endif
    

    
    seed_table  = (p_couple)calloc(total_seeds,sizeof(couple));
    test_alloc(seed_table);
    iterate_and_fill_offsets(seeds_count,gv);
    
    
    
    total_seeds=0;
    ///second loop over fragments  : create the index
    for(unsigned long fragment_id=0;fragment_id<all_predictions.size();fragment_id++){
        
        const char * w = all_predictions[fragment_id  ]->upperCaseSequence.c_str();
#ifdef DEBUG_INDEXING
		printf("indexing in %s\n", all_predictions[fragment_id]->upperCaseSequence.c_str());
#endif
		// read all the seeds present on the fragment
		stop=strlen(w)-gv.size_seeds+1;
		for (i=0;i<stop;i+= gv.index_stride){
            coded_seed=gv.codeSeed(w+i); // init the seed
            hash_fill_kmer_index(seeds_count,&coded_seed,seed_table, fragment_id, i,gv);
            total_seeds++;
		}
	}
    
    
#ifdef DEBUG_INDEXING
    printf("%d seeds\n", total_seeds);
    //check the first sequence indexed:
    
    uint64_t offset_seed;
    uint64_t nb_occurrences;
    const char * w = all_predictions[0  ]->upperCaseSequence.c_str();
    cout<<" checking seeds in "<<w<<endl;
    stop=strlen(w)-gv.size_seeds+1;
    for (i=0;i<stop;i++){
        coded_seed=gv.codeSeed(w+i);
        cout<<"i = "<<i<<" "<<coded_seed<<endl;
        
        if(get_seed_info(seeds_count,&coded_seed,&offset_seed,&nb_occurrences,gv)){
            cout<<"nb_occurrences "<<nb_occurrences<<endl;
            // for each occurrence of this seed on the starter:
            for (int ii=offset_seed; ii<offset_seed+nb_occurrences; ii++) {
                couple * value = &(seed_table[ii]);
                const char * starter = all_predictions[value->a]->upperCaseSequence.c_str();
                cout<<"i="<<i<<" value->b="<<value->b<<" value->a="<<value->a<<endl;
                cout<<"seed"<<w+i<<" starts "<<starter+value->b<<endl;
            }
        }
    }
#endif

    
    
    ///third loop over fragments : for SNPs, store the SNP positions
    for(unsigned long fragment_id=0;fragment_id<all_predictions.size();fragment_id+=2){
        
        if ( all_predictions[fragment_id]->nbOfSnps==0 ) { // This is an indel.
            all_predictions[fragment_id]->SNP_positions = (unsigned int *) malloc (sizeof(unsigned int)); // add a dummy contrained positions
            test_alloc(all_predictions[fragment_id]->SNP_positions);
            all_predictions[fragment_id]->SNP_positions[0] = max(all_predictions[fragment_id  ]->upperCaseSequence.size(), all_predictions[fragment_id+1]->upperCaseSequence.size())+1; // DUMMY SNP
            continue;
        } // the rest applies only for SNPs
        
        const char * seq1 = all_predictions[fragment_id  ]->upperCaseSequence.c_str();
        const char * seq2 = all_predictions[fragment_id+1]->upperCaseSequence.c_str();
        const int size_seq1 = strlen(seq1);
        const int size_seq2 = strlen(seq2);
        const int size_seq = min(size_seq1,size_seq2);
        
//        COMMENTED ON NOV 2017: WITH RAD SEQ DATA, SNP SEQUENCES MAY HAVE DISCTINCT SIZES
//        assert(size_seq == strlen(seq2));
//        if(size_seq != strlen(seq2)){
//            cerr<<"two SNP sequences of distinct sizes. Impossible"<<endl;
//            cerr<<"ID="<<fragment_id<<endl;
//            cerr<<seq1<<endl;
//            cerr<<seq2<<endl;
//            exit(1);
//        }
        
        // compute the number of SNPs:
        int local_number_of_SNPs=0;
        for (i=0; i<size_seq; i++) {
            if (seq1[i]!=seq2[i]) {
                local_number_of_SNPs++;
            }
        }
        all_predictions[fragment_id]->nbOfSnps=local_number_of_SNPs;
        all_predictions[fragment_id]->SNP_positions = (unsigned int *) malloc (sizeof(unsigned int)*(local_number_of_SNPs+1)); // add a dummy SNP
        test_alloc(all_predictions[fragment_id]->SNP_positions);
        
        // we do not fill the all_predictions[fragment_id+1] with the same information
        local_number_of_SNPs=0;
        for (i=0; i<size_seq; i++) {
            if (seq1[i]!=seq2[i]) {
                all_predictions[fragment_id]->SNP_positions[local_number_of_SNPs]=i;
                local_number_of_SNPs++;
            }
        }
        all_predictions[fragment_id]->SNP_positions[local_number_of_SNPs] = max(size_seq1,size_seq2)+1; // DUMMY SNP
	}
    
    
}
