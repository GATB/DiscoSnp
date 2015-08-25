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
 * fragment_info.h
 *
 *  Created on: 15 sept. 2010
 *      Author: ppeterlo
 */

#ifndef FRAGMENT_INFO_H_
#define FRAGMENT_INFO_H_

#include <hash.h>
#include <list.h>
#include <gatb/gatb_core.hpp>

//with CHARQUAL, use an unsigned char to store qual values. To avoid overflow, compute average value 'on the fly'
//#ifndef INTQUAL  // if compilation was NOT done with CFLAGS=-DINTQUAL, use the CHARQUAL option 
#define CHARQUAL 
//#endif

#ifndef CLASSICAL_SPANNING
#define KMER_SPANNING // ask more than only all position covered by at least min_coverage reads. Each kmer spanning each position should be covered by at least min_coverage reads.
#endif

//char prefix(const char *pre, const char *str)
//{
//    return strncmp(pre, str, strlen(pre)) == 0;
//}

class FragmentInfo{
public:
    Sequence sequence;
    string upperCaseSequence;
	// fixed once at the beggining:
    char * SNP_positions; // If the fragment is a SNP, stores the positions of the SNPs in order to avoid to authorize errors at these positions. Coded on char, the SNP positions should not be longer than 255
    int nbOfSnps;                  // if zero: the sequence is generic or an indel. Else, number of predicted SNPs
	bool * read_coherent;          // =for every set of reads, 1 if the fragment is detected as read coherent, else =0
	unsigned char ** local_coverage;           //  number of reads covering this position can be a char, min coverage required is low


    unsigned int * sum_qualities; // sum of the mapped qualities for each read set
    unsigned int * nb_mapped_qualities; // number of quality mapped for each read set. If there is a unique read, this is the number of mapped reads. In case of close SNPs, as a unique read may cover several SNPs, this can be bigger.
    
	int * number_mapped_reads;      //for every set of reads, number of reads starting (REPLACE reads_starting)
    
    
    
    FragmentInfo(Sequence& seq, const int number_of_read_sets){
        sequence=seq;
        upperCaseSequence=getUpperCaseOnly();
        read_coherent =                 (bool*) malloc(sizeof(bool)*number_of_read_sets);                          test_alloc(read_coherent);
		number_mapped_reads =           (int*) malloc(sizeof(int)*number_of_read_sets);                            test_alloc(number_mapped_reads);
		local_coverage =                (unsigned char**) malloc(sizeof(unsigned char*)*number_of_read_sets);      test_alloc(local_coverage);
		sum_qualities =                 (unsigned int*) malloc(sizeof(unsigned int)*number_of_read_sets);          test_alloc(sum_qualities);
        nb_mapped_qualities =           (unsigned int*) malloc(sizeof(unsigned int)*number_of_read_sets);          test_alloc(nb_mapped_qualities);
        
        
        nbOfSnps = 0;
        
        if (strncmp("SNP", sequence.getComment().c_str(), strlen("SNP")) == 0)
        {
//        if (prefix("SNP",sequence.getComment().c_str())) {
            nbOfSnps=1; // We don't know yep how many, at least one.
        }
        
        
		for (int i=0; i<number_of_read_sets; i++)
		{
            read_coherent[i]=false;
            
            
			local_coverage[i] = (unsigned char *) malloc(upperCaseSequence.size()*sizeof(unsigned char));            test_alloc(local_coverage[i]);
			for(int z=0;z<upperCaseSequence.size(); z++) local_coverage[i][z]=(unsigned char)0;
            
            
            
            nb_mapped_qualities[i]=0;
            sum_qualities[i]=0;
			number_mapped_reads[i]=0;
            
            
            
		}
    }
    
//    FragmentInfo(FragmentInfo* ref){
//        sequence=ref->sequence;
//        upperCaseSequence=ref->upperCaseSequence;
//        exit(0);
//    }
    
    ~FragmentInfo(){
        
    }
    
    void set_read_coherent(int read_file_id, GlobalValues gv){
        int i;
        // V1: the whole fragment has to be k_read coherent or V2 where the last k positions have no influence on the coherency of the fragment.
        // V2 is appropriate for the cases where the fragment is the end of a sequence (transcript, chromosome) and thus, no read are "longer" than the sequence:
        //    ----------------- fragment
        //    °°°°°°°°°°°°        read
        //    °°°°°°°°°°°°     read
        //         °°°°°°°°°°°° read
        //         °°°°°°°°°°°° read
#ifdef KMER_SPANNING
        const int stop=upperCaseSequence.size()-gv.minimal_read_overlap;

//        cout<<"YYYY stop "<<stop<<" rfid "<<read_file_id<<endl; //DEB
        if(stop<=0){
            if(local_coverage[read_file_id][0]<gv.min_coverage[read_file_id]) {read_coherent[read_file_id]=false; return;}
            read_coherent[read_file_id]=true; return;
        }

#else
        const int stop=strlen(upperCaseSequence);
#endif
//        for(i=0;i<stop;i++) cout<<i<<"--"<<(unsigned int)local_coverage[read_file_id][i]<< " "<<gv.min_coverage[read_file_id]<<endl; //DEB
        for(i=0;i<stop;i++) if((unsigned int)local_coverage[read_file_id][i]<gv.min_coverage[read_file_id]) {read_coherent[read_file_id]=false; return;}
        read_coherent[read_file_id]=true;
    }

    
    string getUpperCaseOnly(){
        int sizeUpper=0;
        string res="";
        string stringseq=sequence.toString();
        for (int i=0;i<stringseq.size();i++){
            if (stringseq.at(i)>=(unsigned int)'A' && stringseq.at(i)<=(unsigned int)'Z')
                res+=stringseq.at(i);
        }
        return res;
    }
    
};


#endif /* FRAGMENT_INFO_H_ */
