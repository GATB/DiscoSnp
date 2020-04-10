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
 * fragment_index.c
 *
 *  Created on: 16 sept. 2010
 *      Author: ppeterlo
 */

#include<fragment.h>

inline string getUpperCaseOnly(string stringseq){
    string res="";
//    string stringseq=sequence.toString();
    for (unsigned long i=0;i<stringseq.size();i++){
        if (stringseq.at(i)>=(int)'A' && stringseq.at(i)<=(int)'Z')
            res+=stringseq.at(i);
    }
    return res;
}

Fragment::Fragment(Sequence& seq, const int number_of_read_sets){
    sequence=seq;
    upperCaseSequence =             getUpperCaseOnly(seq.toString());
    read_coherent =                 (bool*) malloc(sizeof(bool)*number_of_read_sets);                          test_alloc(read_coherent);
    number_mapped_reads =           (int*) malloc(sizeof(int)*number_of_read_sets);                            test_alloc(number_mapped_reads);
    local_coverage =                (unsigned char*) malloc(sizeof(unsigned char)*upperCaseSequence.size());   test_alloc(local_coverage);
    sum_qualities =                 (unsigned int*) malloc(sizeof(unsigned int)*number_of_read_sets);          test_alloc(sum_qualities);
    nb_mapped_qualities =           (unsigned int*) malloc(sizeof(unsigned int)*number_of_read_sets);          test_alloc(nb_mapped_qualities);
    nbOfSnps = 0;
    
    if (strncmp("SNP", sequence.getComment().c_str(), strlen("SNP")) == 0)
    {
        nbOfSnps=1; // We don't know yep how many, at least one.
    }
    
    
    for (int i=0; i<number_of_read_sets; i++)
    {
        nb_mapped_qualities[i]=0;
        sum_qualities[i]=0;
        number_mapped_reads[i]=0;
    }
}


void Fragment::set_read_coherent(int read_file_id, GlobalValues gv){
    
    int i;
    // V1: the whole fragment has to be k_read coherent or V2 where the last k positions have no influence on the coherency of the fragment.
    // V2 is appropriate for the cases where the fragment is the end of a sequence (transcript, chromosome) and thus, no read are "longer" than the sequence:
    //    ----------------- fragment
    //    °°°°°°°°°°°°        read
    //    °°°°°°°°°°°°     read
    //         °°°°°°°°°°°° read
    //         °°°°°°°°°°°° read
    const int stop=upperCaseSequence.size()-gv.minimal_read_overlap;
    
    //        cout<<"YYYY stop "<<stop<<" rfid "<<read_file_id<<endl; //DEB
    if(stop<=0){
        
        if(local_coverage[0]<gv.min_coverage[read_file_id]) {read_coherent[read_file_id]=false; return;}
        read_coherent[read_file_id]=true; return;
        
    }
 
    
    if (gv.radseq_option){
        //depth should be homogenous along the prediction, we skip the 3 first and last bases in case of short indels in the read, resulting in different length of reads
        
        //we first check the extremities
        for(i=0;i<3;i++){
            if(((unsigned int)local_coverage[i])<gv.min_coverage[read_file_id])
            {read_coherent[read_file_id]=false; return;}
        }
    
        
        for(i=stop-3;i<stop;i++){
            if(((unsigned int)local_coverage[i])<gv.min_coverage[read_file_id])
            {read_coherent[read_file_id]=false; return;}
        }
        
        //then the heart of the prediction
        unsigned int ref_coverage=((unsigned int)local_coverage[2]);
        
        for(i=3;i<(stop-3);i++){
            if(((unsigned int)local_coverage[i])!=ref_coverage)
            {read_coherent[read_file_id]=false; return;}
        }
    }
    
    for(i=0;i<stop;i++){
        if(((unsigned int)local_coverage[i])<gv.min_coverage[read_file_id])
        {read_coherent[read_file_id]=false; return;}
    }
    
    read_coherent[read_file_id]=true;
}

