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

#include<fragment_info.h>


void FragmentInfo::set_read_coherent(int read_file_id, GlobalValues gv){
    
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
        
        if(local_coverage[0]<gv.min_coverage[read_file_id]) {read_coherent[read_file_id]=false; return;}
        read_coherent[read_file_id]=true; return;
        
    }
    
#else
    const int stop=strlen(upperCaseSequence);
#endif
    //        for(i=0;i<stop;i++) cout<<i<<"--"<<(unsigned int)local_coverage[read_file_id][i]<< " "<<gv.min_coverage[read_file_id]<<endl; //DEB
    
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

