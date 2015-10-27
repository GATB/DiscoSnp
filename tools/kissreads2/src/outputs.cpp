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
 * outputs.c
 *
 *  Created on: 27 oct. 2010
 *      Author: ppeterlo
 */

//#define DEBUG_QUALITY


//#if !HAVE_LOG2F
//#define log2f log
//#endif

#include <outputs.h>

///*Calculates the phi coefficient of 2*2 contingency table. Value close to 1 indicates an association between the alleles and the conditions.*/
///*Note that this value is valid if at least 3 out of the 4 values are non 0, or if the .Otherwise it output -1*/
//float phi(int a,int b, int c,int d) {
//    //  int denom=(a+b)*(c+d)*(a+c)*(b+d);
//    //  if (denom==0)
//    //    return 0;
//    //  float Phi = (a*d-b*c)/sqrt(denom);
//    //  return Phi;
//    if((a+b)==0) return 0;
//    if((c+d)==0) return 0;
//    if((a+c)==0) return 0;
//    if((b+d)==0) return 0;
//    // avoid the computation of denom, possibly bigger than an int or an unsigned long long int...
//    return (a*d-b*c)/(sqrt((float)(a+b))*sqrt((float)(c+d))*sqrt((float)(a+c))*sqrt((float)(b+d)));
//}

/*Computes the chi2 value of the matrix 2*number_of_read_sets */
float rank_phi_N (const int *sum_up, const int *sum_lo, const int number_of_read_sets) {
    if (number_of_read_sets==1)
        return 0;
    int i;
    float n=0; for (i=0;i<number_of_read_sets;i++) n+=sum_up[i]+sum_lo[i];
    float all_up=0; for (i=0;i<number_of_read_sets;i++) all_up+=sum_up[i];
    float all_lo=0; for (i=0;i<number_of_read_sets;i++) all_lo+=sum_lo[i];
    float expected;
    
    
    float som=0;
    for(i=0;i<number_of_read_sets;i++){
        // UPPER PATH
        expected=(sum_up[i]+sum_lo[i])*all_up/n;
        if(expected!=0) som+=pow(sum_up[i]-expected,2)/expected;
        // LOWER PATH
        expected=(sum_up[i]+sum_lo[i])*all_lo/n;
        if(expected!=0) som+=pow(sum_lo[i]-expected,2)/expected;
    }
    
    return sqrt(som/n);
}

///*Computes all pairwise phi values for all pairs of conditions and returns the max*/
//float rank_phi(const int *sum_up, const int *sum_lo, const int number_of_read_sets) {
//    float phimax=0;
//    if (number_of_read_sets==1)
//        return 0;
//    else
//    {
//        int i,j;
//        float phicur=0;
//        for (i=0;i<number_of_read_sets;i++)
//            for (j=i+1;j<number_of_read_sets;j++)
//            {
//                phicur=phi(sum_up[i],sum_up[j],sum_lo[i],sum_lo[j]);
//                phimax=MAX(phimax,ABS(phicur));
//            }
//    }
//    return phimax;
//}

/**
 * Computes the log10(Cnk)
 */
const double mylog10choose (const int n,const int k){
    if (n==k){
        return 0;
    }
    double  res=0;
    int i;
    for(i=k+1;i<=n;i++){
        res+=log10(i);
    }
    
    for(i=1;i<=n-k;i++) {
        res-=log10(i);
    }
    
    return res;
}


char * genotype_simple_model(const int c1, const int c2, const float err, const float prior_het){
    
    // LIKELIHOOD
    double lik0 = c1*log10(1-err)+c2*log10(err)+mylog10choose(c1+c2,c1);
    double lik1 = c2*log10(1-err)+c1*log10(err)+mylog10choose(c1+c2,c1);
    double lik2 = (c1+c2)*log10(0.5)+mylog10choose(c1+c2,c1);
    
    // PRIOR HETEROZYGOUS
    lik0+=log10((1-prior_het)/2);
    lik1+=log10((1-prior_het)/2);
    lik2+=log10(prior_het);
    
    // PHRED SCORE
    lik0=floor(-10*lik0);
    lik1=floor(-10*lik1);
    lik2=floor(-10*lik2);
    
    // FORMATING RESULTS
    char * append = (char *)malloc(sizeof(char)*2048); test_alloc(append);
    char geno[4];
    if (lik0<lik1 &&lik0<lik2){
        sprintf(geno, "0/0");
    }
    else
    {
        
        if (lik1<lik0 && lik1<lik2){
            sprintf(geno, "1/1");
        }
        
        else{
            sprintf(geno, "0/1");
        }
    }
    sprintf(append, "%s:%d,%d,%d",geno,(int)lik0,(int)lik2,(int)lik1);
    return append;
}

/**
 * prints a couple using the reads starting position instead of coverage per position
 */
void print_couple_i(ofstream &fasta_out, FragmentIndex & index, int fragment_id, GlobalValues & gv){
	
    
    // on upper path
	int sum_up[gv.number_of_read_sets];
	int avg_up[gv.number_of_read_sets];
    
	// on lower path
	int sum_lo[gv.number_of_read_sets];
	int avg_lo[gv.number_of_read_sets];
    
    
	int read_set_id;
    
    
    
   	//if( gv.qual ){
    
    for(read_set_id=0;read_set_id<gv.number_of_read_sets;read_set_id++){
        avg_up[read_set_id] = 0;
        avg_lo[read_set_id] = 0;
        
        if (index.all_predictions[fragment_id ]->nb_mapped_qualities[read_set_id]>0)
            avg_up[read_set_id] += index.all_predictions[fragment_id  ]->nb_mapped_qualities[read_set_id]==0?0:index.all_predictions[fragment_id  ]->sum_qualities[read_set_id] / index.all_predictions[fragment_id  ]->nb_mapped_qualities[read_set_id];
        if (index.all_predictions[fragment_id+1]->nb_mapped_qualities[read_set_id]>0)
            avg_lo[read_set_id] += index.all_predictions[fragment_id+1]->nb_mapped_qualities[read_set_id]==0?0:index.all_predictions[fragment_id+1]->sum_qualities[read_set_id] / index.all_predictions[fragment_id+1]->nb_mapped_qualities[read_set_id];
    }
    //    }
    
	//	float sum=0;
	for(int read_set_id=0;read_set_id<gv.number_of_read_sets;read_set_id++){
        sum_up[read_set_id]=index.all_predictions[fragment_id]->number_mapped_reads[read_set_id];
        sum_lo[read_set_id]=index.all_predictions[fragment_id+1]->number_mapped_reads[read_set_id];
	}
    const float err = 0.01;
    const float prior_het = 1/(float)3;
    float rank = rank_phi_N(sum_up,sum_lo,gv.number_of_read_sets);
    char genotypes[819200]; genotypes[0]='\0';
    char append[160000];
    
    if(gv.compute_genotypes){
        // CONSTRUCT THE COMMON HEADER COMMENT (Genotypes, Coverages, Qualities, Rank)
        for(read_set_id=0;read_set_id<gv.number_of_read_sets;read_set_id++){
            char * geno_likelihood = genotype_simple_model(sum_up[read_set_id], sum_lo[read_set_id], err, prior_het);
            sprintf(append, "G%d_%s|",read_set_id+1,geno_likelihood);
            free(geno_likelihood);
            strcat(genotypes,append);
        }
    }
    //    cout<<"genotypes"<<genotypes<<endl; //DEB
    // DEAL WITH STANDARD FASTA OR ONE LINE PER COUPLE (STARTING WITH THE RANK)
    char sep;
    if (gv.standard_fasta) {
        
        sep='\n';
    }
    else{
        fasta_out<<setprecision(5)<<rank<<" ";
        //        fprintf(out, "%2f ",rank);
        sep=';';
    }
    
    
    // UPPER PATH
    fasta_out<<">"<<index.all_predictions[fragment_id]->sequence.getComment()<<"|";
    for(read_set_id=0;read_set_id<gv.number_of_read_sets;read_set_id++){
        fasta_out<<"C"<<read_set_id+1<<"_"<<sum_up[read_set_id]<<"|";
    }
    for(read_set_id=0;read_set_id<gv.number_of_read_sets;read_set_id++){
        fasta_out<<"Q"<<read_set_id+1<<"_"<<avg_up[read_set_id]<<"|";
    }
    fasta_out<<genotypes;
    fasta_out<<"rank_"<<setprecision(5)<<rank;
    
    fasta_out<<sep<<index.all_predictions[fragment_id]->sequence.toString()<<sep;
    
    // LOWER PATH
    fasta_out<<">"<<index.all_predictions[fragment_id+1]->sequence.getComment()<<"|";
    for(read_set_id=0;read_set_id<gv.number_of_read_sets;read_set_id++){
        fasta_out<<"C"<<read_set_id+1<<"_"<<sum_lo[read_set_id]<<"|";
    }
    for(read_set_id=0;read_set_id<gv.number_of_read_sets;read_set_id++){
        fasta_out<<"Q"<<read_set_id+1<<"_"<<avg_lo[read_set_id]<<"|";
    }
    fasta_out<<genotypes;
    fasta_out<<"rank_"<<setprecision(5)<<rank;
    
    fasta_out<<sep<<index.all_predictions[fragment_id+1]->sequence.toString()<<endl;
    
    
    
}



/**
 * checks if at least one read set provide read coherency for a path.
 */
inline bool one_coherent(FragmentInfo * fragment, int number_of_read_sets, GlobalValues & gv){
    int i;
    for(i=0;i<number_of_read_sets;i++){
        if(fragment->is_read_coherent(i,gv)) return true;
//        if(fragment->read_coherent[i]) return true;
        
    }
    return false;
}


void print_results_2_paths_per_event(ofstream &coherent_out, ofstream &uncoherent_out,  FragmentIndex &index, GlobalValues & gv){
     index.nb_coherent=0;
     index.nb_uncoherent=0;
    //printf("number ofread sets = %d\n", number_of_read_sets);
    
    //
    //                 C1           C2           C3 ....
    // path1 (i)      [0/1]        [0/1]        [0/1]...
    // path2 (i+1)    [0/1]        [0/1]        [0/1]...
    //
    // event is kept only if each line has at least one "1" per line:
    //
    
    
    
    for(int i=0;i<index.all_predictions.size();i+=2){
        //        cout<<"HEY "<<one_coherent(index.all_predictions[i],gv.number_of_read_sets)<<" -- "<<one_coherent(index.all_predictions[i+1],gv.number_of_read_sets)<<endl; //DEB
        if(one_coherent(index.all_predictions[i],gv.number_of_read_sets,gv) || one_coherent(index.all_predictions[i+1],gv.number_of_read_sets,gv))
        {
            index.nb_coherent++;
            print_couple_i(coherent_out, index, i, gv);
            //                           results_against_set, i, number_of_read_sets, qual, 1, compute_genotype, paired);
        }
        else{
            index.nb_uncoherent++;
            print_couple_i(uncoherent_out, index, i, gv);
            //            , results_against_set, i, number_of_read_sets, qual, 1, compute_genotype, paired);
        }
    }
}





