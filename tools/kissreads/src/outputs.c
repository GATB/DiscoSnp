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
#include <fragment_index.h>
#include <fragment_info.h>
#include <commons.h>
#include <string.h>

//#include <tree.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <assert.h>
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define ABS(a) (((a) < 0) ? -(a) : (a))

//#define DEBUG_QUALITY


#if !HAVE_LOG2F
#define log2f log
#endif

// Operation made on each data set (read_set_id)
//#define op() MIN(corrected_avg_lo[read_set_id],corrected_avg_up[read_set_id]) /	(corrected_avg_lo[read_set_id]+corrected_avg_up[read_set_id]);
#define op() corrected_avg_up[read_set_id] /	(corrected_avg_lo[read_set_id]+corrected_avg_up[read_set_id]);




float rank_max(const float * corrected_avg_up, const float * corrected_avg_lo, const int number_of_read_sets){
	int read_set_id;
	float rank=0;
	float temp;
	//	float sum=0;
	for(read_set_id=0;read_set_id<number_of_read_sets;read_set_id++){
		if(corrected_avg_lo[read_set_id]+corrected_avg_up[read_set_id]){ // not null in the 2 isoforms
			temp = op()
			if(temp > rank)
				rank = temp;
		}
	}
	return rank;
}

float rank_min(const float * corrected_avg_up, const float * corrected_avg_lo, const int number_of_read_sets){
	int read_set_id;
	float rank=INT_MAX;
	float temp;
	for(read_set_id=0;read_set_id<number_of_read_sets;read_set_id++){
		if(corrected_avg_lo[read_set_id]+corrected_avg_up[read_set_id]){ // not null in the 2 isoforms
			temp = op()
			if(temp < rank)
				rank = temp;
		}
	}
	return rank;
}

float rank_sum(const float * corrected_avg_up, const float * corrected_avg_lo, const int number_of_read_sets){
	int read_set_id;
	float rank=0;
	for(read_set_id=0;read_set_id<number_of_read_sets;read_set_id++){
		if(corrected_avg_lo[read_set_id]+corrected_avg_up[read_set_id]){ // not null in the 2 isoforms
			rank += op()
		}
	}
	return rank;
}

float rank_avg(const float * corrected_avg_up, const float * corrected_avg_lo, const int number_of_read_sets){
	return rank_sum(corrected_avg_up,corrected_avg_lo,number_of_read_sets)/(float)number_of_read_sets;
}

float rank_standard_deviation (const float * corrected_avg_up, const float * corrected_avg_lo, const int number_of_read_sets){
	float avg= rank_avg(corrected_avg_up,corrected_avg_lo,number_of_read_sets);
	float temp;
	float sdev=0.0;
	float sum=0.0;
	int read_set_id;
	for(read_set_id=0;read_set_id<number_of_read_sets;read_set_id++){
		if(corrected_avg_lo[read_set_id]+corrected_avg_up[read_set_id]){ // not null in the 2 isoforms
			temp=op()
			sum+=pow(temp-avg,2);
		}
	}
	sdev=sum/(number_of_read_sets-1);
	sdev = sqrt(sdev);
	return sdev;
}

/**
 * TODO: the entropy function raises some problems: any path having zero value creates a -inf value with the log of the fraction...
 */
float rank_entropy(const float * corrected_avg_up, const float * corrected_avg_lo, const int number_of_read_sets){
    int read_set_id;
	float ranksum=rank_sum(corrected_avg_up,corrected_avg_lo,number_of_read_sets);
	float temp;
	float entropy = 0.0;
	float fraction;
    //	printf("hey (sum =%f): ",ranksum);
	for(read_set_id=0;read_set_id<number_of_read_sets;read_set_id++){
        //		printf("(%f/%f) ", corrected_avg_lo[read_set_id], corrected_avg_up[read_set_id]);
		if(corrected_avg_lo[read_set_id]+corrected_avg_up[read_set_id]){ // not null in the 2 isoforms
			temp=op()
			fraction = temp/ranksum;
            //			printf("[fraction=%f] ",fraction);
			entropy -= (fraction * log2f (fraction));
		}
	}
    //	printf(" --> %f\n", entropy);
	return entropy;
}

/*Calculates the phi coefficient of 2*2 contingency table. Value close to 1 indicates an association between the alleles and the conditions.*/
/*Note that this value is valid if at least 3 out of the 4 values are non 0, or if the .Otherwise it output -1*/
float phi(int a,int b, int c,int d) {
    //  int denom=(a+b)*(c+d)*(a+c)*(b+d);
    //  if (denom==0)
    //    return 0;
    //  float Phi = (a*d-b*c)/sqrt(denom);
    //  return Phi;
    if((a+b)==0) return 0;
    if((c+d)==0) return 0;
    if((a+c)==0) return 0;
    if((b+d)==0) return 0;
    // avoid the computation of denom, possibly bigger than an int or an unsigned long long int...
    return (a*d-b*c)/(sqrt((float)(a+b))*sqrt((float)(c+d))*sqrt((float)(a+c))*sqrt((float)(b+d)));
}

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

/*Computes all pairwise phi values for all pairs of conditions and returns the max*/
float rank_phi(const int *sum_up, const int *sum_lo, const int number_of_read_sets) {
    float phimax=0;
    if (number_of_read_sets==1)
        return 0;
    else
    {
        int i,j;
        float phicur=0;
        for (i=0;i<number_of_read_sets;i++)
            for (j=i+1;j<number_of_read_sets;j++)
            {
                phicur=phi(sum_up[i],sum_up[j],sum_lo[i],sum_lo[j]);
                phimax=MAX(phimax,ABS(phicur));
            }
    }
    return phimax;
}

//const char * genotype_simple_model(const int c1, const int c2, const float err){ // WITH PRIOR HET = 0.33
//    const float lik0 = c1*log10(1-err)+c2*log10(err)+log10(Cnp(c1,c1+c2));
//    const float lik1 = c2*log10(1-err)+c1*log10(err)+log10(Cnp(c1,c1+c2));
//    const float lik2 = (c1+c2)*log10(0.5)+log10(Cnp(c1,c1+c2));
//    
//    
////    revoyer les trois avec la plus petite du moins 1à comme reference de genotype.a
//    
//    
//    
//    pow(1-err,c1)*pow(err,c2); // homozygous higher (0/0)
//    const float lik1 = pow(1-err,c2)*pow(err,c1); // homozygous higher (1/1)
//    const float lik2 = pow(0.5,c1+c2);            // heterozygous (0/1)
//    
//    //    printf("%d %d LIKE %f %f %f \n", c1, c2, lik0, lik1, lik2);
//    
//    const float prob0 = lik0*(1-prior_het)/2;
//    const float prob1 = lik1*(1-prior_het)/2;
//    const float prob2 = lik2*prior_het;
//    //    printf("PROB %f %f %f \n", prob0, prob1, prob2);
//    
//    if ( prob0>=prob1 &&  prob0>=prob2) return "0/0";
//    if ( prob1>=prob0 &&  prob1>=prob2) return "1/1";
//    return "0/1";
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
    char * append = (char *)malloc(sizeof(char)*2048);
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
void print_couple_i(char * comment, FILE* out, const p_fragment_info * results_against_set, int cycle_id, int number_of_read_sets, int qual, const char map_snps, const char compute_genotype, const char paired){
	
    
    if (paired) {
        number_of_read_sets/=2;
    }
    // on upper path
	int sum_up[number_of_read_sets];
	int avg_up[number_of_read_sets];
    
	// on lower path
	int sum_lo[number_of_read_sets];
	int avg_lo[number_of_read_sets];
    
    
	int read_set_id;
    
   	if( qual ){
        for(read_set_id=0;read_set_id<number_of_read_sets;read_set_id++){
            avg_up[read_set_id] = 0;
            avg_lo[read_set_id] = 0;
            if(!paired){
                if (results_against_set[cycle_id  ]->nb_mapped_qualities[read_set_id]>0)
                    avg_up[read_set_id] += results_against_set[cycle_id  ]->sum_qualities[read_set_id] / results_against_set[cycle_id  ]->nb_mapped_qualities[read_set_id];
                if (results_against_set[cycle_id+1]->nb_mapped_qualities[read_set_id]>0)
                    avg_lo[read_set_id] += results_against_set[cycle_id+1]->sum_qualities[read_set_id] / results_against_set[cycle_id+1]->nb_mapped_qualities[read_set_id];
            }
            else{ // PAIRED
                
                if(results_against_set[cycle_id  ]->nb_mapped_qualities[read_set_id*2  ] + results_against_set[cycle_id  ]->nb_mapped_qualities[(read_set_id*2)+1] >0)
                    avg_up[read_set_id] += (results_against_set[cycle_id  ]->sum_qualities[read_set_id*2    ] + results_against_set[cycle_id  ]->sum_qualities[(read_set_id*2)+1] )
                    / (results_against_set[cycle_id  ]->nb_mapped_qualities[read_set_id*2  ] + results_against_set[cycle_id  ]->nb_mapped_qualities[(read_set_id*2)+1]) ;
                
                
                if(results_against_set[cycle_id+1]->nb_mapped_qualities[read_set_id*2  ] + results_against_set[cycle_id+1]->nb_mapped_qualities[(read_set_id*2)+1] >0)
                    avg_lo[read_set_id] += (results_against_set[cycle_id+1]->sum_qualities[read_set_id*2    ] + results_against_set[cycle_id+1]->sum_qualities[(read_set_id*2)+1] )
                    / (results_against_set[cycle_id+1]->nb_mapped_qualities[read_set_id*2  ] + results_against_set[cycle_id+1]->nb_mapped_qualities[(read_set_id*2)+1]) ;
                
            }
        }
    }
	
	//	float sum=0;
	for(read_set_id=0;read_set_id<number_of_read_sets;read_set_id++){
        if (paired){
            sum_up[read_set_id]=results_against_set[cycle_id]->number_mapped_reads[read_set_id*2];
            sum_lo[read_set_id]=results_against_set[cycle_id+1]->number_mapped_reads[read_set_id*2];
            sum_up[read_set_id]+=results_against_set[cycle_id]->number_mapped_reads[read_set_id*2+1];
            sum_lo[read_set_id]+=results_against_set[cycle_id+1]->number_mapped_reads[read_set_id*2+1];
        }
        else
        {
            sum_up[read_set_id]=results_against_set[cycle_id]->number_mapped_reads[read_set_id];
            sum_lo[read_set_id]=results_against_set[cycle_id+1]->number_mapped_reads[read_set_id];
        }
	}
    const float err = 0.01;
    const float prior_het = 1/(float)3;
    
    float rank = rank_phi_N(sum_up,sum_lo,number_of_read_sets);
    char genotypes[8192]; genotypes[0]='\0';
    char append[16000];
    
    if(compute_genotype){
        // CONSTRUCT THE COMMON HEADER COMMENT (Genotypes, Coverages, Qualities, Rank)
        for(read_set_id=0;read_set_id<number_of_read_sets;read_set_id++){
            char * geno_likelihood = genotype_simple_model(sum_up[read_set_id], sum_lo[read_set_id], err, prior_het);
            sprintf(append, "G%d_%s|",read_set_id+1,geno_likelihood);
            free(geno_likelihood);
            strcat(genotypes,append);
        }
    }
    
    // DEAL WITH STANDARD FASTA OR ONE LINE PER COUPLE (STARTING WITH THE RANK)
    char sep;
    if (standard_fasta) {
        
        sep='\n';
    }
    else{
        fprintf(out, "%2f ",rank);
        sep=';';
    }
    
    
    // UPPER PATH
    fprintf(out, ">%s%s|",comment,results_against_set[cycle_id]->comment);
    for(read_set_id=0;read_set_id<number_of_read_sets;read_set_id++){
        fprintf(out, "C%d_%d|",read_set_id+1,sum_up[read_set_id]);
    }
    fprintf(out, "%s",genotypes);
    if (qual)
        for(read_set_id=0;read_set_id<number_of_read_sets;read_set_id++){
            fprintf(out, "Q%d_%d|",read_set_id+1,avg_up[read_set_id]);
        }
    fprintf(out, "rank_%.5f",rank);
    
    
    if(map_snps)
        fprintf(out, "%c%s%s%s%c", sep, results_against_set[cycle_id]->left_extension, results_against_set[cycle_id]->w, results_against_set[cycle_id]->right_extension, sep);
    else
        fprintf(out, "%c%s%c", sep, results_against_set[cycle_id]->w, sep);
    
    // LOWER PATH
    fprintf(out, ">%s%s|",comment,results_against_set[cycle_id+1]->comment);
    for(read_set_id=0;read_set_id<number_of_read_sets;read_set_id++){
        fprintf(out, "C%d_%d|",read_set_id+1,sum_lo[read_set_id]);
    }
    fprintf(out, "%s",genotypes);
    if (qual)
        for(read_set_id=0;read_set_id<number_of_read_sets;read_set_id++){
            fprintf(out, "Q%d_%d|",read_set_id+1,avg_lo[read_set_id]);
        }
    fprintf(out, "rank_%.5f",rank);
    if(map_snps)
        fprintf(out, "%c%s%s%s\n", sep,results_against_set[cycle_id+1]->left_extension, results_against_set[cycle_id+1]->w, results_against_set[cycle_id+1]->right_extension);
    else
        fprintf(out, "%c%s\n", sep,results_against_set[cycle_id+1]->w);
}



/**
 * prints a couple using the reads starting position instead of coverage per position
 */
void print_quadruplet_i(FILE* out, const p_fragment_info * results_against_set, int cycle_id, int number_of_read_sets, int qual){
	int cov_1[number_of_read_sets] ; // coverage path 1 au
    int cov_2[number_of_read_sets] ; // coverage path 2 vb
    int cov_3[number_of_read_sets] ; // coverage path 3 av'
    int cov_4[number_of_read_sets] ; // coverage path 4 u'b
    
    
    int qual_1[number_of_read_sets]; // quality path 1
    int qual_2[number_of_read_sets]; // quality path 2
    int qual_3[number_of_read_sets]; // quality path 3
    int qual_4[number_of_read_sets]; // quality path 4
    int read_set_id;
    
    
    
//	if( qual ){// TODO: UNTESTED CODE - APRIL 2013
//        // we are providing results for generic dataset
//        for(read_set_id=0;read_set_id<number_of_read_sets;read_set_id++){
//            qual_1[read_set_id] = 0;
//            qual_2[read_set_id] = 0;
//            qual_3[read_set_id] = 0;
//            qual_4[read_set_id] = 0;
//            for (j=kmer_size-1;j<=strlen(results_against_set[cycle_id]->w)-kmer_size;j++){
//#ifdef CHARQUAL // FIXME: IT SHOULKD BE THE OPOSIT NO ? (PIERRE APRL 2013)
//                if(results_against_set[cycle_id]->read_coherent_positions[read_set_id][j]) qual_1[read_set_id] = qual_1[read_set_id] + results_against_set[cycle_id]->sum_quality_per_position[read_set_id][j];
//                if(results_against_set[cycle_id+1]->read_coherent_positions[read_set_id][j]) qual_2[read_set_id] = qual_2[read_set_id] + results_against_set[cycle_id+1]->sum_quality_per_position[read_set_id][j];
//                if(results_against_set[cycle_id+2]->read_coherent_positions[read_set_id][j]) qual_3[read_set_id] = qual_3[read_set_id] + results_against_set[cycle_id+2]->sum_quality_per_position[read_set_id][j];
//                if(results_against_set[cycle_id+3]->read_coherent_positions[read_set_id][j]) qual_4[read_set_id] = qual_4[read_set_id] + results_against_set[cycle_id+3]->sum_quality_per_position[read_set_id][j];
//#else
//                if(results_against_set[cycle_id]->read_coherent_positions[read_set_id][j]) qual_1[read_set_id] = qual_1[read_set_id] + (results_against_set[cycle_id]->sum_quality_per_position[read_set_id][j] / results_against_set[cycle_id]->read_coherent_positions[read_set_id][j]);
//                if(results_against_set[cycle_id+1]->read_coherent_positions[read_set_id][j]) qual_2[read_set_id] = qual_2[read_set_id] + (results_against_set[cycle_id+1]->sum_quality_per_position[read_set_id][j] / results_against_set[cycle_id+1]->read_coherent_positions[read_set_id][j]);
//                if(results_against_set[cycle_id+2]->read_coherent_positions[read_set_id][j]) qual_3[read_set_id] = qual_3[read_set_id] + (results_against_set[cycle_id+2]->sum_quality_per_position[read_set_id][j] / results_against_set[cycle_id+2]->read_coherent_positions[read_set_id][j]);
//                if(results_against_set[cycle_id+3]->read_coherent_positions[read_set_id][j]) qual_4[read_set_id] = qual_4[read_set_id] + (results_against_set[cycle_id+3]->sum_quality_per_position[read_set_id][j] / results_against_set[cycle_id+3]->read_coherent_positions[read_set_id][j]);
//                
//#endif
//            }
//        }
//    } // END UNTESTED CODE
    
    
    // on upper path
	int sum_up[number_of_read_sets];
    
	// on lower path
	int sum_lo[number_of_read_sets];
    
//    // considering the uncoherent as covered by 0 reads
//	for(read_set_id=0;read_set_id<number_of_read_sets;read_set_id++){
//		cov_1[read_set_id]=results_against_set[cycle_id]->read_coherent[read_set_id]?results_against_set[cycle_id]->number_mapped_reads[read_set_id]:0;
//        cov_2[read_set_id]=results_against_set[cycle_id+1]->read_coherent[read_set_id]?results_against_set[cycle_id+1]->number_mapped_reads[read_set_id]:0;
//        sum_up[read_set_id]=cov_1[read_set_id]+cov_2[read_set_id];
//        cov_3[read_set_id]=results_against_set[cycle_id+2]->read_coherent[read_set_id]?results_against_set[cycle_id+2]->number_mapped_reads[read_set_id]:0;
//        cov_4[read_set_id]=results_against_set[cycle_id+3]->read_coherent[read_set_id]?results_against_set[cycle_id+3]->number_mapped_reads[read_set_id]:0;
//        sum_lo[read_set_id]=cov_3[read_set_id]+cov_4[read_set_id];
//	}
     // not changing the uncoherent
	for(read_set_id=0;read_set_id<number_of_read_sets;read_set_id++){
		cov_1[read_set_id]=results_against_set[cycle_id]->number_mapped_reads[read_set_id];
        cov_2[read_set_id]=results_against_set[cycle_id+1]->number_mapped_reads[read_set_id];
        sum_up[read_set_id]=cov_1[read_set_id]+cov_2[read_set_id];
        cov_3[read_set_id]=results_against_set[cycle_id+2]->number_mapped_reads[read_set_id];
        cov_4[read_set_id]=results_against_set[cycle_id+3]->number_mapped_reads[read_set_id];
        sum_lo[read_set_id]=cov_3[read_set_id]+cov_4[read_set_id];
	}
   

    float rank = rank_phi_N(sum_up,sum_lo,number_of_read_sets);
    
    if (!standard_fasta)
    {
        // PATH1
        fprintf(out, ">%s|",results_against_set[cycle_id]->comment);
        for(read_set_id=0;read_set_id<number_of_read_sets;read_set_id++) fprintf(out, "C%d_%d|",read_set_id+1,cov_1[read_set_id]);
        if (qual) for(read_set_id=0;read_set_id<number_of_read_sets;read_set_id++) fprintf(out, "Q%d_%d|",read_set_id+1,qual_1[read_set_id]);
        fprintf(out, "rank_%.5f",rank);
        fprintf(out, ";%s", results_against_set[cycle_id]->w);
        fprintf(out, ";");
        
        // PATH2
        fprintf(out, ">%s|",results_against_set[cycle_id+1]->comment);
        for(read_set_id=0;read_set_id<number_of_read_sets;read_set_id++) fprintf(out, "C%d_%d|",read_set_id+1,cov_2[read_set_id]);
        if (qual) for(read_set_id=0;read_set_id<number_of_read_sets;read_set_id++) fprintf(out, "Q%d_%d|",read_set_id+1,qual_2[read_set_id]);
        fprintf(out, "rank_%.5f",rank);
        fprintf(out, ";%s", results_against_set[cycle_id+1]->w);
        fprintf(out, ";");
        
        // PATH3
        fprintf(out, ">%s|",results_against_set[cycle_id+2]->comment);
        for(read_set_id=0;read_set_id<number_of_read_sets;read_set_id++) fprintf(out, "C%d_%d|",read_set_id+1,cov_3[read_set_id]);
        if (qual) for(read_set_id=0;read_set_id<number_of_read_sets;read_set_id++) fprintf(out, "Q%d_%d|",read_set_id+1,qual_3[read_set_id]);
        fprintf(out, "rank_%.5f",rank);
        fprintf(out, ";%s", results_against_set[cycle_id+2]->w);
        fprintf(out, ";");
        
        // PATH4
        fprintf(out, ">%s|",results_against_set[cycle_id+3]->comment);
        for(read_set_id=0;read_set_id<number_of_read_sets;read_set_id++) fprintf(out, "C%d_%d|",read_set_id+1,cov_4[read_set_id]);
        if (qual) for(read_set_id=0;read_set_id<number_of_read_sets;read_set_id++) fprintf(out, "Q%d_%d|",read_set_id+1,qual_4[read_set_id]);
        fprintf(out, "rank_%.5f",rank);
        fprintf(out, ";%s", results_against_set[cycle_id+3]->w);
        fprintf(out, "\n");
    }
    else
    {
        // PATH1
        fprintf(out, ">%s|", results_against_set[cycle_id]->comment);
        for(read_set_id=0;read_set_id<number_of_read_sets;read_set_id++) fprintf(out, "C%d_%d|",read_set_id+1,cov_1[read_set_id]);
        if (qual) for(read_set_id=0;read_set_id<number_of_read_sets;read_set_id++) fprintf(out, "Q%d_%d|",read_set_id+1,qual_1[read_set_id]);
        fprintf(out, "rank_%.5f\n",rank);
        fprintf(out, "%s\n", results_against_set[cycle_id]->w);
        
        // PATH2
        fprintf(out, ">%s|", results_against_set[cycle_id+1]->comment);
        for(read_set_id=0;read_set_id<number_of_read_sets;read_set_id++) fprintf(out, "C%d_%d|",read_set_id+1,cov_2[read_set_id]);
        if (qual) for(read_set_id=0;read_set_id<number_of_read_sets;read_set_id++) fprintf(out, "Q%d_%d|",read_set_id+1,qual_2[read_set_id]);
        fprintf(out, "rank_%.5f\n",rank);
        fprintf(out, "%s\n", results_against_set[cycle_id+1]->w);
        
        // PATH3
        fprintf(out, ">%s|", results_against_set[cycle_id+2]->comment);
        for(read_set_id=0;read_set_id<number_of_read_sets;read_set_id++) fprintf(out, "C%d_%d|",read_set_id+1,cov_3[read_set_id]);
        if (qual) for(read_set_id=0;read_set_id<number_of_read_sets;read_set_id++) fprintf(out, "Q%d_%d|",read_set_id+1,qual_3[read_set_id]);
        fprintf(out, "rank_%.5f\n",rank);
        fprintf(out, "%s\n", results_against_set[cycle_id+2]->w);
        
        // PATH4
        fprintf(out, ">%s|", results_against_set[cycle_id+3]->comment);
        for(read_set_id=0;read_set_id<number_of_read_sets;read_set_id++) fprintf(out, "C%d_%d|",read_set_id+1,cov_4[read_set_id]);
        if (qual) for(read_set_id=0;read_set_id<number_of_read_sets;read_set_id++) fprintf(out, "Q%d_%d|",read_set_id+1,qual_4[read_set_id]);
        fprintf(out, "rank_%.5f\n",rank);
        fprintf(out, "%s\n", results_against_set[cycle_id+3]->w);
    }
    
	
}

/**
 * checks if at least one read set provide read coherency for a path.
 */
inline int one_coherent(const p_fragment_info * results_against_set, int cycle_id, int number_of_read_sets){
	int i;
	for(i=0;i<number_of_read_sets;i++){
        if(results_against_set[cycle_id]->read_coherent[i]) return 1;
		
	}
	return 0;
}


void print_results_2_paths_per_event(FILE * coherent_out, FILE * uncoherent_out,  const p_fragment_info * results_against_set, const int number_of_read_sets, int nb_events_per_set, int qual, const char compute_genotype, const char paired){
    int i;
	int nb_read_coherent=0;
	int nb_unread_coherent=0;
	//printf("number ofread sets = %d\n", number_of_read_sets);
    
    //
    //                 C1           C2           C3 ....
    // path1 (i)      [0/1]        [0/1]        [0/1]...
    // path2 (i+1)    [0/1]        [0/1]        [0/1]...
    //
    // event is kept only if each line has at least one "1" per line:
    //
    
	
    
	for(i=0;i<nb_events_per_set*2;i+=2){
		if(one_coherent(results_against_set,i,number_of_read_sets) && one_coherent(results_against_set,i+1,number_of_read_sets))
		{
			nb_read_coherent++;
			print_couple_i("",coherent_out, results_against_set, i, number_of_read_sets, qual, 1, compute_genotype, paired);
		}
		else{
			nb_unread_coherent++;
			print_couple_i("", uncoherent_out, results_against_set, i, number_of_read_sets, qual, 1, compute_genotype, paired);
		}
	}
    
//	printf("Among %d bubbles:\n\t%d read coherent and\n\t%d not read coherent\n",
//           nb_events_per_set, nb_read_coherent, nb_unread_coherent);
	if (!silent) printf("Among %d bubbles: %d are read coherent\n", nb_events_per_set, nb_read_coherent);
}

//#define READ2INV

#ifdef READ2INV

void print_results_invs(FILE * coherent_out, FILE * uncoherent_out,  const p_fragment_info * results_against_set, const int number_of_read_sets, int nb_events_per_set, int qual){
    int i;
	int nb_read_coherent=0;
	int nb_unread_coherent=0;
    //                     C1           C2
    // au  path1 (i)      [0/1]        [0/1]
    // vb  path2 (i+1)    [0/1]        [0/1]
    // av' path2 (i+2)    [0/1]        [0/1]
    // u'b path2 (i+3)    [0/1]        [0/1]
    if(number_of_read_sets!=2){
        fprintf(stderr,"this kind of test is available only on exactly 2 datasets - please use 2 datasets or recompile kissreads after commenting ligne \"#define READ2INV\"\n");
        exit(1);
    }
    printf("\nOUTPUTS ONLY MOTIFS WHERE au-vb is specific to one datasets (non existing in the other) and av'-u'b is specific to the other \n");
    for(i=0;i<nb_events_per_set*4;i+=4){
        char coherent=0;
        if(results_against_set[i]->read_coherent[0] && results_against_set[i+1]->read_coherent[0] && //au and vb coherent in C1
           (!results_against_set[i+2]->read_coherent[0] || !results_against_set[i+3]->read_coherent[0]) && //av' or u'b uncoherent in C1
           (!results_against_set[i]->read_coherent[1] || !results_against_set[i+1]->read_coherent[1]) && //au or vb uncoherent in C2
           results_against_set[i+2]->read_coherent[1] && results_against_set[i+3]->read_coherent[1]) //av' and u'b coherent in C2
            coherent = 1;
        
        if(results_against_set[i]->read_coherent[1] && results_against_set[i+1]->read_coherent[1] && //au and vb coherent in C2
           (!results_against_set[i+2]->read_coherent[1] || !results_against_set[i+3]->read_coherent[1]) && //av' or u'b uncoherent in C2
           (!results_against_set[i]->read_coherent[0] || !results_against_set[i+1]->read_coherent[0]) && //au or vb uncoherent in C1
           results_against_set[i+2]->read_coherent[0] && results_against_set[i+3]->read_coherent[0]) //av' and u'b coherent in C1
            coherent=1;
        
        
        if(coherent){
			nb_read_coherent++;
			print_quadruplet_i(coherent_out, results_against_set, i, number_of_read_sets, qual);
		}
		else{
			nb_unread_coherent++;
			print_quadruplet_i(uncoherent_out, results_against_set, i, number_of_read_sets, qual);
		}
	}
    
	if (!silent)  printf("Among %d inversions:\n\t%d read coherent and\n\t%d not read coherent\n",
           nb_events_per_set, nb_read_coherent, nb_unread_coherent);
    
    
}
#else // !READ2INV

void print_results_invs(FILE * coherent_out, FILE * uncoherent_out,  const p_fragment_info * results_against_set, const int number_of_read_sets, int nb_events_per_set, int qual){
    int i;
	int nb_read_coherent=0;
	int nb_unread_coherent=0;
	//printf("number ofread sets = %d\n", number_of_read_sets);
    
    //
    //                 C1           C2           C3 ....
    // path1 (i)      [0/1]        [0/1]        [0/1]...
    // path2 (i+1)    [0/1]        [0/1]        [0/1]...
    // path2 (i+2)    [0/1]        [0/1]        [0/1]...
    // path2 (i+3)    [0/1]        [0/1]        [0/1]...
    //
    // event is kept only if each line has at least one "1" per line:
    //
    
	
    
	for(i=0;i<nb_events_per_set*4;i+=4){
		if(one_coherent(results_against_set,i,number_of_read_sets) &&
           one_coherent(results_against_set,i+1,number_of_read_sets) &&
           one_coherent(results_against_set,i+2,number_of_read_sets) &&
           one_coherent(results_against_set,i+3,number_of_read_sets)
           )
		{
			nb_read_coherent++;
			print_quadruplet_i(coherent_out, results_against_set, i, number_of_read_sets, qual);
		}
		else{
			nb_unread_coherent++;
			print_quadruplet_i(uncoherent_out, results_against_set, i, number_of_read_sets, qual);
		}
	}
    
	if (!silent) printf("Among %d inversions:\n\t%d read coherent and\n\t%d not read coherent\n",
           nb_events_per_set, nb_read_coherent, nb_unread_coherent);
}
#endif // !READ2INV



/**
 * prints a sequence
 */
void print_sequence_i(FILE* out, const p_fragment_info * results_against_set, int cycle_id, int number_of_read_sets, int qual){
	
	int sum[number_of_read_sets];
	int avg[number_of_read_sets];
    int read_set_id;
    
//	if( qual ){
//        // we are providing results for generic dataset
//        for(read_set_id=0;read_set_id<number_of_read_sets;read_set_id++){
//            avg[read_set_id] = 0;
//            for (j=kmer_size-1;j<=strlen(results_against_set[cycle_id]->w)-kmer_size;j++){
//                if(results_against_set[cycle_id]->read_coherent_positions[read_set_id][j]){
//                    avg[read_set_id] = avg_up[read_set_id] + (results_against_set[cycle_id]->sum_quality_per_position[read_set_id][j] / results_against_set[cycle_id]->read_coherent_positions[read_set_id][j]);
//                }
//            }
//        }
//    }
	
	for(read_set_id=0;read_set_id<number_of_read_sets;read_set_id++){
		sum[read_set_id]=results_against_set[cycle_id]->number_mapped_reads[read_set_id];
	}
    
    fprintf(out, ">%s|", results_against_set[cycle_id]->comment);
    for(read_set_id=0;read_set_id<number_of_read_sets;read_set_id++)
        fprintf(out, "C%d_%d|",read_set_id+1,sum[read_set_id]);
    if (qual)
        for(read_set_id=0;read_set_id<number_of_read_sets;read_set_id++)
            fprintf(out, "Q%d_%d|",read_set_id+1,avg[read_set_id]);
    fprintf(out, "%s%s\n", standard_fasta?"\n":";",results_against_set[cycle_id]->w);
    
    
	
}



void print_generic_results(FILE * coherent_out, FILE * uncoherent_out,  const p_fragment_info * results_against_set, const int number_of_read_sets, int nb_events_per_set, int qual){
    
    int i;
	int nb_read_coherent=0;
	int nb_unread_coherent=0;
    for(i=0;i<nb_events_per_set;i++){
        if(one_coherent(results_against_set,i,number_of_read_sets)){
            nb_read_coherent++;
            print_sequence_i(coherent_out, results_against_set, i, number_of_read_sets, qual);
        }
        else{
            nb_unread_coherent++;
            print_sequence_i(uncoherent_out, results_against_set, i, number_of_read_sets, qual);
        }
    }
    if (!silent)  printf("Among %d sequences:\n\t%d read coherent and\n\t%d not read coherent\n",
           nb_events_per_set, nb_read_coherent, nb_unread_coherent);
}








