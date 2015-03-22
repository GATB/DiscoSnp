/* This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software.  You can  use, 
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info". 
 */

/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file gapless.cpp
 * \brief Module Gapless, responsable de la réalisation de l'alignement sans gap
 * \author Dominique Lavenier
 * \author Damien Fleury
 * \version 5.3
 * \date 28/08/2008
 */

#include <emmintrin.h>
//#include <tmmintrin.h>
#include "gapless.h"
#define _mm_extract_epi32(x, imm) \
_mm_cvtsi128_si32(_mm_srli_si128((x), 4 * (imm)))

#define	Mprint(vec,mess)  printf("VEC %s : %x %x %x %x \n",mess,_mm_extract_epi32(vec,0), _mm_extract_epi32(vec,1), _mm_extract_epi32(vec,2), _mm_extract_epi32(vec,3))

#include "constants.h"
#include "code.h"
#include "Bank.h"
#include "Index.h"
#include "Hit.h"
#include "display.h"
#include "Stat.h"
#include "Doublon.h"
#include "Alignment.h"

 pthread_rwlock_t verrou_out;

using namespace std;





inline  int pre_filtre_align_gless( int left2,int right2, int max_mis, char ** tabprec ,unsigned int fleft1,unsigned int fright1
		       ,unsigned int fleft2,unsigned int fright2)
{
  int tai,errtot=0;
  tai = min(TAI_TAB,left2);
  int decalage =0;
  // A gauche
  if(tai>=1)
    {
      decalage = 2*(4*sizeof(unsigned int) -tai);
      errtot = (int) tabprec[fleft1>>decalage][fleft2>>decalage] ;
    }
  if (errtot> max_mis) return 0;


  // A droite
  tai = min(TAI_TAB,right2);
  if ( tai>=1) 
    {
      decalage = 2*(4*sizeof(unsigned int) -tai);
      errtot += (int) tabprec[fright1>>decalage][fright2>>decalage] ;
    }
  if (errtot> max_mis) return 0;

  return 1;
}



inline unsigned int popcount_mult(unsigned int x)
{
    unsigned int m1 = 0x55555555;
    unsigned int m2 = 0x33333333;
    unsigned int m4 = 0x0f0f0f0f;
    unsigned int h01 = 0x01010101;
    x -= (x >> 1) & m1;               /* put count of each 2 bits into those 2 bits */
    x = (x & m2) + ((x >> 2) & m2);   /* put count of each 4 bits in */ 
    x = (x + (x >> 4)) & m4;          /* put count of each 8 bits in partie droite  4bit piece*/
    return (x * h01) >> 24;           /* returns left 8 bits of x + (x<<8) + ... */ 
}


// x : sequence dentrre, nt : lettre a compter, repetee sur 32 bits
//renvoit tableau de bits avec des 1 la où lettres trouvees, reste a faire poopcount
inline unsigned int count_nt(unsigned int x, unsigned int nt)
{
  unsigned int m= 0x55555555;

  x = ~ (x ^ nt) ;
  x = x & (x >> 1);
  return (x & m);

}




inline int filtre_align_NT_vec_gless(int left2, int right2, int max_mis 
			,unsigned int fleft1,unsigned int fright1
			,unsigned int fleft2,unsigned int fright2)
{
  int min_err=0;
  unsigned int masque ;
  int tairest;
  
  //cote gauche
  tairest =  min(left2,16); //lecture 16 max au total

  masque = ~  (0xffffffff >> tairest >> tairest); 
  
  min_err = tairest - (int)popcount_mult (count_nt(fleft2,fleft1) & masque );

  if (min_err> max_mis) return 0;

  //----------------       cote droit 
  tairest =  min(right2,16);

  masque = ~  (0xffffffff >> tairest >> tairest); 
  
  min_err += tairest -  (int)popcount_mult (count_nt(fright2,fright1) & masque );

  if (min_err> max_mis) return 0;

  return 1;

}

/*
 
inline unsigned int get_code(unsigned int *  vec,int tai_acc,int tai)
{

  int dec = tai_acc/10; //cahque uint sert a lecture de plage de 10 nt
  unsigned int temp = vec[dec];
  tai_acc -= 10*dec; //4*dec
  int    decalage = 2*(4*sizeof(unsigned int) -tai - tai_acc);
  unsigned int masque = (1 << (2*tai)) -1 ;   //masque de 2*tai bit à 1 
  return ((temp>>decalage) &masque);


}

//avec petite table 
int filtre_align_DC_gless(char* seq1, char* seq2, int left2, int right2, int max_mis,char ** tabprec
		    ,unsigned int fleft1,unsigned int fright1
		     ,unsigned int fleft2,unsigned int fright2
		    )
{
  int tai;
  int errtot=0;
  int errd=0;
  int errg=0;
  int tairest ;
  int tailue;
  int decalage =0;

  //
  int nbu;
  // unsigned int * bleft1,* bleft2, *bright1, *bright2 ; //codage en binaire des abords, jusqau bout query , max 64 
  unsigned int * bleft1;
  unsigned int * bleft2;
  unsigned int * bright1;
  unsigned int * bright2;


  //codage cote gauche
      nbu = (left2) / 10;
      nbu += ((left2 -( nbu*10))!=0); //nombre de uint necessaires au stockage chavauchant
      bleft1 = new unsigned int [nbu];
      bleft2 = new unsigned int [nbu];
      bleft1[0]=fleft1; //les 2 premiers sont dans l'index
      bleft2[0]=fleft2;

      tairest = left2 -16; //16 deja lus dans lindex
      tailue=10;//chavauchant donc tailue 10 seulement
      tai = min(16,tairest + 6 );
      nbu=1;
      while (tairest>0)
	{
	  bleft1[nbu]=seq2code_rev(seq1-1-tailue,tai)<<(2*(16-tai));
	  bleft2[nbu]=seq2code_rev(seq2-1-tailue,tai)<<(2*(16-tai));
	  nbu++;
	  
	  tairest -= tai -6  ; // on en a lu tairest -6 nouveaux
	  tailue += tai - 6 ; 
	  tai = min(16,tairest + 6);
	}

      //codage cote droit
      nbu = (right2) / 10;
      nbu += ((right2 -( nbu*10))!=0); 
      bright1 = new unsigned int [nbu];
      bright2 = new unsigned int [nbu];
      bright1[0]=fright1; 
      bright2[0]=fright2;

      tairest = right2 -16; 
      tailue=10;
      tai = min(16,tairest + 6 );
      nbu=1;
      while (tairest>0)
	{
	  bright1[nbu]=seq2code(seq1+SIZE_SEED+tailue,tai)<<(2*(16-tai));
	  bright2[nbu]=seq2code(seq2+SIZE_SEED+tailue,tai)<<(2*(16-tai));
	  nbu++;
	  
	  tairest -= tai - 6 ;
	  tailue += tai - 6 ;
	  tai = min(16,tairest + 6);
	}
  // le codage est fini

   tairest = left2;
   tai = min(TAI_TAB,tairest);

  // A gauche
  if(tai>=1)
    {
      decalage = 2*(4*sizeof(unsigned int) -tai);
      errg = (int) tabprec[fleft1>>decalage][fleft2>>decalage] ;
      errtot += errg;
    }

  if (errtot> max_mis) { delete [] bright1;delete [] bright2;delete [] bleft1;delete [] bleft2;  return errtot;}
  // A droite
  tairest = right2;
  tai = min(TAI_TAB,tairest);

  int tai_acc=0;
  if ( tai>=1) 
    {
       decalage = 2*(4*sizeof(unsigned int) -tai);
       errd += (int) tabprec[fright1>>decalage][fright2>>decalage] ;
      //errd += (int) tabprec[seq2code(seq1 + SIZE_SEED,tai)][seq2code(seq2 + SIZE_SEED,tai)] ;
      errtot += errd;
    }
  if (errtot> max_mis) { delete [] bright1;delete [] bright2;delete [] bleft1;delete [] bleft2;  return errtot;}

	

  /////////// filtre phase 2 : apres les TAI_BIG_TAB premieres bases on doit tester toutes les possibilites	
  //avec boucle partie droite
  tairest-= tai;
  tai_acc+= tai;
  tai = min(TAI_TAB,tairest);
  int errdprec,errgprec;
  while (tai_acc < 64 && tairest>0)
    {
      errdprec = errd;
  
      errd =  errdprec +  (int) tabprec[get_code(bright1,tai_acc,tai)][get_code(bright2,tai_acc,tai)] ;

      errtot = errd +errg;
  if (errtot> max_mis) { delete [] bright1;delete [] bright2;delete [] bleft1;delete [] bleft2;  return errtot ;}

      tairest-= tai;
      tai_acc+= tai;
      tai = min(TAI_TAB,tairest);
	    
    }

  //  puis suite gauche  avec boucle 
  tai_acc = 0;
  tairest = left2;
  tai = min(TAI_TAB,tairest);
  // on reprend la ou on en etait
  tairest-= tai;
  tai_acc+= tai;
  tai = min(TAI_TAB,tairest);


while (tai_acc < 64 && tairest>0)

    {
      errgprec = errg;
      errg =  errgprec + tabprec[get_code(bleft1,tai_acc,tai)][get_code(bleft2,tai_acc,tai)] ;


      errtot = errg+ errd;
      if (errtot> max_mis) { delete [] bright1;delete [] bright2;delete [] bleft1;delete [] bleft2;  return errtot;}

      tairest-= tai;
      tai_acc+= tai;
      tai = min(TAI_TAB,tairest);
	    
    }

delete [] bright1;delete [] bright2;delete [] bleft1;delete [] bleft2; 
  return errtot;

}




*/


/**
 * Fonction permettant de calculer un score d'alignement obtenu entre 2 séquences
 * \param seq1 le pointeur du tableau des données de la première banque
 * \param seq2 le pointeur du tableau des données de la seconde banque
 * \param left le nombre de bases à parcourir à gauche
 * \param right le nombre de bases à parcourir à droite
 * \param num_seed le code de la graine précédente
 * \param max_mis le nombre maximal de mésappariement pour conserver l'alignement
 * \return Le nombre de mésappariements obtenus dans l'alignement
 */
int scoreHit(char* seq1, char* seq2, int left, int right, int num_seed, int max_mis)
{
	char* s1 = seq1-1;
	char* s2 = seq2-1;
	int k = 0;
	int nbmis = 0;

  	/// On parcourt les nucléotides à gauche
	while((nbmis<=max_mis)&&(k<left))
	{
      	/// On vérifie l'appariement
		if(!identNT(*s1,*s2))
		{
			++nbmis;
		}
	
		/// On passe aux caractères suivants
		--s1;
		--s2;
		++k;
	}
	

	s1 = seq1 + SIZE_SEED;
	s2 = seq2 + SIZE_SEED;
	k = 0;

  	/// On parcourt les nucléotides à droite
	while((nbmis<=max_mis)&&(k<right))
	{
      	/// On vérifie l'appariement
		if(!identNT(*s1,*s2))
		{
					++nbmis;
		}
	
		/// On passe aux caractères suivants
		++s1;
		++s2;
		++k; 
	}
  	/// On rend le nombre de mésappariements
	return nbmis;
}


/**
 * Méthode qui à partir des positions d'une graine va effectuer tous les alignements possibles avec les séquences des banques
 * \param ff le fichier de sortie, où seront affichés les alignements obtenus
 * \param BK1 le pointeur de la première banque de séquences
 * \param BK2 le pointeur de la seconde banque de séquences
 * \param I1 le pointeur de l'index de la première banque de séquences
 * \param num_seed l'indice de la graine utilisée pour déterminer les alignements
 * \param idpc le pourcentage de ressemblances minimal entre les 2 séquences pour qu'un alignement soit conservé
 * \return le nombre d'alignements obtenus
 */
/*
int gapless_align(FILE *ff, Bank *BK1, Bank *BK2, Index *I1, int num_seed, int idpc, bool rev_comp,Stat * St)

{
	int  idx1, idx2;
	int  nb_hit;
	int  nb_item;
	int  i;
	int  left1, left2, right1, right2;
	int  max_mis, nb_mis, size2, score;
	int nb_occur = 0;
	int n = 0;
	double e_value = 0;
	int num_sequence,offset_sequence;
	/// Objet Hit permettant de mémoriser un appariement
	Hit* H1;
	H1 = new Hit();
	/// Tableau permettant de stocker les positions des graines
	Hit* H2[NB_MAX_HIT];
	
	char *s1, *s2, *j1, *j2;

	nb_item = 0;
  	/// Pour la seconde banque :
  	
  	/// On récupère l'index de la graine dans le tableau seed
	idx2 = I2->offset_seed[num_seed];

  	/// On récupère le nombre d'occurences de la graine
	nb_occur = I2->nb_seed[num_seed];

  	/// On parcourt toutes les occurences de la graine
	for(n=0; n<nb_occur; ++n)
	{
 		/// On récupère les informations sur la position de la graine
		H2[nb_item] = new Hit(BK2, (I2->seed[idx2]).num_seq, (I2->seed[idx2]).off_seq);
      	/// On passe à l'occurence suivante
		idx2 += 1;
		++nb_item;
	}
  	/// On arête si on n'a aucune occurence
	if (nb_item == 0) return 0;
	nb_hit = 0;

  	/// Pour la première banque :
  	/// On récupère l'index de la graine dans le tableau seed
	idx1 = I1->offset_seed[num_seed];

  	/// On récupère le nombre d'occurences de la graine
	nb_occur = I1->nb_seed[num_seed];
  	/// On parcourt toutes les occurences de la graine
	for(n=0; n<nb_occur; ++n)
	{
      	/// On récupère les informations sur la position de la graine
	        num_sequence = (I1->seed[idx1]).num_seq;
	        offset_sequence =  (I1->seed[idx1]).off_seq ;
		H1->offhit = (BK1->seq[num_sequence] + offset_sequence) ;
		H1->offseq = BK1->seq[num_sequence];
		H1->sizeseq = BK1->size[num_sequence];
		H1->numseq = num_sequence;

		left1     = H1->offhit - H1->offseq; // ou directement (I1->seed[idx1]).off_seq
		right1    = H1->offseq + H1->sizeseq - H1->offhit  - SIZE_SEED; 
		/// On extrait les informations pour la séquence de la première banque
		s1 = (BK1->data) + BK1->seq[H1->numseq];

		for (i=0; i<nb_item; ++i)
		{


	  		/// On extrait les informations pour la séquence de la seconde banque
			left2   = H2[i]->offhit - H2[i]->offseq;
			right2  = H2[i]->offseq + H2[i]->sizeseq - H2[i]->offhit - SIZE_SEED;

	  		/// On vérifie que la seconde séquence (la plus petite) va pouvoir s'aligner entièrement avec la première
	  		/// (Il faut que la seconde séquence soit incluse dans la première)
			if ((left1>=left2)&&(right1>=right2))
			{
				s2 = (BK2->data) + BK2->seq[H2[i]->numseq];
				size2   = H2[i]->sizeseq;
	      		/// On calcule le nombre de mésappariements maximal autorisé afin de satisfaire le pourcentage de différences
	      		max_mis = (size2*(100-idpc)) / 100;
	      		
	      		/// On calcule le score de l'alignement
			nb_mis  = scoreHit(s1 + left1,s2 + left2,left2,right2,num_seed,max_mis);
	      		/// Si le score d'alignement est suffisant, on enregistre l'alignement
				if (nb_mis <=max_mis)
				{
		  			/// On récupère les index des commentaires des séquences
					j1 = (BK1->data) + BK1->com[H1->numseq];
					j2 = (BK2->data) + BK2->com[H2[i]->numseq];
		  			/// Calcul du score
					score = (size2-nb_mis)*MATCH - nb_mis*MISMATCH;
					/// Calcul de la E-value
					e_value = ComputeEvalue(score,BK1->nb_tot_res,BK1->nb_tot_seq,H2[i]->sizeseq,1);
		  			/// Affichage de l'alignement dans le fichier de sortie
					
					/// Verrou pour protéger l'accès au fichier de sortie
		  			pthread_rwlock_wrlock(&verrou_out);
					display_ungap(ff,s1,s2,j1,j2,H1->offhit-left2-H1->offseq,H2[i]->offhit-left2-H2[i]->offseq,size2,H1->sizeseq,e_value,rev_comp);
					pthread_rwlock_unlock(&verrou_out);
					
					++nb_hit;
				}
			}
		}
      	/// On passe à l'occurence suivante pour la première séquence
		idx1 += 1;
	}
	for(i=0;i<nb_item;++i)
	{
		delete H2[i];
	}
	delete H1;
	return nb_hit;
}


*/




#ifdef CPT_FILT

extern long long int cpt_hit;
extern long long int cpt_pref;
extern long long int cpt_ntsse;
extern long long int cpt_fil;
extern long long int cpt_doub;
extern long long int cpt_max;
extern long long int cpt_dint;

#endif





int gapless_align(FILE *ff, Bank *BK1, Bank *BK2, Index * I1, int idpc, bool rev_comp, int num_gaps,char ** tabprec,char ** tabnt,char **tabprec7,Stat * St,Doublon * Doub, int start, int end )
{
  int  idx1;
  int  nb_hit;
  int  i;
  int  left1, left2, right1, right2;
  int  max_mis, nb_mis,  score;
  int nb_occur = 0;
  int n = 0;
  double e_value = 0;

  //fin test sans sse 
  int nquery;
  int size_query;
  int num_sequence,offset_sequence,offseq,sizeseq;
  int offseq2;
  unsigned int fleft=0; // pour filtre
  unsigned int fright=0; 

  unsigned int fleftB,frightB;
   int tair=0;
  nb_hit = 0;
  int num_seed;
  char *s1, *s2, *j1, *j2;
  Alignment* al;

  Seed *  seed_bank;
  al =  new Alignment(BK2->tailleMaxSeq + (int)((BK2->tailleMaxSeq*(100-idpc)) / 100) + 1);
  Hit* H;
  H = new Hit();
  // list<Hit> L;
  int cpt_pos=0;
  int cpt_query;
  unsigned int mc= 0x55555555;
  unsigned int m2 = 0x33333333;
  unsigned int m4 = 0x0f0f0f0f;
  unsigned int h01 = 0x01010101;
	
  __m128i vec_nt;
  __m128i vec_5,vec_3,vec_0f,vec_h1; //3masques pour popcnt
  //init des vecteurs masques 
  vec_5  = _mm_set_epi32(mc, mc, mc, mc);
  vec_3  = _mm_set_epi32(m2, m2, m2, m2);
  vec_0f  = _mm_set_epi32(m4, m4, m4, m4);
  vec_h1  = _mm_set_epi32(h01, h01, h01, h01);
  vec_nt = _mm_set_epi32(0x00000000, 0xaaaaaaaa, 0x55555555, 0xffffffff); // reg 3 2 1 0
	
  //pour chaque query
  for (nquery = start; nquery <= end; nquery++) 
    {
      cpt_query=0;
      // nouvelle liste pour chaque query
      //    L.clear();
      size_query =   BK2->size[nquery] ;  
      num_seed=-1;
      offseq2 = BK2->seq[nquery];
      s2 = (BK2->data) + BK2->seq[nquery];
      ////// printf("traitement query %i  %i--%i  tai %i \n",nquery,start,end,size_query);
      // printf("tai %i  %s\n",size_query,s2);
      fleft=0;fright=0;

      max_mis = (size_query*(100-idpc)) / 100;

      //pour chaque graine dans query
      for(i=0; i<= (size_query - SIZE_SEED); i++)
	{
	  cpt_pos=0;
	  
	  //pour seq B2
	  /// offset_sequence2 =  i ;
	  //	  offhit2 = offseq2 + i ;
	  left2 =  i;
	  right2 =  size_query - i - SIZE_SEED;



	  //	printf("i %i limit %i \n",i,(size_query - SIZE_SEED));	
	  //i-> nquery
	  //j-> i

	  //lecture flanking regions pour filtre
	  if (i==0) // premiere graine, rien a gauche
	    {
	      tair=min(4*sizeof(unsigned int),size_query-i-SIZE_SEED);  
	      fright =  seq2code_special(s2 + i + SIZE_SEED,tair,4*sizeof(unsigned int)) ; // 4 charac par octets4*sizeof(unsigned int)
	    }
	  else if (i<(size_query-SIZE_SEED-4*(int)sizeof(unsigned int)+1)) // il reste des trucs a lire a droite
	    {
	      fleft = (unsigned int)seq2codeLeft(s2+i-1,fleft);
	      fright = (unsigned int) seq2codeRight(s2+i+SIZE_SEED,4*sizeof(unsigned int),fright);  
	    }
	  else
	    {
	      fleft = (unsigned int)seq2codeLeft(s2+i-1,fleft);
	      fright = (fright <<2); 
	    }
	 

	  num_seed = codeSeedRight(&BK2->data[BK2->seq[nquery]+i],num_seed);
	  /// On parcourt toutes les occurences de la graine ds B1
	  //	if(num_seed<0) printf("heho!\n");
	  if(num_seed!=-1)
	    {
	 
	      if(SIZE_SEED>14) 
		{
		  I1->get_seed_info_through_hashTable(num_seed,&nb_occur,&idx1);
		}
	      else
		{
		  nb_occur = I1->nb_seed[num_seed]; 
		  idx1 = I1->offset_seed[num_seed];
		}
	
	      // pour chaque occurence de la graine dans B1  
	      for(n=0; n<nb_occur; n++)
		{
		  seed_bank = &(I1->seed[idx1]);
		  //pour seq B1
		  num_sequence = seed_bank->num_seq;
		  offset_sequence =  seed_bank->off_seq ;
		  offseq = BK1->seq[num_sequence];
	
		  //	  offhit = offseq + offset_sequence ;
		  sizeseq = BK1->size[num_sequence];

		   left1 = offset_sequence;
		   right1 =  sizeseq - offset_sequence - SIZE_SEED;
		   s1 = (BK1->data) + BK1->seq[num_sequence];


		  if ((left1>=left2)&&(right1>=right2))
		    {
		      fleftB=seed_bank->left;
		      frightB=seed_bank->right;


#ifdef CPT_FILT

		      cpt_hit++;
#endif
		               if(pre_filtre_align_gless( left2,right2,max_mis,tabprec,fleftB,frightB,fleft,fright))
			{
#ifdef CPT_FILT

			  cpt_pref++;
#endif

			    if(filtre_align_NT_vec_gless(left2,right2,max_mis 
			  		 ,fleft,fright,fleftB,frightB))
			    {
#ifdef CPT_FILT

			      cpt_ntsse++;
#endif


				{

				  if(! (Doub->getpos(offset_sequence-i)))
				    {
#ifdef CPT_FILT
				      cpt_doub++;
#endif
				      if ((cpt_pos>MAXPOS) && MAXPOS)
					{
#ifdef CPT_FILT
					  cpt_max++;
#endif
					  goto pos_suivante; //MAXPOS==0 nolimit
					}


				      //   if(filtre_align_DC_gless(s1+left1,s2+left2,left2,right2,max_mis,tabprec,
				      //	 fleftB,frightB,fleft,fright ))
					  { //dans ce cas le nombre d'err est inferieur a max_mis


					  //	nb_mis  = filtre_align_DC_gless(s1+left1,s2+left2,left2,right2,max_mis,tabprec,
					  //  fleftB,frightB,fleft,fright );

					  nb_mis  = scoreHit(s1 + left1,s2 + left2,left2,right2,0,max_mis);


	                 if (nb_mis <=max_mis)
	                            {

#ifdef CPT_FILT

					  cpt_dint++;
#endif
					      //enregistrement dans stat
					      St->incr_query();
					      //enregistrement dans doublons
					      Doub->setpos(offset_sequence-i);
					      /// On récupère les index des commentaires des séquences
					      j1 = (BK1->data) + BK1->com[num_sequence];
					      j2 = (BK2->data) + BK2->com[nquery];
					      /// Calcul du score
					      score = (size_query-nb_mis)*MATCH - nb_mis*MISMATCH;
					      /// Calcul de la E-value
					      e_value = ComputeEvalue(score,BK1->nb_tot_res,BK1->nb_tot_seq,size_query,1);
					      /// Affichage de l'alignement dans le fichier de sortie
		    
					      /// Verrou pour protéger l'accès au fichier de sortie
					      pthread_rwlock_wrlock(&verrou_out);
					      // display_withgap(ff,al,j1,j2,e_value);
					      //     display_ungap(ff,s1,s2,j1,j2,H1->offhit-left2-H1->offseq,H2[i]->offhit-left2-H2[i]->offseq,size2,H1->sizeseq,e_value,rev_comp);
					      //A voir
					      display_ungap(ff,s1,s2,j1,j2,offset_sequence-left2,i-left2,size_query,sizeseq,e_value,rev_comp);

					      pthread_rwlock_unlock(&verrou_out);
		    
					      ++nb_hit;
					        if (St->quota_atteint())  goto query_suivante;
				}
					} //fin filtre DC
						cpt_pos++;
				    }
				}
			    }
			}

		      //	} // fin if doub
		    } //fin test s2 incluse ds s1

		  idx1++;
		} //fin boucle occurences dun graine dans B1




	    } // fin test graine valide (passe filtre complexite)

	pos_suivante:;
	} //fin boucle  graines dans query
  
     

      //remise a 0 du tab des doublons
    query_suivante: ;
      Doub->reset();
      St->reset();

    } //fin boucle sur les  query


  delete al;
  return nb_hit;
}





