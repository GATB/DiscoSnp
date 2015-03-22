
/* This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software.  You can  use, 
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info". 
 */

/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file withgap.cpp
 * \brief Module Withgap, responsable de la réalisation de l'alignement avec gap
 * \author Damien Fleury
 * \author Guillaume Rizk
 * \version 6.1
 * \date 21/01/2009
 */
#include "withgap.h"

#include "constants.h"
#include "code.h"
#include "Bank.h"
#include "Index.h"
#include "Hit.h"
#include "display.h"
#include "Alignment.h"
#include "gapless.h"
#include "Stat.h"
#include "Doublon.h"



#if SSE==2
#include <emmintrin.h> 
#include "emul_sse2.h"
#endif

#if SSE==3
#include <emmintrin.h>
#include <tmmintrin.h>
#endif



#include "misc.h"


#define _mm_extract_epi32(x, imm) \
_mm_cvtsi128_si32(_mm_srli_si128((x), 4 * (imm)))

#define	Mprint(vec,mess)  printf("VEC %s : %x %x %x %x \n",mess,_mm_extract_epi32(vec,0), _mm_extract_epi32(vec,1), _mm_extract_epi32(vec,2), _mm_extract_epi32(vec,3))

#include <list>

using namespace std;

/** 
 * Fonction de recherche du meilleur alignement selon le modèle de Wunsch et Needleman
 * \param al l'objet Alignment qui constituera le résultat
 * \param seq1 le pointeur du tableau des caractères de la première banque
 * \param seq2 le pointeur du tableau des caractères de la seconde banque
 * \param size1 le nombre de bases à parcourir dans la première séquence
 * \param size2 le nombre de bases à parcourir dans la seconde séquence
 * \param num_gaps le nombre maximal de gaps autorisés dans la partie de l'alignement
 * \param side indiquant si on fait l'alignement sur le côté gauche ou droit de la graine
 * \param num_seed le code de la graine courante
 * \param max_mis le nombre max de mismatch
 * \return AL_VALID si l'alignement trouvé est correct, AL_ERR si on n'a pas retenu l'alignement
 */
 

int WN_alignment(Alignment* al, char* seq1, char* seq2, int size1, int size2, int num_gaps, int side, int num_seed, int max_mis)
{
  int i,j;

  int binf, bsup;
  int barriere;
  char* s1;
  char* s2;
  /// Les dimensions de la matrice d'alignement
  int dim1 = size1+1;
  int dim2 = size2+1;
  /// La matrice d'alignement, la première séquence constitue la première dimension
  caseMatrix tab[dim1][dim2];


  /// Les tableaux de suivi de la construction de l'alignement
  char mem[size2+num_gaps];
  char mem2[size2+num_gaps];
	
  int cpt = 0;
	
  /// Variables utiles au contrôle de la validité de l'alignement
  int lastCode = num_seed;
  int L = SIZE_SEED;
  int max_gaps = num_gaps;


  if(side == LEFT)
    {

      /// On se place au début des séquences
      s1 = seq1 - size1;
      s2 = seq2 - size2;
    }
  else
    {

      /// Dans ce cas on inverse les séquences
      //  printf("seq1 %s seq2 %s \n",seq1,seq2);

      //   printf("WN align size1 %i  size2 %i ngaps %i side %i numseed %i  max_mis %i \n",size1,size2,num_gaps,side,num_seed,max_mis);
      s1 = new char[size1];
      s2 = new char[size2];
		
      for(i=0; i<size1; ++i)
	{
	  s1[i] = seq1[size1-i-1];
	}
      for(j=0; j<size2; ++j)
	{
	  s2[j] = seq2[size2-j-1];
	}
    }


  
  /// Remplissage de la matrice d'alignement

  //	Initialisation des cases de la matrice
  //trop long, et inutile 
//   for(i = 0; i < dim1; ++i)
//     {
//       for(j = 0; j < dim2; ++j)
// 	{
// 	  tab[i][j].score = 0;
// 	  tab[i][j].path = 99;
// 	}
//     }





  /// On initialise la première case de la matrice
  tab[0][0].score = 0;
  // On initialise les cases utiles de la première colonne
  // Ppour première colonne : gaps gratuits car differentes positions de depart possibles
  for(i=1; i <= min(num_gaps + num_gaps,size1); ++i) //bande de 2*numgaps de large
    {
      tab[i][0].score = 0;
      tab[i][0].path = GA2;
    }
  /// On initialise les cases utiles de la première ligne
  for(j=1; j <= min(num_gaps,size2); ++j)
    {
      tab[0][j].score = j * COUT_GAP;
      tab[0][j].path = GA1;
    }

  /// On remplit le reste de la matrice d'alignement
  /// On ne complète que les cases des diagonales respectant le nombre de gaps autorisés
  for(i = 1; i < dim1; ++i)
    {
      binf = max(i-num_gaps -num_gaps,1); //ok car on commence sur diag dacalee
      bsup = min(i,dim2-1); 
      barriere = 0;
      for(j = binf; j <= bsup; ++j)
	{

	  tab[i][j].score = -10000 ; 
	  tab[i][j].path = MES;
	  /// On vérifie l'appariement des caractères
	  if(identNT(s1[i-1],s2[j-1])
	     && tab[i-1][j-1].score >= -max_mis )
	    {
	      tab[i][j].score =  tab[i-1][j-1].score  + COUT_MATCH;
	      tab[i][j].path = APP;
	      barriere++;
	    }
	  else
	    if (   tab[i-1][j-1].score > -max_mis )

	      {
		tab[i][j].score =  tab[i-1][j-1].score  + COUT_MISMATCH;
		tab[i][j].path = MES;
		barriere++;
	      }
			
	  /// On teste la possibilité d'ajouter un gap
	  if(( tab[i-1][j].score + COUT_GAP > tab[i][j].score) && 
	     ((j-i) < 0) // test depassement bande
	     && tab[i-1][j].score > -max_mis )
	    {
	      tab[i][j].score =  tab[i-1][j].score + COUT_GAP;
	      tab[i][j].path = GA2;
	      barriere++;
	    }
	  if(( tab[i][j-1].score + COUT_GAP > tab[i][j].score) && 
	     ((i-j) < (2*num_gaps)) // depassement bande
	     && tab[i][j-1].score > -max_mis )
	    {
	      tab[i][j].score = tab[i][j-1].score  + COUT_GAP;
	      tab[i][j].path = GA1;
	      barriere++;
	    }
	}
      if(!barriere) { // tous les chemins menent a nombre d'erreurs > max_mis, on peut s'arreter
	if(side == RIGHT)
	  {
	    delete [] s1;
	    delete [] s2;
	  }


	return AL_ERR;

      }
    }


  --i;
  --j;
  /// On parcourt le chemin représentant le meilleur alignement dans la matrice
  /// On commence dans le coin en bas à droite de la matrice, correspondant aux derniers caractères , au niveau de la graine 
  /// On va continuer jusqu'à ce que la séquence la plus courte (la seconde) ait été complètement lue
  while(j>0)
    {

      switch(tab[i][j].path)
	{
	case APP:
	  --i;
	  --j;
	  mem[cpt]= s1[i];
	  mem2[cpt]= s2[j];
		  			
	  break;
	case MES:
	  --i;
	  --j;
	  mem[cpt]= s1[i];
	  mem2[cpt]= s2[j];
	  al->addMis();
	  L = 0;
	  lastCode = 0;
	  break;
	case GA2:
	  // gap sur query
	  if(side == LEFT) al->decStart1(); else al->incEnd1(); //ok
	  --i;
	  mem[cpt]= s1[i];
	  mem2[cpt]= CHAR_GAP;
	  al->addGap();
	  //  al->incAlign();
	  // al->decAlign();			
	  L = 0;
	  lastCode = 0;
				
	  break;
	case GA1:
	  // gap sur bank
	  if(side == LEFT) al->incStart1(); else al->decEnd1(); //ok
	  --j;
	  mem[cpt]= CHAR_GAP;
	  mem2[cpt]= s2[j];
	  al->addGap();
	  //  al->decAlign();
	  // al->incAlign();			
	  L = 0;
	  lastCode = 0;
				
	  break;
	}
		
      /// Si on a déjà utilisé trop de gaps, on arrête la construction de l'alignement
      if((al->getGaps() +al->getMis() )>max_mis)
	{
	  if(side == RIGHT)
	    {
	      delete [] s1;
	      delete [] s2;
	    }
	 
	  return AL_ERR;
	}


      if((al->getGaps() )>max_gaps)
	{
	  if(side == RIGHT)
	    {
	      delete [] s1;
	      delete [] s2;
	    }
	 
	  return AL_ERR;
	}


      ++cpt;
    }


  /// On recopie ce chemin dans l'objet Alignement
  if(side == LEFT)
    {
      while(cpt>0)
	{
	  --cpt;
	  al->addPair(mem[cpt],mem2[cpt]);
	}
    }
  else
    {
      i = 0;
      while(i<cpt)
	{
	  al->addPair(mem[i],mem2[i]);
	  ++i;
	}
      delete [] s1;
      delete [] s2;
    }
  
  return AL_VALID;
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


// x : sequence dentree, nt : lettre a compter, repetee sur 32 bits
//renvoit tableau de bits avec des 1 la où lettres trouvees, reste a faire poopcount
inline unsigned int count_nt(unsigned int x, unsigned int nt)
{
  unsigned int m= 0x55555555;

  x = ~ (x ^ nt) ;
  x = x & (x >> 1);
  return (x & m);

}



 int filtre_align_NT_vec(int left2, int right2, int max_mis 
			  ,int al,int cl,int tl, int gl, int ar,int cr,int tr, int gr 
		  ,unsigned int fleft2,unsigned int fright2)
{
  int a2=0,c2=0,t2=0,g2=0;
  int min_err=0;
  unsigned int masque ;
  int tairest;
  unsigned int ma= 0x00000000;
  unsigned int mt= 0xaaaaaaaa;
  unsigned int mc= 0x55555555;
  
  //cote gauche
  tairest =  min(left2,16); //lecture 16 max au total


  masque = ~  (0xffffffff >> tairest >> tairest); 
  
  // a1 = (int)popcount_mult (count_nt(fleft1,ma) & masque );
  a2 = (int)popcount_mult (count_nt(fleft2,ma) & masque );
  //c1 = (int)popcount_mult (count_nt(fleft1,mc) & masque );
  c2 = (int)popcount_mult (count_nt(fleft2,mc) & masque );
  //t1 = (int)popcount_mult (count_nt(fleft1,mt) & masque );
  t2 = (int)popcount_mult (count_nt(fleft2,mt) & masque );
  //g1 = (int)popcount_mult (count_nt(fleft1,mg) & masque );
  //g2 = (int)popcount_mult (count_nt(fleft2,mg) & masque );
   g2 = tairest - a2-c2-t2;


  min_err = abs (al-a2) + abs (cl-c2) + abs (tl-t2) + abs (gl-g2);
  if ((min_err/2)> max_mis) return 0;


  //----------------       cote droit 
  tairest =  min(right2,16);

  masque = ~  (0xffffffff >> tairest >> tairest); 
  
  // a1 = (int)popcount_mult (count_nt(fright1,ma) & masque );
  a2 = (int)popcount_mult (count_nt(fright2,ma) & masque );
  //c1 = (int)popcount_mult (count_nt(fright1,mc) & masque );
  c2 = (int)popcount_mult (count_nt(fright2,mc) & masque );
  //t1 = (int)popcount_mult (count_nt(fright1,mt) & masque );
  t2 = (int)popcount_mult (count_nt(fright2,mt) & masque );
  //g1 = (int)popcount_mult (count_nt(fright1,mg) & masque );
  //g2 = (int)popcount_mult (count_nt(fright2,mg) & masque );
   g2 = tairest - a2-c2-t2;


  min_err += abs (ar-a2) + abs (cr-c2) + abs (tr-t2) + abs (gr-g2);
  if ((min_err/2)> max_mis) return 0;

  return 1;

}


// premier filtre sur nombre de nt
//renvoit 1 si passe filtre, 0 si l'alignement est éliminé

 int filtre_align_NT(char* seq1, char* seq2, int left1, int left2, int right1, int right2, int max_mis,
		  char ** tabnt
		  ,unsigned int fleft1,unsigned int fright1
		  ,unsigned int fleft2,unsigned int fright2)
{

  // int size =min(min(left1,left2),16);
  // unsigned char nt;
  int a1=0,a2=0,c1=0,c2=0,t1=0,t2=0,g1=0,g2=0;
  //  int a1,a2,c1,c2,t1,t2,g1,g2;

  int min_err=0;
  unsigned int masque ; // masque de 16 bits 65535
  //methode lente pour valider algo
  
  
  int decalage;
  int tairest;
  int tai,tai_acc;
  unsigned int  idx1,idx2;
  int temp_err=0;

  /*
  // faire prefiltre a droite ?
  
   tairest =  min(min(right1,right2),16);
   tai = min (tairest,8);
  if(tai>=1)
    {
   decalage = 2*(16 -tai ); //16 maxi

       idx1 = (fright1 >> decalage);
       idx2 = (fright2 >> decalage);

      a1+=tabnt[0][idx1]; 
      a2+=tabnt[0][idx2];  

      if(abs(a1-a2)> max_mis) return 0; //premier filtrage le plus tot possible avec un seul nt donc pas /2

      c1+=tabnt[1][idx1]; 
      c2+=tabnt[1][idx2];  
      if(abs(c1-c2)> max_mis) return 0; 


  }

  */

  //cote gauche 
  //1er passage
  tairest =  min(min(left1,left2),16);
  tai = min (tairest,8);
  decalage = 2*(16 -tai ); //16 maxi
  //pas besoin de masque au premier passage
  if(tai>=1)
    {
      idx1 = (fleft1 >> decalage);
      idx2 = (fleft2 >> decalage);

	
    //  a1=tabnt[0][idx1]; 
	//  a2=tabnt[0][idx2];  

	//  if(abs(a1-a2)> max_mis) return 0; //premier filtrage le plus tot possible avec un seul nt donc pas /2

	c1=tabnt[1][idx1]; 
	c2=tabnt[1][idx2];  
//	if(abs(c1-c2)> max_mis) return 0; 

	t1=tabnt[2][idx1]; 
	t2=tabnt[2][idx2]; 
//	if(abs(t1-t2)> max_mis) return 0; 

	g1=tabnt[3][idx1]; 
	g2=tabnt[3][idx2]; 

//	if(abs(g1-g2)> max_mis) return 0; 

	a1 = tai - c1 - t1 -g1;
	a2 = tai - c2 - t2 -g2;

	temp_err = abs (a1-a2) + abs (c1-c2) + abs (t1-t2) + abs (g1-g2); //Result intermediaire , est il reutilisable ?

    //  if ((temp_err/2)> max_mis) return 0;


    }
      tairest -= tai;
      if(tairest>0) //Second passage si necessaire 
	{
	  tai_acc = tai;
	  tai = min(8,tairest);
	  decalage = 2*(16 -tai - tai_acc); //16 maxi
	  masque = 65535 >> (2*(8-tai)); 

	  idx1 = (fleft1 >> decalage)&masque;
	  idx2 = (fleft2 >> decalage)&masque;

	  //a1+=tabnt[0][idx1]; 
	  //a2+=tabnt[0][idx2];      
	  c1+=tabnt[1][idx1]; 
	  c2+=tabnt[1][idx2];         	   	       
	  t1+=tabnt[2][idx1]; 
	  t2+=tabnt[2][idx2]; 	       
	  g1+=tabnt[3][idx1]; 
	  g2+=tabnt[3][idx2]; 
	  a1 = tai + tai_acc - c1 - t1 -g1;
	  a2 = tai + tai_acc - c2 - t2 -g2;
	
	  min_err = abs (a1-a2) + abs (c1-c2) + abs (t1-t2) + abs (g1-g2);
	  if ((min_err/2)> max_mis) return 0;

	}else
	{
	  min_err = temp_err;
       if ((min_err/2)> max_mis) return 0;
	}
       




  //cote droit 
  //1er passage



  tairest =  min(min(right1,right2),16);
  tai = min (tairest,8);
  decalage = 2*(16 -tai ); //16 maxi

  if(tai>=1)
    {
      idx1 = (fright1 >> decalage);
      idx2 = (fright2 >> decalage);

	
    //  a1=tabnt[0][idx1]; 
    //  a2=tabnt[0][idx2];  

    //  if(abs(a1-a2)> max_mis) return 0; 

      c1=tabnt[1][idx1]; 
      c2=tabnt[1][idx2];         	   	       
      t1=tabnt[2][idx1]; 
      t2=tabnt[2][idx2]; 	       
      g1=tabnt[3][idx1]; 
      g2=tabnt[3][idx2]; 

      a1 = tai - c1 - t1 -g1;
      a2 = tai - c2 - t2 -g2;
      temp_err = abs (a1-a2) + abs (c1-c2) + abs (t1-t2) + abs (g1-g2); //Result intermediaire , est il reutilisable ?

     // if ((temp_err/2)> max_mis) return 0;


    }
      tairest -= tai;
      if(tairest>0) //Second passage si necessaire 
	{
	  tai_acc = tai;
	  tai = min(8,tairest);
	  decalage = 2*(16 -tai - tai_acc); //16 maxi
	  masque = 65535 >> (2*(8-tai)); 

	  idx1 = (fright1 >> decalage)&masque;
	  idx2 = (fright2 >> decalage)&masque;

//	  a1+=tabnt[0][idx1]; 
//	  a2+=tabnt[0][idx2];      
	  c1+=tabnt[1][idx1]; 
	  c2+=tabnt[1][idx2];         	   	       
	  t1+=tabnt[2][idx1]; 
	  t2+=tabnt[2][idx2]; 	       
	  g1+=tabnt[3][idx1]; 
	  g2+=tabnt[3][idx2]; 
	         a1 = tai + tai_acc - c1 - t1 -g1;
		  a2 = tai + tai_acc - c2 - t2 -g2;
      min_err += abs (a1-a2) + abs (c1-c2) + abs (t1-t2) + abs (g1-g2);

     // if ((min_err/2)> max_mis) return 0;

	}else
	{
	  min_err += temp_err;
	}
       

            if ((min_err/2)> max_mis) return 0;


  return 1;

}







//test pre filtre avec tablee SW precalc sur 1 TAI_TAB seulemnt , petit donc plus rapide ? 
inline  int pre_filtre_align( int left2,int right2, int max_mis, char ** tabprec ,unsigned int fleft1,unsigned int fright1
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
/*
//tai_acc nombre de nt deja lus
//tai nombre de nt a lire
//tabcode 128 bits 64 nt
inline unsigned int get_code(unsigned int * tabcode,int tai_acc,int tai)
{
  char * var = (char*) tabcode;
  int dec = (2*tai_acc)/8; //nombre d'octets à manger
  unsigned int * temp;
  tai_acc -= 4*dec;
  var += dec;
  temp =  (unsigned int *) var;
  int    decalage = 2*(4*sizeof(unsigned int) -tai - tai_acc);
  unsigned int masque = (1 << tai) -1 ;   //masque de 2*tai bit à 1 

  return (((*temp)>>decalage) &masque);

}
*/


  
inline unsigned int get_code(unsigned int *  vec,int tai_acc,int tai)
{

  int dec = tai_acc/10; //cahque uint sert a lecture de plage de 10 nt
  unsigned int temp = vec[dec];
  tai_acc -= 10*dec; //4*dec
  int    decalage = 2*(4*sizeof(unsigned int) -tai - tai_acc);
  unsigned int masque = (1 << (2*tai)) -1 ;   //masque de 2*tai bit à 1 
  return ((temp>>decalage) &masque);

}

int filtre_align_DC(char* seq1, char* seq2, int left2, int right2, int max_mis,char ** tabprec
		    ,unsigned int fleft1,unsigned int fright1
		    ,unsigned int fleft2,unsigned int fright2,int num_gaps
		    )
{
  int err,tai;
  int errtot=0;
  int errd=0;
  int errg=0;
  int tairest ;
  int tailue;
  int decalage =0;

  //
  int nbu;

  //codage cote gauche
      nbu = (left2) / 10;
      nbu += ((left2 -( nbu*10))!=0); //nombre de uint necessaires au stockage chavauchant
      if (!left2) nbu++;
      unsigned int  bleft1[nbu];// vla , pas grand pour seq 500  nbu vaut au plus 50
      unsigned int  bleft2[nbu];// vla

      bleft1[0]=fleft1; //les 2 premiers sont dans l'index
      bleft2[0]=fleft2;

      tairest = left2 - 10; //16 deja lus dans lindex // on en a lu 10  + 6chev
      tailue=10;//chevauchant donc tailue 10 seulement
      tai = min(16,tairest);//Reste au plus tairest a lire
      nbu=1;
      while (tairest>0)
	{
	  bleft1[nbu]=seq2code_rev(seq1-1-tailue,tai)<<(2*(16-tai));
	  bleft2[nbu]=seq2code_rev(seq2-1-tailue,tai)<<(2*(16-tai));
	  nbu++;
	  tairest -= 10  ; //il e nreste 10 de moins
	  tailue += 10; // car on en a lu 10 de plus pour de bon
	  tai = min(16,tairest);
	}
      //codage cote droit
      nbu = (right2) / 10;
      nbu += ((right2 -( nbu*10))!=0); 
      if (!right2) nbu++;
      unsigned int  bright1[nbu];// vla
      unsigned int  bright2[nbu];// vla
      bright1[0]=fright1; 
      bright2[0]=fright2;

      tairest = right2 -10; 
      tailue=10;
      tai = min(16,tairest );
      nbu=1;
      while (tairest>0)
	{
	  bright1[nbu]=seq2code(seq1+SIZE_SEED+tailue,tai)<<(2*(16-tai));
	  bright2[nbu]=seq2code(seq2+SIZE_SEED+tailue,tai)<<(2*(16-tai));
	  nbu++;
	  
	  tairest -= 10;
	  tailue += 10 ;
	  tai = min(16,tairest);
	}
  // le codage est fini

   tairest = left2;
   tai = min(TAI_BIG_TAB,tairest);

  // A gauche
  if(tai>=1)
    {
      decalage = 2*(4*sizeof(unsigned int) -tai);
      errg = (int) tabprec[fleft1>>decalage][fleft2>>decalage] ;
      errtot += errg;
    }

  if (errtot> max_mis) 
    return 0;

  // A droite
  tairest = right2;
  tai = min(TAI_BIG_TAB,tairest);

  int tai_acc=0;
  if ( tai>=1) 
    {
       decalage = 2*(4*sizeof(unsigned int) -tai);
       errd += (int) tabprec[fright1>>decalage][fright2>>decalage] ;
      errtot += errd;
    }
  if (errtot> max_mis) 
    return 0;

	

  /////////// filtre phase 2 : apres les TAI_BIG_TAB premieres bases on doit tester toutes les possibilites	
  //avec boucle partie droite
  tairest-= tai;
  tai_acc+= tai;
  tai = min(TAI_BIG_TAB,tairest);
  int errdprec,errgprec;
  while ( tairest>0) 
    {
      errdprec = errd;
      errd =  errdprec +  (int) tabprec[get_code(bright1,tai_acc,tai)][get_code(bright2,tai_acc,tai)] ;

      for (err=1; err <= num_gaps && err <= tai_acc; err++)
	{

	  errd=min(errd, max(errdprec,err) + (int) tabprec[get_code(bright1,tai_acc,tai)][get_code(bright2,tai_acc-err,tai)] );  
	  errd=min(errd, max(errdprec,err) + (int) tabprec[get_code(bright1,tai_acc-err,tai)][get_code(bright2,tai_acc,tai)] );   
	}


      errtot = errd +errg;
  if (errtot> max_mis) 
    return 0;

      tairest-= tai;
      tai_acc+= tai;
      tai = min(TAI_BIG_TAB,tairest);
	    
    }

  //  puis suite gauche  avec boucle 
  tai_acc = 0;
  tairest = left2;
  tai = min(TAI_BIG_TAB,tairest);
  // on reprend la ou on en etait
  tairest-= tai;
  tai_acc+= tai;
  tai = min(TAI_BIG_TAB,tairest);


while (tairest>0)
    {
      errgprec = errg;
      errg =  errgprec + tabprec[get_code(bleft1,tai_acc,tai)][get_code(bleft2,tai_acc,tai)] ;

      for (err=1; err <= num_gaps && err <= tai_acc ; err++)
	{
	  errg=min(errg, max(errgprec,err) + (int) tabprec[get_code(bleft1,tai_acc,tai)][get_code(bleft2,tai_acc-err,tai)] ); 
	  errg=min(errg, max(errgprec,err) + (int) tabprec[get_code(bleft1,tai_acc-err,tai)][get_code(bleft2,tai_acc,tai)] );
	}
      errtot = errg+ errd;
  if (errtot> max_mis) 
    return 0;

      tairest-= tai;
      tai_acc+= tai;
      tai = min(TAI_BIG_TAB,tairest);	    
    }
  return 1;
}



//renvoit 1 si passe filtre, 0 si l'alignement est éliminé

 int filtre_align( int left2, int right2, int max_mis,
		  char ** tabprec
		  ,unsigned int fleft1,unsigned int fright1
		   ,unsigned int fleft2,unsigned int fright2,int num_gaps)
{
  /// On vérifie le nombre de caractères utilisables dans la première séquence
  // int size = min(left1,left2+max_mis);
  //verif table precalculee : donne nombre erreurs min de l'alignement de 2 seq de taille TAI_TAB
  // si trop d'erreurs on elimine l'alignement
  //
  int err,tai;
  int errtot=0;
  int errd=0;
  int errg=0;
  int tairest = left2;
  tai = min(TAI_TAB,tairest);

  int decalage =0;
  // A gauche
  if(tai>=1)
    {
      decalage = 2*(4*sizeof(unsigned int) -tai);
      errg = (int) tabprec[fleft1>>decalage][fleft2>>decalage] ;

      errtot += errg;
    }
  if (errtot> max_mis) return 0;
  // A droite
  tairest = right2;
  tai = min(TAI_TAB,tairest);

  int tai_acc=0;

  if ( tai>=1) 
    {
      decalage = 2*(4*sizeof(unsigned int) -tai);

      errd += (int) tabprec[fright1>>decalage][fright2>>decalage] ;
      errtot += errd;

    }
  if (errtot> max_mis) return 0;

	

  /////////// filtre phase 2 : apres les TAI_TAB premieres bases on doit tester toutes les possibilites	
  //avec boucle partie droite
  tairest-= tai;
  tai_acc+= tai;
  tai = min(TAI_TAB,tairest); 
  int errdprec,errgprec;
  unsigned char masque ;
  while (tai_acc < 16 && tairest>0)
    {
      decalage = 2*(4*sizeof(unsigned int) -tai - tai_acc);
      errdprec = errd;
      //    masque = 255>>(2*(4-tai));  // 4 en dur ?
   masque = (1 << (2*tai)) -1 ; 
      errd =  errdprec +  (int) tabprec[(fright1>>decalage)&masque][(fright2>>decalage)&masque] ;

      for (err=1; err <= num_gaps; err++) //si err élevé, depassement au premier coup, pas grave mais minore l'erreur
	{
	  errd=min(errd, max(errdprec,err) + (int) tabprec[(fright1>>decalage)&masque][(fright2>>(decalage+2*err))&masque] );
	  errd=min(errd, max(errdprec,err) + (int) tabprec[(fright1>>(decalage+2*err))&masque][fright2>>(decalage)&masque] );	  
	}


      errtot = errd +errg;
      if (errtot> max_mis) return 0;
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


  while (tai_acc < 16 && tairest>0)

    {
      decalage = 2*(4*sizeof(unsigned int) -tai - tai_acc);
      //  masque = 255>>(2*(4-tai)); 
      masque = (1 << (2*tai)) -1 ; 
      errgprec = errg;
      errg =  errgprec +  (int) tabprec[(fleft1>>decalage)&masque][(fleft2>>decalage)&masque] ;

      for (err=1; err <= num_gaps; err++)
	{
	  errg=min(errg, max(errgprec,err) + (int) tabprec[(fleft1>>decalage)&masque][(fleft2>>(decalage+2*err))&masque] );
	  errg=min(errg, max(errgprec,err) + (int) tabprec[(fleft1>>(decalage+2*err))&masque][(fleft2>>decalage)&masque] );
	}
 

      errtot = errg+ errd;
      if (errtot> max_mis) return 0;
      tairest-= tai;
      tai_acc+= tai;
      tai = min(TAI_TAB,tairest);

	    
    }

  return 1;
}




 int filtre_align_NT_vec_sse( int left, int right, int max_mis,unsigned int fleft,unsigned int fright,__m128i nt_query_left,__m128i nt_query_right,
							  __m128i vec_5,__m128i vec_3,
						      __m128i vec_0f,__m128i vec_h1,__m128i vec_nt
					/*	,int * a1,int * a2,int * c1,int * c2,int * t1,int * t2,int * g1,int * g2,int * minerr*/)
{
	int min_err=0;
	unsigned int masque ;

	int tairest;
	

	__m128i vec_cpt1;
//	__m128i vec_nt;
	__m128i temp;
//	__m128i vec_5;
//	__m128i vec_3,vec_0f,vec_h1; //3masques pour popcnt
	//init des vecteurs masques // a faire en dehors ?
//	vec_5  = _mm_set_epi32(mc, mc, mc, mc);
//	vec_3  = _mm_set_epi32(m2, m2, m2, m2);
//	vec_0f  = _mm_set_epi32(m4, m4, m4, m4);
//	vec_h1  = _mm_set_epi32(h01, h01, h01, h01);
//	vec_nt = _mm_set_epi32(0x00000000, 0xaaaaaaaa, 0x55555555, 0xffffffff); // reg 3 2 1 0 


	vec_cpt1 = _mm_set_epi32(fleft, fleft, fleft, fleft);  
	vec_cpt1 = _mm_xor_si128(vec_cpt1,vec_nt ); //xor avec masque lettres
	temp = _mm_srli_epi32(vec_cpt1, 1) ;  //shift de 1 vers droite
	vec_cpt1 = _mm_or_si128(temp,vec_cpt1); //ou
	vec_cpt1 = _mm_andnot_si128(vec_cpt1,vec_5);
	//idem cpt2

	//reste popcnt a faire
	//masque pour selec zone a lire
	tairest =  min(left,16); //lecture 16 max au total
	masque = ~  (0xffffffff >> tairest >> tairest);
	temp  = _mm_set_epi32(masque, masque, masque, masque);
	vec_cpt1 = _mm_and_si128(vec_cpt1,temp);


	temp = _mm_and_si128(vec_cpt1, vec_3);
	vec_cpt1 = _mm_srli_epi32(vec_cpt1, 2) ; 
	vec_cpt1 = _mm_and_si128(vec_cpt1, vec_3);
	vec_cpt1 = _mm_add_epi32(vec_cpt1, temp); 
	temp =  _mm_srli_epi32(vec_cpt1, 4) ; 
	vec_cpt1 = _mm_add_epi32(vec_cpt1, temp);
	vec_cpt1 = _mm_and_si128(vec_cpt1, vec_0f); 
	temp =  _mm_srli_epi32(vec_cpt1, 8) ; 
	vec_cpt1 = _mm_add_epi32(vec_cpt1, temp);
	temp =  _mm_srli_epi32(vec_cpt1, 16) ; 
	vec_cpt1 = _mm_add_epi32(vec_cpt1, temp);

//filtre final

	temp  = _mm_set_epi32(0x3f, 0x3f, 0x3f, 0x3f);
	vec_cpt1 = _mm_and_si128(vec_cpt1, temp);
/*	*a1=_mm_extract_epi32(vec_cpt1,3);
	*t1=_mm_extract_epi32(vec_cpt1,2);	
	*c1=_mm_extract_epi32(vec_cpt1,1);
	*g1=_mm_extract_epi32(vec_cpt1,0);
	*/
	
	vec_cpt1=_mm_sub_epi32(vec_cpt1,nt_query_left);
	vec_cpt1=_mm_abs_epi32(vec_cpt1); //                              // ssse3
	vec_cpt1=_mm_hadd_epi32(vec_cpt1,vec_cpt1);                     // ssse3
	vec_cpt1=_mm_hadd_epi32(vec_cpt1,vec_cpt1);                     //  ssse3
	min_err =_mm_cvtsi128_si32(vec_cpt1);

	if ((min_err/2)> max_mis) return 0;


//puis idem cote droit

/*
	char gr[17];
	code2seq(fright,gr,16);
	gr[17]='\0';
		printf("BANK right %s \n",gr);
*/

		vec_cpt1 = _mm_set_epi32(fright, fright, fright, fright);  
	//	Mprint(vec_cpt1,"vec_cpt1");
		vec_cpt1 = _mm_xor_si128(vec_cpt1,vec_nt ); //xor avec masque lettres
		temp = _mm_srli_epi32(vec_cpt1, 1) ;  //shift de 1 vers droite
		vec_cpt1 = _mm_or_si128(temp,vec_cpt1); //ou
		vec_cpt1 = _mm_andnot_si128(vec_cpt1,vec_5);
		//idem cpt2

		//reste popcnt a faire
		//masque pour selec zone a lire
		tairest =  min(right,16); //lecture 16 max au total
		masque = ~  (0xffffffff >> tairest >> tairest);
		temp  = _mm_set_epi32(masque, masque, masque, masque);
		vec_cpt1 = _mm_and_si128(vec_cpt1,temp);


		temp = _mm_and_si128(vec_cpt1, vec_3);
		vec_cpt1 = _mm_srli_epi32(vec_cpt1, 2) ; 
		vec_cpt1 = _mm_and_si128(vec_cpt1, vec_3);
		vec_cpt1 = _mm_add_epi32(vec_cpt1, temp); 
		temp =  _mm_srli_epi32(vec_cpt1, 4) ; 
		vec_cpt1 = _mm_add_epi32(vec_cpt1, temp);
		vec_cpt1 = _mm_and_si128(vec_cpt1, vec_0f); 
		temp =  _mm_srli_epi32(vec_cpt1, 8) ; 
		vec_cpt1 = _mm_add_epi32(vec_cpt1, temp);
		temp =  _mm_srli_epi32(vec_cpt1, 16) ; 
		vec_cpt1 = _mm_add_epi32(vec_cpt1, temp);

	//filtre final

		temp  = _mm_set_epi32(0x3f, 0x3f, 0x3f, 0x3f);
		vec_cpt1 = _mm_and_si128(vec_cpt1, temp);
	//	Mprint(vec_cpt1,"vec_cpt1");
/*
		*a2=_mm_extract_epi32(vec_cpt1,3);
		*t2=_mm_extract_epi32(vec_cpt1,2);	
		*c2=_mm_extract_epi32(vec_cpt1,1);
		*g2=_mm_extract_epi32(vec_cpt1,0);
*/


		vec_cpt1=_mm_sub_epi32(vec_cpt1,nt_query_right);
		vec_cpt1=_mm_abs_epi32(vec_cpt1); //a lair couteux 2ic en 0,7
		vec_cpt1=_mm_hadd_epi32(vec_cpt1,vec_cpt1);
		vec_cpt1=_mm_hadd_epi32(vec_cpt1,vec_cpt1);
		min_err += _mm_cvtsi128_si32(vec_cpt1);


		
	//	*minerr = min_err;
	if ((min_err/2)> max_mis) return 0;


	return 1;

}


//inutilse maintenant, car filtre TNW ets plus efficace
int filtre_diNT_DEC(int left2,int right2,int max_mis ,unsigned int fleft1,unsigned int fright1 ,unsigned int fleft2, unsigned int fright2)
{

  char tab1[16];
  char tab2[16];
  int tairest;
  int i;
   int min_err;
  unsigned int * zero;
  zero = (unsigned int *) &tab1;
  zero[0]=0;zero[1]=0;zero[2]=0;zero[3]=0;
  zero = (unsigned int *) &tab2;
  zero[0]=0;zero[1]=0;zero[2]=0;zero[3]=0;

  //cote gauche
  tairest =  min(left2,16); 
  for(i=1; i<tairest; i++) // i commence a 1 : si 1 nt : rien a lire
    {
      tab1[(fleft1 & 0xff000000)>> 28 ]++; //poids fort  a lire en premier , le contraire serait plus simple
                                           //mais je sais plus, peut etre sinon prob pour sous table < 4
      fleft1 = fleft1 << 2;

      tab2[(fleft2 & 0xff000000)>> 28]++;
      fleft2 = fleft2 << 2;
    }

  // a faire en sse ?
  min_err = abs (tab1[0]-tab2[0]) + abs (tab1[1]-tab2[1]) + abs (tab1[2]-tab2[2]) + abs (tab1[3]-tab2[3])+
    abs (tab1[4]-tab2[4]) + abs (tab1[5]-tab2[5]) + abs (tab1[6]-tab2[6]) + abs (tab1[7]-tab2[7])+
    abs (tab1[8]-tab2[8]) + abs (tab1[9]-tab2[9]) + abs (tab1[10]-tab2[10]) + abs (tab1[11]-tab2[11])+
    abs (tab1[12]-tab2[12]) + abs (tab1[13]-tab2[13]) + abs (tab1[14]-tab2[14]) + abs (tab1[15]-tab2[15]);

  if ((min_err/4)> max_mis) return 0;


  //remise a zero des compteurs avant cote droit
  zero = (unsigned int *) &tab1;
  zero[0]=0;zero[1]=0;zero[2]=0;zero[3]=0;
  zero = (unsigned int *) &tab2;
  zero[0]=0;zero[1]=0;zero[2]=0;zero[3]=0;
 
  //cote droit
  tairest =  min(right2,16); 
  for(i=1; i<tairest; i++)
    {
      tab1[(fright1& 0xff000000)>> 28 ]++;
      fright1 = fright1 << 2;

      tab2[(fright2 & 0xff000000)>> 28 ]++;
      fright2 = fright2 << 2;
    }
  min_err += abs (tab1[0]-tab2[0]) + abs (tab1[1]-tab2[1]) + abs (tab1[2]-tab2[2]) + abs (tab1[3]-tab2[3])+
    abs (tab1[4]-tab2[4]) + abs (tab1[5]-tab2[5]) + abs (tab1[6]-tab2[6]) + abs (tab1[7]-tab2[7])+
    abs (tab1[8]-tab2[8]) + abs (tab1[9]-tab2[9]) + abs (tab1[10]-tab2[10]) + abs (tab1[11]-tab2[11])+
    abs (tab1[12]-tab2[12]) + abs (tab1[13]-tab2[13]) + abs (tab1[14]-tab2[14]) + abs (tab1[15]-tab2[15]);
  if ((min_err/4)> max_mis) return 0;

  return 1;



}






inline void calc_nt_query(unsigned int fleft, unsigned int fright,__m128i vec_5,__m128i vec_3,__m128i vec_0f,__m128i vec_h1,__m128i vec_nt
		   ,__m128i * pt_nt_query_left,__m128i *  pt_nt_query_right, int left2, int right2)
{

  __m128i temp;
  int tairest;
  unsigned int masque;
  __m128i nt_query_left;
  __m128i nt_query_right;

  nt_query_left = _mm_set_epi32(fleft, fleft, fleft, fleft);  
	
  nt_query_left = _mm_xor_si128(nt_query_left,vec_nt ); //xor avec masque lettres
  temp = _mm_srli_epi32(nt_query_left, 1) ;  //shift de 1 vers droite
  nt_query_left = _mm_or_si128(temp,nt_query_left); //ou
  nt_query_left = _mm_andnot_si128(nt_query_left,vec_5);

  //masque pour selec zone a lire
  tairest =  min(left2,16); //lecture 16 max au total
  masque = ~  (0xffffffff >> tairest >> tairest);
  temp  = _mm_set_epi32(masque, masque, masque, masque);
  nt_query_left = _mm_and_si128(nt_query_left,temp);

  temp = _mm_and_si128(nt_query_left, vec_3);
  nt_query_left = _mm_srli_epi32(nt_query_left, 2) ; 
  nt_query_left = _mm_and_si128(nt_query_left, vec_3);
  nt_query_left = _mm_add_epi32(nt_query_left, temp); 
  temp =  _mm_srli_epi32(nt_query_left, 4) ; 
  nt_query_left = _mm_add_epi32(nt_query_left, temp);
  nt_query_left = _mm_and_si128(nt_query_left, vec_0f); 
  temp =  _mm_srli_epi32(nt_query_left, 8) ; 
  nt_query_left = _mm_add_epi32(nt_query_left, temp);
  temp =  _mm_srli_epi32(nt_query_left, 16) ; 
  nt_query_left = _mm_add_epi32(nt_query_left, temp);
  //filtre final
  temp  = _mm_set_epi32(0x3f, 0x3f, 0x3f, 0x3f);
  nt_query_left = _mm_and_si128(nt_query_left, temp);


  //cote droit


  nt_query_right = _mm_set_epi32(fright, fright, fright, fright);
	  
  nt_query_right = _mm_xor_si128(nt_query_right,vec_nt ); //xor avec masque lettres
  temp = _mm_srli_epi32(nt_query_right, 1) ;  //shift de 1 vers droite
  nt_query_right = _mm_or_si128(temp,nt_query_right); //ou
  nt_query_right = _mm_andnot_si128(nt_query_right,vec_5);

  //masque pour selec zone a lire
  tairest =  min(right2,16); //lecture 16 max au total
  masque = ~  (0xffffffff >> tairest >> tairest);
  temp  = _mm_set_epi32(masque, masque, masque, masque);
  nt_query_right = _mm_and_si128(nt_query_right,temp);

  temp = _mm_and_si128(nt_query_right, vec_3);
  nt_query_right = _mm_srli_epi32(nt_query_right, 2) ; 
  nt_query_right = _mm_and_si128(nt_query_right, vec_3);
  nt_query_right = _mm_add_epi32(nt_query_right, temp); 
  temp =  _mm_srli_epi32(nt_query_right, 4) ; 
  nt_query_right = _mm_add_epi32(nt_query_right, temp);
  nt_query_right = _mm_and_si128(nt_query_right, vec_0f); 
  temp =  _mm_srli_epi32(nt_query_right, 8) ; 
  nt_query_right = _mm_add_epi32(nt_query_right, temp);
  temp =  _mm_srli_epi32(nt_query_right, 16) ; 
  nt_query_right = _mm_add_epi32(nt_query_right, temp);
  //filtre final
  temp  = _mm_set_epi32(0x3f, 0x3f, 0x3f, 0x3f);
  nt_query_right = _mm_and_si128(nt_query_right, temp);

	  

  *pt_nt_query_left = nt_query_left;
  *pt_nt_query_right = nt_query_right;

}
/**
 * Fonction permettant de calculer un score d'alignement obtenu entre 2 séquences
 * \param al l'objet Alignment qui constituera le résultat
 * \param seq1 le pointeur du tableau des caractères de la première banque
 * \param seq2 le pointeur du tableau des caractères de la seconde banque
 * \param left1 le nombre de bases présentes à gauche dans la première séquence
 * \param left2 le nombre de bases présentes à gauche dans la seconde séquence
 * \param right1 le nombre de bases présentes à droite dans la première séquence
 * \param right2 le nombre de bases présentes à droite dans la seconde séquence
 * \param max_mis le nombre maximal de mésappariement pour conserver l'alignement
 * \param num_seed le code de la graine courante
 * \param num_gaps le nombre maximal de gaps autorisés dans les alignements
 * \return Le nombre de mésappariements obtenus dans l'alignement
 */


int withgap_scoreHit(Alignment* al, char* seq1, char* seq2, int left1, int left2, int right1, int right2, int max_mis, int num_seed, int num_gaps)
{	 
  int i = 0;
  int validity;
  int nb_gaps;
  nb_gaps = num_gaps;
  int max_err= max_mis;
  /// On vérifie le nombre de caractères utilisables dans la première séquence
  int size = min(left1,left2+nb_gaps);
	
  // printf("l1 l2 r1 r2  %i %i %i %i    size %i   ngaps %i  taille %i Mo\n", left1,left2,right1,right2,size, num_gaps,size*left2*sizeof(caseMatrix)/1024/1024 );
  //return max_mis+1;

  if(size > 0 && left2 > 0)
    {


      if (left2 > 10000 || size > 10000) 
	// allocation memoire dans pile,et bande uniquement
	validity = WN_very_big_alignment(al, seq1, seq2, size, left2, num_gaps, LEFT, num_seed, max_err); 
      else if (left2 > 1000 || size > 1000) 
	// allocation memoire avec malloc
	validity = WN_big_alignment(al, seq1, seq2, size, left2, num_gaps, LEFT, num_seed, max_err); 
      else
	// allocation memoire dans pile, ne passe pas pour sequence > 1000 si stack 10 Mo
	validity = WN_alignment(al, seq1, seq2, size, left2, num_gaps, LEFT, num_seed, max_err); 

      /// On vérifie la validité de l'alignement

      //  if(validity == AL_DOUBLON)  aldoublon++;
      if(validity != AL_VALID) return max_mis+1;
    }
  /// On ajoute la graine à l'appariement

  while(i<SIZE_SEED)
    {
      al->addPair(seq1[i],seq2[i]);
      ++i;
    }

  char* s1 = seq1 + SIZE_SEED;
  char* s2 = seq2 + SIZE_SEED;

  /// On vérifie le nombre de caractères utilisables dans la première séquence
  size = min(right1,right2+nb_gaps);
	
  if(size > 0 && right2 > 0)
    {

      //  printf("size %i  right2 %i max_err %i  num_seed %i  max_err  %i  right1 %i left1 %i \n", size, right2, max_err,  num_seed, max_err,right1, left1);
      if (right2 > 10000 || size > 10000) 
	// allocation memoire dans pile,et bande uniquement
	validity = WN_very_big_alignment(al, s1, s2, size, right2, num_gaps, RIGHT, num_seed, max_err); 
      else if (right2 > 1000 || size > 1000)
	validity = WN_big_alignment(al, s1, s2, size, right2, num_gaps, RIGHT, num_seed, max_err);
      else
	validity = WN_alignment(al, s1, s2, size, right2, num_gaps, RIGHT, num_seed, max_err);
      
      /// On vérifie la validité de l'alignement
      //  if(validity == AL_DOUBLON)  aldoublon++;
      if(validity != AL_VALID) return max_mis+1;
    }


  al->finish_cigar();	 // indique au string cigar que l'align est fini

  // if( (al->getMis() + al->getGaps()) <= max_mis ) alok++;
  return al->getMis() + al->getGaps();
}


#ifdef CPT_FILT

extern long long int cpt_hit;
extern long long int cpt_pref;
extern long long int cpt_ntsse;
extern long long int cpt_fil;
extern long long int cpt_doub;
extern long long int cpt_max;
extern long long int cpt_dint;
extern long long int cpt_good;

#endif

/**
 * Méthode qui à partir des positions d'une graine va effectuer tous les alignements possibles avec les séquences des banques
 * \param ff le fichier de sortie, où seront affichés les alignements obtenus
 * \param BK1 le pointeur de la première banque de séquences
 * \param BK2 le pointeur de la seconde banque de séquences
 * \param I1 le pointeur de l'index de la première banque de séquences
 * \param num_seed l'indice de la graine utilisée pour déterminer les alignements
 * \param idpc le pourcentage de ressemblances minimal entre les 2 séquences pour qu'un alignement soit conservé
 * \param num_gaps le nombre maximal de gaps autorisés dans les alignements
 * \return le nombre d'alignements obtenus
 */
//traite toutes les query des numeros start a end 
int withgap_align(FILE *ff, Bank *BK1, Bank *BK2, Index * I1, int idpc, bool rev_comp, int num_gaps,char ** tabprec,char ** tabnt,char **tabprec7,Stat * St,Doublon * Doub, int start, int end )
{

  // printf("ENTREE withgap align  maxpos %i nbthreads %i maxhits %i bestal %i  \n",MAXPOS,NBTHREADS,MAXHITS,BESTAL);
  int  nb_hit;
  nb_hit = 0;

  // printf("nummax, %i numalign %i \n",NUMMAX,NUM_ALIGN_COURANT);
  if(NUMMAX && (NUM_ALIGN_COURANT >= NUMMAX)) return nb_hit;

  int  idx1;
  int  i;
  int  left1, left2, right1, right2;
  int  max_mis, nb_mis,  score;
  int nb_occur = 0;
  int n = 0;
   long long absolute_pos;
   /*  int ale, cl, tl, gl, ar,cr,tr, gr ;//pour test sans sse
  unsigned int ma= 0x00000000;
  unsigned int mt= 0xaaaaaaaa;*/
  //fin test sans sse 
  int nquery;
  int size_query;
  int num_sequence,offset_sequence,offseq,sizeseq;
  int offseq2;
  unsigned int fleft=0; // pour filtre
  unsigned int fright=0; 
  char doing_rev_query=0;
  unsigned int fleftB,frightB;
   int tair=0;
  long long  num_seed;
  char *s1, *s2, *s2_rev;
  s2_rev = new char[BK2->tailleMaxSeq];
  Alignment* al;
  Seed *  seed_bank;
  al =  new Alignment(BK2->tailleMaxSeq + (int)((BK2->tailleMaxSeq*(100-idpc)) / 100) + 1);

  int cpt_pos=0;
  int percent_gaps=num_gaps;
  int cpt_pos_amont;
  int cpt_query;
  unsigned int mc= 0x55555555;
  unsigned int m2 = 0x33333333;
  unsigned int m4 = 0x0f0f0f0f;
  unsigned int h01 = 0x01010101;
  unsigned int masque;
	
  int tairest;
  __m128i vec_nt;
  __m128i vec_5,vec_3,vec_0f,vec_h1; //3masques pour popcnt
  //init des vecteurs masques 
  vec_5  = _mm_set_epi32(mc, mc, mc, mc);
  vec_3  = _mm_set_epi32(m2, m2, m2, m2);
  vec_0f  = _mm_set_epi32(m4, m4, m4, m4);
  vec_h1  = _mm_set_epi32(h01, h01, h01, h01);
  vec_nt = _mm_set_epi32(0x00000000, 0xaaaaaaaa, 0x55555555, 0xffffffff); // reg 3 2 1 0  
	
  __m128i nt_query_left;
  __m128i temp;
  __m128i nt_query_right;
  //pour chaque query
  for (nquery = start; nquery <= end; nquery++) 
    {

         if(NUMMAX && (NUM_ALIGN_COURANT >= NUMMAX)) return nb_hit;



        
      cpt_query=0;
      // nouvelle liste pour chaque query
      size_query =   BK2->size[nquery] ;  
      num_seed=-1;
      offseq2 = BK2->seq[nquery];
      s2 = (BK2->data) + BK2->seq[nquery];
      if( doing_rev_query) //revcomp s2
	{
	  for(int p=0;p<=size_query-1;p++)
	    {
	      s2_rev[p] = complNT(s2[size_query-1-p]);
	    }
	  s2 = s2_rev;
	}
      fleft=0;fright=0;

      max_mis = (size_query*(100-idpc)) / 100;
      //si num_gaps non choisi,
      // num_gaps = max_mis ;
      if(NUMGAPS_AUTO) 
	{
	  num_gaps = max_mis ;
	}
      else
	{
	  int max_gaps = (size_query*percent_gaps) / 100;
	  num_gaps = max_gaps;
	}
      //pour chaque graine dans query
      for(i=0; i<= (size_query - SIZE_SEED); i++)
	{
	  cpt_pos = 0;
	  cpt_pos_amont = 0;

	  left2 =  i;
	  right2 =  size_query - i - SIZE_SEED;

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

	  // + calcul nb nt pour filtre nt  //faisable plus rapidement ar decalage, a changer eventuellement
	  //   calc_nt_query(fleft,fright,vec_5,vec_3,vec_0f,vec_h1,vec_nt
	  // ,&nt_query_left,&nt_query_right,left2,right2);
	   
	  nt_query_left = _mm_set_epi32(fleft, fleft, fleft, fleft);  
	  nt_query_left = _mm_xor_si128(nt_query_left,vec_nt ); //xor avec masque lettres
	  temp = _mm_srli_epi32(nt_query_left, 1) ;  //shift de 1 vers droite
	  nt_query_left = _mm_or_si128(temp,nt_query_left); //ou
	  nt_query_left = _mm_andnot_si128(nt_query_left,vec_5);
	  //masque pour selec zone a lire
	  tairest =  min(left2,16); //lecture 16 max au total
	  masque = ~  (0xffffffff >> tairest >> tairest);
	  temp  = _mm_set_epi32(masque, masque, masque, masque);
	  nt_query_left = _mm_and_si128(nt_query_left,temp);
	  temp = _mm_and_si128(nt_query_left, vec_3);
	  nt_query_left = _mm_srli_epi32(nt_query_left, 2) ; 
	  nt_query_left = _mm_and_si128(nt_query_left, vec_3);
	  nt_query_left = _mm_add_epi32(nt_query_left, temp); 
	  temp =  _mm_srli_epi32(nt_query_left, 4) ; 
	  nt_query_left = _mm_add_epi32(nt_query_left, temp);
	  nt_query_left = _mm_and_si128(nt_query_left, vec_0f); 
	  temp =  _mm_srli_epi32(nt_query_left, 8) ; 
	  nt_query_left = _mm_add_epi32(nt_query_left, temp);
	  temp =  _mm_srli_epi32(nt_query_left, 16) ; 
	  nt_query_left = _mm_add_epi32(nt_query_left, temp);
	  //filtre final
	  temp  = _mm_set_epi32(0x3f, 0x3f, 0x3f, 0x3f);
	  nt_query_left = _mm_and_si128(nt_query_left, temp);
	  //cote droit
	  nt_query_right = _mm_set_epi32(fright, fright, fright, fright);
	  nt_query_right = _mm_xor_si128(nt_query_right,vec_nt ); //xor avec masque lettres
	  temp = _mm_srli_epi32(nt_query_right, 1) ;  //shift de 1 vers droite
	  nt_query_right = _mm_or_si128(temp,nt_query_right); //ou
	  nt_query_right = _mm_andnot_si128(nt_query_right,vec_5);
	  //masque pour selec zone a lire
	  tairest =  min(right2,16); //lecture 16 max au total
	  masque = ~  (0xffffffff >> tairest >> tairest);
	  temp  = _mm_set_epi32(masque, masque, masque, masque);
	  nt_query_right = _mm_and_si128(nt_query_right,temp);
	  temp = _mm_and_si128(nt_query_right, vec_3);
	  nt_query_right = _mm_srli_epi32(nt_query_right, 2) ; 
	  nt_query_right = _mm_and_si128(nt_query_right, vec_3);
	  nt_query_right = _mm_add_epi32(nt_query_right, temp); 
	  temp =  _mm_srli_epi32(nt_query_right, 4) ; 
	  nt_query_right = _mm_add_epi32(nt_query_right, temp);
	  nt_query_right = _mm_and_si128(nt_query_right, vec_0f); 
	  temp =  _mm_srli_epi32(nt_query_right, 8) ; 
	  nt_query_right = _mm_add_epi32(nt_query_right, temp);
	  temp =  _mm_srli_epi32(nt_query_right, 16) ; 
	  nt_query_right = _mm_add_epi32(nt_query_right, temp);
	  //filtre final
	  temp  = _mm_set_epi32(0x3f, 0x3f, 0x3f, 0x3f);
	  nt_query_right = _mm_and_si128(nt_query_right, temp);
	  //fin calc nt query
	  /*

	  //calc sans sse
	  tairest =  min(right2,16);
	  masque = ~  (0xffffffff >> tairest >> tairest); 
  
	  ar = (int)popcount_mult (count_nt(fright,ma) & masque );
	  cr = (int)popcount_mult (count_nt(fright,mc) & masque );
	  tr = (int)popcount_mult (count_nt(fright,mt) & masque );
	  gr = tairest - ar-cr-tr;
	  
	  tairest =  min(left2,16);
	  masque = ~  (0xffffffff >> tairest >> tairest); 
	  
	  ale = (int)popcount_mult (count_nt(fleft,ma) & masque );
	  cl = (int)popcount_mult (count_nt(fleft,mc) & masque );
	  tl = (int)popcount_mult (count_nt(fleft,mt) & masque );
	  gl = tairest - ale-cl-tl;

	  */

	  // pourquoi pas s 2 ?
	  ////	  num_seed = codeSeedRight(&BK2->data[BK2->seq[nquery]+i],num_seed);
	  num_seed = codeSeedRight(s2 + i ,num_seed);
	  /// On parcourt toutes les occurences de la graine ds B1
	  if(num_seed!=-1)
	    {
	 
	      // ici acces tablke hachage
	      if(SIZE_SEED>14) 
		{
		  I1->get_seed_info_through_hashTable(num_seed,&nb_occur,&idx1);

		}
	      else
		{
		  nb_occur = I1->nb_seed[num_seed]; //invalid read //probleme avec nummseed peut etre vaut il -1 ?
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

		  absolute_pos = BK1->pos_seq[num_sequence] + offset_sequence - i ;


		  // filtre puis calcul
	    
		   left1 = offset_sequence;
		   right1 =  sizeseq - offset_sequence - SIZE_SEED;
		   s1 = (BK1->data) + BK1->seq[num_sequence];


		   // pour align global-global , faire test egalite.
		  if ((left1>=left2)&&(right1>=right2))
		    {
		      fleftB=seed_bank->left;
		      frightB=seed_bank->right;

		      //filtre maxpos en amont 
		      cpt_pos_amont++;
		      if ((cpt_pos_amont > MAXPOS_AMONT) && MAXPOS)
			goto pos_suivante;  
		      
			    //
			    
#ifdef CPT_FILT

		      cpt_hit++;
#endif

		      if(pre_filtre_align( left2,right2,max_mis,tabprec,fleftB,frightB,fleft,fright)) //TWN(4)
			{
#ifdef CPT_FILT

			  cpt_pref++;
#endif


			  if(filtre_align_NT_vec_sse(left2,right2,max_mis,fleftB, //filtre nt-frequency vectorise
   		  	     frightB,nt_query_left,nt_query_right,
			    		  	     vec_5,vec_3,vec_0f,vec_h1,vec_nt))
			    {
#ifdef CPT_FILT

			      cpt_ntsse++;
#endif

	      

			      if( Doub->getpos2(absolute_pos)) //doublons des mauvais
				goto occ_suivante ;
			      // filtre beaucoup, du coup moins arrive a maxpos
			      // donc sensib augmente car plus de hit parcourus, il faut diminuer les seuils de maxpos en consequence, ou pas 
			      Doub->setpos2(absolute_pos); // setpos du hit, on ne sait pas encore sil est bon ou mauvais
#ifdef CPT_FILT
			      cpt_doub++;
#endif
			  	  


			      //puis filtre diNT --> n apporte rien : filtre a peu pres 2 fois plus que nt sse mais est plus long 
			      //et ne filtre presque rien en plus que filtre_align
			      /*  	  if(filtre_diNT_DEC(left2,right2,max_mis
				  	     ,fleftB,frightB
				  	     ,fleft,fright   
				  	     ))
				  {
					  cpt_dint++;
			      */




			      if(filtre_align(left2,right2,max_mis,tabprec,fleftB,frightB,fleft,fright,num_gaps)) //TWN(16)
				{
#ifdef CPT_FILT

				  cpt_fil++;
#endif




				  if(! (Doub->getpos(absolute_pos))) // test position "fictive", cest a dire si pas de gaps // devenu obsolete il me semble
				    // car si passe pas, signifie  un bon align pos abso a deja ete trouve ici, donc deja filtre par getpos2 plus haut
				    {

				      //    #ifdef CPT_FILT
				      //		      cpt_doub++;
				      //#endif
				      //	{
				      //	      if (cpt_pos>MAXPOS && MAXPOS) goto pos_suivante; //MAXPOS==0 nolimit


					  //	  if(debug) printf("apres maxpos\n");





				      if(filtre_align_DC(s1+left1,s2+left2,left2,right2,max_mis,tabprec7, //TWN(full)
							 fleftB,frightB,fleft,fright,num_gaps ))
					{

					  #ifdef CPT_FILT
					  					  cpt_dint++;
					  #endif

					  al->init(doing_rev_query);
					  /// Enregistrement des indexs de début et de fin des séquences
					  al->setOffsets(offset_sequence-left2,i-left2,size_query,sizeseq,rev_comp);

					  /// On calcule le score de l'alignement
					  nb_mis = max_mis+1;
					  nb_mis = withgap_scoreHit(al,s1+left1,s2+left2,left1,left2,right1,right2,max_mis,num_seed,num_gaps); //WN reel
					  ///	  al->adjust_rev_comp(sizeseq,rev_comp);


					  if (nb_mis<=max_mis )
					    {
					      // il faut tout de meme faire les set pos celui fictif nest pas forcement a 1
					      if( Doub->getpos( BK1->pos_seq[num_sequence] + al->getStart1()-1))//ultime verif non doublon, avec vrai position
						{
						  Doub->setpos(absolute_pos);
						  goto occ_suivante; 
						}
					      
					      Doub->setpos(absolute_pos);
					      
					      //set pos vrai position 
					      Doub->setpos( BK1->pos_seq[num_sequence] + al->getStart1()-1);
					       
					      /// On récupère les index des commentaires des séquences
					      al->j1 = (BK1->data) + BK1->com[num_sequence];
					      al->j2 = (BK2->data) + BK2->com[nquery];
					
					      //on stocke le numero des sequences 
#ifdef NUMSEQSORTIE
					      al->n1=num_sequence; al->n2=nquery;
					      int zz;
					      for(zz=0; zz<BK1->num_part; zz++) al->n1+= BK1->nb_seq[zz];
#endif
					      /// Calcul du score
					      score = (size_query-nb_mis)*MATCH - nb_mis*MISMATCH;
					      /// Calcul de la E-value
					      al->e_value = ComputeEvalue(score,BK1->nb_tot_res,BK1->nb_tot_seq,size_query,1);

					      /// Affichage de l'alignement dans le fichier de sortie
		    
					      if (BESTAL)
						{
                          //  printf("adding al for query %i \n",nquery);
						  St->add_al(*al);
						  //check si quota atteint de scores optimaux dans ce cas
						  // il faut ecrire resu et goto query suivante
						  if ( St->good_enough())
						    {
                              //  printf("good enough \n",nquery);


						      nb_hit+=St->get_nb_al();
						      pthread_rwlock_wrlock(&verrou_out);
						      St->print_als(ff);
						      pthread_rwlock_unlock(&verrou_out);

#ifdef CPT_FILT
						      cpt_good++;
#endif
                                doing_rev_query = 1;
						      goto query_suivante;
						    }
						}
					      else
						{
						  St->incr_query();
						  /// Verrou pour protéger l'accès au fichier de sortie
						  pthread_rwlock_wrlock(&verrou_out);
						  display_withgap(ff,al,al->j1,al->j2,al->e_value);
						  pthread_rwlock_unlock(&verrou_out);
		    
						  ++nb_hit;
						  if (St->quota_atteint())  goto query_suivante;
						}
					    }
					} //fin filtre DC
					  // skip :
				      /////////cpt_pos ++; // ca changerait quoi ici ? 

				            if ((cpt_pos>MAXPOS) && MAXPOS)
					{
					  goto pos_suivante;  
					}
#ifdef CPT_FILT
					    cpt_max++; // ce compteur est approx, ya des hit qui arrivent pas ici  si doublon bon align
#endif
					  cpt_pos ++; // ici ou avant ? 

				      
				    }// fin doub
				   
				} // fin filtre
                              ///////}  // fin filtre di nt
			    } // fin filtre NT sse
			}

		      //	} // fin if doub
		    } //fin test s2 incluse ds s1
		occ_suivante:;
		  idx1++;
		} //fin boucle occurences dun graine dans B1




	    } // fin test graine valide (passe filtre complexite)

	pos_suivante:;
	 
	} //fin boucle  graines dans query
  
	  //on a tout parcouru pour cette query : si bestal il faut sortir les resu
      //si rev comp, sortie que si on a deja fait la reverse
      if((BESTAL && rev_comp==REV_COMP_ALIGN && doing_rev_query)  || (BESTAL && rev_comp==NORMAL_ALIGN) )
	{
        //printf("should output query  %i nbres %i \n",nquery,St->get_nb_al());
	  nb_hit+=St->get_nb_al();
	  pthread_rwlock_wrlock(&verrou_out);
	  St->print_als(ff);
	  pthread_rwlock_unlock(&verrou_out);
	}
      //remise a 0 du tab des doublons
    query_suivante: ;


      if(rev_comp==REV_COMP_ALIGN  )
	{
	  Doub->reset();
	  if(doing_rev_query)  St->reset();
	  doing_rev_query = !doing_rev_query;
	  nquery -= doing_rev_query;
	}
      else
	{
     doing_rev_query = 0;
	  Doub->reset();
	  St->reset();
	}
    } //fin boucle sur les  query


  delete [] s2_rev;

  delete al;
  return nb_hit;
}



/////////////: WN_ alignment grande seq (> 1000)




#define tab(i,j)  tab[(j)+ (i)*dim2]


/** 
 * Fonction de recherche du meilleur alignement selon le modèle de Wunsch et Needleman
 * \param al l'objet Alignment qui constituera le résultat
 * \param seq1 le pointeur du tableau des caractères de la première banque
 * \param seq2 le pointeur du tableau des caractères de la seconde banque
 * \param size1 le nombre de bases à parcourir dans la première séquence
 * \param size2 le nombre de bases à parcourir dans la seconde séquence
 * \param num_gaps le nombre maximal de gaps autorisés dans la partie de l'alignement
 * \param side indiquant si on fait l'alignement sur le côté gauche ou droit de la graine
 * \param num_seed le code de la graine courante
 * \param max_mis le nombre max de mismatch
 * \return AL_VALID si l'alignement trouvé est correct, AL_ERR si on n'a pas retenu l'alignement
 */
 
// avec allocation dans tas au lieu pile
int WN_big_alignment(Alignment* al, char* seq1, char* seq2, int size1, int size2, int num_gaps, int side, int num_seed, int max_mis)
{
  int i,j;

  int binf, bsup;
  int barriere;
  char* s1;
  char* s2;
  /// Les dimensions de la matrice d'alignement
  int dim1 = size1+1;
  int dim2 = size2+1;
  /// La matrice d'alignement, la première séquence constitue la première dimension
  // caseMatrix tab[dim1][dim2];
  caseMatrix * tab = (caseMatrix *) malloc(dim1*dim2*sizeof(caseMatrix));
  //printf(" big align  %i %i  ssizec %li  total %li Mo \n------------------\n",dim1,num_gaps,sizeof(caseMatrix),(dim1*dim2*sizeof(caseMatrix))/1024/1024 );

  /// Les tableaux de suivi de la construction de l'alignement
  char mem[size2+num_gaps];
  char mem2[size2+num_gaps];
	
  int cpt = 0;
	
  /// Variables utiles au contrôle de la validité de l'alignement
  int lastCode = num_seed;
  int L = SIZE_SEED;
  int max_gaps = num_gaps;


  if(side == LEFT)
    {

      /// On se place au début des séquences
      s1 = seq1 - size1;
      s2 = seq2 - size2;
    }
  else
    {

      /// Dans ce cas on inverse les séquences
      //  printf("seq1 %s seq2 %s \n",seq1,seq2);

      //   printf("WN align size1 %i  size2 %i ngaps %i side %i numseed %i  max_mis %i \n",size1,size2,num_gaps,side,num_seed,max_mis);
      s1 = new char[size1];
      s2 = new char[size2];
		
      for(i=0; i<size1; ++i)
	{
	  s1[i] = seq1[size1-i-1];
	}
      for(j=0; j<size2; ++j)
	{
	  s2[j] = seq2[size2-j-1];
	}
    }


  
  /// Remplissage de la matrice d'alignement

  //	Initialisation des cases de la matrice
  //trop long, et inutile 
//   for(i = 0; i < dim1; ++i)
//     {
//       for(j = 0; j < dim2; ++j)
// 	{
// 	  tab[i][j].score = 0;
// 	  tab[i][j].path = 99;
// 	}
//     }





  /// On initialise la première case de la matrice
  tab(0,0).score = 0;
  // On initialise les cases utiles de la première colonne
  // Ppour première colonne : gaps gratuits car differentes positions de depart possibles
  for(i=1; i <= min(num_gaps + num_gaps,size1); ++i) //bande de 2*numgaps de large
    {
      tab(i,0).score = 0;
      tab(i,0).path = GA2;
    }
  /// On initialise les cases utiles de la première ligne
  for(j=1; j <= min(num_gaps,size2); ++j)
    {
      tab(0,j).score = j * COUT_GAP;
      tab(0,j).path = GA1;
    }

  /// On remplit le reste de la matrice d'alignement
  /// On ne complète que les cases des diagonales respectant le nombre de gaps autorisés
  for(i = 1; i < dim1; ++i)
    {
      binf = max(i-num_gaps -num_gaps,1); //ok car on commence sur diag dacalee
      bsup = min(i,dim2-1); 
      barriere = 0;
      for(j = binf; j <= bsup; ++j)
	{

	  tab(i,j).score = -10000 ; 
	  tab(i,j).path = MES;
	  /// On vérifie l'appariement des caractères
	  if(identNT(s1[i-1],s2[j-1])
	     && tab(i-1,j-1).score >= -max_mis )
	    {
	      tab(i,j).score =  tab(i-1,j-1).score  + COUT_MATCH;
	      tab(i,j).path = APP;
	      barriere++;
	    }
	  else
	    if (   tab(i-1,j-1).score > -max_mis )

	      {
		tab(i,j).score =  tab(i-1,j-1).score  + COUT_MISMATCH;
		tab(i,j).path = MES;
		barriere++;
	      }
			
	  /// On teste la possibilité d'ajouter un gap
	  if(( tab(i-1,j).score + COUT_GAP > tab(i,j).score) && 
	     ((j-i) < 0) // test depassement bande
	     && tab(i-1,j).score > -max_mis )
	    {
	      tab(i,j).score =  tab(i-1,j).score + COUT_GAP;
	      tab(i,j).path = GA2;
	      barriere++;
	    }
	  if(( tab(i,j-1).score + COUT_GAP > tab(i,j).score) && 
	     ((i-j) < (2*num_gaps)) // depassement bande
	     && tab(i,j-1).score > -max_mis )
	    {
	      tab(i,j).score = tab(i,j-1).score  + COUT_GAP;
	      tab(i,j).path = GA1;
	      barriere++;
	    }
	}
      if(!barriere) { // tous les chemins menent a nombre d'erreurs > max_mis, on peut s'arreter
	if(side == RIGHT)
	  {
	    delete [] s1;
	    delete [] s2;
	    //   free(tab);
	  }

	free(tab);	
	return AL_ERR;

      }
    }


  --i;
  --j;
  /// On parcourt le chemin représentant le meilleur alignement dans la matrice
  /// On commence dans le coin en bas à droite de la matrice, correspondant aux derniers caractères , au niveau de la graine 
  /// On va continuer jusqu'à ce que la séquence la plus courte (la seconde) ait été complètement lue
  while(j>0)
    {

      switch(tab(i,j).path)
	{
	case APP:
	  --i;
	  --j;
	  mem[cpt]= s1[i];
	  mem2[cpt]= s2[j];
		  			
	  break;
	case MES:
	  --i;
	  --j;
	  mem[cpt]= s1[i];
	  mem2[cpt]= s2[j];
	  al->addMis();
	  L = 0;
	  lastCode = 0;
	  break;
	case GA2:
	  // gap sur query
	  if(side == LEFT) al->decStart1(); else al->incEnd1(); //ok
	  --i;
	  mem[cpt]= s1[i];
	  mem2[cpt]= CHAR_GAP;
	  al->addGap();
	  //  al->incAlign();
	  // al->decAlign();			
	  L = 0;
	  lastCode = 0;
				
	  break;
	case GA1:
	  // gap sur bank
	  if(side == LEFT) al->incStart1(); else al->decEnd1(); //ok
	  --j;
	  mem[cpt]= CHAR_GAP;
	  mem2[cpt]= s2[j];
	  al->addGap();
	  //  al->decAlign();
	  // al->incAlign();			
	  L = 0;
	  lastCode = 0;
				
	  break;
	}
		
      /// Si on a déjà utilisé trop de gaps, on arrête la construction de l'alignement
      if((al->getGaps() +al->getMis() )>max_mis)
	{
	  if(side == RIGHT)
	    {
	      delete [] s1;
	      delete [] s2;
	      //  free(tab);
	    }
	  free(tab);	
	  return AL_ERR;
	}


      if((al->getGaps() )>max_gaps)
	{
	  if(side == RIGHT)
	    {
	      delete [] s1;
	      delete [] s2;
	      // free(tab);
	    }
	  free(tab);	
	  return AL_ERR;
	}


      ++cpt;
    }


  /// On recopie ce chemin dans l'objet Alignement
  if(side == LEFT)
    {
      while(cpt>0)
	{
	  --cpt;
	  al->addPair(mem[cpt],mem2[cpt]);
	}
    }
  else
    {
      i = 0;
      while(i<cpt)
	{
	  al->addPair(mem[i],mem2[i]);
	  ++i;
	}
      delete [] s1;
      delete [] s2;
    }  
   free(tab);	 
  return AL_VALID;
}









// sur une bande : de 0 a 2ngaps , donc 2ng+1 cases,  +1 pour bordure init ( par ex pour  cas i,j = 0,1)   
#define gtab(i,j)  tab[(j-(i)+2*num_gaps)+ (i)*(2*num_gaps+2)]



/** 
 * Fonction de recherche du meilleur alignement selon le modèle de Wunsch et Needleman
 * \param al l'objet Alignment qui constituera le résultat
 * \param seq1 le pointeur du tableau des caractères de la première banque
 * \param seq2 le pointeur du tableau des caractères de la seconde banque
 * \param size1 le nombre de bases à parcourir dans la première séquence
 * \param size2 le nombre de bases à parcourir dans la seconde séquence
 * \param num_gaps le nombre maximal de gaps autorisés dans la partie de l'alignement
 * \param side indiquant si on fait l'alignement sur le côté gauche ou droit de la graine
 * \param num_seed le code de la graine courante
 * \param max_mis le nombre max de mismatch
 * \return AL_VALID si l'alignement trouvé est correct, AL_ERR si on n'a pas retenu l'alignement
 */
 
// avec allocation dans tas au lieu pile, et alloc bande uniquement  **en cours de test***
int WN_very_big_alignment(Alignment* al, char* seq1, char* seq2, int size1, int size2, int num_gaps, int side, int num_seed, int max_mis)
{
  int i,j;
  int binf, bsup;
  int barriere;
  char* s1;
  char* s2;
  /// Les dimensions de la matrice d'alignement
  int dim1 = size1+1;
  int dim2 = size2+1;
  /// La matrice d'alignement, la première séquence constitue la première dimension
  // caseMatrix tab[dim1][dim2];
  // caseMatrix * tab = (caseMatrix *) malloc(dim1*dim2*sizeof(caseMatrix));
  caseMatrix * tab = (caseMatrix *) malloc(dim1*(2*num_gaps +2)*sizeof(caseMatrix)); // alloc bande uniqement
  //  printf("very big align  %i %i  ssizec %li  total %li Mo \n------------------\n",dim1,num_gaps,sizeof(caseMatrix),(dim1*(2*num_gaps +2)*sizeof(caseMatrix))/1024/1024 );
  if(tab==NULL)
    {
      printf("malloc tab failed .. skipping this alignment, query is too large .\n");
      return AL_ERR;
    }
  // printf("size 1 %i  size2 %i \n",size1,size2);
  //printf("SEQ1 \n %s \n\n  sSEQ2\n %s \n\n  ",seq1,seq2);

  /// Les tableaux de suivi de la construction de l'alignement
  char mem[size2+num_gaps];
  char mem2[size2+num_gaps];
	
  int cpt = 0;
	
  /// Variables utiles au contrôle de la validité de l'alignement
  int lastCode = num_seed;
  int L = SIZE_SEED;
  int max_gaps = num_gaps;


  if(side == LEFT)
    {

      /// On se place au début des séquences
      s1 = seq1 - size1;
      s2 = seq2 - size2;
    }
  else
    {

      /// Dans ce cas on inverse les séquences
      //  printf("seq1 %s seq2 %s \n",seq1,seq2);

      //   printf("WN align size1 %i  size2 %i ngaps %i side %i numseed %i  max_mis %i \n",size1,size2,num_gaps,side,num_seed,max_mis);
      s1 = new char[size1];
      s2 = new char[size2];
		
      for(i=0; i<size1; ++i)
	{
	  s1[i] = seq1[size1-i-1];
	}
      for(j=0; j<size2; ++j)
	{
	  s2[j] = seq2[size2-j-1];
	}
    }


  
  /// Remplissage de la matrice d'alignement

  //	Initialisation des cases de la matrice
  //trop long, et inutile 
//   for(i = 0; i < dim1; ++i)
//     {
//       for(j = 0; j < dim2; ++j)
// 	{
// 	  tab[i][j].score = 0;
// 	  tab[i][j].path = 99;
// 	}
//     }





  /// On initialise la première case de la matrice
  gtab(0,0).score = 0;
  // On initialise les cases utiles de la première colonne
  // Ppour première colonne : gaps gratuits car differentes positions de depart possibles
  for(i=1; i <= min(num_gaps + num_gaps,size1); ++i) //bande de 2*numgaps de large
    {
      gtab(i,0).score = 0;
      gtab(i,0).path = GA2;
    }

  // sort de la bande ??
  /// On initialise les cases utiles de la première ligne
  // for(j=1; j <= min(num_gaps,size2); ++j)
//     {
//       gtab(0,j).score = j * COUT_GAP;
//       gtab(0,j).path = GA1;
//     }

//case 0,1
  j=1;i=0;
  gtab(i,j).score = j * COUT_GAP;
  gtab(i,j).path = GA1;

  /// On remplit le reste de la matrice d'alignement
  /// On ne complète que les cases des diagonales respectant le nombre de gaps autorisés
  for(i = 1; i < dim1; ++i)
    {
      binf = max(i-num_gaps -num_gaps,1); //ok car on commence sur diag dacalee
      bsup = min(i,dim2-1); 
      barriere = 0;
      for(j = binf; j <= bsup; ++j)
	{

	  gtab(i,j).score = -10000 ; 
	  gtab(i,j).path = MES;
	  /// On vérifie l'appariement des caractères
	  if(identNT(s1[i-1],s2[j-1])
	     && gtab(i-1,j-1).score >= -max_mis )
	    {
	      gtab(i,j).score =  gtab(i-1,j-1).score  + COUT_MATCH;
	      gtab(i,j).path = APP;
	      barriere++;
	    }
	  else
	    if (   gtab(i-1,j-1).score > -max_mis )

	      {
		gtab(i,j).score =  gtab(i-1,j-1).score  + COUT_MISMATCH;
		gtab(i,j).path = MES;
		barriere++;
	      }
			
	  /// On teste la possibilité d'ajouter un gap
	  if( ((j-i) < 0)  &&// test depassement bande
	      ( gtab(i-1,j).score + COUT_GAP > gtab(i,j).score) 
	      && gtab(i-1,j).score > -max_mis )
	    {
	      gtab(i,j).score =  gtab(i-1,j).score + COUT_GAP;
	      gtab(i,j).path = GA2;
	      barriere++;
	    }
	  if(
	     ((i-j) < (2*num_gaps)) && // depassement bande
	     ( gtab(i,j-1).score + COUT_GAP > gtab(i,j).score) && 
	     gtab(i,j-1).score > -max_mis )
	    {
	      gtab(i,j).score = gtab(i,j-1).score  + COUT_GAP;
	      gtab(i,j).path = GA1;
	      barriere++;
	    }
	}
      if(!barriere) { // tous les chemins menent a nombre d'erreurs > max_mis, on peut s'arreter
	if(side == RIGHT)
	  {
	    delete [] s1;
	    delete [] s2;
	    //   free(tab);
	    // printf("free tab \n");
	  }

	free(tab);// printf("free tab \n");
	return AL_ERR;

      }
    }


  --i;
  --j;
  /// On parcourt le chemin représentant le meilleur alignement dans la matrice
  /// On commence dans le coin en bas à droite de la matrice, correspondant aux derniers caractères , au niveau de la graine 
  /// On va continuer jusqu'à ce que la séquence la plus courte (la seconde) ait été complètement lue
  while(j>0)
    {

      switch(gtab(i,j).path)
	{
	case APP:
	  --i;
	  --j;
	  mem[cpt]= s1[i];
	  mem2[cpt]= s2[j];
		  			
	  break;
	case MES:
	  --i;
	  --j;
	  mem[cpt]= s1[i];
	  mem2[cpt]= s2[j];
	  al->addMis();
	  L = 0;
	  lastCode = 0;
	  break;
	case GA2:
	  // gap sur query
	  if(side == LEFT) al->decStart1(); else al->incEnd1(); //ok
	  --i;
	  mem[cpt]= s1[i];
	  mem2[cpt]= CHAR_GAP;
	  al->addGap();
	  //  al->incAlign();
	  // al->decAlign();			
	  L = 0;
	  lastCode = 0;
				
	  break;
	case GA1:
	  // gap sur bank
	  if(side == LEFT) al->incStart1(); else al->decEnd1(); //ok
	  --j;
	  mem[cpt]= CHAR_GAP;
	  mem2[cpt]= s2[j];
	  al->addGap();
	  //  al->decAlign();
	  // al->incAlign();			
	  L = 0;
	  lastCode = 0;
				
	  break;
	}
		
      /// Si on a déjà utilisé trop de gaps, on arrête la construction de l'alignement
      if((al->getGaps() +al->getMis() )>max_mis)
	{
	  if(side == RIGHT)
	    {
	      delete [] s1;
	      delete [] s2;
  

	    }
	  free(tab);	  //  printf("free tab \n");
	  return AL_ERR;
	}


      if((al->getGaps() )>max_gaps)
	{
	  if(side == RIGHT)
	    {
	      delete [] s1;
	      delete [] s2;
	      // free(tab);	    printf("free tab \n");

	    }
	  free(tab);	  //  printf("free tab \n");
	  return AL_ERR;
	}


      ++cpt;
    }


  /// On recopie ce chemin dans l'objet Alignement
  if(side == LEFT)
    {
      while(cpt>0)
	{
	  --cpt;
	  al->addPair(mem[cpt],mem2[cpt]);
	}
    }
  else
    {
      i = 0;
      while(i<cpt)
	{
	  al->addPair(mem[i],mem2[i]);
	  ++i;
	}
      delete [] s1;
      delete [] s2;
      //  free(tab);	    printf("free tab \n");

    }  
  free(tab);	  //  printf("free tab \n");
  return AL_VALID;
}






