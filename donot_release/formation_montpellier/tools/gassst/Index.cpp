/* This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software.  You can  use, 
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info". 
 */

/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file Index.cpp
 * \brief Classe Index, responsable de l'indexage des banques de séquences
 * \author Dominique Lavenier
 * \author Damien Fleury
 * \author Guillaume Rizk
 * \version 6.1
 * \date 21/01/2009
 */

#include "Index.h"

#include "Seed.h"
#include "Bank.h"
#include "misc.h"
#include "code.h"
#include "constants.h"


/**
 * Constructeur d'Index par défaut
 */
Index::Index()
{
  long long   nombreSeeds = nbSeeds();
  int i;
  /// On alloue de la mémoire pour les tableaux nb_seed et offset_seed de l'index
  if(SIZE_SEED<=14)
    {
      nb_seed = new int[nombreSeeds];
      offset_seed = new int[nombreSeeds];
    }
  else
    {
      storage = new Pool;
      mask1 = (1 << NT_HACH)-1;
      mask2 = mask1  << NT_HACH;
      mask3 = NHACH -1 ;;
      //printf("m1 m2 %llx %llx \n",mask1,mask2);
      dec_haschcode = 2*(SIZE_SEED-(NBITS_HACH/2));
      table_hachage = (cell **) malloc( NHACH * sizeof(cell * ));
      for (i=0; i<NHACH; i++) table_hachage[i]=NULL;
    }
}


/**
 * Destructeur d'Index
 */
Index::~Index()
{
  if(SIZE_SEED<=14)
    {
      delete [] offset_seed;
      delete [] nb_seed;
    }
  else
    {
      free(table_hachage) ;
      delete storage;
    }
}

void Index::free_hashtable()
{
  int i;
   storage->empty_all(); 
   //delete storage;
   //storage = new Pool;
   for (i=0; i<NHACH; ++i) 
     {
       table_hachage[i]=NULL; 
     }



//   int i;
//   cell *  cur,*prec;
//   for (i=0; i<NHACH; ++i) 
//     {
//       cur = table_hachage[i]; 
//       prec=NULL;

//       while(cur!=NULL )
// 	{
// 	  prec = cur; cur = cur->suiv;
// 	  free(prec);
// 	}
//       table_hachage[i]=NULL; 
//     }

}


/**
 * Constructeur d'Index par recopie
 * \param i un objet Index
 */
Index::Index(const Index& i)
{
  *this = i;
}


/**
 * Opérateur d'affectation de Hit
 * \param i un objet Index
 * \return l'objet Index affecté
 */
Index& Index::operator=(const Index& i)
{
  if(this!=&i)
    {
      offset_seed = i.offset_seed;
      nb_seed = i.nb_seed;
      seed = i.seed;
      table_hachage = i.table_hachage;
		
    }
  return *this;
}


/**
 * Méthode permettant d'indexer une banque
 * \param B le pointeur de la banque à indexer
 */
void Index::indexBank(Bank *B,int index_stride)
{
  int i, j, k;
  long long  nombreSeeds = nbSeeds(),x;

  if(SIZE_SEED<=14)
    {
      /// Comptage des graines

      k = 0;
      /// Initialisation de l'index
      for (i=0; i<nombreSeeds; ++i)
	{
	  nb_seed[i]=0;
	}


      /// On parcourt les séquences de la partition
      for (i=0; i<(B->nb_seq[B->num_part]); ++i)
	{
	  x = -1; 
	  /// On parcourt les lettres de la séquences pour compter les graines
	  for (j=0; j<B->size[i]-SIZE_SEED+1; j+=index_stride)
	    {
	      /// Calcul du code de la graine
	      //	x = codeSeedRight(&B->data[B->seq[i]+j],x);
	      x = codeSeedRight_stride(&B->data[B->seq[i]+j],x,index_stride);

	      if (x>=0)
	  	{
		  ++nb_seed[x]; // chercher graine dans table hach puis  inc 
		  ++k;
	  	}
	    }
	}
      /// Fin du comptage des graines

      /// Allocation de la mémoire pour le tableau de graines de l'index
      seed = new Seed[k];

      /// Initialisation des index des graines dans le tableau seed
      offset_seed[0]=0;
      for (i=1; i<nombreSeeds; ++i) 
	{
	  offset_seed[i] = offset_seed[i-1] + nb_seed[i-1];
	  nb_seed[i-1]=0;
	}
      nb_seed[nombreSeeds-1]=0;


      // Etape d'indexation
      unsigned int left=0;
      unsigned int right=0;
      int tair=0;
      /// On parcourt les séquences de la partition
      for (i=0; i<B->nb_seq[B->num_part]; ++i)
	{
	  x = -1;
	  left=0;
	  right=0;
	  /// On parcourt les lettres de la séquence pour compter les graines
	  for (j=0; j<B->size[i]-SIZE_SEED+1; j+=index_stride)
	    {

	      /// Calcul du code de la graine
	      x = codeSeedRight_stride(&B->data[B->seq[i]+j],x,index_stride);

	      if (j==0)
		{
		  tair=min(4*sizeof(unsigned int),B->size[i]-j-SIZE_SEED);

		  right =  seq2code_special(&B->data[B->seq[i]+j] + SIZE_SEED,tair,4*sizeof(unsigned int)) ; // 4 charac par octets4*sizeof(unsigned int)
		}
	      else if (j<(B->size[i]-SIZE_SEED-4*(int)sizeof(unsigned int)+1))
		{

		  left = (unsigned int)seq2codeLeft_stride(&B->data[B->seq[i]+j-1],left,index_stride);
		  right = (unsigned int) seq2codeRight_stride(&B->data[B->seq[i]+j+SIZE_SEED],4*sizeof(unsigned int),right,index_stride);  

		}
	      else // on est au bout on arrete de lire a droite
		//sauf dans le cas de stride > 1
		//il faut que tout soit deja lu ! 
		{
		   
		  left = (unsigned int)seq2codeLeft_stride(&B->data[B->seq[i]+j-1],left,index_stride);
		  //  right = (right <<2);  
		  if (index_stride<2)
		    {
		      right = (right <<2); 
		    }
		  else
		    {
		      tair=min(4*sizeof(unsigned int),B->size[i]-j-SIZE_SEED);
		      right = (unsigned int) seq2codeRight_stride_tail(&B->data[B->seq[i]+j+SIZE_SEED],4*sizeof(unsigned int),right,INDEX_STRIDE,tair);  

		    } 

		}
	  
	      if (x>=0) 
		{

		  k = offset_seed[x] + nb_seed[x];
		  /// On indique le numéro de séquence et l'index de la graine identifiée
		  seed[k] = Seed(i,j,left,right);
		  // seed[k] = Seed(i,j);
		  ++nb_seed[x];
		}
	    

	    }
 


	}
      /// Fin de l'indexation


    }

  else  // grande graine, avec table de hachage
    {

      /// Comptage des graines
      //  char seq[16];

      k = 0;
      int clef;
      cell * newcell, *cur,*prec;
      /// On parcourt les séquences de la partition
      for (i=0; i<(B->nb_seq[B->num_part]); ++i)
	{
	  x = -1;
	  // printf("comptage graines contig %i \n",i);
	  /// On parcourt les lettres de la séquences pour compter les graines
	  for (j=0; j<B->size[i]-SIZE_SEED+1; j+=index_stride)
	    {
	      /// Calcul du code de la graine
	      //	x = codeSeedRight(&B->data[B->seq[i]+j],x);
	      x = codeSeedRight_stride(&B->data[B->seq[i]+j],x,index_stride);
	      //	  code2seq(x,seq,SIZE_SEED);
	      // seq[16]='\0';


	      if (x>=0)
	  	{



		  clef = hashCode(x,NHACH);
		  if (table_hachage[clef]==NULL) //clef contient aucune graine
		    {
		      //    printf("%s vide \n",seq);
		      //     newcell = (cell*) malloc(sizeof(cell)); // faire avec memory pool
		      newcell = storage->allocate_cell_in_pool(); 
		      newcell->nb_seed=1; newcell->suiv=NULL; newcell->graine=x;
		      table_hachage[clef] = newcell;
		    }
		  else
		    {
		      //parcours liste
		      cur = table_hachage[clef]; prec=NULL;
		      while(cur!=NULL &&  cur->graine < x) // la 
			{
			  prec = cur;
			  cur = cur->suiv;
			}
		      if (cur==NULL) //graine non trouvée , insertion a la fin
			{
			  //   printf("%s inser fin \n",seq);

			  //	  newcell = (cell*) malloc(sizeof(cell));// faire avec memory pool
		      newcell = storage->allocate_cell_in_pool(); 

			  newcell->nb_seed=1; newcell->suiv=NULL;newcell->graine=x;
			  prec->suiv = newcell;
			}
		      else if ( cur->graine==x )  (cur->nb_seed)++;  // graine trouvée
		      else // insertion graine  avant cur
			{
			  //    printf("%s inser avant cur \n",seq);

			  //  newcell = (cell*) malloc(sizeof(cell)); // faire avec memory pool
			  newcell = storage->allocate_cell_in_pool(); 

			  newcell->nb_seed=1; newcell->suiv=cur;newcell->graine=x;
			  (prec == NULL ) ? table_hachage[clef]= newcell  :  prec->suiv = newcell;
			}

		    }

		  // ++nb_seed[x]; // chercher graine dans table hach puis  inc  fait
      
		  ++k;
	  	}
	    }
	}
      /// Fin du comptage des graines

      /// Allocation de la mémoire pour le tableau de graines de l'index
      seed = new Seed[k];

      ///printf("init graines offset  \n");

      /// Initialisation des index des graines dans le tableau seed avec parcours de table hachage
      int offset_courant=0;
      //les graines sont triees
      for (i=0; i<NHACH; ++i) 
	{
	  cur = table_hachage[i];

	  while(cur!=NULL )
	    {
	      // code2seq(cur->graine,seq,SIZE_SEED);
	      //  seq[16]='\0';
	      //  printf("%s init \n",seq);

	      cur->offset_seed = offset_courant;
	      offset_courant += cur->nb_seed;
	      cur->nb_seed = 0; // reset pour apres
	      cur = cur->suiv;
	    }
	  
	}

 

      // Etape d'indexation
      unsigned int left=0;
      unsigned int right=0;
      int tair=0;
      /// On parcourt les séquences de la partition
      for (i=0; i<B->nb_seq[B->num_part]; ++i)
	{
	  ///	  printf("index contig %i  \n",i);

	  x = -1;
	  left=0;
	  right=0;
	  /// On parcourt les lettres de la séquence pour compter les graines
	  for (j=0; j<B->size[i]-SIZE_SEED+1; j+=index_stride)
	    {

	      /// Calcul du code de la graine
	      x = codeSeedRight_stride(&B->data[B->seq[i]+j],x,index_stride);

	      if (j==0)
		{
		  tair=min(4*sizeof(unsigned int),B->size[i]-j-SIZE_SEED);

		  right =  seq2code_special(&B->data[B->seq[i]+j] + SIZE_SEED,tair,4*sizeof(unsigned int)) ; // 4 charac par octets4*sizeof(unsigned int)
		}
	      else if (j<(B->size[i]-SIZE_SEED-4*(int)sizeof(unsigned int)+1))
		{

		  left = (unsigned int)seq2codeLeft_stride(&B->data[B->seq[i]+j-1],left,index_stride);
		  right = (unsigned int) seq2codeRight_stride(&B->data[B->seq[i]+j+SIZE_SEED],4*sizeof(unsigned int),right,index_stride);  

		}
	      else // on est au bout on arrete de lire a droite
		//sauf dans le cas de stride > 1
		//il faut que tout soit deja lu ! 
		{
		   
		  left = (unsigned int)seq2codeLeft_stride(&B->data[B->seq[i]+j-1],left,index_stride);
		  //  right = (right <<2);  
		  if (index_stride<2)
		    {
		      right = (right <<2); 
		    }
		  else
		    {
		      tair=min(4*sizeof(unsigned int),B->size[i]-j-SIZE_SEED);
		      right = (unsigned int) seq2codeRight_stride_tail(&B->data[B->seq[i]+j+SIZE_SEED],4*sizeof(unsigned int),right,INDEX_STRIDE,tair);  

		    } 

		}
	  
	      if (x>=0) 
		{

		  /////////////////////////



		  clef = hashCode(x,NHACH);
		  if (table_hachage[clef]==NULL) //clef    non trouvee         , probleme
		    {
		      printf("Probleme critique, clef non trouvee\n");
		    }
		  else
		    {
		      //parcours liste
		      cur = table_hachage[clef]; prec=NULL;
		      while(cur!=NULL &&  cur->graine !=x)
			{
			  prec = cur;
			  cur = cur->suiv;
			}
		      if (cur==NULL) 
			{
			  printf("Probleme critique, clef non trouvee\n");                          
			}
		      else // graine trouvee
			{
			  k = cur->offset_seed + cur->nb_seed;
			  seed[k] = Seed(i,j,left,right);
			  (cur->nb_seed)++; 
			}
		    }

		}
	    

	    }
 


	}
      /// Fin de l'indexation








    }
}


inline int Index::hashCode(long long code, int max)
{
  //return ((code & mask1) |  ((code >> dec_haschcode) & mask2  ));
  //  return (code >> (SIZE_SEED-NT_HACH) );
  //  return (code >> dec_haschcode );
  //return ((code >> (dec_haschcode/2) ) & 268435455LL );
  // return (( ((code & 262143LL) * (code >> 18 )) + code)  & 268435455LL); // best
  // return ( ((code & 262143LL) * (code >> 18 ))  & 268435455LL);

  //  return (((code + (code << 3)) + (code << 12))  & 268435455LL ) ;

  //   long long a  = code;
  //   a = (a+0x7ed55d16) + (a<<12);
  //    a = (a^0xc761c23c) ^ (a>>19);
  //    a = (a+0x165667b1) + (a<<5);
  //    a = (a+0xd3a2646c) ^ (a<<9);
  //    a = (a+0xfd7046c5) + (a<<3);
  //    a = (a^0xb55a4f09) ^ (a>>16);
  //    return a;


  code = (~code) + (code << 18); 
  code = code ^ (code >> 31);
  code = code * 21; 
  code = code ^ (code >> 11);
  code = code + (code << 6);
  code = code ^ (code >> 22);
  return ((int) (code  & mask3) );
  

}





void Index::get_seed_info_through_hashTable(long long seed, int * nb_occur, int * offset_seed )
{

  int clef;
  cell  *cur;
  
  //char seq[19];
  //code2seq(seed,seq,SIZE_SEED);
  //seq[18]='\0';

  clef = hashCode(seed,NHACH);
  // printf(" seed %lli  %s \n",seed,seq);

  *nb_occur = 0;
  // printf("clef  %i  \n",clef);

  //parcours liste
  cur = table_hachage[clef]; 
  while(cur!=NULL &&  cur->graine != seed)
    cur = cur->suiv;
  if (cur!=NULL)
    {
      *nb_occur = cur->nb_seed;
      *offset_seed = cur->offset_seed;
    }
  

}






void Index::printstat()
{
  printf("----------------------Stat Index ---------------------\n");
  //printf("tai cell %lu  cell*: %lu  \n",sizeof(cell),sizeof(cell*));

  int ma=0,mi=99999,cpt=0;
  int i;
  int maxclef=0;
  cell  *cur;
  int distrib_colli[512];
  int nb_cell=0;
  for (i=0;i<512;i++){distrib_colli[i]=0;}

  for (i=0; i<NHACH; ++i) 
    {
      cur = table_hachage[i];
      cpt=0;
      while(cur!=NULL )
	{
	  nb_cell++;
	  cpt++;
	  cur = cur->suiv;
	}
      if(cpt>ma) {ma=cpt;maxclef=i;}
      ma = max(ma,cpt); mi = min(mi,cpt);
      distrib_colli[cpt]++;
    }

  printf("max collisions = %i pour clef %i   nb_cell total %i  %g entree/graine \n",ma,maxclef,nb_cell,(float)(NHACH-distrib_colli[0])/(float)nb_cell);
  for (i=0; i<10; i++) 
    {
      printf("  %9i collisions :  %9i  \n",i,distrib_colli[i]);
    }


  //  cur = table_hachage[maxclef];
  //   char seq[SIZE_SEED+1];
  //   while(cur!=NULL )
  //     {
  //       code2seq(cur->graine,seq,SIZE_SEED);
  //       seq[SIZE_SEED]='\0';
  //       printf("graine %lli seq %s    nbseed %i   offsetseed %i  \n",cur->graine, seq ,cur->nb_seed,cur->offset_seed);
  //       // printf(" seqg %s  code %lli \n",seq, cur->graine );
  //       cur = cur->suiv;  
  //     }
  
}


// void Index::printIndex()
// {
//   int i,len=0,j;
//   int nombreSeeds = nbSeeds();
//   char seq[10];
//   char flank[17];

//   j=0;
//   for (i=0; i < nombreSeeds ; i++)
//     {
//       if(nb_seed [i] !=0)
// 	{
// 	  code2seq(i,seq,SIZE_SEED);
// 	  seq[9]='\0';
// 	  printf("----------------------Graine %i  : %s---------------------\n",i,seq);
// 	  len +=  nb_seed [i];
	  
// 	  for (; j < len ; j++)
// 	    {
// 	      // printf("--------Occurrence  %i -------------\n",j);
// 	      printf("num_seq %i \n",seed[j].num_seq);
// 	      printf("off_seq %i \n",seed[j].off_seq);
// 	      code2seq(seed[j].left,flank,16);
// 	      flank[17]='\0';

// 	      printf("left %u  %s\n",seed[j].left,flank);
// 	      code2seq(seed[j].right,flank,16);
// 	      flank[17]='\0';

// 	      printf("right %u %s \n",seed[j].right,flank);
	      
// 	    }
// 	}
  
//     }


//   // for (i=0; i < len ; i++)
//   //     {
//   //       printf("-----Graine %i -----\n",i);
//   //       printf("num_seq %i \n",seed[i].num_seq);
//   //       printf("off_seq %i \n",seed[i].off_seq);
//   //       printf("left %u \n",seed[i].left);
//   //       printf("right %u \n",seed[i].right);

//   //     }


// }
