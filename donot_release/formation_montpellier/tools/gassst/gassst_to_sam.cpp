
/* This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software.  You can  use,
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 */

/**
 * Logiciel gassst_to_sam (conversion sortie gassst vers sam)
 * \file gassst_to_sam.c
 * \brief Module principal du programme
 * \author Guillaume Rizk
 * \date 01/02/2011
 */





#define min2(a, b) ((a) < (b) ? (a) : (b))




//#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "constants.h"

#include "resu.cpp"
#include <list>

#define NHACHG2S 200000
#define max(a,b) (a > b ? a : b)
#define NA 9999

#define NCHAR 100

#define NHACH_READS 2000000

//typedef unsigned char bool;





typedef struct cell
{
  char * nom;
  char * seq;
  //données du meilleur align 
  char * nom_contig;
  char * cigar;
  unsigned int  pos; // position du meilleur align
  unsigned int nb_gaps;
  unsigned int nb_mis;
 unsigned int score_al;
  unsigned int edit;
  int forward;
  unsigned int edit_second; //edit distance second meilleur
  //  unsigned int nb_gaps_sec; // gaps et misma du second meilleur
  //unsigned int nb_mis_sec;
  unsigned int score_al_sec;
  unsigned int nb_align; //nb align total du read
  unsigned int nb_best;  // nb  align best
  unsigned int tai_read;
  struct cell *suiv;
} cell;



typedef struct ali
{
  char * nom;
  char * seq;
  //données du meilleur align 
  char * nom_contig;
  char * cigar;
  unsigned int  pos; // position du meilleur align
  unsigned int nb_gaps;
  unsigned int nb_mis;
  unsigned int score_al;
  unsigned int edit;
  int forward;
  unsigned int tai_read;
} ali;








//cellule table hachage pour index fichier resultat 
typedef struct hashcell
{
  unsigned int offset; //4
  unsigned int nb; //4
  char * nom;
  struct hashcell *suiv;
} hashcell;



/// tableau index des reads ..
//unsigned int * index_reads;
fpos_t * index_reads;


#define GMATCH  10
#define GMISMATCH 15
#define GGAP 20


inline unsigned char compute_map_qual(cell * align)
{

  //  printf("comp map qual \n ");
  //printf("nb best %i  edit second %i \n ",align->nb_best,align->edit_second );
  //printf("al %i  al sec%i \n ",align->score_al ,align->score_al_sec  );

  if(align->nb_best > 1 ) return 0;
  else
    {

      if (align->edit_second == NA ) return 50; // combien  ?


      unsigned int S1 = align->score_al ;
      unsigned int S2 = align->score_al_sec;

      return (  (unsigned char) ( 250*  ( ((double)(S1-S2)) / ((double)(S1)) ))    );

    }

}






inline unsigned char compute_map_qual_multi(unsigned int nb_best, unsigned int S1, unsigned int S2)
{

  // printf("map qual multi:    nbbest %i  S1 %i  S2 %i    resu %u \n",nb_best,S1,S2,(unsigned char) ( 250*  ( ((double)(S1-S2)) / ((double)(S1)) )));
  //  printf("comp map qual \n ");
  //printf("nb best %i  edit second %i \n ",align->nb_best,align->edit_second );
  //printf("al %i  al sec%i \n ",align->score_al ,align->score_al_sec  );
  if(nb_best > 1 ) return 0;
  else
    {

      // if (nb_align==1) return 50; // combien  ? // cas traite en amont 


      //  unsigned int S1 = align->score_al ;
      //  unsigned int S2 = align->score_al_sec;

      return (  (unsigned char) ( 250*  ( ((double)(S1-S2)) / ((double)(S1)) ))    );

    }

}



//necessite champ nb gaps et nb_mis et tai_read
inline unsigned int  compute_score(unsigned int ngaps, unsigned int nmis, unsigned int nmatch)
{

  return ( GMATCH*nmatch - GMISMATCH * nmis - GGAP *ngaps ) ;

}





inline unsigned short compute_flag(ali * align)
{

  if(align->forward ) return 0;
  else return 16;
}


inline void free_hashcell(hashcell * read)
{
  free(read->nom);
  free(read);
}


inline void free_cell(cell * align)
{
  free (align->seq);
  free (align->nom);
  free (align->nom_contig);
  free (align->cigar);
  free(align);

}

inline void free_ali(ali * align)
{
  free (align->seq);
  free (align->nom);
  free (align->nom_contig);
  free (align->cigar);
  // free(align);

}

int lecture_read_align(	FILE* in,char * nom_read);

int lecture_align(FILE* ,char *,char *,unsigned int *,unsigned int *,unsigned int *,unsigned int *,unsigned int *,char *,int *, char * );
int hashCode(char * nom_seq,int max);


int main(int argc, char** argv)
{

  if (argc!=3 && argc!=4){
    printf(" Usage :\n");
    printf(" %s   <Gassst alignment file>   <output file in sam format>  [nb_aligns_per_read] \n",argv[0]);
    printf(" Default behavior outputs all alignments of a read.\nIf [nb_aligns_per_read] is specified, outputs the \"nb_aligns_per_read\" best alignments of each read  \n");

    return EXIT_FAILURE;
  }


  FILE  *align_file,*sam;
  int i;
  int nb_res_per_read=0; // valeur par defaut : output tous les align

  if (!(align_file = fopen(argv[1], "r")))
    {
      fprintf  (stderr, "cannot open %s\n exit...\n", argv[1]); exit (0);
    }
  
  
  sam = fopen(argv[2], "w");
  
  if (argc==4)   nb_res_per_read = atoi(argv[3]);
  
  if(nb_res_per_read<=0 && argc==4)
    {
      printf("Number of reported alignments per read you specified (%i) makes no sense.\n",nb_res_per_read);
      printf("Switching to defaut behavior: output all alignments of each read");
      nb_res_per_read=0;

    }
  char cigar[1000];
  char nom_read[1000];
  char nom_contig[1000];
  char * seq_read;
  seq_read= (char *) malloc(3000000 * sizeof(char));
  // char seq_read[5000];



  // printf("sizeof pos %lu bytes  %lu  \n",sizeof(fpos_t), sizeof(unsigned int));
  //exit(0);
  //////////SAM header 

  fprintf(sam,"@HD\tVN:0.1.2\tSO:unsorted\n");
  fprintf(sam,"@PG\tID:GASSST\tVN:%s\n",GASSST_VERSION);

  char c; 
  //puis recopie  header  dictionnaire seq
  while(  (c = fgetc(align_file))== '@'  )
    {
      char * line = NULL;size_t len = 0;
      getline(&line, &len, align_file);
      fprintf(sam,"@%s",line);
    }
  
  fseek(align_file,-1,SEEK_CUR);
  
  ///////////

  // Etape 1 : scan du fichier de resultats, compte nb res par reads
  int k ;
  int clef;
  hashcell * newcell, *cur;
  hashcell ** table_hachage_reads;
  unsigned int    pos, gaps,mis,edit_dist,tai_read;
  int forward;

  edit_dist=0;
  table_hachage_reads = (hashcell **) malloc( NHACH_READS * sizeof(hashcell * ));
 for (i=0; i<NHACH_READS; i++) table_hachage_reads[i]=NULL;
     k = 0;

  while( lecture_read_align(align_file, nom_read) )
    {


     clef = hashCode(nom_read,NHACH_READS);
     if (table_hachage_reads[clef]==NULL) // clef contient aucune cellule
       {
	 newcell = (hashcell*) malloc(sizeof(hashcell));
	 newcell->nb = 1;
	 newcell->nom =  (char *)  malloc((strlen(nom_read)+1)*sizeof(char));
	 strcpy(newcell->nom, nom_read);
	 newcell->suiv=NULL; //this was missing !
	 table_hachage_reads[clef] = newcell;
       }
     else
       {
	 //parcours liste
	 cur = table_hachage_reads[clef];
	 while(cur!=NULL && strcmp(cur->nom,nom_read))
	   cur = cur->suiv;

	  if (cur==NULL) //read non trouvé
	    {
	      newcell = (hashcell*) malloc(sizeof(hashcell));
	      newcell->nb = 1;
	      newcell->nom =  (char *)  malloc((strlen(nom_read)+1)*sizeof(char));
	      strcpy(newcell->nom, nom_read);

	      newcell->suiv = table_hachage_reads[clef]; // push front 
	      table_hachage_reads[clef] = newcell;
	    }
	  else //read trouvé
	    {
	      (cur->nb)++;
	    }


       }

     k++;

    }


  // fin etape 1 ,comptage du nb de res par reads



  // etape 2 , calcul des index des reads dans la table index_reads

     /// Allocation de la mémoire pour le tableau de graines de l'index
  index_reads = (fpos_t *)   malloc(k*sizeof(fpos_t));

      unsigned int offset_courant=0;
      for (i=0; i<NHACH_READS; ++i) 
	{
	  cur = table_hachage_reads[i];
	  while(cur!=NULL )
	    {
	      cur->offset= offset_courant;
	      offset_courant += cur->nb;
	      cur->nb = 0; // reset pour apres
	      cur = cur->suiv;
	    }
	  
	}

 
      //etape 3, re-scan du fichier resultat, et remplissage tableau index_reads

  fseek(align_file,0,SEEK_SET);
  //passe le header du fichier ...
  while(  (c = fgetc(align_file))== '@'  )
    {
      char * line = NULL;size_t len = 0;
      getline(&line, &len, align_file);
    }  
  fseek(align_file,-1,SEEK_CUR);
  
  fpos_t pos_courante,pos_suivante;
  fgetpos(align_file,&pos_courante); // pointe vers premiere ligne 

  while( lecture_read_align(align_file, nom_read) )
    {
      fgetpos(align_file,&pos_suivante); // pointe vers la ligne suivante

      clef = hashCode(nom_read,NHACH_READS);
      if (table_hachage_reads[clef]==NULL) //clef    non trouvee, probleme
	{
	  printf("Probleme critique, clef non trouvee\n");
	}
      else
	{
	  //parcours liste
	  cur = table_hachage_reads[clef];
	  while(cur!=NULL && strcmp(cur->nom,nom_read))
	    cur = cur->suiv;
	  if (cur==NULL) //read non trouvé
	    {
	      printf("Probleme critique, clef non trouvee\n");
	    }
	  else //read trouve dans la hashtable
	    {
	      k = cur->offset + cur->nb;
	      index_reads[k]=pos_courante;
	      // seed[k] = Seed(i,j,left,right);
	      (cur->nb)++; 
	    }


	}
      pos_courante = pos_suivante;
    }


  /// derniere etape: scan de la hashtable pour conversion en SAM:
  fpos_t pos_res;
  unsigned int offset_read;
  unsigned int nb_res;
  unsigned int j;
    hashcell * cell_used;

    ali * tab_ali;
    std::list<resu> listali;
    
    ali * newali, *ali_courant, *ali_suivant;
    std::list<resu>::iterator p;
    
 unsigned int nb_best;
 unsigned int nb_align;
 int nb_sam;



 /*
 //////////////avec respect ordre des query : 
 //repositionnement dans align_file :
  fseek(align_file,0,SEEK_SET);
  //passe le header du fichier ...
  while(  (c = fgetc(align_file))== '@'  )
    {
      char * line = NULL;size_t len = 0;
      getline(&line, &len, align_file);
    }  
  fseek(align_file,-1,SEEK_CUR);
  //////////////
  while( lecture_read_align(align_file, nom_read) )
    {
 clef = hashCode(nom_read,NHACH_READS);
      if (table_hachage_reads[clef]==NULL) //clef    non trouvee, probleme
	{
	  printf("Probleme critique, clef non trouvee\n");
	}
      else
	{
	  //parcours liste
	  cur = table_hachage_reads[clef];
	  while(cur!=NULL && strcmp(cur->nom,nom_read))
	    cur = cur->suiv;
	  if (cur==NULL) //read non trouvé
	    {
	      printf("Probleme critique, clef non trouvee\n");
	    }
	  else //read trouve dans la hashtable
	    {

	    ///traitement pour un read ......


	    //printf("%s     :   nb %u  offset %u \n",cur->nom,cur->nb, cur->offset);
	    nb_res = cur->nb; offset_read = cur->offset;
	    tab_ali = (ali * ) malloc(nb_res * sizeof(ali)); // stockage pour les  nb_res alignements du read

	    nb_align= nb_res;
	    //remplissage en memoire
	    for (j=0; j<nb_res; j++)
	      {
		pos_res = index_reads[offset_read+j];
		fsetpos(align_file,&pos_res);// placement au bon endroit dans le fichier 
		// puis lecture de la ligne correspondante
		lecture_align(align_file, nom_read,nom_contig, &gaps,&mis,&edit_dist,&pos,&tai_read,seq_read,&forward, cigar);

		//stockage des donnees en memoire : 
		newali = tab_ali +j;
		unsigned int nscore = compute_score(gaps,mis,tai_read-edit_dist);
		
		
		//newali = (cell*) malloc(sizeof(cell));
		newali->nom =  (char *)  malloc((strlen(nom_read)+1)*sizeof(char));
		strcpy(newali->nom,nom_read);
		newali->tai_read=tai_read;
		newali->seq = (char *)  malloc((tai_read+1)*sizeof(char));
		strcpy(newali->seq,seq_read);
		//données bestal
		newali->nom_contig =  (char *)  malloc((strlen(nom_contig)+1)*sizeof(char));
		newali->cigar =  (char *)  malloc((strlen(cigar)+1)*sizeof(char));
		strcpy(newali->nom_contig,nom_contig);
		strcpy(newali->cigar,cigar);
		newali->pos=pos; newali->nb_gaps=gaps; newali->nb_mis=mis;  newali->edit=edit_dist;
		newali->score_al=nscore;
		newali->forward = forward; 
		
		listali.push_back(resu(nscore,j));
      		//printf("--------- Resultat  %i --------------\n",j+1);
		//		printf("%s pos  :%u\n",nom_read,pos);
		//printf("-------------------------------------\n");
	      }

	    // tri des n align du read  :
	    listali.sort();

	    //premier parcours pour nb_best
	   p = listali.begin();
	   nb_best=0;
	   unsigned int prev_score;
	   prev_score= p->score_al;

	   while(p != listali.end() && (prev_score == p->score_al) ) {
	     
	     nb_best++;
	     prev_score = p->score_al;
	     p++;
	   }

	   // nb_best est calculé 
	   nb_sam = min2(nb_res,(unsigned int)nb_res_per_read);
	   if (nb_res_per_read==0) nb_sam= nb_res;

	   p = listali.begin();
	   int cpt=0;
	   while(p != listali.end() && (cpt<nb_sam) ) {
	     ali_courant = &(tab_ali[p->index]); 
	    

	     fprintf(sam,"%s\t",ali_courant->nom);
	     fprintf(sam,"%u\t",compute_flag(ali_courant)); 
	     fprintf(sam,"%s\t",ali_courant->nom_contig);
	     fprintf(sam,"%u\t",ali_courant->pos);
	     if(cpt!=0) // tous les align subopt : mapqual = 0
	       fprintf(sam,"%u\t",0); //mapping quality
	     else
	       {
		 // si unique : 50 
		 if (nb_align==1) 
		   fprintf(sam,"%u\t",50); //mapping quality 50 si best et seul align
		 else
		   {
		     p++;
		     ali_suivant = &(tab_ali[p->index]); 
		     p--;
		     fprintf(sam,"%u\t",compute_map_qual_multi(nb_best,ali_courant->score_al,ali_suivant->score_al)); //mapping quality

		   }
	       }
	     fprintf(sam,"%s\t",ali_courant->cigar); // CIGAR
	     fprintf(sam,"*\t"); // MRNM
	     fprintf(sam,"0\t"); // MPOS
	     fprintf(sam,"0\t"); // ISIZE
	     fprintf(sam,"%s\t",ali_courant->seq); // SEQ
	     fprintf(sam,"*\t"); // QUAL read 
	     fprintf(sam,"NM:i:%u\t",ali_courant->edit); //NM 
	     fprintf(sam,"NH:i:%u\t",nb_align); //NH
	     fprintf(sam,"HI:i:%i\t",cpt+1); //HI
	     fprintf(sam,"IH:i:%i\t",nb_sam); //IH
	     fprintf(sam,"X0:i:%u\t",nb_best); //X0
	     fprintf(sam,"X1:i:%u",nb_align - nb_best); //X1
	     fprintf(sam,"\n"); // fin align


	     prev_score = p->score_al;
	     p++;cpt++;
	   }


	   // clear memoire
	   listali.clear();
	   for (j=0; j<nb_res; j++)
	     {
	       free_ali( tab_ali+j );//
	     }
	    free(tab_ali);
	    
	    cell_used = cur;
	    cur = cur->suiv;
	    free_hashcell(cell_used);
	    ///fin traitement du read


	    }


	}
    }

  /////////// fin nouvel ordre
  */
  

 for (i=0; i<NHACH_READS; i++) //du coup ordre aleatoire des resu , provient hashcode
      {
	cur = table_hachage_reads[i];
	while(cur!=NULL )
	  {
	    ///traitement pour un read ......
	    //	printf("----------------------------------------------------------------------------\n");
	    //printf("----------------------------------------------------------------------------\n");


	    //printf("%s     :   nb %u  offset %u \n",cur->nom,cur->nb, cur->offset);
	    nb_res = cur->nb; offset_read = cur->offset;
	    tab_ali = (ali * ) malloc(nb_res * sizeof(ali)); // stockage pour les  nb_res alignements du read

	    nb_align= nb_res;
	    //remplissage en memoire
	    for (j=0; j<nb_res; j++)
	      {
		pos_res = index_reads[offset_read+j];
		fsetpos(align_file,&pos_res);// placement au bon endroit dans le fichier 
		// puis lecture de la ligne correspondante
		lecture_align(align_file, nom_read,nom_contig, &gaps,&mis,&edit_dist,&pos,&tai_read,seq_read,&forward, cigar);

		//stockage des donnees en memoire : 
		newali = tab_ali +j;
		unsigned int nscore = compute_score(gaps,mis,tai_read-edit_dist);
		
		
		//newali = (cell*) malloc(sizeof(cell));
		newali->nom =  (char *)  malloc((strlen(nom_read)+1)*sizeof(char));
		strcpy(newali->nom,nom_read);
		newali->tai_read=tai_read;
		newali->seq = (char *)  malloc((tai_read+1)*sizeof(char));
		strcpy(newali->seq,seq_read);
		//données bestal
		newali->nom_contig =  (char *)  malloc((strlen(nom_contig)+1)*sizeof(char));
		newali->cigar =  (char *)  malloc((strlen(cigar)+1)*sizeof(char));
		strcpy(newali->nom_contig,nom_contig);
		strcpy(newali->cigar,cigar);
		newali->pos=pos; newali->nb_gaps=gaps; newali->nb_mis=mis;  newali->edit=edit_dist;
		newali->score_al=nscore;
		newali->forward = forward; 
		
		listali.push_back(resu(nscore,j));
      		//printf("--------- Resultat  %i --------------\n",j+1);
		//		printf("%s pos  :%u\n",nom_read,pos);
		//printf("-------------------------------------\n");
	      }

	    // tri des n align du read  :
	    listali.sort();

	    //premier parcours pour nb_best
	   p = listali.begin();
	   nb_best=0;
	   unsigned int prev_score;
	   prev_score= p->score_al;

	   while(p != listali.end() && (prev_score == p->score_al) ) {
	     
	     nb_best++;
	     prev_score = p->score_al;
	     p++;
	   }

	   // nb_best est calculé 
	   nb_sam = min2(nb_res,(unsigned int)nb_res_per_read);
	   if (nb_res_per_read==0) nb_sam= nb_res;

	   p = listali.begin();
	   int cpt=0;
	   while(p != listali.end() && (cpt<nb_sam) ) {
	     ali_courant = &(tab_ali[p->index]); 
	    

	     fprintf(sam,"%s\t",ali_courant->nom);
	     fprintf(sam,"%u\t",compute_flag(ali_courant)); 
	     fprintf(sam,"%s\t",ali_courant->nom_contig);
	     fprintf(sam,"%u\t",ali_courant->pos);
	     if(cpt!=0) // tous les align subopt : mapqual = 0
	       fprintf(sam,"%u\t",0); //mapping quality
	     else
	       {
		 // si unique : 50 
		 if (nb_align==1) 
		   fprintf(sam,"%u\t",50); //mapping quality 50 si best et seul align
		 else
		   {
		     p++;
		     ali_suivant = &(tab_ali[p->index]); 
		     p--;
		     fprintf(sam,"%u\t",compute_map_qual_multi(nb_best,ali_courant->score_al,ali_suivant->score_al)); //mapping quality

		   }
	       }
	     fprintf(sam,"%s\t",ali_courant->cigar); // CIGAR
	     fprintf(sam,"*\t"); // MRNM
	     fprintf(sam,"0\t"); // MPOS
	     fprintf(sam,"0\t"); // ISIZE
	     fprintf(sam,"%s\t",ali_courant->seq); // SEQ
	     fprintf(sam,"*\t"); // QUAL read 
	     fprintf(sam,"NM:i:%u\t",ali_courant->edit); //NM 
	     fprintf(sam,"NH:i:%u\t",nb_align); //NH
	     fprintf(sam,"HI:i:%i\t",cpt+1); //HI
	     fprintf(sam,"IH:i:%i\t",nb_sam); //IH
	     fprintf(sam,"X0:i:%u\t",nb_best); //X0
	     fprintf(sam,"X1:i:%u",nb_align - nb_best); //X1
	     fprintf(sam,"\n"); // fin align


	     prev_score = p->score_al;
	     p++;cpt++;
	   }


	   // clear memoire
	   listali.clear();
	   for (j=0; j<nb_res; j++)
	     {
	       free_ali( tab_ali+j );//
	     }
	    free(tab_ali);
	    
	    cell_used = cur;
	    cur = cur->suiv;
	    free_hashcell(cell_used);

	    
	  } //fin while cur

      }




  ///////////////


  /*
  cell ** tab_hach = (cell **) malloc( NHACHG2S * sizeof(cell * ));
  for (i=0; i<NHACHG2S; i++) tab_hach[i]=NULL;


  //  printf("sizeof cell %i octets \n",sizeof(cell));

  int clef; 
  int cpt_debug=0;
  cell * cur;
  unsigned int    pos, gaps,mis,edit_dist,tai_read;
  int forward;
  //   printf("Début construction table hachage ...\n");

    while( lecture_align(align_file, nom_read,nom_contig, &gaps,&mis,&edit_dist,&pos,&tai_read,seq_read,&forward, cigar) )
    {
      //printf("%i  ",cpt_debug);
      // printf(" %s %s  gaps %u mis %u edit %u  pos %u  tairead %u \n",nom_read,nom_contig,gaps,mis,edit_dist,pos,tai_read);

	      unsigned int nscore = compute_score(gaps,mis,tai_read-edit_dist);


      cell * newcell;
      clef = hashCode (nom_read,NHACHG2S);
      if (tab_hach[clef]==NULL) //clef contient aucune cellule
	{

	  newcell = (cell*) malloc(sizeof(cell));
	  newcell->nom =  (char *)  malloc((strlen(nom_read)+1)*sizeof(char));
	  strcpy(newcell->nom,nom_read);
	  newcell->tai_read=tai_read;
	  newcell->seq = (char *)  malloc((tai_read+1)*sizeof(char));
	  strcpy(newcell->seq,seq_read);
	  //données bestal
	  newcell->nom_contig =  (char *)  malloc((strlen(nom_contig)+1)*sizeof(char));
	  newcell->cigar =  (char *)  malloc((strlen(cigar)+1)*sizeof(char));

	  strcpy(newcell->nom_contig,nom_contig);
	  strcpy(newcell->cigar,cigar);
	  newcell->pos=pos; newcell->nb_gaps=gaps; newcell->nb_mis=mis;  newcell->edit=edit_dist;
	  newcell->score_al=nscore;
	  newcell->forward = forward; 


	  //
	  newcell->nb_align=1; newcell->nb_best=1;
	  newcell->edit_second = NA;  newcell->score_al_sec = NA; // not applicable

	  newcell->suiv=NULL;
	  tab_hach[clef]=newcell;
	}
      else
	{
	  //parcours liste
	  cur = tab_hach[clef];
	  while(cur!=NULL && strcmp(cur->nom,nom_read))
	    cur = cur->suiv;
	  if (cur==NULL) //read non trouvé
	    {

	      //creation nouvelle cellule
	      newcell = (cell*) malloc(sizeof(cell));
	      newcell->nom =  (char *)  malloc((strlen(nom_read)+1)*sizeof(char));

	      strcpy(newcell->nom,nom_read);
	      newcell->tai_read=tai_read;
	      newcell->seq = (char *)  malloc((tai_read+1)*sizeof(char));
	      strcpy(newcell->seq,seq_read);
	      //données bestal
	      newcell->nom_contig =  (char *)  malloc((strlen(nom_contig)+1)*sizeof(char));
	      newcell->cigar =  (char *)  malloc((strlen(cigar)+1)*sizeof(char));

	      strcpy(newcell->nom_contig,nom_contig);
	      strcpy(newcell->cigar,cigar);
	      newcell->pos=pos; newcell->nb_gaps=gaps; newcell->nb_mis=mis;  newcell->edit=edit_dist;
	      newcell->score_al = nscore;
	      newcell->forward = forward; 
	      //
	      newcell->nb_align=1; newcell->nb_best=1;
	      newcell->edit_second = NA; // NA
	      
	      newcell->suiv = tab_hach[clef]; // push front 
	      tab_hach[clef] = newcell;
  
	    }
	  else //read trouvé
	    {
	      (cur->nb_align)++;
	      //  if (cur->edit > edit_dist) // le nouveau est meilleur
	      if (cur->score_al < nscore) // le nouveau est meilleur
		{
		  cur->nb_best=1;

		  cur->nom_contig =  (char *)  realloc(cur->nom_contig,(strlen(nom_contig)+1)*sizeof(char));
		  cur->cigar =  (char *)  realloc(cur->cigar,(strlen(cigar)+1)*sizeof(char));
		  cur->seq =  (char *)  realloc(cur->seq,(tai_read+1)*sizeof(char));

		  strcpy(cur->nom_contig,nom_contig);
		  strcpy(cur->cigar,cigar);
		  strcpy(cur->seq,seq_read); // recopie read aussi car il peut changer si passage rev /forward ...

		  cur->edit_second = cur->edit; cur->score_al_sec = cur->score_al;   //cur->nb_gaps_sec=cur->nb_gaps; cur->nb_mis_sec=cur->nb_mis;  //le précedent devient second meilleur
		  cur->pos=pos; cur->nb_gaps=gaps; cur->nb_mis=mis;  cur->edit=edit_dist;
		  cur->score_al = nscore;

		  cur->forward = forward;


		}
	      //  else if (cur->edit == edit_dist) 
	      else if (cur->score_al == nscore) 
		{
		  //second align de mem score, map qaulity à 0, on garde le premier (arbitraire)
		  (cur->nb_best)++;
		}
	      //     else if (cur->edit < edit_dist) 
	      else if (cur->score_al > nscore) 

		// moins  bon que precedent meilleur , on teste si devient meilleur second
		{
		  //  if (edit_dist<cur->edit_second ) { cur->edit_second =edit_dist; cur->score_al_sec = nscore; }
		  if (nscore > cur->score_al_sec ) { cur->edit_second =edit_dist; cur->score_al_sec = nscore; }
		}

	    }

	}
      cpt_debug++;
    }


 


    ///////////////////// etape 2, parcours table hachage
    
    cell * cell_used;
    for (i=0; i<NHACHG2S; i++)
      {
	cur = tab_hach[i];
	
	while(cur!=NULL )
	  {
	    //print align	  
     

	    fprintf(sam,"%s\t",cur->nom);
	    fprintf(sam,"%u\t",compute_flag(cur)); 
	    fprintf(sam,"%s\t",cur->nom_contig);
	    fprintf(sam,"%u\t",cur->pos);
	    fprintf(sam,"%u\t",compute_map_qual(cur)); //mapping quality
	    fprintf(sam,"%s\t",cur->cigar); // CIGAR
	    fprintf(sam,"*\t"); // MRNM
	    fprintf(sam,"0\t"); // MPOS
	    fprintf(sam,"0\t"); // ISIZE
	    fprintf(sam,"%s\t",cur->seq); // SEQ
	    fprintf(sam,"*\t"); // QUAL read 
	    fprintf(sam,"NM:i:%u\t",cur->edit); //NM 
	    fprintf(sam,"NH:i:%u\t",cur->nb_align); //NH
	    fprintf(sam,"HI:i:1\t"); //HI
	    fprintf(sam,"X0:i:%u\t",cur->nb_best); //X0
	    fprintf(sam,"X1:i:%u",cur->nb_align - cur->nb_best); //X1
	    fprintf(sam,"\n"); // fin align
	  
	    cell_used = cur;
	    cur = cur->suiv;
	    free_cell(cell_used);
	  }
  
      }

  */

  fclose(align_file);
  fclose(sam);
  free(index_reads);
  free(table_hachage_reads);
  free(seq_read);
  exit(0);
}









int lecture_read_align(
		FILE* in,
		char * nom_read
		)
{

  // printf("lecture align \n");
  //necessite in deja ouvert
  char c;
  int res;
    {
      if ((c = fgetc(in)) != EOF) 
	{
	  fseek(in,-1,SEEK_CUR);

	  res=fscanf(in,"%s %*[^\n]\n",nom_read);
	  if(res==EOF) return 0;
	  // printf("res %i  %i \n",res,EOF);
	  return 1;
	}
      else
	{
	  return 0;
	}
    }
    
}











int lecture_align(
		FILE* in,
		char * nom_read,
		char * nom_contig,
		unsigned int * gaps,
		unsigned int * mis,
		unsigned int * edit_dist, // edit_dist
		unsigned int * pos,
		unsigned int * tai_read,
		char * seq_read,
		int * forward,
		char * cigar

		)
{

  // printf("lecture align \n");
  //necessite in deja ouvert
  char c;
  int res;
    {
      if ((c = fgetc(in)) != EOF) 
	{
	  fseek(in,-1,SEEK_CUR);

	  res=fscanf(in,"%s %s %i %u %s %u %u %s %u",nom_read,nom_contig,forward,pos,cigar,mis,gaps,seq_read,tai_read);
	  if(res==EOF) return 0;
	  // printf("res %i  %i \n",res,EOF);
	  *forward = !(*forward) ; // reverse dans fichier
	  *edit_dist= ((*mis+*gaps)); // = edit distance

	  return 1;
	}
      else
	{
	  return 0;
	}
    }
    
}


int hashCode(char * nom_seq,int max) {
  int h = 0 ;
  unsigned int k;
  for ( k = 0 ; k < strlen(nom_seq); k++)
    h = 31 * h + nom_seq[k] ;
  return (abs(h) % max);
}

