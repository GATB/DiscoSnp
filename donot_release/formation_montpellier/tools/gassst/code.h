#ifndef CODE_H
#define CODE_H

/* This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software.  You can  use, 
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info". 
 */

/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file code.h
 * \brief Module Code, responsable du codage des séquences de nucléotides
 * \author Dominique Lavenier
 * \author Damien Fleury
 * \author Guillaume Rizk
 * \version 6.1
 * \date 21/01/2009
 */

#include <cctype>
#include "Seed.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "misc.h"
#include <list>
#include <stdio.h> 


/*A*/ int min3(int a, int b, int c);



/**
 * Calcul du nombre total de graines
 * \return le nombre de graines
 */
/*A*/ long long  nbSeeds();


/**
 * Méthode permettan d'obtenir le code d'un caractère
 * (A,a --> 0 ; C,c --> 1 ; G,g --> 3 ; T,t --> 2)
 * \param c, le caractère dont on veut connaitre le code
 * \return le code du caractère en paramètre
 */
/*A*/ int codeNT(char c);




//obtenir char du code
/*A*/ char NTcode(int c);


/**
 * Méthode de calcul du code d'une séquence
 * \param seq, un pointeur de caractères vers une séquence
 * \param tai, la taille de la sequence
 * \return le code de la séquence en paramètre
 */
/*A*/ int seq2code(char *seq,int tai);

//pour premier appel quand on est sur premiere pos de la seq : remplissage de voisinage droit
//non changé par INDEX_STRIDE
/*A*/ int seq2code_special(char *seq,int tailue,int tai_zone_ecrite);
/**
 * Méthode de calcul du code d'une séquence
 * \param seq, un pointeur de caractères vers la fin de la seq
 * seq a lire de la fin vers le debut (renversement)
 * \param tai, la taille de la sequence
 * \return le code de la séquence en paramètre
 */
/*A*/ int seq2code_rev(char *seq,int tai);




//obtenir sequence correspondante a code
// avec codage suivant :
//A->0  C->1 G->3 T->2
//code      A  A  C  T       :   00  00  01  10 
//  bits poids faible :  derniere lettre     

/*A*/  void code2seq(unsigned int code, char * seq, int taille);








/**
 * Méthode qui teste l'égalité entre 2 caractères a et b
 * \param a un caractère
 * \param b un caractère
 * \return 1 si a et b sont identiques, 0 sinon
 */
/*A*/ int identNT(char a, char b);
/*
 inline  int identNT(char a, char b)
{
	switch(a)
	{
		case 'a':
			return (b=='a' || b=='A');
		case 'c':
			return (b=='c' || b=='C');
		case 't':
			return (b=='t' || b=='T');
		case 'g':
			return (b=='g' || b=='G');
		case 'A':
			return (b=='a' || b=='A');
		case 'C':
			return (b=='c' || b=='C');
		case 'T':
			return (b=='t' || b=='T');
		case 'G':
			return (b=='g' || b=='G');
		default:
			return 0;
	}
}
*/
/**
 * Méthode permettant d'obtenir la base complémentaire
 * \param c un caractère correspondant à une base
 * \return le caractère correspondant à la base complémentaire
 */
/*A*/ char complNT(char c);

//version avec gaps
/*A*/ char complNTG(char c);


/**
 * Méthode de calcul du code d'une séquence
 * \param seq, un pointeur de caractères vers une séquence
 * \return le code de la séquence en paramètre
 */
/*A*/ long long  codeSeed(char *seq);


/**
 * Méthode permettant de calculer le code d'une graine à partir de celui de la graine précédente.
 * \param seq, le pointeur de la séquence
 * \param val_seed, le code de la graine précédente (donc correspondant à un décalage de une lettre vers la gauche)
 */
/*A*/ long long  codeSeedRight(char *seq, long long  val_seed);

// avec decalage de stride par rapport a valeur lue precedemment 

long long  codeSeedRight_stride(char *seq,  long long  val_seed,int stride);







//pointeur seq vers debut sequence right
// val_seq valeur de seq precedente (une lettre vers la gauche)
//tai = tai zone lecture
/*A*/ unsigned int seq2codeRight(char *seq, int tai, unsigned int val_seq);


unsigned int seq2codeRight_stride(char *seq, int tai, unsigned int val_seq, int stride);



/*
unsigned int seq2codeRight_stride(char *seq, int tai, unsigned int val_seq, int stride,int tair)
{
  int i;
  unsigned int temp = val_seq;
  //il n en reste que tair a droite
  //or on en a deja lu 16 - stride (mangés a gauche)
  // il n en reste que tair -16 + stride a lire au plus 
  int alire= min(stride, tair - 4*sizeof(unsigned int) +stride);
  for( i=0; i<stride; i++)
    {
      temp = (temp << 2 ) ;
      if (i<alire) temp +=  codeNT( seq[tai-i] ) ;

    }
  return temp;
}
*/

//pointeur seq vers fin seq left
// val_seq valeur de seq precedente (une lettre vers la gauche)
// bool ==1  si seq precedente non disponible pour le calcul
/*A*/ unsigned int seq2codeLeft(char *seq, unsigned int val_seq);

// pourqoui -tai + 1 ?


unsigned int seq2codeLeft_stride(char *seq, unsigned int val_seq,int stride);





unsigned int seq2codeRight_stride_tail(char *seq, int tai, unsigned int val_seq, int stride,int tair);




/**
 * Méthode qui rend la majuscule correspondant à un caractère
 * \param c, un caractère
 * \return le caractère majuscule correspondant
 */
/*A*/ char majuscule(char c);


/**
 * Méthode qui rend la minuscule correspondant à un caractère
 * \param c, un caractère
 * \return le caractère minuscule correspondant
 */
/*A*/ char minuscule(char c);


/** 
 * Méthode qui calcule la valeur du ScoreBit d'un alignement en fonction de son score
 * \param s le score de l'alignement
 * \return la valeur ScoreBit de l'alignement
 */
double ScoreBit(double s);

/**
 * Méthode calculant la E-value des banques pour un certain score
 * \param s un entier représentant un score
 * \param n1 la taille de la première banque
 * \param i1 le nombre de sequences dans la première banque
 * \param n2 la taille de la ceconde banque
 * \param i2 le nombre de séquences dans la seconde banque
 */
double ComputeEvalue(int s, long long int n1, int i1, int n2, int i2);



//renvoit score alignement des 2 seq de longueur ncar
int mini_align(char *s1, char* s2,int ncar);



/**
 * remplissage de la table de precalcul des alignements avec gap
 * le tableau table_prec doit etre deja alloué
 * \param table_prec, le tableau bidimensionnel
 * \param ncar, nombre de bases precalculees
 * \return le tableau rempli en parametre
 */
/*void compute_table(char ** table_prec, int ncar)
{
  long int dimt = (long int)pow(4,ncar);
  //  char* s1 = malloc(ncar);
  // char* s2 = malloc(ncar);
  char * s1 = new char [ncar];
  char * s2 = new char [ncar];

  int cs1,cs2;
  for (cs1=0; cs1< dimt; cs1++)
    {
      code2seq(cs1,s1,ncar);
      for (cs2=0; cs2< dimt; cs2++)
	{
	  code2seq(cs2,s2,ncar);
	  table_prec[cs2][cs1]= (char) mini_align(s1,s2,ncar);
	  // printf("%i",(int) table_prec[cs1][cs2]);

	}    
      //  printf("\n");
    }
  
  delete [] s1;
  delete [] s2;
}
*/

void compute_table(char ** table_prec, int ncar);




void cpt_nt(char *s1, char *a, char *c, char *t, char*g,int ncar);


void compute_table_nt(char ** table_nt, int ncar);


inline int cpt_err(char * s1, char * s2, int ncar);


void compute_table_gless(char ** table_prec, int ncar);



int find_seed_len(int len, int nerr, int neval,int nquantile);











#endif
