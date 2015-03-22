/* This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software.  You can  use, 
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info". 
 */

/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file display.cpp
 * \brief Module Display, responsable de l'enregistrement du résultat dans un fichier de sortie
 * \author Dominique Lavenier
 * \author Damien Fleury
 * \version 5.2
 * \date 28/08/2008
 */
 
#include "display.h"

#include "constants.h"
#include "code.h"
#include "misc.h"
#include "Alignment.h"

int OUTPUT_FORMAT;

/**
 * Méthode de vérification du format de sortie
 * \param k le type de format demandé
 */
void checkOutputFormat(int k)
{
  if (k==STD_OUTPUT_FORMAT) return;
  if (k==M8_OUTPUT_FORMAT) return;
  if (k==SAM_READY_FORMAT) return;

  fprintf (stderr,"wrong output format, should be %i, %i or %i",STD_OUTPUT_FORMAT,M8_OUTPUT_FORMAT,SAM_READY_FORMAT);
  exit(0);
  // ExitError("wrong output format");
}

/**
 * Fonction permettant d'enregistrer un alignement dans le fichier de sortie
 * \param ff le pointeur vers l'objet FILE du fichier de sortie
 * \param s1 le pointeur vers la première séquence
 * \param s2 le pointeur vers la seconde séquence
 * \param c1 le pointeur vers le début du commentaire de la séquence de la première banque
 * \param c2 le pointeur vers le début du commentaire de la séquence de la seconde banque
 * \param start1 la position de départ de l'alignement pour la première séquence
 * \param start2 la position de départ de l'alignement pour la seconde séquence
 * \param len la longueur de l'alignement (longueur de la plus courte séquence)
 * \param l_seq1 la longueur de la première séquence, la plus longue
 * \param eval la valeur de la E-value
 * \param rev_comp un booléen indiquant si on a un alignement sur la séquence normale (false) ou inversée et complémentée (true)
 */
void display_ungap(FILE *ff, char* s1, char* s2, char* c1, char* c2, int start1, int start2, int len, int l_seq1, double eval, bool rev_comp)
{
  int i,j,k;
  int tailleColonneM8 = 30;
  double f;
  
  /// Calcul des index réels de début et de fin de l'alignement
  int deb1,deb2,stop1,stop2;
  deb2 = start2 + 1;
  stop2 = start2 + len;
  
  /// On calcule les bornes de l'alignement selon le type d'alignement
  switch(rev_comp)
  {
  	case NORMAL_ALIGN:
		/// Cas d'un alignement dans la séquence normale
		deb1 = start1 + 1;
		stop1 = start1 + len;
		break;
	case REV_COMP_ALIGN:
		///Cas d'un alignement dans la séquence inversée et complémentée
		deb1 = l_seq1 - start1;
		stop1 = l_seq1 - start1 - len + 1;
		break;
	default:
		deb1 = start1 + 1;
		stop1 = start1 + len ;
  }
  
  /// Cas du format de sortie standard
  if (OUTPUT_FORMAT == STD_OUTPUT_FORMAT)
  {
  	  /// Adresse de départ de la première séquence
      fprintf (ff,"BANK  %10d  ",deb1);
      /// Impression de la première séquence
      for (i=0; i<len; ++i)
      {
      	fprintf (ff,"%c",s1[start1+i]);
      }
      /// Adresse de fin de la première séquence
      fprintf (ff,"  %10d\t\t",stop1);
      i=1;
      /// Affichage du nom de la première séquence
      while ((c1[i]!='>')&&(c1[i]!='\0'))
      {
      	fprintf(ff,"%c",c1[i]);
      	++i;
      }
      fprintf (ff,"\n");
      fprintf (ff,"                  ");
      k=0;
      /// Impression des barres pour marquer les appariements entre les deux séquences
      for (i=0,j=0; i<len; ++i, ++j)
      {
		if (identNT(s1[start1+i],s2[start2+j])==1) fprintf (ff,"|"); 
		else { fprintf (ff," "); ++k; }
      }
      /// Affichage du nombre de mésappariments et de la e-value
      fprintf (ff,"   \t\t\t# mismatche(s): %d    e-value: %2.0e\n",k,eval);
      /// Adresse de départ de la seconde séquence
      fprintf (ff,"QUERY %10d  ",deb2);
      /// Impression de la seconde séquence
      for (j=0; j<len; ++j) fprintf (ff,"%c",s2[start2+j]);
      /// Adresse de fin de la seconde séquence
      fprintf (ff,"  %10d\t\t",stop2);
      i=1;
      /// Affichage du nom de la seconde séquence
      while ((c2[i]!='>')&&(c2[i]!='\0'))
      {
      	fprintf(ff,"%c",c2[i]);
      	++i;
      }
      fprintf (ff,"\n\n");
      return;
  }
  /// Cas du format de sortie m8
  if (OUTPUT_FORMAT == M8_OUTPUT_FORMAT)
  {
  	i = 1;
  	/// Affichage du nom de la première séquence
    while ((c2[i]!='>') && (c2[i]!='\0') /*&& (i<tailleColonneM8)*/)
	{
		fprintf (ff,"%c",c2[i]);
	  	++i;
	}
	while(i<tailleColonneM8)
	{
		fprintf (ff," ");
	  	++i;
	}
	fprintf (ff,"\t");
    i = 1;
    /// Affichage du nom de la seconde séquence
    while ((c1[i]!='>') && (c1[i]!='\0') /*&& (i<tailleColonneM8)*/)
	{
		fprintf (ff,"%c",c1[i]);
		++i;
	}
	while(i<tailleColonneM8)
	{
		fprintf (ff," ");
	  	++i;
	}
	fprintf (ff,"\t");

    k=0;
    /// Calcul du pourcentage d'identité
    for(i=0,j=0; i<len; ++i,++j)
    {
    	if(identNT(s1[start1+i],s2[start2+j])==1) ++k;
    }
    f = k*100;
    f = f/len;
    
    fprintf (ff,"%5.2f\t",f);			/// Pourcentage d'identité
    fprintf (ff,"%6d\t",len);			/// Taille de l'alignement
    fprintf (ff,"%6d\t",len-k);			/// Nombre de mismatch
    fprintf (ff,"%6d\t",0);				/// Nombre de gaps (pour être conforme au format Blast, ici 0)
    fprintf (ff,"%6d\t",deb2);			/// Debut alignement séquence 2
    fprintf (ff,"%6d\t",stop2);			/// Fin alignement séquence 2
    fprintf (ff,"%6d\t",deb1);			/// Debut alignement séquence 1
    fprintf (ff,"%6d\t",stop1);			/// Fin alignement séquence 1
    fprintf (ff,"%6.0e\t",eval);		/// E-value
    f = ScoreBit(k*MATCH);
    fprintf (ff,"%6.1f",f);				/// Bit score
    fprintf (ff,"\n");
    return;
  }
}


/**
 * Fonction permettant d'enregistrer un alignement avec gap(s) dans le fichier de sortie
 * \param ff le pointeur vers l'objet FILE du fichier de sortie
 * \param al l'objet Alignment contenant les informations sur l'alignement
 * \param c1 le pointeur vers le début du commentaire de la séquence de la première banque
 * \param c2 le pointeur vers le début du commentaire de la séquence de la seconde banque
 * \param eval la valeur de la E-value
 */
void display_withgap(FILE *ff, Alignment* al, char* c1, char* c2, double eval)
{
  int i,j,n;
  double k;
  int tailleColonneM8 = 30;

  double f;
  int len = al->getLength();
  /// Cas du format de sortie standard
  NUM_ALIGN_COURANT ++;
 al->apply_rev_comp(OUTPUT_FORMAT);
  if (OUTPUT_FORMAT == STD_OUTPUT_FORMAT)
  {
      
      /// Adresse de départ de la première séquence
      fprintf (ff,"BANK  %10d  ",al->getStart1());
      /// Impression de la première séquence
      fprintf (ff,"%s",al->getSeq1());
      /// Adresse de fin de la première séquence
      fprintf (ff,"  %10d\t\t",al->getEnd1());
      i=1;
#ifdef NUMSEQSORTIE
	fprintf(ff,"CONTIG:%u",al->n1);
#else
      /// Affichage du nom de la première séquence
      while ((c1[i]!='>')&&(c1[i]!='\0'))
      {
      	fprintf(ff,"%c",c1[i]);
      	++i;
      }
#endif

      fprintf (ff,"\n");
      fprintf (ff,"                  ");
      
      /// Impression des barres pour marquer les appariements entre les deux séquences
      for (i=0,j=0; i<len; ++i, ++j)
      {
		if (identNT(al->getChar1(i),al->getChar2(j)) == 1) fprintf (ff,"|"); 
		else
		  fprintf (ff," ");
      }
      
      /// Affichage du nombre de mésappariments et de la e-value
      n = al->getMis() ;
      //+ al->getGaps();
      fprintf (ff,"   \t\t\t# mismatche(s): %d    gap(s): %d    e-value: %2.0e\n",n,al->getGaps(),eval);
      /// Adresse de départ de la seconde séquence
      fprintf (ff,"QUERY %10d  ",al->getStart2());
      /// Impression de la seconde séquence
      fprintf (ff,"%s",al->getSeq2());
      /// Adresse de fin de la seconde séquence
      fprintf (ff,"  %10d\t\t",al->getEnd2());
      i=1;

      /// Affichage du nom de la seconde séquence
      while ((c2[i]!='>')&&(c2[i]!='\0'))
      {
      	fprintf(ff,"%c",c2[i]);
      	++i;
      }
      fprintf (ff,"\n\n");
      return;
  }
  
  /// Cas du format de sortie m8
  if (OUTPUT_FORMAT == M8_OUTPUT_FORMAT)
  {
  	i = 1;
  	/// Affichage du nom de la première séquence
    while ((c2[i]!='>') && (c2[i]!='\0') /*&& (c2[i]!=' ')*/ /* && (i<tailleColonneM8)*/)
	{
        if(c2[i]!=' ')
        {
            fprintf (ff,"%c",c2[i]);
        }
        else
            fprintf (ff,"_");	  	++i;
	}
	while(i<tailleColonneM8)
	{
		fprintf (ff," ");
	  	++i;
	}
	fprintf (ff,"\t");
    i = 1;
    /// Affichage du nom de la seconde séquence
    while ((c1[i]!='>') && (c1[i]!='\0') /*&& (c1[i]!=' ')*/  /*&& (i<tailleColonneM8)*/)
	{
        if(c1[i]!=' ')
        {
            fprintf (ff,"%c",c1[i]);
        }
        else
            fprintf (ff,"_");		++i;
	}
	while(i<tailleColonneM8)
	{
		fprintf (ff," ");
	  	++i;
	}
	fprintf (ff,"\t");
	
	/// Calcul du pourcentage d'identité
    k = (double) ((al->getLength() - al->getMis() - al->getGaps() )*100) / (double) al->getLength();

    fprintf (ff,"%5.2f\t",k);					/// Pourcentage d'identité
    fprintf (ff,"%6d\t",al->getLength());		/// Taille de l'alignement
    fprintf (ff,"%6d\t",al->getMis());			/// Nombre de mismatch
    fprintf (ff,"%6d\t",al->getGaps());			/// Nombre de gaps
    fprintf (ff,"%6d\t",al->getStart2());		/// Debut alignement séquence 2
    fprintf (ff,"%6d\t",al->getEnd2());			/// Fin alignement séquence 2
    fprintf (ff,"%6d\t",al->getStart1());		/// Debut alignement séquence 1
    fprintf (ff,"%6d\t",al->getEnd1());			/// Fin alignement séquence 1
    fprintf (ff,"%6.0e\t",eval);				/// E-value
    f = ScoreBit(k*MATCH);
    fprintf (ff,"%6.1f",f);						/// Bit score
    fprintf (ff,"\t%s",al->getCigar());						/// cigar
    fprintf (ff,"\n");
    return;
  }


 if (OUTPUT_FORMAT == SAM_READY_FORMAT)
  {
  	i = 1;
  	/// Affichage du nom de la première séquence, read
    while ((c2[i]!='>') && (c2[i]!='\0') /*&& (c2[i]!=' ')*/)
	{
        if(c2[i]!=' ')
        {
            fprintf (ff,"%c",c2[i]);
        }
        else
            fprintf (ff,"_");
        
	  	++i;
	}
	fprintf (ff,"\t");
    i = 1;
    /// Affichage du nom de la seconde séquence , reference
    //fprintf (ff,"%s\t",&c1[1]);

     while ((c1[i]!='>') && (c1[i]!='\0') && (c1[i]!=' '))
 	{
 		fprintf (ff,"%c",c1[i]);
 		++i;
 	}
	
 	fprintf (ff,"\t");
	
    fprintf (ff,"%i\t",al->rev_comp);			/// reverse
    fprintf (ff,"%6d\t",al->getStart1());		///  position left sur ref forward
    fprintf (ff,"%s\t",al->getCigar());			/// cigar
    fprintf (ff,"%6d\t",al->getMis());			/// Nombre de mismatch
    fprintf (ff,"%6d\t",al->getGaps());			/// Nombre de gaps
    fprintf (ff,"%s\t",al->getSeq2());			/// seq 2
    fprintf (ff,"%6d",al->getLength());		/// Taille de l'alignement



    fprintf (ff,"\n");
    return;
  }


}

