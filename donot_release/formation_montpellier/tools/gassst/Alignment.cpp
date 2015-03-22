/* This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software.  You can  use, 
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info". 
 */

/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file Alignment.cpp
 * \brief Class Alignment, définissant un alignement entre deux séquences
 * \author Damien Fleury
 * \author Guillaume Rizk
 * \date 15/07/10
 */


#include "Alignment.h"

#include <iostream>
#include "constants.h"
#include "misc.h"
#include "code.h"
#include "display.h"
#include "withgap.h"
#include <string.h>


int extended_cigar;

/**
 * Le consructeur par défaut d'Alignment
 */
Alignment::Alignment() : length(0), nb_mismatches(0), nb_gaps(0),etat_cigar(0),cpt_cigar(0), cigar_index(0)
{
	sequence1 = NULL;
	sequence2 = NULL;
	cigar = NULL;
}

/**
 * Un consructeur d'Alignment
 * \param size la taille de l'alignement
 */
Alignment::Alignment(int size) : length(0), nb_mismatches(0), nb_gaps(0),etat_cigar(0),cpt_cigar(0), cigar_index(0),e_value(0)
{
	sequence1 = new char[size+1];
	sequence2 = new char[size+1];
	cigar = new char[4*size+1];
}

/**
 * La fonction d'initialisation d'un Alignment
 */
void Alignment::init()
{
	length = 0;
	nb_mismatches = 0;
	nb_gaps = 0;
	start1 = 0;
	end1 = 0;
	start2 = 0;
	end2 = 0;
	e_value = 0;
	j1 = NULL;
	j2 = NULL;
	rev_comp = 0;
	sequence1[length] = '\0';
	sequence2[length] = '\0';
	cigar[length] = '\0';
	cpt_cigar = 0;etat_cigar=0; cigar_index=0;

}


/**
 * La fonction d'initialisation d'un Alignment
 */
void Alignment::init(char rever)
{
	length = 0;
	nb_mismatches = 0;
	nb_gaps = 0;
	start1 = 0;
	end1 = 0;
	start2 = 0;
	end2 = 0;
	e_value = 0;
	j1 = NULL;
	j2 = NULL;
	rev_comp = 0;
	sequence1[length] = '\0';
	sequence2[length] = '\0';
	rev_comp=rever;
	cigar[length] = '\0';
	cpt_cigar = 0;etat_cigar=0;cigar_index=0;
}


/**
 * Le consructeur par recopie d'Alignment
 * \param al un objet Alignment
 */
Alignment::Alignment(const Alignment& al)
{
	*this = al;
}

/**
 * L'opérateur d'affectation d'Alignment
 */
Alignment& Alignment::operator=(const Alignment& al)
{
	if(this != &al)
	{
		length = al.length;
		nb_mismatches = al.nb_mismatches;
		nb_gaps = al.nb_gaps;
		sequence1 = new char[length+20];
		sequence2 = new char[length+20];
		strncpy(sequence1, al.sequence1, length); 
		strncpy(sequence2, al.sequence2, length);
		sequence1[length] = '\0';
		sequence2[length] = '\0';

		int lencigar = strlen(al.cigar) +1; // +1 pour '\0'
		cigar = new char[lencigar];
		strncpy(cigar, al.cigar, lencigar);

		start1 = al.start1;
		start2 = al.start2;
		end1 = al.end1;
		end2 = al.end2;

		e_value = al.e_value;
		j1 = al.j1;
		j2 = al.j2;
		n1 = al.n1;
		n2 = al.n2;
		rev_comp = al.rev_comp;
		cpt_cigar = al.cpt_cigar;
		etat_cigar = al.etat_cigar;
		cigar_index = al.cigar_index;

	}
	return *this;
}

/**
 * Le destructeur par défaut d'Alignment
 */
Alignment::~Alignment()
{
	delete [] sequence1;
	delete [] sequence2;
	delete [] cigar;
}

/**
 * Méthode d'ajout d'un couple d'un caractère aux séquences
 * \param c1 le caractère à ajouter à la première séquence
 * \param c2 le caractère à ajouter à la seconde séquence
 */
void Alignment::addPair(char c1, char c2)
{

	if(OUTPUT_FORMAT==SAM_READY_FORMAT || OUTPUT_FORMAT== M8_OUTPUT_FORMAT)
	  {
	    //ajout  calcul string cigar
	    char  etat_courant;
	    if(c1==CHAR_GAP) etat_courant='I';
	    else if(c2==CHAR_GAP) etat_courant='D';
	    else
        {
            if(extended_cigar)
            {
                if(c1!=c2)
                    etat_courant=c2;
                else
                    etat_courant='M';
            }else
                etat_courant='M';
        }

	    // cas initialisation
	    if(etat_cigar == 0) {cpt_cigar=1; etat_cigar=etat_courant;}
	    else if (etat_cigar == etat_courant) cpt_cigar++;
	    else // chgt d'etat
	      {
		cigar_index += sprintf(cigar + cigar_index ,"%i%c",cpt_cigar,etat_cigar);
		cpt_cigar=1; etat_cigar = etat_courant;
	      }
	
	    if(c2!=CHAR_GAP)
	      {
		sequence2[length] = c2;
		length++;
		sequence2[length] = '\0';
	      }
	    return;
	  }

	sequence1[length] = c1;
	sequence2[length] = c2;
	++length;
	sequence1[length] = '\0';
	sequence2[length] = '\0';
	
}

/**
 * Méthode permettant d'initialiser les adresses de début et de fin de l'alignement dans les séquences
 * \param deb1 l'indice de départ du fragment de la première séquence
 * \param deb2 l'indice de départ du fragment de la seconde séquence
 * \param len la longueur de l'alignement
 * \param sizeSeq1 la longueur de la plus longue séquence
 * \rev_comp un booléen indiquant si la première séquence est inversée et complémentée ou non
 */
void Alignment::setOffsets(int deb1, int deb2, int len, int sizeSeq1, bool rev_comp)
{
	/// Calcul des index réels de début et de fin de l'alignement
	start2 = deb2 + 1;
	end2 = deb2 + len;

	start1 = deb1 + 1;
	end1 = deb1 + len;
		
}

void Alignment::adjust_rev_comp( int sizeSeq1, bool rev_comp)
{
  if(rev_comp==REV_COMP_ALIGN)
    {
      start1 = sizeSeq1 - (start1-1);
      end1 =  sizeSeq1 - (end1-1);
    }
}


//suivant format de sortie, pas la meme chose a faire
// pour seq reverse

// a la base si align sur strand reverse, le read a ete rev complemente
// donc on re rev comp ref et read pour STD_OUTPUT_FORMAT ou M8  car on veut ref reverse et read normal
// pour sam cest le contraire, on affiche toujours align sur forward ref,
// et on affiche le complem du read ( donc rien a faire pour read,ni pour cigar, bien calcule)
void Alignment::apply_rev_comp( int output_format)
{

  if(rev_comp)
    {
      if(output_format==STD_OUTPUT_FORMAT || output_format==M8_OUTPUT_FORMAT)
	{
	  char * tmp; 
	  tmp = new char [length]; 
	  for(int p=0;p<length;p++)
	    {
	      tmp[p] = sequence1[p];
	    }
	  for(int p=0;p<=length-1;p++)
	    {
	      sequence1[p] = complNTG(tmp[length-1-p]);
	    }

	  for(int p=0;p<length;p++)
	    {
	      tmp[p] = sequence2[p];
	    }
	  for(int p=0;p<=length-1;p++)
	    {
	      sequence2[p] = complNTG(tmp[length-1-p]);
	    }

	  int stemp=start1; start1=end1;  end1=stemp;
	  delete [] tmp;
	}
      else if (output_format==SAM_READY_FORMAT)
	{
	  // rien a faire 
	  ;

	}
    }
}


/**
 * Méthode permettant d'obtenir le caractère à l'indice en paramètre dans la première séquence
 * \param indice l'indice du caractère dans le tableau séquence1
 * \return le caractère à l'indice voulu dans la première séquence
 */
char Alignment::getChar1(int indice)
{
	if(indice < length && indice >= 0)
	{
		return sequence1[indice];
	}
	else return '#';
}

/**
 * Méthode permettant d'obtenir le caractère à l'indice en paramètre dans la seconde séquence
 * \param indice l'indice du caractère dans le tableau séquence2
 * \return le caractère à l'indice voulu dans la seconde séquence
 */
char Alignment::getChar2(int indice)
{
	if(indice < length && indice >= 0)
	{
		return sequence2[indice];
	}
	else return '#';
}
		
/**
 * Méthode d'affichage d'un Aligment
 */
void Alignment::affichage()
{
	printf("seq1 :  ");
	printf("%s ",getSeq1());
	printf("\nseq2 :  ");
	printf("%s ",getSeq2());
	printf("\n");
}








		
/**
 * Méthode de fin de  calcul du cigar
 */
void Alignment::finish_cigar()
{

  if(etat_cigar==0) return;
  sprintf(cigar+cigar_index,"%i%c",cpt_cigar,etat_cigar);
  
}
