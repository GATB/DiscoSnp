/* This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software.  You can  use, 
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info". 
 */

/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file Hit.cpp
 * \brief Classe Hit, définissant un alignement entre 2 séquences
 * \author Dominique Lavenier
 * \author Damien Fleury
 * \version 5.2
 * \date 28/08/2008
 */
 
#include "Hit.h"

#include "Bank.h"


/**
 * Constructeur par défaut
 */
Hit::Hit() : offhit(0), offseq(0), sizeseq(0), numseq(0)
{
}


/**
 * Constructeur de Hit
 * \param BK un pointeur vers la banque où est créé l'alignement 
 * \param num_sequence le numéro de la séquence où est situé l'alignement
 * \param offset_sequence la position de l'alignement dans la séquence
 */
Hit::Hit(Bank *BK, int num_sequence, int offset_sequence) : offhit(BK->seq[num_sequence] + offset_sequence), offseq(BK->seq[num_sequence]), sizeseq(BK->size[num_sequence]),numseq(num_sequence)
{
}


/**
 * Constructeur de Hit par recopie
 * \param h un objet Hit
 */
Hit::Hit(const Hit& h) : offhit(h.offhit),	offseq(h.offseq), sizeseq(h.sizeseq), numseq(h.numseq)
{
}


/**
 * Opérateur d'affectation de Hit
 * \param h un objet Hit
 * \return l'objet Hit affecté
 */
Hit& Hit::operator=(const Hit& h)
{
	if(this!=&h)
	{
		offhit = h.offhit;
		offseq = h.offseq;
		sizeseq = h.sizeseq;
		numseq = h.numseq;
	}
	return *this;
}


/**
 * Destructeur de Hit
 */
Hit::~Hit()
{
}



//operateur  <  ordonne suivant position dans le genome 
int Hit::operator<(const Hit &h) const
{

   if( this->offhit < h.offhit ) return 1;
   return 0;


}
