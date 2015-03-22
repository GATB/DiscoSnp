/* This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software.  You can  use, 
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info". 
 */

/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file filter.h
 * \brief Module Filter, définit le filtre Low Complexity qui détecte les zones non pertinentes des séquences afin de ne pas les indexer
 * \author Dominique Lavenier
 * \author Damien Fleury
 * \version 5.2
 * \date 28/08/2008
 */
 
#include "filter.h"

#include "code.h"

/**
 * Méthode permettant d'appliquer le filtre Low Complexity
 * Les zones inintéressantes des séquences sont mises en lettres minuscules
 * \param data le tableau des données de la banque de séquences
 * \param id l'index du début de la séquence
 * \param lenseq la longueur de la séquence
 * \return 1 si l'opération s'est déroulée correctement
 */
int filterLowComplexity (char* data, int lenseq)
{
	int DUSTSCORE[64];
	int i, j, k, s, m;

	for (i=0; i<64; i++) DUSTSCORE[i]=0;

	for (j=2; j<WS-1 && j<lenseq; ++j)
	{
		k = codeNT(data[j-2])*16 + codeNT(data[j-1])*4 + codeNT(data[j]);
		++DUSTSCORE[k];
	}

	s=0;
	for (i=0; i<64; ++i)
	{
		m = DUSTSCORE[i];
		s = s + m*(m-1)/2;
	}

	for (j=0; j<lenseq-2; ++j)
	{
		if (j-WS-2>=0)
		{
			k = codeNT(data[j-WS-2])*16 + codeNT(data[j-WS-1])*4 + codeNT(data[j-WS]);
			m = DUSTSCORE[k];
			s = s - m + 1;
			--DUSTSCORE[k];
		}
		if (j+WS<lenseq)
		{
			k = codeNT(data[j+WS-2])*16 + codeNT(data[j+WS-1])*4 + codeNT(data[j+WS]);
			m = DUSTSCORE[k];
			s = s + m;
			++DUSTSCORE[k];
		}
		if (s>=20) data[j]=minuscule(data[j]);
		else data[j]=majuscule(data[j]);
	}
	return 1;
}
