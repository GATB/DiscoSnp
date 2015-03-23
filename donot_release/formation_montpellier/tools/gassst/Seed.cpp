/* This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software.  You can  use, 
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info". 
 */
/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file Seed.cpp
 * \brief Classe Seed, définissant une graine
 * \author Dominique Lavenier
 * \author Damien Fleury
 * \author Guillaume Rizk
 * \version 6.1
 * \date 19/12/2008
 */
 
 #include "Seed.h"


/**
 * Constructeur par défaut
 */
Seed::Seed() : num_seq(-1), off_seq(-1), left(0), right (0)
{
}

/**
 * Constructeur de Seed
 * \param num le numéro de la séquence
 * \param offset, l'index dans séquence
 */
Seed::Seed(int num, int offset, unsigned int seqleft, unsigned int seqright)
{
	num_seq = num;
	off_seq = offset;
	left = seqleft;
	right = seqright;
}

/**
 * Constructeur par recopie
 * \param s un objet Seed
 */
Seed::Seed(const Seed& s)
{
	*this = s;
}

/**
 * Destructeur
 */
Seed::~Seed()
{
}

/**
 * Opérateur d'affectation
 * \param s un objet Seed
 * \return l'objet Seed affecté
 */
Seed& Seed::operator=(const Seed& s)
{
	if(this!=&s)
	{
		num_seq = s.num_seq;
		off_seq = s.off_seq;
		left = s.left;
		right = s.right;
	}
	return *this;
}
