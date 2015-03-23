#ifndef SEED_H
#define SEED_H

/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file Seed.h
 * \brief Classe Seed, définissant une graine
 * \author Dominique Lavenier
 * \author Damien Fleury
 * \author Guillaume Rizk
 * \version 6.1
 * \date 19/12/2008
 */

/**
 * Le nombre de graines différentes
 */
extern long long  NB_DIFF_SEED;

/**
 * La taille de la graine
 */
extern int SIZE_SEED;

/**
 * \class Seed, Une graine correspond à une petite séquence d'ADN
 * \brief Cette class définit une graine, c'est à dire une courte séquence qui permettra la recherche des alignements
 */
class Seed{
	public:
	
	/**
	 * Numéro de la séquence
	 */
	int num_seq;
	/**
	 * Index de la graine dans la séquence
	 */
	int off_seq;
	/**
	 *  code séquence des 16 carac a gauche de la graine
	 */
	unsigned int left;

	/**
	 *  code séquence des 16 carac a droite de la graine
	 */
	unsigned int right;
	
	/**
	 * Constructeur par défaut
	 */
	Seed();
	
	/**
	 * Constructeur de Seed
	 * \param num le numéro de la séquence
	 * \param offset, l'index dans séquence
	 */
	Seed(int num, int offset, unsigned int seqleft, unsigned int seqright);
	/**
	 * Constructeur par recopie
	 * \param s un objet Seed
	 */
	Seed(const Seed& s);
	
	/**
	 * Destructeur
	 */
	~Seed();
	
	/**
	 * Opérateur d'affectation
	 * \param s un objet Seed
	 * \return l'objet Seed affecté
	 */
	Seed& operator=(const Seed& s);
};

#endif
