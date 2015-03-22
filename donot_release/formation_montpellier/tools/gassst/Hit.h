#ifndef HIT_H
#define HIT_H

/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file Hit.h
 * \brief Classe Hit, définissant un alignement entre 2 séquences
 * \author Dominique Lavenier
 * \author Damien Fleury
 * \version 5.2
 * \date 28/08/2008
 */

class Bank;

/**
 * \class Hit, Une position d'une graine dans une séquence
 * \brief Cette classe correspond à une position d'une graine dans une séquence d'une banque
 */
class Hit{
	 public:

	/**
	 * Index de la position de la graine
	 */
	int offhit;
 	/**
	 * Index de la séquence
	 */
	int offseq;
	/**
	 * Taille de la séquence
	 */
	int sizeseq;
	/** 
	 * Numéro de la séquence
	 */
	int numseq;
	
	//	unsigned char seqleft;
	//	unsigned char seqright;
	/**
	 * Constructeur par défaut
	 */
	Hit(); 
	
	/**
	 * Constructeur de Hit
	 * \param BK un pointeur vers la banque où est créé l'alignement
	 * \param num_sequence le numéro de la séquence où est situé l'alignement
	 * \param offset_sequence la position de l'alignement dans la séquence
	 */
	Hit(Bank *BK, int num_sequence, int offset_sequence);
	
	/**
	 * Constructeur de Hit par recopie
	 * \param h un objet Hit
	 */
	Hit(const Hit& h);
	
	/**
	 * Destructeur de Hit
	 */
	~Hit();
	
	/**
	 * Opérateur d'affectation de Hit
	 * \param h un objet Hit
	 * \return l'objet Hit affecté
	 */
	Hit& operator=(const Hit& h);

	//operateur ordre utilisé pour tri
	int operator<(const Hit &h) const;

};

#endif
