#ifndef resu_H
#define resu_H

/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file resu.h
 * \author Guillaume Rizk
 * \date 01/02/2011
 */


/**
 * \class resu, Une graine correspond à une petite séquence d'ADN
 * \brief Cette class définit une graine, c'est à dire une courte séquence qui permettra la recherche des alignements
 */
class resu{
 public:
  
  
  unsigned int score_al;
  
  unsigned int index;




	
	/**
	 * Constructeur par défaut
	 */
	resu();
	
	/**
	 * Constructeur de resu
	 * \param num le numéro de la séquence
	 * \param offset, l'index dans séquence
	 */
	resu(  unsigned int score_in,
	       unsigned int index_in);

	
	resu(const resu &);
	~resu(){};
	resu &operator=(const resu &rhs);
	int operator==(const resu &rhs) const;
	int operator<(const resu &rhs) const;


};

#endif
