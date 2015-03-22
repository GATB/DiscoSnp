#ifndef ALIGNMENT_H
#define ALIGNMENT_H

/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file Alignment.h
 * \brief Class Alignment, définissant un alignement entre deux séquences
 * \author Damien Fleury
 * \version 5.2
 * \date 28/08/2008
 */


/**
 * \class Alignment, Un alignement entre deux séquences
 * \brief Un alignement contient toutes les informations identifiant l'alignement et permettant de l'enregistrer
 * Cette classe est principalement utilisée pour le cas des alignements avec gap
 */

extern int extended_cigar;

class Alignment
{
private:
	/**
	 * La première séquence de l'alignement, c'est-à-dire celle de la base de données
	 * La chaîne est terminée par un \0
	 */
	char* sequence1;
	
	/**
	 * La seconde séquence de l'alignement, c'est-à-dire celle de la base de requêtes
	 * La chaîne est terminée par un \0
	 */
	char* sequence2;

	/**
	 * La longueur de l'alignement en nombre de caractères
	 */
	int length;
	
	/**
	 * Le nombre de mésappariements dans l'alignement
	 */
	int nb_mismatches;
	
	/**
	 * Le nombre de gaps dans l'alignement
	 */
	int nb_gaps;
	
	/**
	 * L'indice de départ de l'alignement dans la première séquence
	 */
	int start1;
	
	/**
	 * L'indice de départ de l'alignement dans la seconde séquence
	 */
	int start2;
	
	/**
	 * L'indice de fin de l'alignement dans la première séquence
	 */
	int end1;
	
	/**
	 * L'indice de fin de l'alignement dans la seconde séquence
	 */
	int end2;
	

	//2variables utiles pour le calcul du cigar
	/**
	 * operation cigar en cours : 'M'  'I' ou 'D'
	 */
	char etat_cigar;

	/**
	 *  longueur de l'operation en cours du cigar
	 */
	unsigned int cpt_cigar; 

	/**
	 *  nombre de char deja imprimes du cigar
	 */
	unsigned int cigar_index;

public:
	
	/**
	 * indique si alignement sur strand -, dans ce cas il faudra reverse la sortie
	 */
	char rev_comp;

	//score evalue de l alignement
	double e_value ;

	//index des commentaires des sequences
	char * j1;

	char * j2;

	//le numero des sequences,relatif a la partition
	unsigned int  n1,n2;

	//description de l'align en cigar, (pour format sam)
	char * cigar;


	void finish_cigar();


	/**
	 * Le consructeur par défaut d'Alignment
	 */
	Alignment();
	
	/**
	 * Un consructeur d'Alignment
	 * \param size la taille de l'alignement
	 */
	Alignment(int size);
	
	/**
	 * La fonction d'initialisation d'un Alignment
	 * \param size la taille de l'alignement
	 */
	void init();
	
	/**
	 * La fonction d'initialisation d'un Alignment
	 * \param rever si align sur strand -
	 */
	void init(char rever);
	/**
	 * Le consructeur par recopie d'Alignment
	 * \param al un objet Alignment
	 */
	Alignment(const Alignment& al);
	
	/**
	 * L'opérateur d'affectation d'Alignment
	 * \param al un objet Alignment
	 * \return l'objet Alignment affecté
	 */
	Alignment& operator=(const Alignment& al);
	
	/**
	 * Le destructeur par défaut d'Alignment
	 */
	virtual ~Alignment();
	
	/**
	 * Méthode d'obtention de la longueur de l'alignement
	 * \return la longueur de l'objet Alignment appelant
	 */
	inline int getLength();
	
	/**
	 * Méthode d'obtention du nombre de mésappariements
	 * \return le nombre de mésappariements de l'objet Alignment appelant
	 */
	inline int getMis();
	
	/**
	 * Méthode d'obtention du nombre de gaps
	 * \return le nombre de gaps de l'objet Alignment appelant
	 */
	inline int getGaps();
	
	/**
	 * Méthode d'incrémentation du nombre de mésappariements
	 */
	inline void addMis();
	
	/**
	 * Méthode d'incrémentation du nombre de gaps
	 */
	inline void addGap();

	/**
	 * Méthode permettant d'initialiser les adresses de début et de fin de l'alignement dans les séquences
	 * \param deb1 l'indice de départ du fragment de la première séquence
	 * \param deb2 l'indice de départ du fragment de la seconde séquence
	 * \param len la longueur de l'alignement
	 * \param sizeSeq1 la longueur de la plus longue séquence
	 * \rev_comp un booléen indiquant si la première séquence est inversée et complémentée ou non
	 */
	void setOffsets(int deb1, int deb2, int len, int sizeSeq1, bool rev_comp);
	
	//complement the start position if rev _comp
	void adjust_rev_comp( int sizeSeq1, bool rev_comp);



	//complement the alignment for correct output
	void apply_rev_comp(int output_format);
	/**
	 * Fonction utilisée pour incrémenter la valeur de l'indice de départ dans la première séquence
	 * On peut ainsi tenir compte des gaps utilisés
	 */
	inline void incAlign();
	
	/**
	 * Fonction utilisée pour décrémenter la valeur de l'indice de départ dans la première séquence
	 * On peut ainsi tenir compte des gaps utilisés
	 */
	inline void decAlign();
	
	inline void incEnd1();
	inline void decStart1();

	inline void decEnd1();
	inline void incStart1();
	/**
	 * Méthode d'accès à l'attribut start1
	 * \return la valeur de start1
	 */
	inline int getStart1();
	
	/**
	 * Méthode d'accès à l'attribut end1
	 * \return la valeur de end1
	 */
	inline int getEnd1();
	
	/**
	 * Méthode d'accès à l'attribut start2
	 * \return la valeur de start2
	 */
	inline int getStart2();
	
	/**
	 * Méthode d'accès à l'attribut end2
	 * \return la valeur de end2
	 */
	inline int getEnd2();
	
	/**
	 * Méthode d'ajout d'un couple d'un caractère aux séquences
	 * \param c1 le caractère à ajouter à la première séquence
	 * \param c2 le caractère à ajouter à la seconde séquence
	 */
	 void addPair(char c1, char c2);

	/**
	 * Méthode permettant de récupérer la première séquence
	 * \return la première séquence de l'Alignment
	 */
	inline char* getSeq1();
	
	/**
	 * Méthode permettant de récupérer la seconde séquence
	 * \return la seconde séquence de l'Alignment
	 */
	inline char* getSeq2();
	 
	/**
	 * Méthode permettant de récupérer le string cigar
	 * \return le string cigar
	 */
	inline char* getCigar();
	/**
	 * Méthode permettant d'obtenir le caractère à l'indice en paramètre dans la première séquence
	 * \param indice l'indice du caractère dans le tableau séquence1
	 * \return le caractère à l'indice voulu dans la première séquence
	 */
	char getChar1(int indice);
	
	/**
	 * Méthode permettant d'obtenir le caractère à l'indice en paramètre dans la seconde séquence
	 * \param indice l'indice du caractère dans le tableau séquence2
	 * \return le caractère à l'indice voulu dans la seconde séquence
	 */
	char getChar2(int indice);	
	
	/**
	 * Méthode d'affichage d'un Aligment
	 */
	void affichage();
};


/**
 * Méthode d'obtention de la longueur de l'alignement
 * \return la longueur de l'objet Alignment apelant
 */
int Alignment::getLength()
{
	return length;
}

/**
 * Méthode d'obtention du nombre de mésappariements
 * \return le nombre de mésappariements de l'objet Alignment appelant
 */
int Alignment::getMis()
{
	return nb_mismatches;
}

/**
 * Méthode d'obtention du nombre de gaps
 * \return le nombre de gaps de l'objet Alignment appelant
 */
int Alignment::getGaps()
{
	return nb_gaps;
}

/**
 * Méthode d'incrémentation du nombre de mésappariements
 */
void Alignment::addMis()
{
	++nb_mismatches;
}

/**
 * Méthode d'incrémentation du nombre de gaps
 */
void Alignment::addGap()
{
	++nb_gaps;
}

/**
 * Fonction utilisée pour incrémenter la valeur de l'indice de départ dans la première séquence
 * On peut ainsi tenir compte des gaps utilisés
*/
void Alignment::incAlign()
{
	++start1;
}

   void Alignment::incEnd1()
{
	++end1;

}
   void Alignment::decStart1()
{
	--start1;

}

 void Alignment::decEnd1()
{
	--end1;

}
   void Alignment::incStart1()
{
	++start1;

}
/**
 * Fonction utilisée pour décrémenter la valeur de l'indice de départ dans la première séquence
 * On peut ainsi tenir compte des gaps utilisés
 */
void Alignment::decAlign()
{
	--start1;
}

/**
 * Méthode d'accès à l'attribut start1
 * \return la valeur de start1
 */
inline int Alignment::getStart1()
{
	return start1;
}

/**
 * Méthode d'accès à l'attribut end1
 * \return la valeur de end1
 */
inline int Alignment::getEnd1()
{
	return end1;
}

/**
 * Méthode d'accès à l'attribut start2
 * \return la valeur de start2
 */
inline int Alignment::getStart2()
{
	return start2;
}
	
/**
 * Méthode d'accès à l'attribut end2
* \return la valeur de end2
 */
inline int Alignment::getEnd2()
{
	return end2;
}

/**
 * Méthode permettant de récupérer la première séquence
 * \return la première séquence de l'Alignment
 */
inline char* Alignment::getSeq1()
{
	return sequence1;
}

/**
 * Méthode permettant de récupérer la seconde séquence
 * \return la seconde séquence de l'Alignment
 */
inline char* Alignment::getSeq2()
{
	return sequence2;
}


/**
 * Méthode permettant de récupérer le string cigar
 * \return le string cigar
 */
inline char* Alignment::getCigar()
{
  return cigar;
}


#endif
