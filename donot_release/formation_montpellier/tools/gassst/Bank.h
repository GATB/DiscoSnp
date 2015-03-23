#ifndef BANK_H
#define BANK_H

/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file Bank.h
 * \brief Classe Bank, responsable de la récupération des données dans les banques de séquences
 * \author Dominique Lavenier
 * \author Damien Fleury
 * \version 5.2
 * \date 28/08/2008
 */

#include <iostream>

/// Nombre maximal de partitions
#define MAX_PART 1024
/// Taille d'une ligne d'un fichier au format fasta
#define SIZE_LINE 1024
/// TAille du nom des banques de séquences
#define TAILLE_NOM 1024

/**
 * \class Bank, Une banque de séquences d'ADN
 * \brief Cette classe définit une banque de séquences d'ADN, ainsi que les méthodes permettant d'en récupérer les informations
 */
class Bank{

	public:
	/**
	 * Nom de la banque
	 */
	char fileBank[TAILLE_NOM];
	/**
	 * Nombre de séquences total
	 */
	int  nb_tot_seq;
	/**
	 * Nombre de résidus total
	 */
	long long  int  nb_tot_res;
	/**
	 * Nombre de partitions
	 */
	int  nb_part;
	/**
	 * Numéro de la partition en cours de traitement
	 */
	int  num_part;
	/**
	 * Numéro de la partition suivante
	 */
	int  next_part;
	/**
	 * Adresses de départ des partitions dans le fichier
	 */
	long* start_offset;
	/**
	 * Adresses de fin des partitions dans le fichier
	 */
	long* stop_offset;
	/**
	 * Nombre de séquences dans chaque partition
	 */
	int* nb_seq;
	/**
	 * Image mémoire de la partie de banque en cours de traitement
	 */
	char* data;
	/**
	 * Tableau des index des sequences , en int donc taille max partition = 2G
	 */
	int* seq;
	/**
	 * Tableau des index des commentaires
	 */
	int* com;
	/**
	 * Tableau des tailles des séquences
	 */
	long long* size;
	
	/**
	 * Tableau des positions globales deu début des séquences
	 */
	long long* pos_seq;
	/**
	 * Longueur maximale des séquences de la partition indexée de la banque
	 */
	int tailleMaxSeq;

	int tSeq;


	/**
	 * Constructeur de banque de séquences
	 * \param fname, un pointeur de caractères contenant le nom de la banque de séquences
	 * \param size_max, la taille maximale de la banque de séquences
	 * \param FILE, fichier de sortie, pour mettre header avec noms des contig
	 * \param bankref indique qu on lit la banque de reference
	 */
	Bank(char *fname,long long size_max, FILE *ff, char bankref);
	
	/**
	 * Constructeur de banque par recopie
	 * \param bk, une banque de séquences
	 */
	Bank(const Bank& bk);
	
	/**
	 * Opérateur d'affectation
	 * \param bk une banque de séquences
	 * \return l'objet Bank affecté
	 */
	Bank& operator=(const Bank& bk);
	
	/**
	 * Destructeur de Bank
	 */
	~Bank();
	
	/**
	* Méthode d'affichage d'une banque de séquences
	*/
	void writeInfoBank();
	
	/**
	* Méthode de réinitialisation des partitions d'une banque de séquences
	*/
	void resetBank();
	
	/**
	 * Méthode permettant d'indexer une partition de la banque
	 * \param lx, un booléen qui indique si on utilise le filtre Low Complexity
	 * \return 1 si une partition a été indexée, 0 si toutes les partitions sont déjà été indexées
	 */
	int readBank(bool lx);
	/**
	 * Méthode permettant de faire le "reverse complement" des séquences de la partition courante
	 * de la banque
	 * Toutes les séquences de la partition courante de la banque seront inversées et leurs bases
	 * seront remplacées par les bases complémentaires
	 */
	void reverseComplement();
};

#endif
