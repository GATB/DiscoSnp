#ifndef WITHGAP_H
#define WITHGAP_H

/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file withgap.h
 * \brief Module Withgap, responsable de la réalisation de l'alignement avec gap
 * \author Damien Fleury
 * \author Guillaume Rizk
 * \version 6.1
 * \date 21/01/2009
 */

#include <pthread.h>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include <sys/sysinfo.h>
#include <sys/time.h>

class Bank;
class Index;
class Hit;
class Alignment;
class Stat;
class Doublon;
#define LEFT 0
#define RIGHT 1

#define AL_VALID 1
#define AL_DOUBLON 2
#define AL_ERR 0

/**
 * La structure caseMatrix permet une mémorisation du chemin utilisé pour déterminer
 * le meilleur alignement dans la matrice de l'algorithme
 */
typedef struct
{
	/**
	 * Le score correspondant à la case
	 */
	int score;
	/**
	 * L'origine de la valeur du score
	 */
	int path;
  	/**
	 * Le nombre de gap
	 */
  //	int ngaps;
	/**
	 * Le nombre de mismatch
	 */
  //	int nmis;
}caseMatrix;


/// Les différentes valeurs pour l'attribut path de la structure caseMatrix
#define APP 0
#define MES 1
#define GA1 2
#define GA2 3


/**
 * Fonction de recherche du meilleur alignement selon le modèle de Wunsch et Needleman
 * \param al l'objet Alignment qui constituera le résultat
 * \param seq1 le pointeur du tableau des caractères de la première banque
 * \param seq2 le pointeur du tableau des caractères de la seconde banque
 * \param size le nombre de bases à parcourir
 * \param num_gaps le nombre maximal de gaps autorisés dans la partie de l'alignement
 * \param side indiquant si on fait l'alignement sur le côté gauche ou droit de la graine
 */
//void WN_alignment(Alignment* al, char* seq1, char* seq2, int size, int num_gaps, int side);
int WN_alignment(Alignment* al, char* seq1, char* seq2, int size1, int size2, int num_gaps, int side, int num_seed, int max_mis);

int WN_big_alignment(Alignment* al, char* seq1, char* seq2, int size1, int size2, int num_gaps, int side, int num_seed, int max_mis);

int WN_very_big_alignment(Alignment* al, char* seq1, char* seq2, int size1, int size2, int num_gaps, int side, int num_seed, int max_mis);

/**
 * Fonction permettant de calculer un score d'alignement obtenu entre 2 séquences
 * \param al l'objet Alignment qui constituera le résultat
 * \param seq1 le pointeur du tableau des caractères de la première banque
 * \param seq2 le pointeur du tableau des caractères de la seconde banque
 * \param left1 le nombre de bases présentes à gauche dans la première séquence
 * \param left2 le nombre de bases présentes à gauche dans la seconde séquence
 * \param right1 le nombre de bases présentes à droite dans la première séquence
 * \param right2 le nombre de bases présentes à droite dans la seconde séquence
 * \param max_mis le nombre maximal de mésappariement pour conserver l'alignement
 * \param num_gaps le nombre maximal de gaps autorisés dans les alignements
 * \return Le nombre de mésappariements obtenus dans l'alignement
 */
int withgap_scoreHit(Alignment* al, char* seq1, char* seq2, int left1, int left2, int right1, int right2, int max_mis, int num_gaps);


/**
 * Méthode qui à partir des positions d'une graine va effectuer tous les alignements possibles avec les séquences des banques
 * \param ff le fichier de sortie, où seront affichés les alignements obtenus
 * \param BK1 le pointeur de la première banque de séquences
 * \param BK2 le pointeur de la seconde banque de séquences
 * \param I1 le pointeur de l'index de la première banque de séquences
 * \param start num premiere query a faire
 * \param end num derniere query a faire
 * \param idpc le pourcentage de ressemblances minimal entre les 2 séquences pour qu'un alignement soit conservé
 * \param num_gaps le nombre maximal de gaps autorisés dans les alignements
 * \return le nombre d'alignements obtenus
 */
int withgap_align(FILE *ff, Bank *BK1, Bank *BK2, Index * I1, int idpc, bool rev_comp, int num_gaps, char ** tabprec, char ** tabnt, char **tabprec7, Stat * St,Doublon * Doub,int start,int end);



#endif
