#ifndef GAPLESS_H
#define GAPLESS_H

/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file gapless.h
 * \brief Module Gapless, responsable de la réalisation de l'alignement sans gap
 * \author Dominique Lavenier
 * \author Damien Fleury
 * \version 5.2
 * \date 28/08/2008
 */


#include <pthread.h>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include <sys/sysinfo.h>
#include <sys/time.h>


/// Nombre maximal d'occurences des graines
#define NB_MAX_HIT 100000


class Bank;
class Index;
class Hit;
class Stat;
/// Verrou pour les entrées sorties fichiers de threads
extern pthread_rwlock_t verrou_out;
class Doublon;

/**
 * Fonction permettant de calculer un score d'alignement obtenu entre 2 séquences
 * \param seq1 le pointeur du tableau des données de la première banque
 * \param seq2 le pointeur du tableau des données de la seconde banque
 * \param left le nombre de bases à parcourir à gauche
 * \param right le nombre de bases à parcourir à droite
 * \param num_seed le code de la graine précédente
 * \param max_mis le nombre maximal de mésappariement pour conserver l'alignement
 * \return Le nombre de mésappariements obtenus dans l'alignement
 */
int scoreHit(char* seq1, char* seq2, int left, int right, int num_seed, int max_mis);

/**
 * Méthode qui à partir des positions d'une graine va effectuer tous les alignements possibles avec les séquences des banques
 * \param ff le fichier de sortie, où seront affichés les alignements obtenus
 * \param BK1 le pointeur de la première banque de séquences
 * \param BK2 le pointeur de la seconde banque de séquences
 * \param I1 le pointeur de l'index de la première banque de séquences
 * \param num_seed l'indice de la graine utilisée pour déterminer les alignements
 * \param idpc le pourcentage de ressemblances minimal entre les 2 séquences pour qu'un alignement soit conservé
 * \return le nombre d'alignements obtenus
 */
//int gapless_align(FILE *ff, Bank *BK1, Bank *BK2, Index * I1, Index * I2, int num_seed, int idpc, bool rev_comp, Stat * St);


int gapless_align(FILE *ff, Bank *BK1, Bank *BK2, Index * I1, int idpc, bool rev_comp, int num_gaps,char ** tabprec,char ** tabnt,char **tabprec7,Stat * St,Doublon * Doub, int start, int end );

#endif
