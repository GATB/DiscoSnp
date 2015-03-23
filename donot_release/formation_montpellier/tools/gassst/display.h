#ifndef DISPLAY_H
#define DISPLAY_H

/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file display.h
 * \brief Module Display, responsable de l'enregistrement du résultat dans un fichier de sortie
 * \author Dominique Lavenier
 * \author Damien Fleury
 * \version 5.2
 * \date 28/08/2008
 */

#include <iostream>

class Alignment;

/// Variable contenant le type de format de sortie du programme
extern int OUTPUT_FORMAT;


/**
 * Méthode de vérification du format de sortie
 * \param k le type de format demandé
 */
void checkOutputFormat(int k);


/**
 * Fonction permettant d'enregistrer un alignement dans le fichier de sortie
 * \param ff le pointeur vers l'objet FILE du fichier de sortie
 * \param s1 le pointeur vers la première séquence
 * \param s2 le pointeur vers la seconde séquence
 * \param c1 le pointeur vers le début du commentaire de la séquence de la première banque
 * \param c2 le pointeur vers le début du commentaire de la séquence de la seconde banque
 * \param start1 la position de départ de l'alignement pour la première séquence
 * \param start2 la position de départ de l'alignement pour la seconde séquence
 * \param len la longueur de l'alignement (longueur de la plus courte séquence)
 * \param l_seq1 la longueur de la première séquence, la plus longue
 * \param eval la valeur de la E-value
 * \param rev_comp un booléen indiquant si on a un alignement sur la séquence normale (false) ou inversée et complémentée (true)
 */
void display_ungap(FILE *ff, char* s1, char* s2, char* c1, char* c2, int start1, int start2, int len, int l_seq1, double eval, bool rev_comp);


/**
 * Fonction permettant d'enregistrer un alignement avec gap(s) dans le fichier de sortie
 * \param ff le pointeur vers l'objet FILE du fichier de sortie
 * \param al l'objet Alignment contenant les informations sur l'alignement
 * \param c1 le pointeur vers le début du commentaire de la séquence de la première banque
 * \param c2 le pointeur vers le début du commentaire de la séquence de la seconde banque
 * \param eval la valeur de la E-value
 */
void display_withgap(FILE *ff, Alignment* al, char* c1, char* c2, double eval);

#endif
