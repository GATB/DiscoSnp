#ifndef MISC_H
#define MISC_H

/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file misc.h
 * \brief Module Misc, contenant des fonctions utilitaires, et de gestion des options et des erreurs
 * \author Dominique Lavenier
 * \author Damien Fleury
 * \version 5.2
 * \date 28/08/2008
 */


#include <iostream>
#include <stdlib.h>
#include <string.h>


/// Variables globales utiles à la gestion des options
/// Elles permettent la mémorisation des options spécifiées
extern char KEY0[64];
extern char KEY1[64];
extern char PARAM0[64][64];
extern char PARAM1[64][64];
extern int  NB_KEY1;
extern int  NB_KEY0;
extern char PROG_NAME[1024];
extern int NBTHREADS;
extern int MAXHITS;
extern int MAXPOS;
extern int MAXPOS_AMONT;
extern int BESTAL;
extern int SLEVEL;
extern int BITSTAT;
extern int NUMMAX;
extern int NUM_ALIGN_COURANT;
extern int NUMGAPS_AUTO;
extern int GLOBAL_GLOBAL;


/*
 * Fontion minimum de deux valeurs
 * \param a un entier
 * \param b un entier
 * \return la valeur minimale entre a et b
 */
inline int min(int a, int b)
{
	return a < b ? a : b;
}

/*
 * Fontion maximum de deux valeurs
 * \param a un entier
 * \param b un entier
 * \return la valeur maximale entre a et b
 */
inline int max(int a, int b)
{
	return a > b ? a : b; 
}

/**
 * Fonction d'impession d'un message d'erreur et de fermeture du programme
 * \param msg la chaîne du message d'erreur
 */
void ExitError(char *msg);


/**
 * Fonction d'affichage d'un erreur dans les paramètres du programme et de fermeture du programme
 */
void SyntaxError();


/**
 * Méthode de génération des différentes options du programme
 * \param key le caractère permettant d'identifier l'option
 * \param strg_in le nom de l'option
 * \param opt indique si l'option est obligatoire, 1 si oui, 0 sinon
 */
void option(char key, char *strg_in, int opt);


/**
 * Méthode permettant de vérifier les options spécifiées lors de l'exécution du programme et d'en récupérer les paramètres
 * \param key le caractère permettant d'identifier l'option
 * \param strg_out la valeur utilisée pour l'option
 * \param argc le nombre de paramètres utilisés lors de l'appel au programme
 * \param argv un pointeur vers le tableau des paramètres du programme
 * \param opt indique si l'option est obligatoire, 1 si oui, 0 sinon
 * \return 1 si l'option a été spécifiée, 0 si l'option est facultative et non utilisée
 */
int getoption(char key, char *strg_out, int argc, char *argv[], int opt);


/**
 * Méthode de vérification des options
 * \param argc le nombre de paramètres utilisés lors de l'appel au programme
 * \param argv un pointeur vers le tableau des paramètres du programme
 */
void checkoption(int argc, char *argv[]);

#endif
