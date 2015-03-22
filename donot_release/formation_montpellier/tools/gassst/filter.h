#ifndef FILTER_H
#define FILTER_H

/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file filter.h
 * \brief Module Filter, définit le filtre Low Complexity qui détecte les zones non pertinentes des séquences afin de ne pas les indexer
 * \author Dominique Lavenier
 * \author Damien Fleury
 * \version 5.2
 * \date 28/08/2008
 */

/// Fenêtre du filtre Low complexity
#define WS 12


using namespace std;


/**
 * Méthode permettant d'appliquer le filtre Low Complexity
 * \param data le tableau des données de la banque de séquences
 * \param id l'index du début de la séquence
 * \param lenseq la longueur de la séquence
 * \return 1 si l'opération s'est déroulée correctement
 */
int filterLowComplexity (char* data, int lenseq);

#endif

