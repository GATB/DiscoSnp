#ifndef CONSTANTS_H
#define CONSTANTS_H

/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file constants.h
 * \brief Module Constants, contenant toutes les constantes utilisées par le logiciel
 * Ces constantes peuvent être manipulées pour paramètrer le fonctionnement de Gassst
 * \author Damien Fleury
 * \author Guillaume Rizk
 * \version 5.2
 * \date 28/08/2008
 */

#define GASSST_VERSION "1.29"

//pour compteurs sur les differents filtres

//2 compteurs interessants, ne marchent qu avec un seul thread
//#define CPT_FILT
//#define T_CLOCK

//pour sortie du numero du contig a la place de son nom 
//#define NUMSEQSORTIE

// Module Gassst
////////////////////////////////////////////////////////////////////////////////////////
/// Constante permettant de définir si on active ou non le calcul du temps d'exécution
/// (A ajouter pour intégrer le calcul du temps)
#define EXEC_TIME

/// Taille de graine par défaut
#define DEFAULT_SIZE_SEED 12
#define INDEX_STRIDE 1

//inutile maintenant, car soit automatique -> max gaps = max erreurs
// soit défini par l utilisateur
#define DEFAULT_GAP_ALLOWED 1



/// Activation/désactivation du filtre Low Complexity
#define DUST_OFF false
#define DUST_ON true

/// Activation/désactivation de la recherche reverse complement
#define REV_OFF false
#define REV_ON true

/// Le nombre de threads utilisés dans le programme par défaut 
#define DEFAULT_NBTHREADS 1


/// Le nombre de hits renvoyés par reads (0 = nolimit)
#define DEFAULT_MAXHITS 100

#define DEFAULT_MAXPOS 50
#define DEFAULT_SLEVEL 2

#define NBITS_HACH 24
#define NT_HACH 14 // inutilise avec nouvelle hashfonction
#define NHACH (1<<NBITS_HACH)
//marge utilisee pour CPT MAXPOS_AMONT
// on autorise au plus 20 * le nombre moyen de hits  par graines
#define MARGE 20 
#define MARGE_HAUTE 200 

/// La partition des graines attribuée à un thread
#define SIZE_PARTITION_THREAD 10000

//par defaut ne renvoit pas forcement les meilleurs 
#define DEFAULT_BESTAL 1

////////////////////////////////////////////////////////////////////////////////////////



// Module code
////////////////////////////////////////////////////////////////////////////////////////
/// Valeurs nécessaires au calcul de la e-value
#define CONST_H 1.31
#define CONST_K 0.711
#define LAMBDA  1.37
// Taille de la table precalculee utilisee (la table fait 2^(4*TAI_TAB) octets)
#define TAI_TAB 4
//pour second filtre plus long, plus puissant
#define TAI_BIG_TAB 5

////////////////////////////////////////////////////////////////////////////////////////



// Modules code et gapless
////////////////////////////////////////////////////////////////////////////////////////
/// Score obtenu pour un appariement
#define MATCH 1
/// Score obtenu pour un mésappariement
#define MISMATCH 3
////////////////////////////////////////////////////////////////////////////////////////



// Module withgap
////////////////////////////////////////////////////////////////////////////////////////
/// Les coûts des différentes possibilités dans la matrice de calcul de l'alignement optimal avec gaps
#define COUT_MATCH 0 
#define COUT_MISMATCH -1
#define COUT_GAP -1
 



/// Le caractère représentant un gap dans un alignement
#define CHAR_GAP '-'
////////////////////////////////////////////////////////////////////////////////////////



// Module display
////////////////////////////////////////////////////////////////////////////////////////
/// Les  types de formats de sortie possibles
#define STD_OUTPUT_FORMAT 0
#define M8_OUTPUT_FORMAT 8
// type de format facilement convertible en sam, une ligne par align, avec cigar
#define SAM_READY_FORMAT 1


/// Paramètres de la recherche d'alignement : normal ou reverse complement
#define NORMAL_ALIGN false
#define REV_COMP_ALIGN true
////////////////////////////////////////////////////////////////////////////////////////
 
#endif
