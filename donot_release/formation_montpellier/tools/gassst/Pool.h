#ifndef POOL_H
#define POOL_H

/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file Pool.h
 * \brief Classe Pool, utilisee pour allouer/desallouer rapidement memoire necessaire a liste chainee pour index avec graine  >14
 * \author Guillaume Rizk
 * \date 09/06/2010
 */
#include <stdlib.h>
#include "Seed.h"

typedef struct cell
{
  long long graine; //8
  int offset_seed; //4
  int nb_seed; //4
  struct cell *suiv; // 8
} cell;




#define  TAI_POOL 10000000 // represente 228 Mo  ( 24 octet par cell)
#define  N_POOL   1000 // soit 10 G cells max 
/**
 * \class Pool, 
 * \brief Cette class définit une pool memoire pour allocation rapide de la table de hachage utilisee quand seed >14
 */
class Pool{
	public:
	
	/**
	 * table de cell, pour usage courant, 
	 */
	cell * pool_courante;
	/**
	 * stockage de tous les pointeurs  pool
	 */
	cell ** tab_pool;
	/**
	 *  nombre de piscines remplies
	 */
	unsigned int n_pools;

	/**
	 *  niveau de remplissage de la piscine courante 
	 */
	unsigned int n_cells;



	/**
	 * Constructeur par défaut
	 */
	Pool();
	
	/**
	 * alloue une cellule dans la piscine
	 */
	cell *  allocate_cell_in_pool();


	/**
	 * vide toutes piscines
	 * (garde juste une pool vide)
	 */
	void  empty_all();

	/**
	 * Destructeur
	 */
	~Pool();

};

#endif
