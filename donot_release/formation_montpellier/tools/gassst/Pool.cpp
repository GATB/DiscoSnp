/* This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software.  You can  use, 
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info". 
 */
/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file Pool.h
 * \brief Classe Pool  utilisee pour allouer/desallouer rapidement memoire necessaire a liste chainee pour index avec graine  >14
 * \author Guillaume Rizk
 * \date 09/06/2010
 */


 #include "Pool.h"


/**
 * Constructeur par défaut
 */
Pool::Pool()
{
  n_pools = 0; n_cells=0;
  //allocation table de pool :
  tab_pool = (cell**)  malloc(N_POOL*sizeof(cell *) );

  //allocation de la premiere pool : 
  pool_courante =(cell*)  malloc(TAI_POOL*sizeof(cell) );
  tab_pool[n_pools] = pool_courante;
  n_pools++;
}



/**
 * Destructeur
 */
Pool::~Pool()
{

 unsigned  int i;

  for(i=0;i<n_pools;i++)
    {
      free( tab_pool[i] );
    }

  free(tab_pool);
}


void Pool::empty_all()
{

 unsigned  int i;

  for(i=1;i<n_pools;i++) // garde la premiere pool pour usage futur 
    {
      free( tab_pool[i] );
    }

  //on repasse sur premiere pool 
  pool_courante = tab_pool[0];
  n_cells=0;
  n_pools=1;

}






cell *  Pool::allocate_cell_in_pool()
{

  // ncells = nb de cells deja utilisees
  if (n_cells <TAI_POOL) 
    {
      n_cells ++;
      return (pool_courante + n_cells -1 );
    }
  else // la piscine est pleine, on en alloue une nouvelle
    {
      pool_courante =(cell*)  malloc(TAI_POOL*sizeof(cell) );
      tab_pool[n_pools] = pool_courante;
      n_pools++;
      n_cells = 1;
      return (pool_courante);

    }

}


