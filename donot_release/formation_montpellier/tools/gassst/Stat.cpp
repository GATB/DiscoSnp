/* This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software.  You can  use, 
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info". 
 */

/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file Stat.cpp
 * \brief Classe Stat, comptage des align trouvés par sequence query
 * \author Guillaume Rizk
 * \version 6.1
 * \date 20/02/2009
 */

#include "Stat.h"
#include "misc.h"
#include "constants.h"
#include "display.h"


//construct
Stat::Stat(int nb_reads) : nb_align(0), best_score(9999) , worse_score(9999)
    {

    }



//destruct
Stat::~Stat()
{
  
}


//faire le mutex dedans ou dehors ? 
//pour le moment mutexé du dehors
  void Stat::print_als(FILE * ff)
{

  // int tai = listal.size();
 // Alignment * ali;
  list<Alignment>::iterator it = listal.begin();
  while(it != listal.end())
    {
      //  ali =  &(*it_list);
       //  display_withgap(ff,ali,ali->j1,ali->j2,ali->e_value);
      //  display_withgap(ff,(*it_list),(*it_list).j1,(*it_list).j2,(*it_list).e_value);
       

  display_withgap(ff,&(*it),it->j1,it->j2,it->e_value);
		it++;

    }
  
      //avec vector
      // for(int i=0; i<tai; i++ ) // on ecrit chque align
	   //{
      //display_withgap(ff,&(listal[i]),(listal[i]).j1,(listal[i]).j2,(listal[i]).e_value);
      //}

}
