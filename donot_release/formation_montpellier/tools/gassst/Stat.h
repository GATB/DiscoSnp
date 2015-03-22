#ifndef STAT_H
#define STAT_H


/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file Stat.h
 * \brief Classe Stat, comptage des align trouvés par sequence query
 * \ et gestion meilleur align
 * \author Guillaume Rizk
 * \version 6.1
 * \date 20/02/2009
 */
#include "misc.h"
#include "Alignment.h"
#include <list>
using namespace std;


class Stat{

  // public:
 private :

  int nb_align; // du coup par partition seulement

  int best_score;

  int worse_score; //plus mauvais score de la liste
 public :

  std::list < Alignment >  listal;

  //construct
  Stat(int nb_reads) ;


  ~Stat();

  inline int quota_atteint ();

  inline void incr_query();

  inline void add_al(Alignment al);

  inline void reset();

  inline int good_enough();

  inline int get_nb_al();

  void print_als(FILE * ff);

};


  inline int Stat::get_nb_al()
{
  return nb_align;
};



 inline int Stat::good_enough()
{
  return (quota_atteint() && best_score==0 && worse_score==0  ) ;
}

//pour le moment stockage de liste d align de scores egaux , les meilleurs trouves
inline void Stat::add_al(Alignment al)
{
  int insere = 0;
  int score_courant;
  int score = al.getMis() + al.getGaps() ;
  /*
  if(score < best_score) // si on trouve meilleur ,efface  la liste
    {
      best_score = score;
      listal.clear();
      nb_align = 1;
      listal.push_back(al);
     
    }
  else if (score == best_score && !quota_atteint() )//si aussi bon, on lajoute a la liste
    {
      listal.push_back(al);
      nb_align++;
    }
  */

//version avec liste dalign de score differents
//insertion si meilleur que le plus mauvais ou si il reste de la place
	if(score < worse_score || (!quota_atteint ()) ) 
	{
		list <Alignment>::iterator it_list = listal.begin();
		while(it_list != listal.end())
		{
			score_courant = (*it_list).getMis() + (*it_list).getGaps() ;
			if(score_courant <= score)
			{
				listal.insert(it_list,al);
				it_list=listal.end();
				insere = 1;
			}
			else	++it_list;
		}
		
		
  		if(!insere) {listal.push_back(al); }//au bout, cest le meilleur
		nb_align++;
         
       //si il y en a de trop on vire le plus mauvais 
        if(MAXHITS  && (nb_align>MAXHITS)) {listal.pop_front(); nb_align--;}

       // ali =  &(* listal.begin() );
       // worse_score=  ali->getMis() + ali->getGaps();
         worse_score=  listal.front().getMis() + listal.front().getGaps();

       //   ali =  &(* listal.end() );
       //best_score=  ali->getMis() + ali->getGaps();
        best_score=  listal.back().getMis() + listal.back().getGaps();
	
    }
  
  //sinon , on ne fait rien
 
}



 inline int Stat::quota_atteint ()
  {
    return ( MAXHITS  && (nb_align>=MAXHITS));
  }


inline  void Stat::incr_query()
  {
    nb_align++;
  }
  


inline  void Stat::reset()
  {
    nb_align = 0;
    if(BESTAL)
      {
	best_score = 9999;
	worse_score = 9999;
	listal.clear();
      }
  }
  

#endif
