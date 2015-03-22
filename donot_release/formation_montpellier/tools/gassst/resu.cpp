/* This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software.  You can  use, 
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info". 
 */
/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file resu.cpp
 * \author Guillaume Rizk
 * \date 01/02/2011
 */
 
 #include "resu.h"


/**
 * Constructeur par défaut
 */
resu::resu() : score_al(0), index(0)
{
}

/**
 * Constructeur de resu
 */
resu::resu(
	   unsigned int score_in,
	   unsigned int index_in
	   )
{
  score_al = score_in;
  index= index_in;
}

resu::resu(const resu &copyin)   // Copy constructor to handle pass by value.
{                            
  score_al = copyin.score_al;
  index = copyin.index;
}
	 

	 
resu& resu::operator=(const resu &rhs)
{
  this->score_al = rhs.score_al;
  this->index = rhs.index;
  return *this;
}

int resu::operator==(const resu &rhs) const
{
  if( this->index != rhs.index) return 0;
  if( this->score_al != rhs.score_al) return 0;
  return 1;
}

// tri dans ordre decroissant
int resu::operator<(const resu &rhs) const
{
  return( this->score_al > rhs.score_al);
}
