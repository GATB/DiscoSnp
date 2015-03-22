/* This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software.  You can  use, 
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info". 
 */

/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file Doublon.cpp
 * \brief Classe Stat, comptage des align trouvés par sequence query
 * \author Guillaume Rizk
 * \version 6.1
 * \date 09/06/2009
 */

#include "Doublon.h"
#include "misc.h"
#include "constants.h"


void  Doublon::setpos (unsigned int pos)
{
  /////////   pos = pos >> 3; //  par zone de 8 nt
  //puis stockage sur 1 bit  de cette zone de 8 nt  : 
  pos = pos & masque; // prend le modulo 67108864
   unsigned int index =  pos >> 3  ; //   / 8
  //  unsigned int reste =  pos & 7;  // modulo 8 
  //128 : 10000000   puis decalage du 1 de (pos mod 8) cases
  unsigned char masque_res = 128 >> (pos & 7) ; 
  
  listpos.push_back(index);
  
  align[index] =  align[index] | masque_res;
  




}
void  Doublon::setpos2 (unsigned int pos)
{
  /////////   pos = pos >> 3; //  par zone de 8 nt
  //puis stockage sur 1 bit  de cette zone de 8 nt  : 
  pos = pos & masque; // prend le modulo 67108864
   unsigned int index =  pos >> 3  ; //   / 8
  //  unsigned int reste =  pos & 7;  // modulo 8 
  //128 : 10000000   puis decalage du 1 de (pos mod 8) cases
  unsigned char masque_res = 128 >> (pos & 7) ; 
  
  listpos2.push_back(index);
  
  align2[index] =  align2[index] | masque_res;
  




}

/*
 unsigned char  Doublon::getpos (unsigned int pos)
{

  unsigned int index =  pos >> 3  ;
  unsigned char masque = 128 >> (pos & 7) ; 
  unsigned char resu ;

  resu = align[index];

  return(resu & masque);

}
*/
void  Doublon::reset ()
{
  int tai = listpos.size();

   for(int i=0; i<tai; i++ )
    {
      align[ listpos[i] ] = 0; 
      // printf("%i ;"lispos[i]);
    }
  
  listpos.clear();

  tai = listpos2.size();

   for(int i=0; i<tai; i++ )
    {
      align2[ listpos2[i] ] = 0; 
      // printf("%i ;"lispos[i]);
    }
  
  listpos2.clear();
  //  for(int i=0; i<length; i++) align[i]=0; // a faire avec vector
}
