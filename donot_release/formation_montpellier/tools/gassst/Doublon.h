#ifndef DOUBLON_H
#define DOUBLON_H


//#define NB_BITS 67108864 // soit 2^26 
//#define NB_BITS 536870912 //soit 2^29
//#define NB_BITS 1048576 //20
#define NBITS 20 //soit 2^29

/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file Doublon.h
 * \brief Classe Stat, comptage des align trouvés par sequence query
 * \author Guillaume Rizk
 * \version 6.1
 * \date 09/06/2009
 */
#include "misc.h"
#include <vector>

class Doublon{

 public:
  //tableau de bits disant si align deja trouve a cette position ou non
  char * align;
  char * align2;

  //vecteur des positions trouvees
  //est ce plus rapide avec des char * ? pointant vers case du tab align ?
  std::vector < unsigned int>  listpos;
  std::vector < unsigned int>  listpos2;

  int length;
  unsigned int masque; 
  
  //constructor
 Doublon(int taille) /*: length(NB_BITS/8)*/  // independant taille donnee
  {
    length = (1 << BITSTAT) / 8 ;

    //  align = new char [(taille/8) +100]; // a ajuster
    //for(int i=0; i<length; i++) align[i]=0;
    //align2 = new char [(taille/8) +100]; // a ajuster
    //for(int i=0; i<length; i++) align2[i]=0;

    align = new char [length]; 
    for(int i=0; i<length; i++) align[i]=0;
    align2 = new char [length]; 
    for(int i=0; i<length; i++) align2[i]=0;

    masque = (1<<BITSTAT)-1; //ca fait 26 bits a 1
    
    ///  printf("Taille tableau doublon : %i Ko\n",(length)/1024);
  };
  
  //destructor
  ~Doublon()
    {
      delete [] align;
      listpos.clear();
      delete [] align2;
      listpos2.clear();
    };
  

  //met a 1 la posi, indique quon a trouve un align a cette position
  void  setpos (unsigned int pos);
  void  setpos2 (unsigned int pos);

  // est ce qu on a deja trouve un align a cette position
  // 0 non
  // !=0 oui
 inline unsigned char  getpos (unsigned int pos);
 inline unsigned char  getpos2 (unsigned int pos);

  //remise a 0 du tableau align
  void reset();

};



 unsigned char  Doublon::getpos (unsigned int pos)
{
///////  pos = pos >> 3; // numero de la zone de 8nt correspondante
  //puis recup bit correspondant
  pos = pos & masque; // prend le modulo 67108864
  unsigned int index =  pos >> 3  ;
  unsigned char masque_res = 128 >> (pos & 7) ; 
  unsigned char resu ;

  resu = align[index];

   return(resu & masque_res);
  // return resu;

}

 unsigned char  Doublon::getpos2 (unsigned int pos)
{
///////  pos = pos >> 3; // numero de la zone de 8nt correspondante
  //puis recup bit correspondant
  pos = pos & masque; // prend le modulo 67108864
  unsigned int index =  pos >> 3  ;
  unsigned char masque_res = 128 >> (pos & 7) ; 
  unsigned char resu ;

  resu = align2[index];

   return(resu & masque_res);
  // return resu;

}

#endif
