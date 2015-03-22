/* This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software.  You can  use, 
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info". 
 */
/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file misc.cpp
 * \brief Module Misc, contenant des fonctions de gestion des options et des erreurs
 * \author Dominique Lavenier
 * \author Damien Fleury
 * \version 5.2
 * \date 28/08/2008
 */
 
#include "misc.h"
#include "constants.h"
#include <stdio.h> 

char KEY0[64];
char KEY1[64];
char PARAM0[64][64];
char PARAM1[64][64];
char PROG_NAME[1024];

int  NB_KEY1 = 0;
int  NB_KEY0 = 0;
int NBTHREADS;
int MAXHITS;
int MAXPOS;
int MAXPOS_AMONT;
int BESTAL;
int SLEVEL;
int BITSTAT;
int NUMGAPS_AUTO;
int NUMMAX;
int NUM_ALIGN_COURANT;
/**
 * Fonction d'impession d'un message d'erreur
 * \param msg la chaîne du message d'erreur
 */
void ExitError(char *msg)
{
	fprintf (stderr,"%s\n",msg);
	exit(0);
}


/**
 * Fonction d'affichage d'un erreur dans les paramètres du programme
 */
void SyntaxError()
{
	int i;
	fprintf(stderr,"Gassst version %s\n",GASSST_VERSION);

	fprintf (stderr,"\nSyntax Error\n  Syntax: %s",PROG_NAME);
	for (i=0; i<NB_KEY1; i++)
		fprintf (stderr," -%c <%s>",KEY1[i],PARAM1[i]);
	for (i=0; i<NB_KEY0; i++)
		fprintf (stderr," [-%c <%s>]",KEY0[i],PARAM0[i]);
	fprintf(stderr,"\n");
	ExitError(" ");
}


/**
 * Méthode de génération des différentes options du programme
 * \param key le caractère permettant d'identifier l'option
 * \param strg_in le nom de l'option
 * \param opt indique si l'option est obligatoire, 1 si oui, 0 sinon
 */
void option(char key, char *strg_in, int opt)
{
	if (opt==1)
	{
		KEY1[NB_KEY1]=key;
		strcpy(PARAM1[NB_KEY1],strg_in);
		NB_KEY1++;
	}
	else 
	{
		KEY0[NB_KEY0]=key;
		strcpy(PARAM0[NB_KEY0],strg_in);
		NB_KEY0++;
	}
}


/**
 * Méthode permettant de vérifier les options spécifiées lors de l'exécution du programme et d'en récupérer les paramètres
 * \param key le caractère permettant d'identifier l'option
 * \param strg_out la valeur utilisée pour l'option
 * \param argc le nombre de paramètres utilisés lors de l'appel au programme
 * \param argv un pointeur vers le tableau des paramètres du programme
 * \param opt indique si l'option est obligatoire, 1 si oui, 0 sinon
 * \return 1 si l'option a été spécifiée, 0 si l'option est facultative et non utilisée
 */
int getoption(char key, char *strg_out, int argc, char *argv[], int opt)
{
	int i;

	strcpy(PROG_NAME,argv[0]);

	if ((argc%2)==0) SyntaxError(); 

	for (i=1; i<argc; i=i+2)
	{
		if (argv[i][0]!='-') SyntaxError();
		else
		{
			if (argv[i][1] == key)
			{
				strcpy(strg_out,argv[i+1]); 
				return 1;
			}
		}
	}
	if (opt==1) 
	{
		SyntaxError();
	}
	return 0;
}


/**
 * Méthode de vérification des options
 * \param argc le nombre de paramètres utilisés lors de l'appel au programme
 * \param argv un pointeur vers le tableau des paramètres du programme
 */
void checkoption(int argc, char *argv[])
{
	int i,j,k;
	char msg[128];
	for (i=1; i<argc; i=i+2)
	{
		k=0;
		for (j=0; j<NB_KEY1; j++)
		{
			if (KEY1[j]==argv[i][1]) k=1;
		}
		for (j=0; j<NB_KEY0; j++)
		{
			if (KEY0[j]==argv[i][1]) k=1;
		}
		if (k==0)
		{
			sprintf (msg,"option -%c unknown",argv[i][1]); 
			ExitError(msg); 
		}
	}
}
