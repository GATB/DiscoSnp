/* This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software.  You can  use, 
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info". 
 */

/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file gassst.cpp
 * \brief Module principal du programme
 * \author Dominique Lavenier
 * \author Damien Fleury
 * \author Guillaume Rizk
 * \version 6.1
 * \date 19/12/2008
 */

long long  NB_DIFF_SEED;
int SIZE_SEED;

#include <pthread.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#ifndef OSX
#include <sys/sysinfo.h>
#endif
#include <sys/time.h>



#include "constants.h"
#include "code.h"
#include "Pool.h"
#include "filter.h"
#include "Bank.h"
#include "Index.h"
#include "Seed.h"
#include "Hit.h"
#include "misc.h"
#include "display.h"
#include "gapless.h"
#include "Server.h"
#include "withgap.h"
#include "Stat.h"
#include "Doublon.h"
//#include <valgrind/callgrind.h>

#ifdef T_CLOCK
clock_t t1,t2,t3,t4;
clock_t tindex=0,tcalc=0;
#endif


#ifdef CPT_FILT

long long int cpt_hit=0;
long long int cpt_pref=0;
long long int cpt_ntsse=0;
long long int cpt_fil=0;
long long int cpt_doub=0;
long long int cpt_max=0;
long long int cpt_dint=0;
long long int cpt_good=0;

#endif
//int  alok=0;
//int  aldoublon=0;


using namespace std;

/**
 * Structure définissant les arguments de la fonction appelée par un thread
 */
typedef struct
{
	FILE * file;
	Bank * B1;
	Bank * B2;
	Index * I1;
        Stat * St;
        Doublon * Doub; // chaque thread a le sien
	int percent_ident;
	bool type_align;
	Server * serveur;
	int num_gaps;
        char ** tabprec;
        char ** tabnt;
        char ** tabprec7;
        char ** tabprec_gless;

  int thread_id;
} thread_args;

/**
 * Méthode permettant d'afficher le temps d'exécution
 * \param temps le nombre total de secondes correspondant au temps d'exécution
 */
void afficherTemps(int temps)
{
	int secondes, minutes, heures;
	secondes = temps%60;
	minutes = temps/60;
	heures = minutes/60;
	minutes = minutes%60;
	
	cout << "Execution time = ";

	if(heures>0)
	{
		cout << heures << " hour(s) ";
	}
	if(minutes>0)
	{
		cout << minutes << " minute(s) ";
	}
	cout << secondes << " second(s)";
	cout << endl << endl;
}


/**
 * Méthode associée à un thread
 * \param un pointeur vers une structure d'arguments thread_args
 * \return un pointeur void correspondant au nombre d'alignements trouvés.
 */
void * f_thread(void * arguments)
{
	task* t = new task();
	thread_args *arg = (thread_args*) arguments;
	long nb_hit = 0;
	
	int iter = arg->serveur->give_task(t);

// 	if(arg->num_gaps == 0)
// 	{
// 		/// Cas d'un alignement sans gap
// 		while(iter)
// 		  {
// 		    /// On traite chaque graine
// 		    //	for (num_seed=(t->start); num_seed<=(t->end); ++num_seed)
// 		    //{
// 		    /// On recherche des alignements en utilisant le module gapless
// 		    //	nb_hit = nb_hit + gapless_align(arg->file,arg->B1,arg->B2,arg->I1,arg->I2,num_seed,
// 		    //					arg->percent_ident,arg->type_align,arg->St);
// 		    //
// 		    nb_hit = nb_hit + gapless_align(arg->file,arg->B1,arg->B2,arg->I1,
// 						    arg->percent_ident,arg->type_align,arg->num_gaps,arg->tabprec_gless,arg->tabnt,arg->tabprec7,arg->St,arg->Doub,t->start,t->end);



// 		    // }
// 			iter = arg->serveur->give_task(t);
// 		}
// 	}
// 	else
	{

	  //	  CALLGRIND_START_INSTRUMENTATION

		/// Cas d'un alignement avec gap(s)
		while(iter)
		{
#ifdef T_CLOCK
		  t3=clock();
#endif

		    //  printf("thread %i start %i end %i \n",arg->thread_id,t->start,t->end);
		  if(NUMMAX && (NUM_ALIGN_COURANT >= NUMMAX)) return (void *) nb_hit;

		    nb_hit = nb_hit + withgap_align(arg->file,arg->B1,arg->B2,arg->I1,
						    arg->percent_ident,arg->type_align,arg->num_gaps,arg->tabprec,arg->tabnt,arg->tabprec7,arg->St,arg->Doub,t->start,t->end);

// 			/// On traite chaque graine
// 			for (num_seed=(t->start); num_seed<=(t->end); ++num_seed)
// 			{
// 				/// On recherche des alignements en utilisant le module withgaps
// 				nb_hit = nb_hit + withgap_align(arg->file,arg->B1,arg->B2,arg->I1,arg->I2,num_seed,
// 								arg->percent_ident,arg->type_align,arg->num_gaps,arg->tabprec,arg->tabnt,arg->St);
// 			}
			iter = arg->serveur->give_task(t);
#ifdef T_CLOCK
		    t4=clock();
		    tcalc+=t4-t3;

#endif
		    //    printf("trhead %i  deb %i fin %i\n",arg->thread_id,t->start,t->end);

		}
		//	CALLGRIND_STOP_INSTRUMENTATION
	}
	delete t;
	return (void *) nb_hit;
}


/**
 * Méthode permettant de calculer les alignements pour une partition de la première banque donnée.
 * On va donc parcourir toutes les partitions de la seconde banque
 * \param B1 le pointeur de la première banque de séquences
 * \param B2 le pointeur de la seconde banque de séquences
 * \param I1 le pointeur de l'index de la première banque de séquences
 * \param file le fichier de sortie, où seront affichés les alignements obtenus
 * \param dust_switch l'état d'activation du filtre Low complexity
 * \param percent_ident le pourcentage de ressemblances minimal entre les 2 séquences pour qu'un alignement soit conserver
 * \param type_align le type des alignements, soit dans la séquence normale ou dans la séquence inverse complémentée
 * \param num_gaps le nombre maximal de gaps autorisés dans les alignements
 * \return un entier long corresspondant au nombre d'alignements trouvés
 */
long align_loop(Bank * B1, Bank * B2, Index * I1, FILE * file, bool dust_switch, int percent_ident, bool type_align, int num_gaps, char** tabprec, char** tabnt,
		char** tabprec7, char** tabprec_gless, Stat ** St, Doublon ** Doub)
{
  //// Server * serveur = new Server(NB_DIFF_SEED, SIZE_PARTITION_THREAD); // achanger nquery, nquery par thread
	
	/// Le tableau des threads
	//new
	pthread_t *tab_threads;
	tab_threads = new pthread_t [NBTHREADS];

	int num_thread;
	long nb_hit = 0;
	void * tmp = 0;
	
	thread_args * tab_arg = new thread_args[NBTHREADS];
	
	//thread_args * tab_arg = malloc (sizeof(thread_args)*NBTHREADS);


	/// Initialisation du verrou d'accès au fichier de sortie
	pthread_rwlock_init(&verrou_out,NULL);
	


	while (B2->readBank(dust_switch)) //normalement un seul passage pour LQ
	{

	  
	        Server * serveur = new Server(B2->nb_seq[B2->num_part], SIZE_PARTITION_THREAD);

		/// Construction de la structure thread_args contenant les arguments des threads

		//meme objet stat pour tous,
		//mais objet doub differents
		//changement <- objet stat indiviuel, sans tab car LQ <- a faire
		for(num_thread=0; num_thread<NBTHREADS; ++num_thread)
		  {
		    tab_arg[num_thread].file = file;
		    tab_arg[num_thread].B1 = B1;
		    tab_arg[num_thread].B2 = B2;
		    tab_arg[num_thread].I1 = I1;
		    tab_arg[num_thread].St = St[num_thread];
		    tab_arg[num_thread].percent_ident = percent_ident;
		    tab_arg[num_thread].type_align = type_align;
		    tab_arg[num_thread].serveur = serveur;
		    tab_arg[num_thread].num_gaps = num_gaps;
		    tab_arg[num_thread].tabprec= tabprec;
		    tab_arg[num_thread].tabprec7= tabprec7;
		    tab_arg[num_thread].tabprec_gless= tabprec_gless;
		    tab_arg[num_thread].tabnt= tabnt;
		    tab_arg[num_thread].Doub = Doub[num_thread];
		    tab_arg[num_thread].thread_id=num_thread;
		  }
		cout << "Partition bank: " << B1->num_part <<", partition request: " << B2->num_part << endl;
		/// Initialisation de la seconde banque
		
		////
		pthread_attr_t attr;

		pthread_attr_init(&attr);
#ifdef DEBUG
		size_t stacksize;
		pthread_attr_getstacksize (&attr, &stacksize);
		printf("Default stack size = %li\n", stacksize);
#endif
		pthread_attr_setstacksize (&attr, 10485760);
#ifdef DEBUG
		pthread_attr_getstacksize (&attr, &stacksize);
		printf("new stack size = %li\n", stacksize);
#endif


		/// Calcul de l'alignement
		/// Création des threads
		for(num_thread=0; num_thread<NBTHREADS; ++num_thread)
		{
			pthread_create(&tab_threads[num_thread],&attr,f_thread, &tab_arg[num_thread]);
		}

		/// Récupération des résultats des threads
		for(num_thread=0; num_thread<NBTHREADS; ++num_thread)
		{
			pthread_join(tab_threads[num_thread], &tmp);
			//nb_hit += (long) f_thread(&arg);
			nb_hit += (long)tmp;
		}


		/// On réinitialise l'attribut tailleMaxSeq
		B2->tailleMaxSeq = 0;

		delete serveur;
	}

	delete [] tab_threads;
	delete [] tab_arg;

	return nb_hit;
}


/**
 * Méthode responsable du calcul des alignements selon certains paramètres
 * \param B1 le pointeur de la première banque de séquences
 * \param B2 le pointeur de la seconde banque de séquences
 * \param I1 le pointeur de l'index de la première banque de séquences
 * \param file le fichier de sortie, où seront affichés les alignements obtenus
 * \param fileName le nom du fichier de sortie
 * \param rev_comp_state un booléen qui indique si on utilise l'option reverse complement
 * \param dust_switch l'état d'activation du filtre Low complexity
 * \param percent_ident le pourcentage de ressemblances minimal entre les 2 séquences pour qu'un alignement soit conservé
 * \param num_gaps le nombre maximal de gaps autorisés dans les alignements
 */
void alignment_search(Bank * B1, Bank * B2, Index * I1, FILE * file, char *fileName, bool rev_comp_state, bool dust_switch, int percent_ident, int num_gaps, Stat ** St)
{
  long nb_hit = 0 ;

   //Objet pour eviter les doublons
  Doublon * Doub[NBTHREADS];
  /// Lecture et indexation des banques par morceaux

 
  //precalcul de la table d alignement
  long int dimt = (long int) pow(4,TAI_TAB);
  long int dimbigt = (long int) pow(4,TAI_BIG_TAB);
  int i;
  char ** table_prec;
  char ** table_prec_gless;

  char ** table_prec7;
  char ** table_nt;

  unsigned int cs1,cs2;
  //allocation de la table
  
  table_prec = (char **) malloc(dimt*sizeof(char *));
  table_prec_gless = (char **) malloc(dimt*sizeof(char *));
  for (i=0; i< dimt; i++)
    {
      table_prec[i] = (char *) malloc(dimt*sizeof(char));
      table_prec_gless[i] = (char *) malloc(dimt*sizeof(char));
    }



table_prec7 = (char **) malloc(dimbigt*sizeof(char *));
  for (i=0; i< dimbigt; i++)
    {
      table_prec7[i] = (char *) malloc(dimbigt*sizeof(char));
    }

  cs1=0;cs2=0;
  //alloc table nt
  long int dim_tabnt= (long int) pow(4,8);

  table_nt = (char **) malloc(4*sizeof(char *));
  for (i=0; i< 4; i++)
    {
      table_nt[i] = (char *) malloc(dim_tabnt*sizeof(char));
    }


  //puis calcul
  compute_table_gless(table_prec_gless,TAI_TAB);
  compute_table(table_prec,TAI_TAB);
  compute_table(table_prec7,TAI_BIG_TAB);
  compute_table_nt(table_nt,8); // table precalc des occur des nt dans seq de 8 bases
  

  B1->resetBank();
  while (B1->readBank(dust_switch))
    {
      /// Indexation de la première banque
  


#ifdef T_CLOCK
      t1=clock();
#endif
      I1->indexBank(B1,INDEX_STRIDE);
      //  I1->printIndex();
#ifdef T_CLOCK
      
      t2=clock();
      tindex+=t2-t1;
#endif

      for (i=0; i< NBTHREADS; i++)
	{
	     Doub[i] = new Doublon(B1->tailleMaxSeq); // sur 1 seq seulement
	  //  Doub[i] = new Doublon(B1->nb_tot_res); // sur tout genome
	}
      B2->resetBank();
		
	
      if(!rev_comp_state)
	cout << endl << endl << "    [[ Normal research ]]" << endl;
      else
	cout << endl << endl << "    [[ Normal + Reverse complement research ]]" << endl;
  

	
      nb_hit += align_loop(B1,B2,I1,file,dust_switch,percent_ident,rev_comp_state,num_gaps,table_prec,table_nt,table_prec7,table_prec_gless,St,Doub);
   
      /////   delete [] I1->seed;
      /////if(SIZE_SEED>14)  I1->free_hashtable();
      
      /// On fait la recherche "reverse complement" si elle est activée 

      //faite en faisant reverse de query      
/*
      if(rev_comp_state)
	{
	  cout << "    [[ Reverse complement research ]]" << endl;
	

		
	  //	  B1->reverseComplement();
	  /// B2->reverseComplement();
#ifdef T_CLOCK
	  
	  t1=clock();
#endif
	  
	  /////	  I1->indexBank(B1,INDEX_STRIDE);
#ifdef T_CLOCK
	  
	  t2=clock();
	  tindex+=t2-t1;
#endif

	  B2->resetBank();
	  nb_hit2 += align_loop(B1,B2,I1,file,dust_switch,percent_ident,REV_COMP_ALIGN,num_gaps,table_prec,table_nt,table_prec7,table_prec_gless,St,Doub);

	 ///////// delete [] I1->seed;

	  /////////// if(SIZE_SEED>14)  I1->free_hashtable();


	}
*/
      
        delete [] I1->seed;
	if(SIZE_SEED>14)  I1->free_hashtable();

  for (i=0; i< NBTHREADS; i++)
     {
      delete Doub[i] ;
    }
     // delete Doub;
    }


  // if(T_CLOCK){
#ifdef T_CLOCK

  printf("-------------------------------------------------------\n");
  printf("----------------------- BILAN -------------------------\n");
  printf("-------------------------------------------------------\n");

  printf("Indexation  : %.4lf s \n", (tindex)/(double)CLOCKS_PER_SEC);
  printf("Calcul  : %.4lf s \n \n", (tcalc)/(double)CLOCKS_PER_SEC);
  //  printf("ALok %i ALdoublon %i \n \n", alok,aldoublon);
#endif

#ifdef CPT_FILT

    printf("HIT    :  %12lld\n",cpt_hit);
    printf("PREF   :  %12lld\n",cpt_pref);
    printf("NTSSE  :  %12lld\n",cpt_ntsse);
    printf("DOUB   :  %12lld\n",cpt_doub);
    printf("FILTRE :  %12lld\n",cpt_fil);
    printf("MAX    :  %12lld\n",cpt_max);
    printf("DC     :  %12lld\n",cpt_dint);

    printf("\nGood enough     :  %12lld\n",cpt_good);

#endif

 


  /// On affiche les résultats trouvés
  if(!rev_comp_state)
     cout << endl << "Forward strand search only:  " << nb_hit << " Alignments found --> " << fileName << endl;
  else
    cout << endl << "Forward and Reverse search :  " << nb_hit << " Alignments found --> " << fileName << endl;
  
  cout << endl;



  for (i=0; i< dimt; i++) 	free(table_prec[i]);
  free(table_prec);



  for (i=0; i< dimt; i++) 	free(table_prec_gless[i]);
  free(table_prec_gless);



  for (i=0; i< dimbigt; i++) 	free(table_prec7[i]);
  free(table_prec7);

  for (i=0; i< 4; i++) 	free(table_nt[i]);
  free(table_nt);

}



/**
 * Méthode main du programme
 */
int main(int argc,char * argv[]){
	
  /// Définition des banques de séquences
  Bank * B1;
  Bank * B2;
  /// Définition des index associés aux banques
  Index * I1;

 
  /// Taille des partitions
  long long  size_partition;

  /// Etat du filtre Low Complexity
  int dust_switch = DUST_OFF;
  
  /// Etat d'activation de la recherche reverse complement
  int Rev_switch = REV_ON;
  
  /// Nombre de gaps
  int num_gaps = DEFAULT_GAP_ALLOWED;

  /// Fichiers de sorties
  FILE * fileRES;
  
  int percent_ident;
  char fileNameBANK[1024];
  char fileNameRQT[1024];
  char fileNameRES[1024];
  char tmp[1024];
  int i;
  int Tpart = 0;
  long long totalram;
#ifndef OSX
  struct sysinfo s_info;
  int error;
#endif

#ifdef EXEC_TIME
  /// Variables de calcul du temps
  int temps_debut, temps_fin;
  struct timeval tv;
  int temps_execution;

  /// mémorisation du temps de départ
  gettimeofday(&tv, 0);
  temps_debut = tv.tv_sec;
#endif

  /// Initialisations
  NB_KEY1 = 0;
  NB_KEY0 = 0;
  // SIZE_SEED = DEFAULT_SIZE_SEED;
  SIZE_SEED = 0;
  NBTHREADS= DEFAULT_NBTHREADS;
  MAXHITS = DEFAULT_MAXHITS;
  MAXPOS = DEFAULT_MAXPOS;
  BESTAL = DEFAULT_BESTAL;
  SLEVEL = DEFAULT_SLEVEL;
  //GLOBAL_GLOBAL =
  BITSTAT = 22;
  NUMGAPS_AUTO = 1;
  NUMMAX = 0;
  NUM_ALIGN_COURANT = 0 ; 
  /// Gestion des options du programme

  /// Définition des options
  option('d',"bank_file",1);
  option('i',"query_file",1);
  option('o',"output_file",1);
  option('p',"identity_percentage",1);
  option('w',"size_seed",0);
  option('m',"output_format",0);
  option('l',"complexity_filter",0);
  option('t',"size_partition",0);
  option('r',"reverse_complement",0);
  option('e',"extended cigar",0);

  option('g',"gaps_percentage",0);
  option('n',"thread_number",0);
  option('h',"max_hits_per_read",0);
  // option('s',"max_hits_per_seed",0);
  option('s',"sensitivity level",0);
  option('b',"output best alignments",0);

  option('z',"obscure option : bitstat",0);

  option('y',"option rayan",0);



  /// option taille du tableau des doublons
  if (getoption('y',tmp,argc,argv,0)==1)
  {
  	sscanf(tmp,"%d",&NUMMAX);
  }


    extended_cigar= 0;
    
    if (getoption('e',tmp,argc,argv,0)==1)
    {
        sscanf(tmp,"%d",&extended_cigar);
    }



  /// option taille du tableau des doublons
  if (getoption('z',tmp,argc,argv,0)==1)
  {
  	sscanf(tmp,"%d",&BITSTAT);
  }



  /// Option filtre Low Complexity
  if (getoption('b',tmp,argc,argv,0)==1)
  {
  	sscanf(tmp,"%d",&BESTAL);
  	if(BESTAL !=0 && BESTAL !=1) ExitError("best alignments must be 0 or 1");
  }


  /// Récupération du nom du fichier de sortie
  getoption('o',fileNameRES,argc,argv,1);
  if ((fileRES=fopen(fileNameRES,"w"))==NULL) ExitError("cannot open output file");

  /// Option pourcentage d'identité
  getoption('p',tmp,argc,argv,1);
  sscanf(tmp,"%d",&percent_ident);
  if ((percent_ident<0)||(percent_ident>100)) ExitError("wrong identity percentage, identity persentage must be in [0 ; 100]"); 

  /// Option format de sortie
   OUTPUT_FORMAT = SAM_READY_FORMAT;
  //OUTPUT_FORMAT =  STD_OUTPUT_FORMAT; 
  if (getoption('m',tmp,argc,argv,0)==1)
  {
	  sscanf(tmp,"%d",&OUTPUT_FORMAT);
	  checkOutputFormat(OUTPUT_FORMAT);
  }

    


  /// Option taille des graines
  if (getoption('w',tmp,argc,argv,0)==1)
  {
	  sscanf(tmp,"%d",&SIZE_SEED);
	  if ((SIZE_SEED<6)||(SIZE_SEED>20)) ExitError("size seed must be in [6 ; 20]");
  } 

  /// Option nombre de threads
  if (getoption('n',tmp,argc,argv,0))
  {
	  sscanf(tmp,"%d",&NBTHREADS);
	  if (NBTHREADS<1) ExitError("number of threads must be 1 minimum");
  }

  /// Option nombre max d align par read
  if (getoption('h',tmp,argc,argv,0))
  {
	  sscanf(tmp,"%d",&MAXHITS);
	  if (MAXHITS<0) ExitError("parameter h number of alignments per query sequence must be postive, (0 = no limit)");
  }

  if (MAXHITS==0) BESTAL=0; // si on les renvoit tous , best na plsu de sens, desactive pour enlever gestion liste
/// Option nombre max d align par read
  if (getoption('s',tmp,argc,argv,0))
  {
    //	  sscanf(tmp,"%d",&MAXPOS);
    //	  if (MAXPOS<0) ExitError("parameter h number of alignments per query sequence must be postive, (0 = no limit)");
    	  sscanf(tmp,"%d",&SLEVEL);
	  if (SLEVEL<0 || SLEVEL > 5) ExitError("parameter s sensibility level must be in [0-5] ,0 is fastest, 5 is best sensibility (default = 2)");

  }


  /// Option filtre Low Complexity
  if (getoption('l',tmp,argc,argv,0)==1)
  {
  	sscanf(tmp,"%d",&dust_switch);
  	if(dust_switch !=0 && dust_switch !=1) ExitError("complexity_filter must be 0 or 1");
  }

  /// Option taille de la partition // en Mo,  un octet=un nt (ou un char de nom seq)
  if (getoption('t',tmp,argc,argv,0)==1)
  {
	  Tpart = 1;
	  sscanf(tmp,"%llu",&size_partition);
	  if ((size_partition<1)||(size_partition>1600)) ExitError("size partition must be in [1 ; 1 600]");
	  /// On transforme la taille en octets
	  size_partition = size_partition * 1048576LL;
  }

  /// Si la taille des partitions n'est pas définie, on la définit à 2% de la taille de la mémoire principale (RAM)
  // car l'index de la banque prend beaucoup plus de place que la banque (environ 16 fois plus)
#ifndef OSX

  error = sysinfo(&s_info);
  if(error) ExitError("Memory Error");
  totalram = s_info.totalram * s_info.mem_unit;
#else
    totalram =  2048* 1048576LL;  //20148
#endif
  if(!Tpart)
    {
      // premiere valeur par defaut, corrigee ensuite si trop gros
      size_partition = (2 * totalram) /100LL;
      size_partition =  size_partition <  (1600LL*1048576LL) ?      size_partition :   (1600LL*1048576LL);
      // printf("Total memory : %lli Mo   Size partition auto : %lli Mo  \n",totalram/1048576LL,size_partition/1048576LL);
    }
 
  /// Option reverse complement, on récupère le nom du fichier de sortie
  if (getoption('r',tmp,argc,argv,0)==1)
  {
  	sscanf(tmp,"%d",&Rev_switch);
  	if(Rev_switch !=0 && Rev_switch !=1) ExitError("reverse_complement must be 0 or 1");
  }

    
    
  /// Option nombre de gaps
  if (getoption('g',tmp,argc,argv,0)==1)
    {
      sscanf(tmp,"%d",&num_gaps);

      if ((num_gaps<0)||(num_gaps>100)) ExitError("wrong gap percentage, gap percentage must be in [0 ; 100]"); 

      if ((100-percent_ident) < num_gaps )
	{
	  num_gaps = 100-percent_ident;
	  printf("\nWarning! max gap percentage adjusted to %i %%   because identity percentage is set to %i %%  \n \n",num_gaps,percent_ident);
	}
      //	if(num_gaps <0 ) ExitError("number of gaps must be >= 0 ");
      NUMGAPS_AUTO=0;
    }



  /// Définition de l'objet contenant statistiques des align trouvés
  Stat * St[NBTHREADS]; //cela marche til ?

  /// Récupération des noms des fichiers des banques
  getoption('d',fileNameBANK,argc,argv,1);
  getoption('i',fileNameRQT,argc,argv,1);
  checkoption(argc,argv);

  /// Fin de gestion des options


  printf("Running Gassst version %s\n",GASSST_VERSION);




  B2 = new Bank(fileNameRQT,2000000000,fileRES,0); //nolimit car pas dindex de la banque 2 (en octets)


  if(SIZE_SEED==0)
    {
      int qlen= B2->tSeq; int nerr= (qlen*(100-percent_ident)) / 100;
      //peut etre changer seuil de maniere coordonnee a option s
      int w_estimate = find_seed_len(qlen,nerr,100,7);
      w_estimate = min(w_estimate,20); w_estimate = max(w_estimate,12);

      cout << "Seed length not set, automatically assigning value according to read length and error rate  ("<< nerr<<"/"<<qlen<<")  . . . .   ---->     " << w_estimate << " \n";;

      SIZE_SEED = w_estimate;
    }
  /// Calcul du nombre de graines différentes
  NB_DIFF_SEED = 1;
  for (i=0; i<SIZE_SEED; ++i) NB_DIFF_SEED = NB_DIFF_SEED * 4;


 long long isize=0;
 if (SIZE_SEED>14) 
   isize = size_partition *(sizeof(cell) +sizeof(Seed)) + NHACH * sizeof(cell **);
 else 
   isize = size_partition * sizeof(Seed) + (NB_DIFF_SEED*sizeof(int)*2);


 if(isize>totalram)
   {
     if(Tpart)
       {
	 cout << "This computer has " << totalram/1048576LL  << " Mo memory, \n";
	 cout << "but the specified partition size  (" << size_partition/1048576LL << ") requires "<< isize/1048576LL << " Mo memory for its "<< SIZE_SEED << "-nt-long seed index. \n";
       }

     size_partition = (2 * totalram) /100;
     size_partition =  size_partition <  (1600LL*1048576LL) ?      size_partition :   (1600LL*1048576LL);
     if (SIZE_SEED>14) 
       {

	 size_partition = ((long long ) (totalram * 0.95f ) - NHACH * sizeof(cell **))  
	   / (sizeof(cell) +sizeof(Seed));

	 isize = size_partition *(sizeof(cell) +sizeof(Seed)) + NHACH * sizeof(cell **);

       }
     else 
       {
	 size_partition = ((long long ) (totalram * 0.95f ) -  (NB_DIFF_SEED*sizeof(int)*2)) / sizeof(Seed) ;

	 isize = size_partition * sizeof(Seed) + (NB_DIFF_SEED*sizeof(int)*2);
       }
     
     size_partition =  size_partition <  (1600LL*1048576LL) ?      size_partition :   (1600LL*1048576LL);

     if(Tpart)
	 cout << "Resuming with partition size adjusted to  " << size_partition/1048576LL << " Mo, requiring " << isize/1048576LL << " Mo memory.\n";
     else
	 cout << "Maximum estimated Index size = "<< isize/1048576LL <<" Mo,  for a partition size  automatically set to "<< size_partition/1048576LL << " Mo, and Seed length="<<SIZE_SEED << "nt \n";
   }
 else
   {
     if(Tpart)
       cout << "Maximum estimated Index size = "<< isize/1048576LL <<" Mo,  for a partition size t="<< size_partition/1048576LL << " Mo, and Seed length="<<SIZE_SEED << "nt \n";
     else
       cout << "Maximum estimated Index size = "<< isize/1048576LL <<" Mo,  for a partition size  automatically set to "<< size_partition/1048576LL << " Mo, and Seed length="<<SIZE_SEED << "nt \n";
   }

 /// On affiche la taille maximale de partitions
  cout << "\nMaximum partition size = "<< size_partition <<" bytes\n\n";


  /// I2 = new Index(); //devenu inutile

  /// Initialisation des banques
  B1 = new Bank(fileNameBANK,size_partition,fileRES,1);
  B1->writeInfoBank();
  B2->writeInfoBank();

  


  /// Initialisation des index
  //il faut connaitre taille graine ici
  I1 = new Index();

  for (i=0; i< NBTHREADS; i++)
    {
      St[i] = new Stat(B2->nb_tot_seq);
    }
  cout << endl;
  cout << B1->fileBank << " --> " << B1->nb_part << " partitions" << endl;
  cout << B2->fileBank << " --> " << B2->nb_part << " partitions" << endl;

  //MAXPOS= MARGE*  (B1->nb_tot_res/NB_DIFF_SEED); //*15 de marge 


  //calcul nombre moyen de hits par graine
  //le compteur maxpos est apres la plupart des filtres, donc ca fait une marge correspondant au taux de filtrage
  
  float maxpos_f =  ( min(size_partition,B1->nb_tot_res)   / ((float)NB_DIFF_SEED)); 
  //marge_maxpos_amont : ce  compteur est avant les filtres, on rajoute un facteur de marge
  int marge_maxpos_amont;
  if(SIZE_SEED<=16) marge_maxpos_amont = MARGE;  else  marge_maxpos_amont = MARGE_HAUTE;
  float maxpos_amont_f =   marge_maxpos_amont  * ( min(size_partition,B1->nb_tot_res)   /((float)NB_DIFF_SEED));
  switch (SLEVEL)
    {
    case  0 : 
      MAXPOS = max(1, (int) round(maxpos_f /6));
      MAXPOS_AMONT = max(1,(int) round(maxpos_amont_f /3));
      break;
    case  1 :
      MAXPOS = max(1,(int) round(maxpos_f /3));
      MAXPOS_AMONT =  max(1,(int) round(maxpos_amont_f *1));
      break;
    case  2 :
      MAXPOS = max(2,(int) round(maxpos_f * 1)) ;
      MAXPOS_AMONT = max(2,(int) round(maxpos_amont_f *2));
      break;
    case  3 :
      MAXPOS = max(6,(int) round(maxpos_f * 2)) ;
      MAXPOS_AMONT = max(6,(int) round(maxpos_amont_f *4));
      break;
    case  4 :
      MAXPOS = max(13,(int) round(maxpos_f * 4)) ;
      MAXPOS_AMONT = max(13,(int) round(maxpos_amont_f *8));
      break;
    case  5 :
      MAXPOS = 0;//0 means no limit
      MAXPOS_AMONT = 0;
      break; 
    }



  // printf("MAXPOS auto : %i %i %i\n",B1->nb_tot_res,NB_DIFF_SEED,B1->nb_tot_res/NB_DIFF_SEED);
  //: printf("MAXPOS auto : %i MAXPOS amont : %i  Slevel %i\n",MAXPOS,MAXPOS_AMONT,SLEVEL);

  alignment_search(B1, B2, I1, fileRES, fileNameRES, Rev_switch, dust_switch, percent_ident,num_gaps,St);
  
  fclose(fileRES);
  
#ifdef EXEC_TIME
  /// Calcul et affichage du temps d'execution
  gettimeofday(&tv, 0);
  temps_fin = tv.tv_sec;
  temps_execution = temps_fin - temps_debut;
  afficherTemps(temps_execution);
#endif

  /// Libération de la mémoire
  delete I1;
  delete B1;
  delete B2;
  for (i=0; i< NBTHREADS; i++)
     {
      delete St[i] ;
    }

  return 0;
 
}
