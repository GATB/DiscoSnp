#ifndef SERVER_H
#define SERVER_H

/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file Server.h
 * \brief Class Server, responsable du partage des tâches entre les threads
 * \author Damien Fleury
 * \version 5.2
 * \date 28/08/2008
 */


#include <pthread.h>

/**
 * Structure définissant une portion de graines à traiter 
 */
typedef struct
{
	/**
	 * Le numéro de la première graine de la portion
	 */
	int start;
	
	/**
	 * Le numéro de la dernière graine de la portion
	 */
	int end;
} task;

/**
 * \class Server, Un objet permettant de répartir le travail de recherche d'alignements entre différents threads
 * \brief Un objet Server gère l'ensemble des graines à traiter lors de la recherche, les threads vont 
 * demander une tâche à l'objet Server qui va leur attribuer une partie de l'ensemble des graines.
 * Le Server va donc répartir le travail le plus équitablement possible entre les différents threads 
 */
class Server
{
private:
	/**
	 * Le numéro de la dernière graine à traiter
	 */
	int last_seed;
	
	/**
	 * La taille d'un ensemble de graines attribué lors d'une demande par un thread
	 */
	int size_task;
	
	/**
	 * Un indice permettant de répertorier les graines qui n'ont pas encore été attribuées à un thread
	 */
	int index;
	
	/**
	 * Le nombre de sous-ensembles de graines à attribuer aux threads
	 */
	int nb_part;
	
	/**
	 * Un verrou empêchant les accès concurentiels aux variables du Server par les threads
	 */
	pthread_rwlock_t verrou_serveur;

public:
	/**
	 * Le constructeur par défaut de Server
	 */
	Server();
	
	/**
	 * Le destructeur par défaut de Server
	 */
	virtual ~Server();
	
	/**
	 * Le constructeur par recopie de Server
	 * \param s un objet Server
	 */
	Server(const Server& s);
	
	/**
	 * L'opérateur d'affectation de Server
	 * \param s un objet Server
	 * \return l'objet server
	 */
	Server& operator=(const Server& s);
	
	/**
	 * Constructeur de Server
	 * `param nbSeeds un entier correspondant au nombre total de graines à traiter
	 * \param nb_partitions un entier correspondant au nombre de partitions à faire avec l'ensemble des graines
	 */
	Server(int nbSeeds, int nb_partitions);
	
	/**
	 * Méthode permettant d'obtenir un sous-ensemble de graines à traiter
	 * \param t un pointeur vers une structure task, correspondant à une partition de l'ensemble des graines
	 * \return 1 si il reste des graines à traiter, on a dans ce cas affecté en conséquence les champs de la structure task pointée par t, 0 si toutes les graines ont été attribuées, les champs de la structure pointée par t sont mis à -1
	 */
	int give_task(task* t);
	
	/**
	 * Méthode permettant de réinitialiser l'index du server
	 */
	void reset();
};

#endif
