#ifndef FILTER_H
#define FILTER_H
//Copyright inria / irisa (2013)
//
//
//raluca.uricaru@gmail.com
//pierre.peterlongo@inria.fr
//
//This software is a computer program whose purpose is to call SNPs from NGS reads.
//
//This software is governed by the CeCILL license under French law and
//abiding by the rules of distribution of free software.  You can  use,
//modify and/ or redistribute the software under the terms of the CeCILL
//license as circulated by CEA, CNRS and INRIA at the following URL
//"http://www.cecill.info".
//
//As a counterpart to the access to the source code and  rights to copy,
//modify and redistribute granted by the license, users are provided only
//with a limited warranty  and the software's author,  the holder of the
//economic rights,  and the successive licensors  have only  limited
//liability.
//
//In this respect, the user's attention is drawn to the risks associated
//with loading,  using,  modifying and/or developing or reproducing the
//software by the user in light of its specific status of free software,
//that may mean  that it is complicated to manipulate,  and  that  also
//therefore means  that it is reserved for developers  and  experienced
//professionals having in-depth computer knowledge. Users are therefore
//encouraged to load and test the software's suitability as regards their
//requirements in conditions enabling the security of their systems and/or
//data to be ensured and,  more generally, to use and operate it in the
//same conditions as regards security.
//
//The fact that you are presently reading this means that you have had
//knowledge of the CeCILL license and that you accept its terms.

/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file filter.h
 * \brief Module Filter, definit le filtre Low Complexity qui detecte les zones non pertinentes des sequences afin de ne pas les indexer
 * \author Dominique Lavenier
 * \author Damien Fleury
 * \version 5.2
 * \date 28/08/2008
 */

/// Fenêtre du filtre Low complexity
/// #define WS 12


using namespace std;


/**
 * Methode permettant d'appliquer le filtre Low Complexity
 * \param data le tableau des donnees de la banque de sequences
 * \param id l'index du debut de la sequence
 * \param lenseq la longueur de la sequence
 * \return 1 si l'operation s'est deroulee correctement
 */
int filterLowComplexity (char* data, int lenseq, int threshold);
int filterLowComplexity2Paths (char* seq1, char *seq2, int lenseq, int threshold);

#endif

