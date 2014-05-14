/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/
#ifndef FILTER_H
#define FILTER_H

/********************************************************************************/

/** Methode permettant d'appliquer le filtre Low Complexity
 * \param data le tableau des donnees de la banque de sequences
 * \param id l'index du debut de la sequence
 * \param lenseq la longueur de la sequence
 * \return 1 si l'operation s'est deroulee correctement
 */
int filterLowComplexity       (char* data, int lenseq, int threshold);
int filterLowComplexity2Paths (char* seq1, char *seq2, int lenseq, int threshold);

/********************************************************************************/

#endif

