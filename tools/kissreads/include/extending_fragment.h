/*****************************************************************************
 *   DiscoMore: discovering polymorphism from raw unassembled NGS reads
 *   A tool from the GATB (Genome Assembly Tool Box)
 *   Copyright (C) 2014  INRIA
 *   Authors: P.Peterlongo, E.Drezen
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
/*
 * extending_fragment.h
 *
 *  Created on: 15 sept. 2010
 *      Author: ppeterlo
 */

#ifndef EXTENDING_FRAGMENT_H_
#define EXTENDING_FRAGMENT_H_

#include <fragment_info.h>

// conserve a link between a fragment that has been extended and its owner (the initial fragment)
typedef struct{
	char * extension;
	p_fragment_info corresponding_fragment;
} extension_fragment, * p_extension_fragment;

#endif /* EXTENDING_FRAGMENT_H_ */
