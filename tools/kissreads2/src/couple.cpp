/*****************************************************************************
 *   discoSnp++: discovering polymorphism from raw unassembled NGS reads
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
 * couple.c
 *
 *  Created on: 25 oct. 2010
 *      Author: ppeterlo
 */
#include <couple.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <commons.h>
p_couple create_couple(uint64_t a,uint64_t b){
	p_couple cpl = (p_couple)malloc(sizeof(couple));
	test_alloc(cpl);
	cpl->a=a;
	cpl->b=b;
	return cpl;
}


void free_couple (const void * v_c){
	couple * c = (couple *) v_c;
	if(c != NULL){
		free(c);
	}
}
