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
 * coherent_algorithm.c
 *
 *  Created on: 16 sept. 2010
 *      Author: pierre peterlongo
 */

#include<coherence_algorithm.h>
#include<commons.h>
#include <string.h>
#include <stdlib.h>


/**
 *  | pwi (may be negative)
 *     --------------  fragment
 *  *************      read
 *
 *  Tests if the overlapping part between read and fragment do not have more than <code>subst_allowed</code> substitions
 *  returns 1 if true between read and fragment, 0 else
 */
char read_coherent_generic(const int pwi, const char * fragment, const char * read, const int subst_allowed){
    
//    printf("read coherent, subst_allowed = %d\n", subst_allowed);
	int substitution_seen=0; // number of seen substitutions for now

	int pos_on_read, pos_on_fragment; // where to start



	/*
	 *  | pwi (negative)
	 *     --------------  fragment
	 *  *************      read
	 *     | we start here
	 */
	if(pwi<0) {
		pos_on_fragment=0;
		pos_on_read=-pwi;
	}

	/*
	 *        | pwi (positive)
	 *     --------------  fragment
	 *        *************      read
	 *        | we start here
	 */
	else{
		pos_on_fragment=pwi;
		pos_on_read=0;

	}


	// walk the read and the fragment together, detecting substitutions.
	// stop if the number of substitution is too high
	while(fragment[pos_on_fragment]!='\0' && read[pos_on_read]!='\0'){
		//if(fragment[pos_on_fragment]!=read[pos_on_read]) && fragment[pos_on_fragment]!='*') {// one subsitution
		if(fragment[pos_on_fragment]!=read[pos_on_read] &&
           fragment[pos_on_fragment]!='*' &&
           fragment[pos_on_fragment]!='?' &&
           fragment[pos_on_fragment]!='N'){ // one subsitution
			substitution_seen++;
			if(substitution_seen>subst_allowed) break; // too much subsitutions
		}
		pos_on_fragment++;
		pos_on_read++;
	}
//    printf("%d subst seen\n", substitution_seen);
	if(substitution_seen<=subst_allowed) return 1;
	return 0;
}

/**
 *  | pwi (may be negative)
 *     --------------  fragment
 *  *************      read
 *
 *  Tests if the overlapping part between read and fragment do not have more than <code>subst_allowed</code> substitions
 * In case of SNPs, we need to avoid any substitution on the central fragment position (the one containing the SNP)
 * Thus in this function, we return 0 if any substitution occurs on this central position, whatever the number of substitution_seen
 *  returns 1 if true between read and fragment, 0 else
 */
char read_coherent_SNP(const int pwi, const char * fragment, const char * read, const int subst_allowed, const char * SNP_positions){
	int substitution_seen=0; // number of seen substitutions for now
    
	int pos_on_read, pos_on_fragment; // where to start
    
    
    
	/*
	 *  | pwi (negative)
	 *     --------------  fragment
	 *  *************      read
	 *     | we start here
	 */
	if(pwi<0) {
		pos_on_fragment=0;
		pos_on_read=-pwi;
	}
    
	/*
	 *        | pwi (positive)
	 *     --------------  fragment
	 *        *************      read
	 *        | we start here
	 */
	else{
		pos_on_fragment=pwi;
		pos_on_read=0;
	}
    
    //const int snp_pos = strlen(fragment)/2;
    char snp_pos = SNP_positions[0];
    int id_array_SNP_position=0;
    
	// walk the read and the fragment together, detecting substitutions.
	// stop if the number of substitution is too high
	while(fragment[pos_on_fragment]!='\0' && read[pos_on_read]!='\0'){
        if (pos_on_fragment>snp_pos)
        {
            id_array_SNP_position++;
            snp_pos = SNP_positions[id_array_SNP_position];
        }
		if (fragment[pos_on_fragment]!=read[pos_on_read] &&
            fragment[pos_on_fragment]!='*' &&
            fragment[pos_on_fragment]!='?' &&
            fragment[pos_on_fragment]!='N'){ // one subsitution
			substitution_seen++;
			if(substitution_seen>subst_allowed) return 0; // too much subsitutions
            if(pos_on_fragment==snp_pos) {
                return 0; // substition should not be on the snp
            }
		}
		pos_on_fragment++;
		pos_on_read++;
	}
	return 1;
}
