/* vonNeumann, a minOver sampler of solutions to a von Neumann problem
*
* Copyright (C) 2015 Francesco A. Massucci
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "substring.h"

/* A function to copy a string up to a certain index (included) to another string */
char *get_substring (char *full_string, char *cutoff){
    
    char *sub_copy, *dummy1 = full_string, *dummy2;
    
    /* allocate memory for the new string */
    sub_copy = (char *) malloc ( (int) (cutoff - full_string + 1) * sizeof(char) );
    
    dummy2 = sub_copy;
    
    /* loop over the original string up to the cutoff value (included) */
    while ( dummy1 < cutoff + 1) {
        
        *dummy2 = *dummy1;
        
        dummy1++;
        
        dummy2++;
    }
    
    /* the last element is set to null */
    *(sub_copy + (int) (cutoff - full_string) ) = '\0';
    
    return sub_copy;
    
}

char *get_substring_no_copy (char *full_string, char *cutoff){

    *cutoff = '\0';

    return full_string;
}
