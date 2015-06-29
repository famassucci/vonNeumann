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

#include "locked_r.h"

/* get the number of locked reactions from the optional string -L */
int get_n_locked (char *locked){
    
    int n_l = 0, locked_length = strlen( locked );
    char *l_index = locked;
    
    /* locked reactions are supplied as ''' index: value ''' */
    /* the number of such reaction is retrieved counting the number of ":" in the optional -L argument */
    while ( (int) (l_index - locked) < locked_length && strchr(l_index, ':') != NULL )  {
        
        n_l++;
        
        l_index = strchr(l_index, ':') + 1;
        
    }
    
    return n_l;
    
}

/* A function to fix the locked reactions by reading the optional argument -L */
int fix_locked (int n_locked, double *s, double **s_locked, double *lock_v, char *locked){
    
    char *l_index1 = locked, *l_index2, *l_index3, *copy1, *copy2;
    
    int locked_string_length = strlen( locked ), which, n_zeros=0;
    
    double l_value, **dummy_s = s_locked;
    
    
    /* if there is only one locked reaction */
    if ( n_locked == 1 ){
        
        /* cut the string at ":" separating the index of reaction to lock from the value to assign */
        l_index1 = strchr(locked, ':');
        
        /* copy the index bit of the string to a new string */
        copy1 = get_substring (locked, l_index1);
        
        /* get the index from the substring */
        which = atoll (copy1) - 1;
        
        /* get the value to assign from the rest of the string */
        l_value = atof (l_index1 + 1);
        
        /* if value = 0 flag that there is a null reaction*/
        if (l_value == 0.) n_zeros += 1;
        
        /* assign the pointer to the fixed reaction */
        *dummy_s = s + which;
        
        /* assign the value to the array of fixed values lock_v*/
        *lock_v = l_value;
        
        /* free the substring */
        free ( copy1 );
        
    }
    
    /* otherwise */
    else if ( n_locked > 1 ){
        
        /* for each fixed reaction (found through the ":")*/
        while ( (int) (l_index1 - locked) < locked_string_length + 1 && strchr(l_index1, ':') != NULL )  {
            
            /* find the substring containing the index of the reaction */
            l_index2 = strchr(l_index1, ':');
            
            /* find the substring containing the value to assign */
            l_index3 = strchr(l_index1, ',');
            
            /* if it is the last reaction there are no more "," - the value substring lasts till the end */
            if ( l_index3 == NULL ) l_index3 = locked + locked_string_length;
            
            /* copy the index substring to a new string */
            copy1 = get_substring (l_index1, l_index2);
            
            /* get the numerical index of the index */
            which = atoi (copy1) - 1;
            
            /* copy the value substring to a new string */
            copy2 = get_substring (l_index2 + 1 , l_index3);
            
            /* get the numerical value of the flux value */
            l_value = atof (copy2);
            
            /* if value = 0 flag that there is a null reaction*/
            if (l_value == 0.) n_zeros += 1;
            
            /* assign a pointer to the locked reaction */
            *dummy_s = s + which;
            
            /* assign the value to the array */
            *lock_v = l_value;
            
            /* free the first substring */
            free ( copy1 );
            
            /* free the second substring */
            free ( copy2 );
            
            /* get to next pointer to locked reactions */
            dummy_s++;
            
            /* get to next array element of locked values*/
            lock_v++;
            
            /* iterate over the whole string up to the next "," */
            l_index1 = l_index3 + 1;
            
        }
        
    }
    
    /* return the number of zero reactions*/
    /* useful to check the problem is solvable, i.e. no metabolite is consumed only */
    return n_zeros;
    
}

/* A function to assign pointers to the null reactions (i.e. reactions locked to 0) */
/* useful because these reactions can be removed from the system */
double  **assign_null_reactions (double **s_locked, double *lock_v, int n_null, int n_locked){
    
    double **s_null = (double **) malloc( n_null * sizeof (double *) ), **dummy_s, **dummy_s_null, *dummy_v;
    
    dummy_s = s_locked;
    
    dummy_s_null = s_null;
    
    /* loop over array of locked value */
    for (dummy_v = lock_v; dummy_v < lock_v + n_locked; dummy_v++) {
        
        /* if a reaction is locked to 0 then */
        if ( *dummy_v == 0. ){
            
            /* store the associated pointer into array s_null */
            *dummy_s_null = *dummy_s;
            
            /* iterate over the pointer to null reactions */
            dummy_s_null++;
        }
        
        /* iterate over the pointer to locked reactions */
        dummy_s++;
    }
    
    /* return the array of pointers to null reactions */
    return s_null;
}

/* A function to add reactions that may eventually be forced to zero to render the system solvable */
int update_null_reactions (double **s_locked, double *lock_v, int n_null, int n_null_final, int n_locked, double **s_null) {
    
    double **dummy_s, **dummy_s2, *dummy_v;
    
    dummy_s2 = s_locked + n_locked;
    
    dummy_v = lock_v + n_locked;
    
    /* loop over the array of pointers to null reactions, from the last "known" value */
    for (dummy_s = s_null + n_null; dummy_s < s_null + n_null_final; dummy_s++){
        
        /* assign the new pointer */
        *dummy_s2 = *dummy_s;
        
        /* assign the zero value to the value array */
        *dummy_v = 0.;
        
        /* iterate over the array of locked reactions */
        dummy_s2++;
        
        /* iterate over the array of values */
        dummy_v++;
        
    }
    
    /* return the final number of locked reactions */
    return n_locked - n_null + n_null_final;
    
}
