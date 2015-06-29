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

#include "remove_r.h"

/* A function to 'remove' reactions from a metabolite i/o list */
void remove_reaction(double **react_list, double **to_remove, int n_reacs, double *coeff_list){
    
    double *tmp, c_tmp;
    
    /* the reaction is not really removed, but pushed to the end of the array */
    tmp = *to_remove;
    
    /* the array location to_remove will now point to reaction pointed by the last guy in the array **react_list */
    *to_remove = *(react_list + n_reacs - 1);
    
    /* the last element of the array points to the reaction originally pointed by to_remove */
    *(react_list + n_reacs - 1) = tmp;
    
    /* stoichiometric coefficients must be swapped too */
    c_tmp = *(coeff_list + (int)( to_remove - react_list) );
    
    /* assign array element associated to to_remove to last array value */
    *(coeff_list + (int)( to_remove - react_list) ) = *(coeff_list + n_reacs - 1 );
    
    /* end swapping by associating to last array elemtnt the coefficient previously associated with to_remove*/
    *(coeff_list + n_reacs - 1 ) = c_tmp;
    
    
}
