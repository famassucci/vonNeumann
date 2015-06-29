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

#include "cascades.h"

/* A function to check whether metabolites are only consumed */
/* Also, null reactions are removed from the system */
int check_cascades (metabolite *metabs, int Nmet, double ***s_zeros, int n_zeros, double *s, FILE *outfile){
    
    metabolite *dummy_m;
    double **s_locked = *s_zeros, **dummy_s1, **dummy_s2;
    int n_zeros_new=0;
    
    for (dummy_s1 = s_locked; dummy_s1 < s_locked + n_zeros; dummy_s1++){
        
        fprintf(outfile, "Reaction %d is zero\n", (int)(*dummy_s1 - s) + 1);
    }
    
    /* loop over null reactions and remove them from the system */
    for (dummy_s1 = s_locked; dummy_s1 < s_locked + n_zeros; dummy_s1++){
        
        /* loop over metabolites to check whther they participate to null reaction dummy_s1 */
        for (dummy_m = metabs; dummy_m < metabs + Nmet; dummy_m++){
            
            /* initalise pointer to dummy_m input reactions array */
            dummy_s2 = dummy_m -> input.react;
            
            /* loop over input reactions for metabolite dummy_m to check whether dummy_s1 is present */
            while ( dummy_s2 < dummy_m -> input.react + dummy_m -> input.n_react && *dummy_s2 != *dummy_s1) dummy_s2++;
            
            /* if while loop halted within the range dummy_m -> input.n_react, then dummy_s1 is present */
            if ( (int) (dummy_s2 - dummy_m -> input.react ) < dummy_m -> input.n_react ) {
                
                fprintf(outfile, "Removing reaction %d from inputs of metabolite %d... \n", (int) (*dummy_s2 - s) + 1, (int) (dummy_m - metabs) + 1);
                
                /* remove reaction dummy_s1 from input reactions of metabolite dummy_m */
                remove_reaction(dummy_m -> input.react, dummy_s2, dummy_m -> input.n_react, dummy_m -> input.coeff);
                
                /* reduce the number of reactions in which dummy_m appears as an input */
                dummy_m -> input.n_react -= 1;
            }
            
            /* initalise pointer to dummy_m output reactions array */
            dummy_s2 = dummy_m -> output.react;
            
            /* loop over output reactions for metabolite dummy_m to check whether dummy_s1 is present */
            while ( dummy_s2 < dummy_m -> output.react + dummy_m -> output.n_react && *dummy_s2 != *dummy_s1) dummy_s2++;
            
            /* if while loop halted within the range dummy_m -> output.n_react, then dummy_s1 is present */
            if ( (int) (dummy_s2 - dummy_m -> output.react ) < dummy_m -> output.n_react ) {
                
                fprintf(outfile, "Removing reaction %d from outputs of metabolite %d... \n", (int) (*dummy_s2 - s) + 1, (int) (dummy_m - metabs) + 1);
                
                /* remove reaction dummy_s1 from output reactions of metabolite dummy_m */
                remove_reaction(dummy_m -> output.react, dummy_s2, dummy_m -> output.n_react, dummy_m -> output.coeff);
                
                /* reduce the number of reactions in which dummy_m appears as an output */
                dummy_m -> output.n_react -= 1;
                
            }
            
        }
        
    }
    
    /* loop over all metabolites and check whether they are only consumed */
    for (dummy_m = metabs; dummy_m < metabs + Nmet; dummy_m++){
        
        /* check if metabolite is only consumed */
        if ( dummy_m -> output.n_react == 0 && dummy_m -> input.n_react != 0){
            /* if so, set to zero all reactions that consume it */
            
            fprintf(outfile, "Metabolite %d is now only consumed...\n", (int) (dummy_m - metabs) + 1);
            
            /* realloc the array of locked reactions to add all reactions that consume metabolite dummy_m */
            s_locked = (double **) realloc (s_locked, (n_zeros + n_zeros_new + dummy_m -> input.n_react) * sizeof (double *) );
            
            /* initialise the pointer to the old last element of s_locked */
            dummy_s1 = s_locked + n_zeros + n_zeros_new;
            
            /* loop over reactions that consume metabolite dummy_m */
            for (dummy_s2 = dummy_m -> input.react; dummy_s2 < dummy_m -> input.react + dummy_m -> input.n_react; dummy_s2++){
                
                fprintf(outfile, "Setting reaction %d to zero.....\n", (int) (*dummy_s2 - s) + 1);
                
                /* add reaction to s_locked */
                *dummy_s1 = *dummy_s2;
                
                /* iterate over s_locked */
                dummy_s1++;
                
            }
            
            /* the number of null reactions is increased */
            n_zeros_new += dummy_m -> input.n_react;
            
        }
    }
    
    
    
    fprintf(outfile, "\n\n");
    
    /* if new reactions are set to zero, check recursively feasibility of the system */
    if ( n_zeros_new > 0 )  return check_cascades (metabs, Nmet, s_zeros, n_zeros + n_zeros_new, s, outfile);
    
    /* otherwise return the final number of null reactions */
    else return n_zeros;
    
}
