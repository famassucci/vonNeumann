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

#include "minover.h"

/*run the minover algorithm for fixed rho value*/
int minover (metabolite *metabs, double *s, double **locked, int n_locked, double *lock_value, double rho, int Nmet, int max_step, double eta, int Nreac){
    
    /* pointers to iterate over the metabolites*/
    metabolite *dummy_m, *dummy0;
    
    double c, cmu0, **s_d, *coeff, check_sign;
    
    int step=0;
    
    /* iterate the algorithm until all constraints are satisfied*/
    do{
        
        
        /*loop over metabolites*/
        for (dummy_m = metabs; dummy_m < metabs + Nmet; dummy_m++){
            
            c=0;
            
            /* loop over (input) reactions attached to the metabolite dummy_m */
            coeff = dummy_m -> input.coeff;
            
            for (s_d = dummy_m -> input.react; s_d < dummy_m -> input.react + dummy_m -> input.n_react; s_d++){
                
                c -= **s_d * (*coeff) * rho;
                
                coeff++;
            }
            
            /* loop over (output) reactions attached to the metabolite dummy_m */
            coeff = dummy_m -> output.coeff;
            
            for (s_d = dummy_m -> output.react; s_d < dummy_m -> output.react + dummy_m -> output.n_react; s_d++){
                
                c += **s_d * (*coeff);
                
                coeff++;
            }
            
            /* store the minimum constraint */
            if( dummy_m == metabs || c < cmu0) {
                
                cmu0 = c;
                
                dummy0 = dummy_m;
            }
            
        }
        
        /* if some constraint is unsatisfied, update fluxes */
        if (cmu0 < 0 ){
            
            /*if (rho > 0.986) fprintf(stderr, "\r mu0 %d cmu0 %g \r", (int)(dummy0 - metabs) +1, cmu0);*/
            /* update (input) reactions attached to metabolite dummy0 */
            coeff = dummy0 -> input.coeff;
            
            for (s_d = dummy0 -> input.react; s_d < dummy0 -> input.react + dummy0 -> input.n_react; s_d++){
                
                /* update the flux according to minover rule. A factor eta is added to ease convergence */
                **s_d -=  (*coeff) * rho * eta;
                
                /* if negative, set to 0 */
                check_sign = sign(**s_d);
                
                **s_d *= (1. + check_sign )/ 2.;
                
                coeff++;
            }
            
            /* update (output) reactions attached to metabolite dummy0 */
            coeff = dummy0 -> output.coeff;
            
            for (s_d = dummy0 -> output.react; s_d < dummy0 -> output.react + dummy0 -> output.n_react; s_d++){
                
                /* update the flux according to minover rule. A factor eta is added to ease convergence */
                **s_d +=  (*coeff) * eta;
                
                coeff++;
            }
            
        }
        
        /* keep the locked reaction fixed*/
        coeff = lock_value;
        for (s_d = locked; s_d < locked + n_locked; s_d++) {
            
            **s_d = *coeff;
            
            coeff++;
        }
        
        normalise_fluxes (s, Nreac, locked, n_locked, lock_value);
        
        step++;
        
    } while (cmu0 < 0 && step < max_step );
    
    return step;
    
}
