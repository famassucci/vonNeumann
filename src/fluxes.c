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

#include "fluxes.h"

/* A function to back up flux values that work for a given value of rho */
void backup_fluxes ( double *s, double *s_backup, int Nreact){
    
    double *dummy1, *dummy2;
    
    dummy2 = s;
    
    /* loop over flux values*/
    for (dummy1 = s_backup; dummy1 < s_backup + Nreact; dummy1++) {
        
        /* and swap with backup */
        *dummy1 = *dummy2;
        
        dummy2++;
    }
}

/* A function to restore flux values that worked for a smaller value of rho */
void restore_backup ( double *s, double *s_backup, int Nreact){
    
    double *dummy1, *dummy2;
    
    dummy2 = s_backup;
    
    /* loop over flux values*/
    for (dummy1 = s; dummy1 < s + Nreact; dummy1++) {
        
        /* and swap with backup */
        *dummy1 = *dummy2;
        
        dummy2++;
    }
}

/* Initialise the flux values. Fluxes are normalised so that their sum = N reactions */
void initialise_fluxes (double *s, int Nreact, double **locked, int n_locked, double *lock_value){
    
    double *dummy, **dummy_l, *lock_v, Z = 0., Z_l = 0.;
    
    for (dummy = s; dummy < s + Nreact; dummy++) {
        
        /* the initial values are absolute values of a gaussian random v.*/
        /* *dummy = fabs( 1.e-6 + 1.e-4 * gaussdev() );*/
        *dummy = fabs(  gaussdev() );
        
        /* update normalisation coefficient */
        Z += *dummy;
        
    }
    
    for (dummy = s; dummy < s + Nreact; dummy++) {
        
        /* normalise reaction, resting the contribution of locked reactions */
        *dummy *= ( (double) Nreact - Z_l)/Z;

        
    }
    
    /* keep the locked reaction fixed*/
    lock_v = lock_value;
    
    /* re-loop over locked reactions to assign the value */
    /* (this should be faster than putting an if statement in the above loop) */
    for (dummy_l = locked; dummy_l < locked + n_locked; dummy_l++) {
        
        /* assign the right value */
        **dummy_l = *lock_v;
        
        
        lock_v++;
    }
}

/* A function to re-normalise fluxes after each MinOver run */
void normalise_fluxes (double *s, int Nreact, double **locked, int n_locked, double *lock_value){
    
    double *dummy, **dummy_l, Z = 0., Z_l = 0.;
    
    /* compute the normalisation coefficient*/
    for (dummy = s; dummy < s + Nreact; dummy++) Z += *dummy;
    
    /* rest the locked values */
    for (dummy_l = locked; dummy_l < locked + n_locked; dummy_l++) {
        
        Z -= **dummy_l;
        
        Z_l += **dummy_l;
    }
    
    for (dummy = s; dummy < s + Nreact; dummy++) {
        
        /* normalise reaction, resting the contribution of locked reactions */
        *dummy *= ( (double) Nreact - Z_l)/Z;
        
    }
    
    /* keep the locked reaction fixed*/
    dummy = lock_value;
    
    /* re-loop over locked reactions to assign the value */
    /* (this should be faster than putting an if statement in the above loop) */
    for (dummy_l = locked; dummy_l < locked + n_locked; dummy_l++) {
        
        /* assign the right value */
        **dummy_l = *dummy;
        
        
        dummy++;
    }

}
