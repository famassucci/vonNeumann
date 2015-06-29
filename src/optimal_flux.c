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

#include "optimal_flux.h"

/* compute the fluxes up to max rho (or min step), starting from an initial rho value*/
double optimal_flux (metabolite *metabs, double *s, double **locked, int n_locked, double *lock_value, double *s_backup, int Nmet, int Nreac, int max_step_init, double step_init, double step_min, double rho_min, double rho_max, double eta){
    
    if (rho_min > rho_max) {
        
        fprintf(stderr, "RHO min cannot be larger than RHO max \n");
        
        exit (EXIT_FAILURE);
        
    }
    
    if (step_init < step_min) {
        
        fprintf(stderr, "Initial step size cannot be smaller than final step size\n");
        
        exit (EXIT_FAILURE);
        
    }
    
    double rho = rho_min, step = step_init, eta_factor = 10.;
    int n_step, max_step = max_step_init;
    
    
    /* initialise fluxes */
    initialise_fluxes (s, Nreac, locked, n_locked, lock_value);
    
    
    /* and store values */
    backup_fluxes (s, s_backup, Nreac);
    
    /* starting from an initial rho value, satisfy constraints for a given rho and increase rho up to rho_max */
    /* if the step becomes too small (i.e. convergence very slow) exit from loop and return last succesful flux array */
    while (rho < rho_max && step > step_min){
        
        /* run minover at given rho */
        n_step = minover (metabs, s, locked, n_locked, lock_value, rho, Nmet, max_step, eta*eta_factor, Nreac);
        
        /* minover returns the number of steps to reach convergence*/
        /* if n steps > max step -> no convergence, restore last succesful value and reduce rho */
        if (n_step >= max_step) {
            
            /* restore the last succesful array */
            restore_backup ( s, s_backup, Nreac);
            
            /*decrease rho to last succesful value */
            rho -= step;
            
            /* decrease the step value */
            step /= 1.5;
            
            /* decrease the eta factor */
            eta_factor /= 1.2;
            
            max_step *= 1.5;
            
        }
        
        else {
            
            /* if successful, normalise the fluxes */
            normalise_fluxes (s, Nreac, locked, n_locked, lock_value);
            
            
            /* and then backup the flux values */
            backup_fluxes (s, s_backup, Nreac);
            
            /* if approaching non - convergence, reduce the step size */
            if ( n_step > 2 * max_step / 3) {
                
                step /= 1.2;
                
                /* decrease the eta factor */
                eta_factor /= 1.1;
                
                max_step *= 1.2;
            }
        }
        
        if (rho == rho_min) eta_factor /= 15.;
        
        /* increase rho value */
        rho += step;
    }
    
    return rho -step;
}

/* A verbose version of the function above: the only difference is to print rho values while sampling */
/* compute the fluxes up to max rho (or min step), starting from an initial rho value*/
double optimal_flux_verbose (metabolite *metabs, double *s, double **locked, int n_locked, double *lock_value, double *s_backup, int Nmet, int Nreac, int max_step_init, double step_init, double step_min, double rho_min, double rho_max, double eta, FILE *log_file){
    
    if (rho_min > rho_max) {
        
        fprintf(stderr, "RHO min cannot be larger than RHO max \n");
        
        exit (EXIT_FAILURE);
        
    }
    
    if (step_init < step_min) {
        
        fprintf(stderr, "Initial step size cannot be smaller than final step size\n");
        
        exit (EXIT_FAILURE);
        
    }
    
    double rho = rho_min, step = step_init, eta_factor = 10.;
    int n_step, max_step = max_step_init;
    
    
    /* initialise fluxes */
    initialise_fluxes (s, Nreac, locked, n_locked, lock_value);
    
    
    /* and store values */
    backup_fluxes (s, s_backup, Nreac);
    
    /* starting from an initial rho value, satisfy constraints for a given rho and increase rho up to rho_max */
    /* if the step becomes too small (i.e. convergence very slow) exit from loop and return last succesful flux array */
    while (rho < rho_max && step > step_min){
        
        /* run minover at given rho */
        n_step = minover (metabs, s, locked, n_locked, lock_value, rho, Nmet, max_step, eta*eta_factor, Nreac);
        
        /* minover returns the number of steps to reach convergence*/
        /* if n steps > max step -> no convergence, restore last succesful value and reduce rho */
        if (n_step >= max_step) {
            
            /* restore the last succesful array */
            restore_backup ( s, s_backup, Nreac);
            
            /*decrease rho to last succesful value */
            rho -= step;
            
            /* decrease the step value */
            step /= 1.5;
            
            /* decrease the eta factor */
            eta_factor /= 1.2;
            
            max_step *= 1.5;
            
        }
        
        else {
            
            /* if successful, normalise the fluxes */
            normalise_fluxes (s, Nreac, locked, n_locked, lock_value);
            
            
            /* and then backup the flux values */
            backup_fluxes (s, s_backup, Nreac);
            
            /* if approaching non - convergence, reduce the step size */
            if ( n_step > 2 * max_step / 3) {
                
                step /= 1.2;
                
                /* decrease the eta factor */
                eta_factor /= 1.1;
                
                max_step *= 1.2;
            }
        }
        
        if (rho == rho_min) eta_factor /= 15.;
        
        /* increase rho value */
        rho += step;
        
        
        
        fprintf(log_file, "\r step %g rho %g ", step, rho);
    }
    
    fprintf(log_file, "\n");
    return rho -step;
}
