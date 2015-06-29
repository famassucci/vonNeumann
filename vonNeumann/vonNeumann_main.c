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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <unistd.h>

#include "vonNeumann.h"


#ifndef LOG_FILE
#define LOG_FILE "von_Neumann.log"
#endif

#ifndef N_SOL
#define N_SOL 1
#endif

#ifndef STEP_INIT
#define STEP_INIT 1.e-4
#endif

#ifndef STEP_MIN
#define STEP_MIN 1.e-6
#endif

#ifndef N_STEP_MAX
#define N_STEP_MAX 1e6
#endif

#ifndef RHO_INIT
#define RHO_INIT 0.95
#endif

#ifndef RHO_MAX
#define RHO_MAX 0.999
#endif

#ifndef ETA
#define ETA 0.0001
#endif

void print_usage () {
    printf ("\n");
    printf ("NAME\n");
    printf ("\tvonNeumann -- a minOver sampler of solutions to a von Neumann problem\n\n");
    printf ("SYNOPSIS\n");
    printf ("\tvonNeumann [options] file\n\n");
    printf ("DESCRIPTION\n");
    printf ("\tvonNeumann reads an input file and seeks s solutions to the (von Neumann) problem:\n");
    printf ("\t\t s (a - rho b) >= 0 \n");
    printf ("\t up to a maximal rho. Input files may be an adjacency list, a stoichiometric matrix, or a reaction list.\n\n");
}

void print_help () {
    
    printf ("\n");
    printf ("\tThe following options are available:\n");
    printf ("\t-e [ETA] Specify the factor eta for the update step. Default ETA=%g.\n", ETA);
    printf ("\t-h: print this help and exit.\n");
    printf ("\t-L \"...\" Lock reactions. A comma separated list of reaction indices : lock values parameters must be provided in apices.\n");
    printf ("\t-M [MAX_STEP] Fix the maximum number of steps of minOver algorithm. Default MAX_STEP=%g.\n", N_STEP_MAX);
    printf ("\t-n [N_SOL] Specify the number of solutions. Default N_SOL=%d.\n", N_SOL);
    printf ("\t-o [FILE] Specify the output file. Default stdout.\n");
    printf ("\t-r [RHO_INIT] Specify the initial rho value. Default RHO_INIT=%g.\n", RHO_INIT);
    printf ("\t-R [RHO_MAX] Specify maximum rho value. Default RHO_MAX=%g.\n", RHO_MAX);
    printf ("\t-S [INIT_STEP_SIZE] Specify the size of the initial step used to update fluxes. Default INIT_STEP_SIZE=%g.\n", STEP_INIT);
    printf ("\t-S [MIN_STEP_SIZE] Specify the minimum step size that can be handled by minOver. Default MIN_STEP_SIZE=%g.\n", STEP_MIN);
    printf ("\t-v Verbose. Print logfile to stderr instead of %s.\n\n", LOG_FILE);
    
    printf ("Note:\n");
    printf ("\t -- Output values are normalised to the number of reactions.\n");
}

int main (int argc, char *argv[] ){
    
    int c, vflag = 0, Lflag = 0, nflag = 0, Sflag = 0, sflag = 0, Mflag = 0, rflag = 0, Rflag = 0, eflag = 0, oflag = 0;
    
    char *LOCKED;
    
    int Nreact, Nmetabs, n_locked = 0, n_null=0, n_null_final;
    
    int sol, n_sol = N_SOL, n_step_max = N_STEP_MAX;
    
    double step_init = STEP_INIT, step_min = STEP_MIN, rho_init = RHO_INIT, rho_max = RHO_MAX, eta = ETA;
    
    double *s, *s_backup, **s_locked = (double **)NULL, *lock_v = (double *)NULL, **s_null = (double **)NULL, rho;
    
    FILE *log_file, *out_file = stdout;
    
    
    /* parse command line options */
    while ((c = getopt (argc, argv, "vhL:n:S:s:M:r:R:e:o:")) != -1) {
        switch (c) {
            
                /* help flag */
            case 'h' :
                print_usage ();
                print_help ();
                exit (0);
                break;
                
                /* eta flag -- fix step factor size */
            case 'e':
                
                eflag = 1;
                
                eta = atof ( optarg );
                
                break;

                /* lock flag, fix locked reactions */
            case 'L' :
                Lflag = 1;
                
                LOCKED = (char *) malloc( strlen(optarg)*sizeof(char) );
                strcpy(LOCKED , optarg);
                n_locked = get_n_locked (LOCKED);
                
                break;
                
                /* step max flag, fix the maximum number of steps for minover */
            case 'M':
                
                Mflag = 1;
                
                n_step_max = atoi ( optarg );
                
                break;

                /* nsol flag, fix the number of solutions */
            case 'n':
                
                nflag = 1;
                
                n_sol = atoi ( optarg );
                
                break;
            
                /* output flag, specify the output file */
            case 'o':
                
                oflag = 1;
                
                out_file = fopen( optarg, "w");
                
                break;
            
                /* rho init flag, fix the initial value of rho */
            case 'r':
                
                rflag = 1;
                
                rho_init = atof ( optarg );
                
                break;
                
                /* rho max flag, fix the maximum value of rho */
            case 'R':
                
                Rflag = 1;
                
                rho_max = atof ( optarg );
                
                break;
                
                /* step init flag, fix the initial step size */
            case 'S':
                
                Sflag = 1;
                
                step_init = atof ( optarg );
                
                break;
                
                /* step min flag, fix the minimum step size */
            case 's':
                
                sflag = 1;
                
                step_min = atof ( optarg );

                break;
                
                /* verbose flag, print log to stderr */
            case 'v' :
                vflag = 1;
                
                log_file = stderr;
                
                break;
                
            default :
                print_usage ();
                exit (EXIT_FAILURE);
        }
    };
    
    /* check that we have the right number of arguments */
    if (optind+1!=argc) {
        fprintf (stderr, "Incorrect usage...\n");
        print_usage ();
        exit (EXIT_FAILURE);
    }

    
    /* if not verbose, open the log file */
    if (vflag == 0) log_file = fopen(LOG_FILE, "w");
    
    /* open the input file and store it into the file_wrappwer structure */
   file_wrapper *input_data = handle_input_file (argv[ argc - 1], 10, log_file);
    
    /* retrieve the number of reactions and metabolites from the file_wrapper struct */
    Nreact = input_data -> Nreact;
    
    Nmetabs = input_data -> Nmet;
    
    /* keep track of everything in the log file */
    fprintf(log_file, "The system has %d metabolites and %d Reactions\n", Nmetabs, Nreact);
    
    /* allocate reactions */
    s = (double *) malloc ( Nreact * sizeof(double) );
    
    /* allocate space to backup reactions */
    s_backup = (double *) malloc ( Nreact * sizeof(double) );
    
    /* allocate space for metabolite structure */
    metabolite *metabs;
    
    metabs = (metabolite *) malloc( Nmetabs * sizeof (metabolite) );
    
    alloc_system (input_data, metabs, Nmetabs, s, log_file);
    
    /* free the file wrapper structure */
    file_wrapper_free(&input_data);
    
    /* if locking some reactions */
    if ( n_locked > 0 ){
        
        /* create an array of pointers associated to the locked reactions */
        s_locked = (double **) malloc ( n_locked*sizeof(double*) );
        
        /* and the lock value */
        lock_v = (double *) malloc ( n_locked*sizeof(double) );
        
        /* get also the number of zero reactions (they can be effectively removed from the system) */
        n_null = fix_locked (n_locked, s, s_locked, lock_v, LOCKED);
        
        /* if there are zero reactions, remove them */
        if ( n_null > 0) s_null = assign_null_reactions (s_locked, lock_v, n_null, n_locked);
    }
 
    /* check feasibility of the system, i.e. whether there are metabolites that are only consumed */
    n_null_final = check_cascades (metabs, Nmetabs, &s_null, n_null, s, log_file);
    
    /* if to make the system feasible, some reactions have been forced to zero... */
    if ( n_null_final != n_null){
        
        /* realloc the null reactions and values */
        s_locked = (double **) realloc (s_locked, n_locked + (n_null_final - n_null) * sizeof(double*) );
        
        lock_v = (double *) realloc (lock_v, n_locked + (n_null_final - n_null) * sizeof(double) );
        
        /* keep track of all locked reactions */
        n_locked = update_null_reactions (s_locked, lock_v, n_null, n_null_final, n_locked, s_null);
        
    }
    
    /* keep track of everything in the log file */
    fprintf(log_file, "\n\nThe system has %d locked reactions (%d of them null)\n", n_locked, n_null_final);

    /* initialise the random number generator */
    srand48( time (NULL) );
    
    /* for non verbose output */
    /* for n_sol times */
    for (sol = 0; sol < n_sol * (1-vflag); sol++) {
        
        /* sample reactions up to the maximum rho */
        rho = optimal_flux (metabs, s, s_locked, n_locked, lock_v , s_backup, Nmetabs, Nreact, n_step_max, step_init, step_min, rho_init, rho_max, eta);
    
        /* keep track to the max rho (may be smaller than rho_max if the step decreases too much) */
        fprintf(log_file, "Solution %d, rho = %g\n", sol+1, rho);
        
        /* print the fluxes sampled at max rho */
        print_fluxes (s, Nreact, out_file);
    }
    
    /* for verbose output */
    for (sol = 0; sol < n_sol*vflag; sol++) {
        
        /* sample reactions up to the maximum rho */
        rho = optimal_flux_verbose (metabs, s, s_locked, n_locked, lock_v , s_backup, Nmetabs, Nreact, n_step_max, step_init, step_min, rho_init, rho_max, eta, log_file);
        
        /* keep track to the max rho (may be smaller than rho_max if the step decreases too much) */
        fprintf(log_file, "Solution %d, rho = %g\n", sol+1, rho);
        
        /* print the fluxes sampled at max rho */
        print_fluxes (s, Nreact, out_file);
    }
    
    /* close output files */
    if (log_file != stderr) fclose(log_file);
    
    if (out_file != stdout) fclose(out_file);
    
    /* free the system */
    free (s);
    
    free (s_backup);
    
    metabolite_free (&metabs, Nmetabs);
        
    return 0;

}
