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

#include "vN_io.h"

void print_matrix (metabolite *metabs, int Nmet, double *s, int Nreac, FILE *outfile){
    
    metabolite *dummy;
    double *dummy_s0, **dummy_s, *c;
    int done;
    
    for ( dummy = metabs; dummy < metabs + Nmet; dummy ++){
        
        fprintf(outfile, "%d: ", (int)(dummy - metabs) + 1);
        
        for ( dummy_s0 = s; dummy_s0 < s + Nreac; dummy_s0++){
            
            done = 0;
            
            c = dummy -> input.coeff;
            
            for( dummy_s = dummy -> input.react; dummy_s < dummy -> input.react + dummy -> input.n_react; dummy_s++){
                
                if ( *dummy_s == dummy_s0) {
                    
                    fprintf(outfile, "%g ", -(*c) );
                    
                    done = 1;
                }
                
                c++;
            }
            
            c = dummy -> output.coeff;
            
            for( dummy_s = dummy -> output.react; dummy_s < dummy -> output.react + dummy -> output.n_react; dummy_s++){
                
                if ( *dummy_s == dummy_s0) {
                    
                    fprintf(outfile, "%g ", (*c) );
                    
                    done = 1;
                }
                
                c++;
            }
            
            if (done == 0) fprintf(outfile, "0 ");
        }
        
        fprintf(outfile, "\n");
    }
    
}

/* A function to re-print the adjacency list of the system, mainly for debug purposes */
void print_adj_list (metabolite *metabs, int Nmet, double *s, FILE *outfile){
    
    metabolite *dummy;
    double **s_dummy, *coeff;
    
    /* Loop over emtabolites */
    for (dummy = metabs; dummy < metabs + Nmet; dummy++){

        fprintf(outfile, "%d: ", (int)(dummy - metabs) + 1);
        
        /* initialise coeff to the first element of input.coeff array */
        coeff = dummy -> input.coeff;
        
        /* loop over input reactions */
        for (s_dummy = dummy -> input.react; s_dummy < dummy -> input.react + dummy -> input.n_react; s_dummy++){
            
            /* print reaction index & coefficient */
            fprintf(outfile, "%d %g ", (int) (*s_dummy - s) + 1, -(*coeff));
            
            /* iterate the coefficient along the input.coeff array */
            coeff++;
        }
        
        
        /* initialise coeff to the first element of output.coeff array */
        coeff = dummy -> output.coeff;
        
        /* loop over output reactions */
        for (s_dummy = dummy -> output.react; s_dummy < dummy -> output.react + dummy -> output.n_react; s_dummy++){
            
            /* print reaction index & coefficient */
            fprintf(outfile, "%d %g ", (int) (*s_dummy - s) + 1, *coeff);
            
            /* iterate the coefficient along the output.coeff array */
            coeff++;
        }
        
        fprintf(outfile, "\n");
    }
}

/* A function to re-print the system, mainly for debug purposes */
void print_system (metabolite *metabs, int Nmet, double *s, FILE *outfile){
    
    metabolite *dummy;
    double **s_dummy, *coeff;
    
    /* Loop over emtabolites */
    for (dummy = metabs; dummy < metabs + Nmet; dummy++){
        
        fprintf(outfile, "Metabolite %d enters:\n", (int)(dummy - metabs) + 1);
        
        fprintf(outfile, "As input in %d reactions:\n", dummy -> input.n_react);
        
        /* initialise coeff to the first element of input.coeff array */
        coeff = dummy -> input.coeff;
        
        /* loop over input reactions */
        for (s_dummy = dummy -> input.react; s_dummy < dummy -> input.react + dummy -> input.n_react; s_dummy++){
            
            /* print reaction index & coefficient */
            fprintf(outfile, "\t %d with coefficient %g\n", (int) (*s_dummy - s) + 1, *coeff);
            
            /* iterate the coefficient along the input.coeff array */
            coeff++;
        }
        
        fprintf(outfile, "As output in %d reactions:\n", dummy -> output.n_react);
        
        /* initialise coeff to the first element of output.coeff array */
        coeff = dummy -> output.coeff;
        
        /* loop over output reactions */
        for (s_dummy = dummy -> output.react; s_dummy < dummy -> output.react + dummy -> output.n_react; s_dummy++){
            
            /* print reaction index & coefficient */
            fprintf(outfile, "\t %d with coefficient %g\n", (int) (*s_dummy - s) + 1, *coeff);
            
            /* iterate the coefficient along the output.coeff array */
            coeff++;
        }
        
        fprintf(outfile, "\n\n");
    }
}

/* A function to re-print the reactions, mainly for debug purposes */
void print_reactions (double *s, int Nreac, metabolite *metabs, int Nmet, FILE *out_file){
    
    double *dummy_s, **io_r, *c;
    metabolite *dummy_m;
    int n_input, n_output;
    
    /* loop over reactions */
    for (dummy_s = s; dummy_s < s + Nreac; dummy_s++){
        
        /* initialise the number of metabolite of input & output */
        n_input = 0;
        n_output = 0;
        
        /* print reaction index */
        fprintf(out_file, "%d: ", (int)(dummy_s - s) + 1);
        
        /* loop over metabolites & look for reaction dummy_s */
        for (dummy_m = metabs; dummy_m < metabs + Nmet; dummy_m++){
            
            /* initialise c to the first element of input.coeff array */
            c = dummy_m -> input.coeff;
            
            /* loop over input reactions of metabolite dummy_m */
            for (io_r = dummy_m -> input.react; io_r < dummy_m -> input.react + dummy_m -> input.n_react; io_r++){
                
                /* if dummy_s is in the reaction list then */
                if ( *io_r == dummy_s){
                    
                    /* if it is not the first input metabolite, we want to print a "+" */
                    if (n_input > 0 ) fprintf(out_file, " + ");
                    
                    /* if coefficient is 1 we don't want to print it out */
                    /* just print the metabolite index */
                    if (*c == 1.) fprintf(out_file, "%d ", (int)(dummy_m - metabs) + 1 );
                    
                    /* otherwise, print the coefficient & metabolite index */
                    else fprintf(out_file, "(%g) %d ", *c, (int)(dummy_m - metabs) + 1 );
                    
                    /* update the number of input metabolites */
                    n_input++;
                }
                
                /* iterate the coefficient along the input.coeff array */
                c++;
                
            }
        }
        
        /* get to the other side of reaction */
        fprintf(out_file, " --> ");
        
        /* loop over metabolites & look for reaction dummy_s */
        for (dummy_m = metabs; dummy_m < metabs + Nmet; dummy_m++){
            
            /* initialise c to the first element of output.coeff array */
            c = dummy_m -> output.coeff;
            
            /* loop over output reactions of metabolite dummy_m */
            for (io_r = dummy_m -> output.react; io_r < dummy_m -> output.react + dummy_m -> output.n_react; io_r++){
                
                /* if dummy_s is in the reaction list then */
                if ( *io_r == dummy_s){
                    
                    /* if it is not the first output metabolite, we want to print a "+" */
                    if (n_output > 0 ) fprintf(out_file, " + ");
                    
                    /* if coefficient is 1 we don't want to print it out */
                    /* just print the metabolite index */
                    if (*c == 1.) fprintf(out_file, "%d + ", (int)(dummy_m - metabs) + 1);
                    
                    /* otherwise, print the coefficient & metabolite index */
                    else fprintf(out_file, "(%g) %d + ", *c, (int)(dummy_m - metabs) + 1);
                    
                    /* update the number of input metabolites */
                    n_input++;
                    
                }
                
                /* iterate the coefficient along the output.coeff array */
                c++;
                
            }
        }
        
        fprintf(out_file, "\n");
    }
    
}


/* A function to print the von Neumann constraints, mainly to check everithing is >= 0. */
void print_constraints (double *s, metabolite *metabs, int Nmet, double rho, FILE *outfile){
    
    metabolite *dummy_m;
    double c, *coeff, **s_d;
    
    /*loop over metabolites*/
    for (dummy_m = metabs; dummy_m < metabs + Nmet; dummy_m++){
        
        c=0;
        
        /* loop over (input) reactions attached to the metabolite dummy_m */
        coeff = dummy_m -> input.coeff;
        
        for (s_d = dummy_m -> input.react; s_d < dummy_m -> input.react + dummy_m -> input.n_react; s_d++){
            
            /* calculate negative contribution to constraint */
            c -= **s_d * (*coeff) * rho;
            
            coeff++;
        }
        
        /* loop over (output) reactions attached to the metabolite dummy_m */
        coeff = dummy_m -> output.coeff;
        
        for (s_d = dummy_m -> output.react; s_d < dummy_m -> output.react + dummy_m -> output.n_react; s_d++){
            
            /* calculate positive contribution to constraint (rho does not affect it) */
            c += **s_d * (*coeff);
            
            coeff++;
        }
        
        fprintf (outfile, "Metabolite %d c %g\n", (int) (dummy_m - metabs), c);
        
    }
}

/* A function to print the flux values */
void print_fluxes (double *s, int Nreac, FILE *outfile){
    
    double *dummy_s, Z=0.;
    
    /* run over fluxes */
    for (dummy_s = s; dummy_s < s + Nreac; dummy_s++) {
        
        fprintf(outfile, "%g ", *dummy_s);
        
        Z += *dummy_s;
    }
    fprintf(outfile, "\n");
    
    /* UNCOMMENT LINE TO ENSURE NORMALISATION IS ATTAINED fprintf(outfile, "\n##\t Z = %g \n", Z); */
    
}
