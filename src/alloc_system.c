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

#include "alloc_system.h"

/* A function to allocate memory for the problem reading the input from a reaction list */
/* all relevant stoichiometric information has been previously stored in a "metabolite_parse" structure*/
void alloc_from_filetype2 (metabolite *system, double *s, metabolite_parse *data, int Nmet, FILE *outfile){
    
    metabolite *dummy;
    metabolite_parse *parser;
    int *which, cnt;
    double *coeff;
    
    dummy = system;
    
    /* loop over the stored data*/
    for ( parser = data; parser < data + Nmet; parser++) {
        
        /* get the # of reactions where dummy metabolite appears as an input */
        dummy -> input.n_react = parser -> input.n_react;
        
        /* allocate memory for reaction/coefficients pointers */
        dummy -> input.react = (double **) malloc( dummy -> input.n_react * sizeof ( double *) );
        
        dummy -> input.coeff = (double *) malloc( dummy -> input.n_react * sizeof ( double ) );
        
        /* for each reaction, get memory position and coefficient */
        coeff = parser -> input.coeff;
        
        /* loop over the "metabolite_parse" to get the reactions attached to metabolite dummy */
        cnt = 0;
        for (which=parser -> input.which_r; which < parser -> input.which_r + parser -> input.n_react; which++) {
            
            
            /* assign pointer to reaction */
            *( dummy -> input.react + cnt ) = s + (*which);
            
            /* assign coefficient value */
            *( dummy -> input.coeff + cnt ) = *coeff;
            
            coeff++;
            cnt++;
        }
        
        /* get the # of reactions where dummy metabolite appears as an ouptut */
        dummy -> output.n_react = parser -> output.n_react;
        
        /* allocate memory for reaction/coefficients pointers */
        dummy -> output.react = (double **) malloc( dummy -> output.n_react * sizeof ( double *) );
        
        dummy -> output.coeff = (double *) malloc( dummy -> output.n_react * sizeof ( double ) );
        
        /* for each reaction, get memory position and coefficient */
        coeff = parser -> output.coeff;
        
        /* loop over the "metabolite_parse" to get the reactions attached to metabolite dummy */
        cnt = 0;
        for (which=parser -> output.which_r; which < parser -> output.which_r + parser -> output.n_react; which++) {
            
            
            /* assign pointer to reaction */
            *( dummy -> output.react + cnt ) = s + (*which);
            
            /* assign coefficient value */
            *( dummy -> output.coeff + cnt ) = *coeff;
            
            coeff++;
            cnt++;
        }

        /* increase dummy to next metabolite */
        dummy++;
    }
    
}

/* A function to allocate memory for the problem reading the input from a stoichiometric matrix */
/* all relevant stoichiometric information is stored in a string named "file_content" */
void alloc_from_filetype1(char *file_content, metabolite *mets, int Nmet, double *s, FILE *outfile){
  
    char *s1 = file_content, *colon, *s2, *s3;
    int n;
    double **react_in, **react_out, *coeff_in, *coeff_out, c;
    metabolite *dummy = mets;
    
    /* separate the file_content string per lines */
    while ( strchr(s1, '\n') != NULL ){
        
        /* get a pointer to the end of line */
        s2 = strchr(s1, '\n');
        
        /* only read if it is not a comment */
        if ( *s1 != '#'){
            
            /* a descriptor of the line may be present with a ":" sign to separate */
            colon = strchr(s1, ':');
            
            /* if there is a ":", ignore whatever comes before it */
            if (colon != NULL ){
                
                s1 = colon + 1;
                
            }
            
            /* get rid of trailing space */
            while (*s1 == ' ') s1++;
            
            /* store a pointer to the beginning of line */
            colon = s1;
            
            /* keep track of everything in the log file*/
            fprintf(outfile, "Allocating metabolite %d... ", (int)(dummy - mets) + 1); fflush(outfile);
            
            /* set a null sign at the end of line to stop reading there */
            *s2 = '\0';
            
            /* initialise the number of i/o reactions for metabolite dummy */
            dummy -> input.n_react = 0;
            
            dummy -> input.n_react = 0;
            
            /* for each column of the line */
            while ( strchr(s1, ' ') != NULL ){
                
                /* get to the next space */
                s3 = strchr(s1, ' ');
                
                /* stop the string there */
                *s3 = '\0';
                
                /* get the stoichiometric coefficient */
                c = atof(s1);
                
                /* if negative, increase the # of reaction where dummy appears as input */
                if (c < 0. ) dummy -> input.n_react++;
                
                /* if positive, increase the # of reaction where dummy appears as output */
                if (c > 0. ) dummy -> output.n_react++;
                
                /* undo the null sign */
                *s3 = ' ';
                
                /* keep reading */
                s1 = s3 + 1;
            }
            
            /* for the last column */
            /* get the stoichiometric coefficient */
            c = atof(s1);
            
            /* if negative, increase the # of reaction where dummy appears as input */
            if (c < 0. ) dummy -> input.n_react++;
            
            /* if positive, increase the # of reaction where dummy appears as output */
            if (c > 0. ) dummy -> output.n_react++;
            
            /* allocate memory for reaction/coefficients pointers */
            dummy -> input.react = (double **) malloc( dummy -> input.n_react * sizeof ( double *) );
            
            dummy -> input.coeff = (double *) malloc( dummy -> input.n_react * sizeof ( double ) );
            
            /* allocate memory for reaction/coefficients pointers */
            dummy -> output.react = (double **) malloc( dummy -> output.n_react * sizeof ( double *) );
            
            dummy -> output.coeff = (double *) malloc( dummy -> output.n_react * sizeof ( double ) );
            
            /* get back to the beginning of line and assign now reactions & coefficients to dummy */
            s1 = colon;
            
            /* an index to keep track of the matrix column i.e. reaction index*/
            n = 0;
            
            /* initialise pointers to reactions and coefficients */
            react_in = dummy -> input.react;
            
            react_out = dummy -> output.react;
            
            coeff_in = dummy -> input.coeff;
            
            coeff_out = dummy -> output.coeff;
            
            /* read again per column */
            while ( strchr(s1, ' ') != NULL ){
                
                /* get to the next space */
                s3 = strchr(s1, ' ');
                
                /* stop the string there */
                *s3 = '\0';
                
                /* get the stoichiometric coefficient */
                c = atof(s1);
                
                /* if negative */
                if (c < 0. ) {
                    
                    /* associate reaction input pointer to the corresponding reaction (via column n) */
                    *react_in = s + n;
                    
                    /* associate coefficient to the right one (coefficients are all positive in minOver) */
                    *coeff_in = -c;
                    
                    /* increase pointer to input reactions */
                    react_in++;
                    
                    /* increase pointer to input coefficients */
                    coeff_in++;
                    
                }
                
                /* if positive */
                if (c > 0. ) {
                    
                    /* associate reaction output pointer to the corresponding reaction (via column n) */
                    *react_out = s + n;
                    
                    /* associate coefficient to the right one */
                    *coeff_out = c;
                    
                    /* increase pointer to output reactions */
                    react_out++;
                    
                    /* increase pointer to output coefficients */
                    coeff_out++;
                    
                }
                
                /* undo the null sign */
                *s3 = ' ';
                
                /* keep reading */
                s1 = s3 + 1;
                
                /* increase the column (reaction) index */
                n++;
            }
            
            /* for the last column */
            /* get the stoichiometric coefficient */
            c = atof(s1);
            
            /* if negative */
            if (c < 0. ) {
                
                /* associate inputs as above */
                *react_in = s + n;
                
                *coeff_in = -c;
                
                react_in++;
                
                coeff_in++;
                
            }
            
            /* if positive */
            if (c > 0. ) {
                
                /* associate outputs as above */
                *react_out = s + n;
                
                *coeff_out = c;
                
                react_out++;
                
                coeff_out++;
                
            }
            
            /* undo the null sign at the end of line */
            *s2 = '\n';
            
            /* keep track of everything in the log file */
            fprintf(outfile, "Done.\n"); fflush(outfile);
            
            /* increase pointer to metabolites */
            dummy++;

        }
        
        /* get to the next line */
        s1 = s2 +1;
        
    }
    
}

/* A function to allocate memory for the problem reading the input from an adjacency list */
/* all relevant stoichiometric information is stored in a string named "file_content" */
void alloc_from_filetype0(char *file_content, metabolite *mets, int Nmet, double *s, FILE *outfile){
    
    char *s1 = file_content, *colon, *s2, *s3;
    int n, which;
    double **react_in, **react_out, *coeff_in, *coeff_out, c;
    metabolite *dummy = mets;
    
    /* separate the file_content string per lines */
    while ( strchr(s1, '\n') != NULL ){
        
        /* get a pointer to the end of line */
        s2 = strchr(s1, '\n');
        
        /* only read if it is not a comment */
        if ( *s1 != '#'){
            
            /* a descriptor of the line may be present with a ":" sign to separate */
            colon = strchr(s1, ':');
            
            /* if there is a ":", ignore whatever comes before it */
            if (colon != NULL ){
                
                s1 = colon + 1;
                
            }
            
            /* get rid of trailing space */
            while (*s1 == ' ') s1++;
            
            /* store a pointer to the beginning of line */
            colon = s1;
            
            /* keep track of everything in the log file*/
            fprintf(outfile, "Allocating metabolite %d... ", (int)(dummy - mets) + 1); fflush(outfile);
            
            /* set a null sign at the end of line to stop reading there */
            *s2 = '\0';
            
            /* initialise the number of i/o reactions for metabolite dummy */
            dummy -> input.n_react = 0;
            
            dummy -> output.n_react = 0;
            
            /* initialise the column index */
            n = 0;
            
            /* for each column of the line */
            while ( strchr(s1, ' ') != NULL ){
                
                /* get to the next space */
                s3 = strchr(s1, ' ');
                
                /* stop the string there */
                *s3 = '\0';
                
                /* get the number */
                c = atof(s1);
                
                /* if negative and column odd, increase the # of reaction where dummy appears as input */
                if (c < 0. && n%2 != 0) dummy -> input.n_react++;
                
                /* if positive and column odd, increase the # of reaction where dummy appears as output */
                if (c > 0. && n%2 != 0) dummy -> output.n_react++;
                
                /* undo the null sign */
                *s3 = ' ';
                
                /* keep reading */
                s1 = s3 + 1;
                
                /* increase the column index */
                n++;
            }
            
            /* for the last column */
            /* surely it cannot be a reaction index, so get the stoichiometric coefficient */
            c = atof(s1);
            
            /* if negative, increase the # of reaction where dummy appears as input */
            if (c < 0. ) dummy -> input.n_react++;
            
            /* if positive, increase the # of reaction where dummy appears as output */
            if (c > 0. ) dummy -> output.n_react++;
            
            /* allocate memory for reaction/coefficients pointers */
            dummy -> input.react = (double **) malloc( dummy -> input.n_react * sizeof ( double *) );
            
            dummy -> input.coeff = (double *) malloc( dummy -> input.n_react * sizeof ( double ) );
            
            /* allocate memory for reaction/coefficients pointers */
            dummy -> output.react = (double **) malloc( dummy -> output.n_react * sizeof ( double *) );
            
            dummy -> output.coeff = (double *) malloc( dummy -> output.n_react * sizeof ( double ) );
            
            
            /* get back to the beginning of line and assign now reactions & coefficients to dummy */
            
            s1 = colon;
            
            /* reinitialise the column index*/
            n = 0;
            
            /* initialise pointers to reactions and coefficients */
            react_in = dummy -> input.react;
            
            react_out = dummy -> output.react;
            
            coeff_in = dummy -> input.coeff;
            
            coeff_out = dummy -> output.coeff;
            
            /* read again per column */
            while ( strchr(s1, ' ') != NULL ){
                
                /* get to the next space */
                s3 = strchr(s1, ' ');
                
                /* stop the string there */
                *s3 = '\0';
                
                /* if even column, get reaction index*/
                if (n %2 == 0) which = atoi(s1);
                
                /* get a potential stoichiometric coefficient */
                c = atof(s1);
                
                /* if negative and column odd */
                if (c < 0. && n %2 != 0) {
                    
                    /* associate reaction output pointer to the corresponding reaction (via index which) */
                    *react_in = s + which - 1;
                    
                    /* associate coefficient */
                    *coeff_in = -c;
                    
                    /* increase pointer to input reactions */
                    react_in++;
                    
                    /* increase pointer to input coefficients */
                    coeff_in++;
                    
                }
                
                /* if positive and column odd */
                if (c > 0. && n %2 != 0) {
                    
                    /* associate reaction output pointer to the corresponding reaction (via index which) */
                    *react_out = s + which - 1;
                    
                    /* associate coefficient */
                    *coeff_out = c;
                    
                    /* increase pointer to output reactions */
                    react_out++;
                    
                    /* increase pointer to output coefficients */
                    coeff_out++;
                    
                }
                
                /* undo the null sign */
                *s3 = ' ';
                
                /* keep reading */
                s1 = s3 + 1;
                
                /* increase the column index */
                n++;
            }
            
            /* for the last column */
            /* surely it cannot be a reaction index, so get the stoichiometric coefficient */
            c = atof(s1);
            
            /* if negative */
            if (c < 0. ) {
                
                /* associate inputs as above */
                *react_in = s + which - 1;
                
                *coeff_in = -c;
                
                react_in++;
                
                coeff_in++;
                
            }
            
            /* if positive */
            if (c > 0. ) {
                
                /* associate outputs as above */
                *react_out = s + which - 1;
                
                *coeff_out = c;
                
                react_out++;
                
                coeff_out++;
                
            }
            
            /* undo the null sign at the end of line */
            *s2 = '\n';
            
            /* keep track of everything in the log file */
            fprintf(outfile, "Done.\n"); fflush(outfile);
            
            /* increase pointer to metabolites */
            dummy++;
        }
        
        /* get to the next line */
        s1 = s2 +1;
    }

    return;
}

/* A function to allocate memory for the problem once the input file has been read */
/* the input data is stored in a structure called "file_wrapper" */
void alloc_system (file_wrapper *input_file, metabolite *mets, int Nmet, double *s, FILE *outfile){
    
    /* if input data is an adjacency list */
    if ( input_file -> filetype == 0 ) {
        
        /* record it on the log file */
        fprintf(outfile, "Allocating from adjacency list\n");
        
        /* use the corresponding function to allocate memory */
        alloc_from_filetype0( input_file -> file_content, mets, Nmet, s, outfile);
        
    }
    
    /* if input data is a stoichiometric matrix */
    else if ( input_file -> filetype == 1 ) {
        
        /* record it on the log file */
        fprintf(outfile, "Allocating from stoichiometric matrix\n");
        
        /* use the corresponding function to allocate memory */
        alloc_from_filetype1( input_file -> file_content, mets, Nmet, s, outfile);
    }
    
    /* if input data is a reaction list */
    else if ( input_file -> filetype == 2 ) {
        
        /* record it on the log file */
        fprintf(outfile, "Allocating from reaction list\n");
        
        /* use the corresponding function to allocate memory */
        alloc_from_filetype2 (mets, s, input_file -> parser, Nmet, outfile);
        
    }
    
    /* else... */
    else{
        
        /* the input file is unkown, flag it out */
        fprintf(outfile," UNKNWON INPUT FILE TYPE\n");
        
        fprintf(stderr," UNKNWON INPUT FILE TYPE\n");
        
        exit (EXIT_FAILURE);
    }
    
}

/* *** THIS FUNCTION IS NOT USED ANYMORE.... *** */
/* A function to build the system to work with */
void build_network (FILE *in_file, metabolite *metabs, double *s, int Nmet, FILE *log_file){
    
    metabolite *dummy;
    int i, cnt;
    double coeff;
    
    /* loop over meatbolites appearing as input */
    for (dummy = metabs; dummy < metabs + Nmet; dummy++){
        
        fprintf(log_file, "Reading metabolite %d as input... ", (int)(dummy - metabs) ); fflush(log_file);
        
        /* get # reacs. where metabolit *dummy appears as input */
        fscanf(in_file,"%d", &(dummy -> input.n_react) );
        
        /* allocate memory for reaction/coefficients pointers */
        dummy -> input.react = (double **) malloc( dummy -> input.n_react * sizeof ( double *) );
        
        dummy -> input.coeff = (double *) malloc( dummy -> input.n_react * sizeof ( double ) );
        
        /* for each reaction, get memory position and coefficient */
        for (cnt=0; cnt < dummy -> input.n_react; cnt++) {
            
            /* get reaction index & coefficient from file */
            fscanf(in_file,"%d %lf", &i, &coeff);
            
            /* assign pointer to reaction */
            *( dummy -> input.react + cnt ) = s + (i - 1);
            
            /* assign coefficient value */
            *( dummy -> input.coeff + cnt ) = coeff;
        }
        
        fprintf(log_file, "Done.\n"); fflush(log_file);
    }
    
    /* loop over over meatbolites appearing as output */
    for (dummy = metabs; dummy < metabs + Nmet; dummy++){
        
        fprintf(log_file, "Reading metabolite %d as output... ", (int)(dummy - metabs) ); fflush(log_file);
        
        /* get # reacs. where metabolit *dummy appears as output */
        fscanf(in_file,"%d", &(dummy -> output.n_react) );
        
        /* allocate memory for reaction/coefficients pointers */
        dummy -> output.react = (double **) malloc( dummy -> output.n_react * sizeof ( double *) );
        
        dummy -> output.coeff = (double *) malloc( dummy -> output.n_react * sizeof ( double ) );
        
        /* for each reaction, get memory position and coefficient */
        for (cnt=0; cnt < dummy -> output.n_react; cnt++) {
            
            /* get reaction index & coefficient from file */
            fscanf(in_file,"%d %lf", &i, &coeff);
            
            /* assign pointer to reaction */
            *( dummy -> output.react + cnt ) = s + (i - 1);
            
            /* assign coefficient value */
            *( dummy -> output.coeff + cnt ) = coeff;
        }
        
        fprintf(log_file, "Done.\n"); fflush(log_file);
    }
    
    fprintf(log_file, "\n\n");
    
}
