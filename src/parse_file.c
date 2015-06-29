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

#include "parse_file.h"

/********************************************************************************************/
/*                                                                                          */
/*          THE FUNCTIONS BELOW ARE USED TO DEAL WITH REACTION LIST FILES ONLY              */
/*                                                                                          */
/********************************************************************************************/

/* A function to free the parser structure */
/* this structure is used to store metabolites when reading a reaction list */
void parser_free ( metabolite_parse **met, int Nmet){
    
    metabolite_parse *dummy;
    
    /* loop over parsed metabolites */
    for (dummy = *met; dummy < (*met) + Nmet; dummy++) {
        
        /* free the name string */
        free( dummy -> name);
        
        /* for the inputs */
        /* free the index of reactions they're attached to */
        free( dummy -> input.which_r);
        
        /* free the coefficients they are attached to */
        free( dummy -> input.coeff);

        /* for the outputs */
        /* free the index of reactions they're attached to */
        free( dummy -> output.which_r);
        
        /* free the coefficients they are attached to */
        free( dummy -> output.coeff);
    }
    
    free(*met);
    
    return;
}


/* A function to check wheter a metabolite has been already encountered */
metabolite_parse *find_metabolite (metabolite_parse *raw_metabs, int n_mets, char *met_name){
    
    metabolite_parse *raw = raw_metabs;
    
    /* loop over parsed metabolites */
    while ( raw < raw_metabs + n_mets && strcmp(met_name, raw -> name) != 0 )  raw++;
    
    /* return a pointer: */
    /* if closer to raw_metabs than n_mets, the matbolite is already in the list */
    return raw;
}

/* A function to append a new metabolite to the metabolite_parse structure */
void append_new_met (char *met_name, metabolite_parse **raw_metabs, int n_mets, int *n_allowed) {
    
    char *dummy1, *dummy2;
    int name_l = strlen(met_name);
    
    /* if the # of metabolites is too large, realloc the structure */
    if ( n_mets >= *n_allowed ) {
        
        /* double the number of allowed metabolites */
        *n_allowed *= 2;
        
        /* realloc */
        *raw_metabs = (metabolite_parse *) realloc(*raw_metabs, *n_allowed * sizeof( metabolite_parse ) );
        
    }
    
    /* allocate a char array for the name string */
    ( (*raw_metabs) + n_mets )  -> name = (char *) malloc( (name_l + 1) * sizeof(char) );
    
    /* associate the name to the newly added metabolite */
    /* initialise a pointer metabolite_parse name string */
    dummy1 = ((*raw_metabs) + n_mets)  -> name;
    
    /* initialise a pointer to the name */
    dummy2 = met_name;
    
    /* copy one string onto the other */
    while ( *dummy2 != '\0'){
        
        *dummy1 = *dummy2;
        
        dummy1++;
        
        dummy2++;
    }
    
    /* finalise the copied string */
    *dummy1 = '\0';
    
    /* initialise the reaction attached to the newly found metabolite */
    ((*raw_metabs) + n_mets) -> input.n_react = 0;
    
    ((*raw_metabs) + n_mets) -> output.n_react = 0;
    
    /* initialise the number of reactions that can be associated to the metabolite */
    ((*raw_metabs) + n_mets) -> input.n_allowed = 10;
    
    ((*raw_metabs) + n_mets) -> output.n_allowed = 10;
    
    /* allocate memory for the reactions and the coefficients */
    ((*raw_metabs) + n_mets) -> input.which_r = (int *) malloc( ( ( (*raw_metabs) + n_mets ) -> input.n_allowed ) * sizeof( int) );
    
    ((*raw_metabs) + n_mets) -> output.which_r = (int *) malloc( ( ( (*raw_metabs) + n_mets ) -> output.n_allowed ) * sizeof( int) );
    
    ((*raw_metabs) + n_mets) -> input.coeff = (double *) malloc( ( ( (*raw_metabs) + n_mets ) -> input.n_allowed ) * sizeof( double ) );
    
    ((*raw_metabs) + n_mets) -> output.coeff = (double *) malloc( ( ( (*raw_metabs) + n_mets ) -> output.n_allowed ) * sizeof( double) );
    
    
}

/* A function to update the i/o coefficients & reactions */
/* when a metabolite is found in some new reaction       */
void update_io (adjacency_parse *io, int which, double coeff){
    
    /* if exceeding the number of allowed reactions, realloc the arrays */
    if ( io -> n_react >= io -> n_allowed ){
        
        /* double the number of allowed reactions */
        io -> n_allowed *= 2;
        
        /* realloc */
        io -> which_r = (int *) realloc (io -> which_r, (io -> n_allowed) * sizeof ( int) );
        
        io -> coeff = (double *) realloc (io -> coeff, (io -> n_allowed) * sizeof ( double) );
    }
    
    /* associate the reaction index to the reactions array */
    *(io -> which_r + io -> n_react) = which;
    
    /* associate the coefficient to the coefficient array */
    *(io -> coeff + io -> n_react) = coeff;
    
    /* increase the number of associated reactions */
    io -> n_react++;
    
}

/* A function to update the metabolite_parse structure */
/* when a metabolite is found in some new reaction     */
void update_metabolite (metabolite_parse *met, int react, double coeff, int io) {
    
    /* if io negative update the input information */
    if (io == -1) update_io ( &(met->input), react, coeff);
    
    /* if io negative update the output information */
    else update_io ( &(met->output), react, coeff);
}

/* A function to parse a string containing information on a metabolite */
void single_met_parse (char *parse1, char *parse2, int *n_allowed, metabolite_parse **raw_metabs, int *n_mets, int which_r, int io){
    
    char *parse3, *coeff_met_raw, *coeff_met, *met_name;
    double coeff;
    metabolite_parse *new_met;
    
    /* get the substring containing the metabolite name and the stoichiometric coefficient */
    coeff_met_raw = get_substring (parse1, parse2);
    
    coeff_met = coeff_met_raw;
    
    /* get rid of the trailing space */
    while ( *coeff_met == ' ') coeff_met ++;
    
    /* get rid of space at the end of string */
    while ( * ( coeff_met_raw + strlen(coeff_met_raw) - 1 ) == ' ') * ( coeff_met_raw + strlen(coeff_met_raw) - 1 ) = '\0';
    
    /* find a space in the middle of the string */
    parse3 = strchr(coeff_met, ' ');
    
    /* if there is a space, the string is "coeff_space_metabolite" */
    if (parse3 != NULL){
        
        /* the metabolite name is right after the space */
        met_name = parse3+1;
        
        /* get rid of trailing space */
        while ( *met_name == ' ') met_name ++;
        
        /* get rid of space at the end of name */
        while ( * ( met_name + strlen(met_name) - 1 ) == ' ') * ( met_name + strlen(met_name) - 1 ) = '\0';
        
        /* check if newly encountered metabolite */
        new_met = find_metabolite (*raw_metabs, *n_mets, met_name);
        
        /* if new metabolite, append it to the list and increase number of metabolites */
        if ( (int) (new_met - *raw_metabs) >= *n_mets) {
            
            /* append the metabolite */
            append_new_met (met_name, raw_metabs, *n_mets, n_allowed);
            
            /* get a pointer to the appended metabolite */
            new_met = (*raw_metabs) + (*n_mets);
            
            /* increase the number of metabolites */
            (*n_mets)++;
        }
        
        /* cut the string at the space position */
        *parse3 = '\0';
        
        /* get the coefficient value */
        coeff = atof(coeff_met);
        
        /* undo the truncation */
        *parse3 = ' ';
        
        /* update metabolite information */
        update_metabolite (new_met, which_r, coeff, io);
        
        
    }
    
    /* if there is no space, there is no coefficient, i.e. coefficient = 1 */
    else{
        
        /* get the metabolite's name from the string */
        met_name = coeff_met;
        
        /* get rid of trailing space */
        while ( *met_name == ' ') met_name ++;
        
        /* get rid of space at the end of name */
        while ( * ( met_name + strlen(met_name) - 1 ) == ' ') * ( met_name + strlen(met_name) - 1 ) = '\0';
        
        /* sinks/sources don't have products/substrates: add a metabolite only if its name is longer than 0 */
        if ( strlen(met_name) > 0 ){
            /* if so... */
            /* check whether it is a new metabolite*/
            new_met = find_metabolite (*raw_metabs, *n_mets, met_name);
            
            /* if new metabolite, append it to the list and increase number of metabolites */
            if ( (int) (new_met - *raw_metabs) >= *n_mets) {
                
                /* append the metabolite */
                append_new_met (met_name, raw_metabs, *n_mets, n_allowed);
                
                /* get a pointer to the appended metabolite */
                new_met = (*raw_metabs) + (*n_mets);
                
                /* increase the number of metabolites */
                (*n_mets)++;
            }
            
            /* coefficient can be 1 only */
            coeff = 1.;
            
            /* update metabolite information */
            update_metabolite (new_met, which_r, coeff, io);
        }
    }
    
    /* free the substring containing the metabolite name and the stoichiometric coefficient */
    free( coeff_met_raw );
}

/* A function to parse the substrates and products of a reaction, which separates metabolites via the "+" sign */
void parse_subs_prods ( char *react, int *n_allowed, metabolite_parse **raw_metabs, int *n_all_mets, int which_r, int io){
    
    char *parse1 = react, *parse2;
    
    /* if there is no "+" sign, there's one metabolite only */
    if ( strchr(parse1, '+') == NULL ) {
        
        /* get to the end of the string */
        parse2 = parse1 + strlen(parse1);
        
        /* parse the metabolite */
        single_met_parse (parse1, parse2, n_allowed, raw_metabs, n_all_mets, which_r, io);
        
        
    }
    
    /* otherwise, separate metabolites by the "+" */
    else {
        
        /* for all substrates (or products) */
        while ( strchr(parse1, '+') != NULL ){
            
            /* get the substring within 2 consecutive "+" */
            parse2 = strchr(parse1, '+');
            
            /* parse the metabolite */
            single_met_parse (parse1, parse2, n_allowed, raw_metabs, n_all_mets, which_r, io);
            
            /* get to the next metabolite */
            parse1 = parse2 +1;
            
            
        }
        
        /* get to the last metabolite */
        parse2 = parse1 + strlen(parse1);
        
        /* parse the metabolite */
        single_met_parse (parse1, parse2, n_allowed, raw_metabs, n_all_mets, which_r, io);
        
    }
    
    return;
}

/* A function to parse a reaction string */
void parse_react_line (char *line, int *n_allowed, metabolite_parse **raw_metabs, int *n_all_mets, int which_r){
    
    /* arrow1 is a pointer to the end of substrates, arrow2 to the beginning of products */
    char *subs, *prods, *colon = strchr(line, ':'), *arrow1 = strchr(line, '-'), *arrow2 = strchr(line, '>');
    
    /* ignore whatever comes before the ":" */
    if (colon == NULL ) colon = line;
    
    /* remove trailing space from substrates */
    if ( *( colon + 1) == ' ') colon +=1;
    
    /* remove trailing space from products */
    if ( *( arrow2 + 1) == ' ') arrow2 +=1;
    
    /* get string of substrates */
    subs = get_substring (colon + 1, arrow1);
    
    /* get string of proucts */
    prods = get_substring (arrow2 + 1, line + strlen( line ) );
    
    /* parse substrates */
    parse_subs_prods (subs, n_allowed, raw_metabs, n_all_mets, which_r, -1);
    
    /* parse products */
    parse_subs_prods (prods, n_allowed, raw_metabs, n_all_mets, which_r, +1);
    
    /* free substrates string */
    free (subs);
    
    /* free products string */
    free (prods);
    
    return;
    
}

/***********************************************************************************/
/*                                                                                 */
/*          END OF FUNCTIONS RELATED WITH REACTION LIST FILES ONLY                 */
/*                                                                                 */
/***********************************************************************************/

/* A function to get the (char) file size */
long get_file_size (FILE *in_stream){
    
    long input_file_size;
    
    /* get to the end of file */
    fseek(in_stream, 0, SEEK_END);
    
    /* get the file size at the end */
    input_file_size = ftell(in_stream);
    
    /* get back to the beginning of file */
    rewind(in_stream);
    
    /* return the chars in file */
    return input_file_size;
    
}

/* A file to get the number of columns in a line                     */
/* also, the max for even columns is stored in max_r if exceeding it */
int get_ncolumns ( char *line, int *max_r) {
    
    int n = 0, m;
    char *colon = strchr(line, ':'), *r, *parse1 = line, *parse2;
    
    /* ignore whatever comes before the ":" */
    if (colon != NULL) parse1 = colon + 1;
    
    /* remove trailing space */
    while ( *parse1 == ' ') parse1++;
    
    /* for all columns */
    while ( strchr(parse1, ' ') != NULL ){
        
        /* if column even */
        if (n%2 == 0){
            
            parse2 = strchr(parse1, ' ');
            
            /* get the column value */
            r = get_substring (parse1, parse2);
            
            /* turn it to int */
            m = atoi( r );
            
            /* if greater than max_r, store it */
            if (m > *max_r) *max_r = m;
        }
        
        /* increase the column counter */
        n++;
        
        /* get to next column */
        parse1 = strchr(parse1, ' ') + 1;
        
    }
    
    /* return the number of columns */
    return n;
    
}

/* A function that guessues the type of the input file and that stores its content in a file_wrapper struct */
int guess_file_type (char *file_content, metabolite_parse **raw_metabs, int initial_n, int *n_metabs, int *n_reacs){
    
    char *s1 = file_content, *s2, *line;
    int file_type0 = 0, file_type1 = 1, file_type2 = 0, n_columns0, n_columnsvar, max_react = -1, n_allowed = initial_n;
    
    int n_metabs0 = 0, n_metabs2 = 0, n_react0 = -1, n_react1= 0, n_react2 = 0;
    
    /* for all file lines */
    while ( strchr(s1, '\n') != NULL ){
        
        /* get to the end of line */
        s2 = strchr(s1, '\n');
        
        /* copy the line content */
        line = get_substring (s1, s2);
        
        /* only evaluate if not a comment */
        if ( *line != '#'){
            
            /* if the ">" sign in line, the file is probably a reaction list (i.e. featuring "-->")*/
            if ( strchr(s1, '>') != NULL ){
                
                /* filetype 2, reaction list */
                file_type2 = 1;
                
                /* parse the line as a reaction */
                parse_react_line (line, &n_allowed, raw_metabs, &n_metabs2, n_react2);
                
                /* increase the number of reactions for filetype 2*/
                n_react2 += 1;
            }
            
            else{
                
                /* otherwise, it may either be a matrix, or an adjacency list */
                if (s1 == file_content) {
                    
                    /* get the number of columns of the first line */
                    n_columns0 = get_ncolumns (line, &max_react);
                    
                    /* if a matrix, the number of columns equals the number of reactions */
                    n_react1 = n_columns0;
                }
                
                /* get the number of columns for every line */
                n_columnsvar = get_ncolumns (line, &n_react0);
                
                /* if the number of columns varies, then it cannot be a matrix */
                if (n_columnsvar != n_columns0) file_type1 = 0;
                
                /* each line of the matrix or the adjacency list is a metabolite, increase the number of metabolites */
                n_metabs0++;
            }
            
        }
        
        /* get to the next line */
        s1 = s2 + 1;
    }
    
    /* if not a reaction list, nor a matrix, then assume it is an adjacency list */
    if ( file_type1 == 0 && file_type2 == 0) {
        
        /* flag it is a filetype 0 */
        file_type0 = 1;
        
        /* record the number of reactions according to the adjacency list standard */
        *n_reacs = n_react0;
        
        /* record the number of metabolites according to the adjacency list standard */
        *n_metabs = n_metabs0;
        
    }
    
    /* otherwise, if it cannot be a reaction list */
    else if (file_type1 == 1 && file_type2 == 0){
        
        /* record the number of reactions according to the matrix standard */
        *n_reacs = n_react1;
        
        /* record the number of metabolites according to the matrix standard */
        *n_metabs = n_metabs0;
        
    }
    
    /* otherwise, if reaction list */
    else if (file_type2 == 1 ){
        
        /* record the number of reactions according to the reaction list standard */
        *n_reacs = n_react2;
        
        /* record the number of metabolites according to the reaction list standard */
        *n_metabs = n_metabs2;
        
    }
    
    /* otherwise unknown */
    else {
        fprintf(stderr," UNKNWON INPUT FILE TYPE\n");
        
        exit (EXIT_FAILURE);
    }
    
    /* return the filetype */
    return file_type1 + file_type2;
}










