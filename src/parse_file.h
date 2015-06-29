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

#ifndef __PARSE_FILE_H__
#define __PARSE_FILE_H__

#include "substring.h"

/* a structure to store I/O metabolite/reactions adjacency lists*/
typedef struct{
    
    /* number of reactions attached */
    int n_react;
    
    /* number of reactions that cab potentially be stored */
    int n_allowed;
    
    /* an array of pointers to the value of the reaction */
    int *which_r;
    
    /* an array of stoichiometric coefficients */
    double *coeff;
    
}adjacency_parse;

/* a metabolite structure */
typedef struct{
    
    /*metabolite name*/
    char *name;
    
    /* input reactions*/
    adjacency_parse input;
    
    /*output reactions*/
    adjacency_parse output;
}metabolite_parse;

void parser_free ( metabolite_parse **, int);

metabolite_parse *find_metabolite (metabolite_parse *, int, char *);

void append_new_met (char *, metabolite_parse **, int, int *);

void update_io (adjacency_parse *, int, double);

void update_metabolite (metabolite_parse *, int, double, int);

void single_met_parse (char *, char *, int *, metabolite_parse **, int *, int, int);

void parse_subs_prods ( char *, int *, metabolite_parse **, int *, int, int);

void parse_react_line (char *, int *, metabolite_parse **, int *, int);

long get_file_size (FILE *);

int get_ncolumns ( char *line, int *max_r);

int guess_file_type (char *, metabolite_parse **, int, int *, int *);

#endif
