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

#include "file_wrapper.h"

void file_wrapper_free (file_wrapper **file_data){

    free( (*file_data) -> file_content);

    if ( (*file_data) -> filetype == 2)
        
        parser_free ( &((*file_data) -> parser), (*file_data) -> Nmet );
}

file_wrapper *handle_input_file (char *filename, int initial_n, FILE *log_file){
    
    FILE *in_stream = fopen(filename, "rb");
    
    long file_size = get_file_size (in_stream);
    
    
    file_wrapper *input_data = (file_wrapper*) malloc(1*sizeof(file_wrapper));
    
    input_data -> Nmet = 0;
    input_data -> Nreact = 0;
    
    input_data -> parser = (metabolite_parse *) malloc( initial_n * sizeof( metabolite_parse ) );
    
    input_data -> file_content = (char *) malloc ( file_size * (sizeof(char) ) );
    
    fread(input_data -> file_content, sizeof(char), file_size, in_stream);
    
    fclose(in_stream);
    
    input_data -> filetype = guess_file_type (input_data -> file_content, &(input_data -> parser), initial_n, &(input_data -> Nmet), &(input_data -> Nreact));
    
    
    fprintf(log_file, "Data from filetype %d\n", input_data -> filetype);
    
    
    return input_data;
    
}

