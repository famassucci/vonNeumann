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

#ifndef __VN_IO_H__
#define __VN_IO_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "metabolites.h"

void print_matrix (metabolite *, int, double *, int, FILE *);

void print_adj_list (metabolite *, int, double *, FILE *);

void print_system (metabolite *, int, double *, FILE *);

void print_reactions (double *, int, metabolite *, int, FILE *);

void print_constraints (double *, metabolite *, int, double, FILE *);

void print_fluxes (double *, int, FILE *);

#endif
