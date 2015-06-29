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

#ifndef __ALLOC_SYSTEM_H__
#define __ALLOC_SYSTEM_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "metabolites.h"

#include "file_wrapper.h"

void alloc_from_filetype2 (metabolite *, double *, metabolite_parse *, int, FILE *);

void alloc_from_filetype1(char *, metabolite *, int, double *, FILE *);

void alloc_from_filetype0(char *, metabolite *, int, double *, FILE *);

void alloc_system (file_wrapper *, metabolite *, int, double *, FILE *);

void build_network (FILE *, metabolite *, double *, int, FILE *);

#endif

