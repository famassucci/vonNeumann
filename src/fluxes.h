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

#ifndef __FLUXES_H__
#define __FLUXES_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "gauss.h"

void backup_fluxes ( double *, double *, int);

void restore_backup ( double *, double *, int);

void initialise_fluxes (double *, int, double **, int, double *);

void normalise_fluxes (double *, int, double **, int, double *);

#endif
