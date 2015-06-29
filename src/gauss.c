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

#include "gauss.h"

/* A function to get a normally distributed variable */
double gaussdev(){
    
    double   g1, g2, j, x1, x2, r_s;
        
        do {
            /* pick two random numbers from -1 to 1*/
            x1 = 2.0 * drand48() - 1.0;
            
            x2 = 2.0 * drand48() - 1.0;
            
            r_s = x1 * x1 + x2 * x2;
        }
        /* while they do not fall in the unit circle */
        while (r_s >= 1.0 || r_s == 0.0);

         /* get two normally distributed variables */

        j = sqrt(-2.0 * log(r_s) / r_s);
        
        g1 = x1 * j;

        g2 = x2 * j;
    
        return g1;
    
}
