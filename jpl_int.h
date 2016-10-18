/* jpl_int.cpp: internal definitions for JPL ephemeris functions

Copyright (C) 2011, Project Pluto

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA.    */

            /* A JPL binary ephemeris header contains five doubles and */
            /* (up to) 41 int32_t integers,  so:                          */
#define JPL_HEADER_SIZE (5 * sizeof( double) + 41 * sizeof( int32_t))

#pragma pack(1)

struct interpolation_info
   {
   double pc[18],vc[18], twot;
   int np, nv;
   };

struct jpl_eph_data {
   double ephem_start, ephem_end, ephem_step;
   int32_t ncon;
   double au;
   double emrat;
   int32_t ipt[13][3];
   int32_t ephemeris_version;
               /* This is the end of the file header.  Following are */
               /* items computed within my code.                     */
   int32_t kernel_size, recsize, ncoeff;
   int32_t swap_bytes;
   int32_t curr_cache_loc;
   double pvsun[6];
   double pvsun_t;
   double *cache;
   struct interpolation_info iinfo;
   FILE *ifile;
   };
#pragma pack()


