/* moid4.cpp: computes MOID (Minimum Orbital Intersection Distance)
between two orbits.

Copyright (C) 2010, Project Pluto

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

#ifdef TEST_VERSION
#include <stdio.h>
#endif

#include <string.h>
#include <math.h>
#include "watdefs.h"
#include "comets.h"

#define PI 3.141592653589793238462643383279502884197399375105820

double find_moid( const ELEMENTS *elem1, const ELEMENTS *elem2); /* moid4.c */
void setup_planet_elem( ELEMENTS *elem, const int planet_idx,
                                          const double t_cen);   /* moid4.c */

static void fill_matrix( double mat[3][3], const ELEMENTS *elem)
{
   memcpy( mat[0], elem->perih_vec, 3 * sizeof( double));
   memcpy( mat[1], elem->sideways, 3 * sizeof( double));
         /* mat[2] is the cross-product of mat[0] & mat[1]: */
   mat[2][0] = mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1];
   mat[2][1] = mat[0][2] * mat[1][0] - mat[0][0] * mat[1][2];
   mat[2][2] = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
}

static void compute_posn_and_derivative( const ELEMENTS *elem,
            const double true_anom, const double matrix[3][3],
            double *posn, double *vel)
{
   const double cos_true_anom = cos( true_anom);
   const double sin_true_anom = sin( true_anom);
   const double denom = 1. + elem->ecc * cos_true_anom;
   const double true_r = elem->q * (1. + elem->ecc) / denom;
   const double x = true_r * cos_true_anom;
   const double y = true_r * sin_true_anom;
   const double dx_dtheta = -y / denom;
   const double dy_dtheta = (x + elem->ecc / denom) / denom;
   int i;

   for( i = 0; i < 3; i++)
      posn[i] = x * matrix[i][0] + y * matrix[i][1];
   if( vel)
      for( i = 0; i < 3; i++)
         vel[i] = dx_dtheta * matrix[i][0] + dy_dtheta * matrix[i][1];
}

#define dot_prod( a, b) (a[0] * b[0] + a[1] * b[1] + a[2] * b[2])

static void compute_improvement( const double *delta, const double *v1,
               const double *v2, double *d1, double *d2)
{
   const double b = 2. * dot_prod( delta, v1);
   const double c = 2. * dot_prod( delta, v2);
   const double d = 2. * dot_prod( v1, v2);
   const double e = dot_prod( v1, v1);
   const double f = dot_prod( v2, v2);

   *d1 = (d * c - 2. * f * b) / (4. * e * f - d * d);
// *d2 = (d * b - 2. * e * c) / (4. * e * f - d * d);
   *d2 = -(b + 2. * e * *d1) / d;
}

#ifdef TEST_VERSION
static double true_anomaly_to_eccentric( const double true_anom,
                                         const double ecc)
{
   const double r = (1. - ecc * ecc) / (1. + ecc * cos( true_anom));
   const double x = r * cos( true_anom) + ecc;
   const double y = r * sin( true_anom) / sqrt( 1. - ecc * ecc);
   const double ecc_anom = PI + atan2( -y, -x);

   return( ecc_anom);
}
#endif

#define N_STEPS 360

/* A simple,  but effective,  MOID-finder.  It starts by computing 3x3
orthonormal matrices mat1 and mat2 that correspond to the base vectors of
the two orbits (i.e.,  each has a vector pointing toward perihelion;
one pointing 90 degrees "ahead" in the orbit;  and one perpendicular to
the plane of the orbit.)  Multiplying one by the inverse of the other
gives the transformation matrix from the first orbit to the second.
Then,  we can work as if one orbit is in the xy plane with perihelion
toward the positive x-axis.   */

double find_moid( const ELEMENTS *elem1, const ELEMENTS *elem2)
{
   double mat1[3][3], mat2[3][3], xform_matrix[3][3];
   const double identity_matrix[3][3] = {
            { 1., 0., 0.},
            { 0., 1., 0.},
            { 0., 0., 1.} };
   double least_dist_squared = 10000.;
   int i, j;

   fill_matrix( mat1, elem1);
   fill_matrix( mat2, elem2);
   for( i = 0; i < 3; i++)
      for( j = 0; j < 3; j++)
         xform_matrix[i][j] = dot_prod( mat1[i], mat2[j]);

   for( i = 0; i < N_STEPS; i++)
      {
      double vect1[3], vect2[3], dist_squared = 0.;
      double deriv1[3], deriv2[3];
      double true_anomaly2 = 2. * PI * (double)i / (double)N_STEPS;
      double true_anomaly1, delta_true1, delta_true2;
      int loop_count = 0, solution_found = 0;

      do
         {
         compute_posn_and_derivative( elem2, true_anomaly2, xform_matrix,
                        vect2, deriv2);
         if( !loop_count)
            true_anomaly1 = atan2( vect2[1], vect2[0]);
         compute_posn_and_derivative( elem1, true_anomaly1, identity_matrix,
                        vect1, deriv1);
         for( j = 0; j < 3; j++)
            vect1[j] -= vect2[j];
         compute_improvement( vect1, deriv1, deriv2, &delta_true1, &delta_true2);
         true_anomaly1 += delta_true1;
         true_anomaly2 -= delta_true2;
         if( fabs( delta_true1) < 2. * PI / N_STEPS)
            if( fabs( delta_true2) < 2. * PI / N_STEPS)
               {
               for( j = 0; j < 3; j++)
                  vect1[j] += delta_true1 * deriv1[j] + delta_true2 * deriv2[j];
               solution_found = 1;
               }
         loop_count++;
         }
         while( solution_found && loop_count < 3);
      dist_squared = dot_prod( vect1, vect1);
      if( dist_squared < least_dist_squared)
         least_dist_squared = dist_squared;
#ifdef TEST_VERSION
      printf( "%3d %c%8.6lf%8.2lf%8.2lf%8.2lf%8.2lf%15lf%15lf\n", i,
                     (solution_found ? '*' : ' '),
                     sqrt( dot_prod( vect1, vect1)),
                     true_anomaly1 * 180. / PI,
                     true_anomaly2 * 180. / PI,
                     true_anomaly_to_eccentric( true_anomaly1, elem1->ecc) * 180. / PI,
                     true_anomaly_to_eccentric( true_anomaly2, elem2->ecc) * 180. / PI,
                     dot_prod( vect1, deriv1),
                     dot_prod( vect1, deriv2));
//    printf( "%3d%15lf%15lf%15lf%15lf%15lf\n", i, x, y,
//                      vect[0], vect[1], vect[2]);
#endif
      }
   return( sqrt( least_dist_squared));
}

#define GAUSS_K .01720209895
#define SOLAR_GM (GAUSS_K * GAUSS_K)

void setup_planet_elem( ELEMENTS *elem, const int planet_idx,
                                             const double t_cen)
{
/* Taken straight from http://ssd.jpl.nasa.gov/elem_planets.html */
/* Gives a, ecc, incl, Omega=asc node, omega=arg per, L=longit */

static const double planet_elem[9 * 6] = {
/* Merc */  0.38709893, .20563069,  7.00487,  48.33167,  77.45645, 252.25084,
/* Venu */  0.72333199, .00677323,  3.39471,  76.68069, 131.53298, 181.97973,
/* Eart */  1.00000011, .01671022,  0.00005, -11.26064, 102.94719, 100.46435,
/* Mars */  1.52366231, .09341233,  1.85061,  49.57854, 336.04084, 355.45332,
/* Jupi */  5.20336301, .04839266,  1.30530, 100.55615,  14.75385,  34.40438,
/* Satu */  9.53707032, .05415060,  2.48446, 113.71504,  92.43194,  49.94432,
/* Uran */ 19.19126393, .04716771,  0.76986,  74.22988, 170.96424, 313.23218,
/* Nept */ 30.06896348, .00858587,  1.76917, 131.72169,  44.97135, 304.88003,
/* Plut */ 39.48168677, .24880766, 17.14175, 110.30347, 224.06676, 238.92881};

static const double planet_elem_rate[9 * 6] = {
/* Merc */  0.00000066,  0.00002527, -23.51,   -446.30,   573.57, 538101628.29,
/* Venu */  0.00000092, -0.00004938,  -2.86,   -996.89,  -108.80, 210664136.06,
/* Eart */ -0.00000005, -0.00003804, -46.94, -18228.25,  1198.28, 129597740.63,
/* Mars */ -0.00007221,  0.00011902, -25.47,  -1020.19,  1560.78,  68905103.78,
/* Jupi */  0.00060737, -0.00012880,  -4.15,   1217.17,   839.93,  10925078.35,
/* Satu */ -0.00301530, -0.00036762,   6.11,  -1591.05, -1948.89,   4401052.95,
/* Uran */  0.00152025, -0.00019150,  -2.09,  -1681.40,  1312.56,   1542547.79,
/* Nept */ -0.00125196,  0.0000251,   -3.64,   -151.25,  -844.43,    786449.21,
/* Plut */ -0.00076912,  0.00006465,  11.07,    -37.33,  -132.25,    522747.90};
   const double *pdata = planet_elem + (planet_idx - 1) * 6;
   const double *rate_data = planet_elem_rate + (planet_idx - 1) * 6;
   double elem_array[6];
   int i;

   for( i = 0; i < 6; i++)
      if( i < 2)
         elem_array[i] = pdata[i] + rate_data[i] * t_cen;
      else
         elem_array[i] = (pdata[i] + rate_data[i] * t_cen / 3600.) * PI / 180.;
   memset( elem, 0, sizeof( ELEMENTS));
   elem->ecc = elem_array[1];
   elem->q = (1. - elem->ecc) * elem_array[0];
   elem->incl = elem_array[2];
   elem->asc_node = elem_array[3];
   elem->arg_per = elem_array[4] - elem_array[3];
      /* l = (100.46435 + (129597740.63 / 3600.) * t_cen) * PI / 180.; */
   derive_quantities( elem, SOLAR_GM);
}

#ifdef TEST_VERSION
#include <stdlib.h>

static void show_elements( const ELEMENTS *elem)
{
   printf( "q=%8.5lf e=%8.6lf i=%8.4lf asc_node=%8.4lf arg_per=%8.4lf\n",
               elem->q, elem->ecc,
               elem->incl * 180. / PI,
               elem->asc_node * 180. / PI,
               elem->arg_per * 180. / PI);
}

int main( const int argc, const char **argv)
{
   ELEMENTS elem, earth_elem;
   double t_cen = 0.06;

   memset( &elem, 0, sizeof( ELEMENTS));
   if( argc == 6)
      setup_planet_elem( &earth_elem, 3, t_cen);
   else
      {
      memset( &earth_elem, 0, sizeof( ELEMENTS));
      sscanf( argv[6], "%lf,%lf,%lf,%lf,%lf",
                  &earth_elem.q,
                  &earth_elem.ecc,
                  &earth_elem.incl,
                  &earth_elem.asc_node,
                  &earth_elem.arg_per);
      earth_elem.incl *= PI / 180.;
      earth_elem.asc_node *= PI / 180.;
      earth_elem.arg_per *= PI / 180.;
      }
   elem.q = atof( argv[1]);
   elem.ecc = atof( argv[2]);
   if( elem.q < 0.)             /* actually the semimajor axis was given; */
      elem.q *= elem.ecc - 1.;  /* cvt it to a perihelion distance */
   elem.incl = atof( argv[3]) * PI / 180.;
   elem.asc_node = atof( argv[4]) * PI / 180.;
   elem.arg_per = atof( argv[5]) * PI / 180.;
   derive_quantities( &elem, SOLAR_GM);
   derive_quantities( &earth_elem, SOLAR_GM);
   printf( "MOID = %lf\n", find_moid( &earth_elem, &elem));
   show_elements( &elem);
   if( argc != 6)
      show_elements( &earth_elem);
   return( 0);
}
#endif
