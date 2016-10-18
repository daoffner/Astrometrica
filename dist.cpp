#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*   The following code implements three methods for computing the
distance between two lat/lon points on the earth.  The first formula
uses a method due to H. Andoyer,  from _Annuaire du Bureau des
Longitudes pour 1950_, p. 145.  It takes the earth's flattening into
account,  resulting in a relative error of the order of the square
of the earth's flattening (i.e.,  about one part in 90000).  I got
the formula "second hand" from Jean Meeus' _Astronomical
Algorithms_,  pp 80-82.

   The second method just uses plain spherical trigonometry,  as if
the earth really was round.  This can cause an error of about .5
percent.  As written,  it also loses precision for very short
distances and those near the antipodes,  where cos_d is almost
exactly 1 or -1.  In some cases,  roundoff might push you beyond
those limits, causing a domain error.

   The third method is due to Thaddeus Vincenty,  and is documented at
http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf .  An implementation in
JavaScript is available at
http://www.movable-type.co.uk/scripts/latlong-vincenty.html , and at
http://www.ga.gov.au/geodesy/datums/vincenty_inverse.jsp .  The method
is supposed to be accurate to within .5 millimeters.  I've made some
modifications to the basic method,  as described below.
*/

double earth_dist( const double lat1, const double lon1,
                   const double lat2, const double lon2,
                   const double flattening)
{
   const double f = (lat1 + lat2) / 2.;
   const double g = (lat1 - lat2) / 2.;
   const double lambda = (lon1 - lon2) / 2.;
   const double sin_g = sin( g), cos_g = cos( g);
   const double sin_f = sin( f), cos_f = cos( f);
   const double sin_lambda = sin( lambda), cos_lambda = cos( lambda);
   const double sin_g2 = sin_g * sin_g, cos_g2 = cos_g * cos_g;
   const double sin_f2 = sin_f * sin_f, cos_f2 = cos_f * cos_f;
   const double sin_lambda2 = sin_lambda * sin_lambda;
   const double cos_lambda2 = cos_lambda * cos_lambda;
   const double s = sin_g2 * cos_lambda2 + cos_f2 * sin_lambda2;
   const double c = cos_g2 * cos_lambda2 + sin_f2 * sin_lambda2;
   const double omega = atan( sqrt( s / c));
   const double r = sqrt( s * c) / omega;
   const double d = 2. * omega;
   const double h1 = (3. * r - 1.) / (2. * c);
   const double h2 = (3. * r + 1.) / (2. * s);
   const double rval = d * (1. + flattening *
                       (h1 * sin_f2 * cos_g2 - h2 * cos_f2 * sin_g2));

   return( rval);
}

/* This implementation for the spherical earth is mathematically exact,
but because the acos( ) function loses precision for cos_d near to
+/- 1,  the precision of the distance is apt to be poor for points
that are either close together (cos_d near 1) or near the antipodes
(cos_d near -1).  If you want to evade this problem,  and/or want to
compute the bearing,  too,  see dist_pa.cpp.  (Contrary to popular
belief,  switching to haversines will _not_ cure all ills,  though
it helps.)                  */

double spherical_earth_dist( const double lat1, const double lon1,
                             const double lat2, const double lon2)
{
   const double cos_d = sin( lat1) * sin( lat2) +
                                cos( lat1) * cos( lat2) * cos( lon1 - lon2);

   return( acos( cos_d));
}

/* Some modifications were made in this implementation of the
Vincenty method:

  -- The "reduced latitudes",  u1 and u2,  are calculated slightly
differently to avoid domain errors at the poles:  with the original
formulation,  if a latitude happens to be _exactly_ +/-90 degrees,
you got a domain error in the tangent function.  If it rounded off to
a hair greater than +/-90 degrees, you got a reduced latitude at the
opposite pole.

  -- It only computes the forward azimuth.  That can be done almost
trivially using values computed inside the loop.  The backward
azimuth would have taken more work (not a lot more,  admittedly).
If you really want a backward azimuth, you can always call the
function again with points 1 and 2 reversed.  Or use:
   back_azimuth = atan2( cos_u1 * sin_lambda,
             cos_u1 * sin_u2 * cos_lambda - sin_u1 * cos_u2);

  -- The core of the code is basically a root-finding algorithm,
trying to find the value of lambda for which delta_lambda (in the
original formulation,  lambda - lambdaprime) is zero.  If we are
far from the antipodes,  the root-finding proceeds briskly with
the original method (and slightly faster with the secant method
implemented below).  Closer to the antipodes,  though,  delta_lambda
becomes a much more complicated function,  and finding the actual
root requires some more care.  The original method could fail to
converge,  or enter a loop between two values.  Also, near the
antipodes,  there can be two solutions,  one leading to a
maximum-distance geodesic and the other to a minimum-distance
geodesic.  We want the former (corresponding to lambda < 180),
but the original code could get the latter (lambda > 180).

   To fix this,  I added some bracketing code:  we know that
delta_lambda is positive for lambda=0 and negative for lambda=180,
so that's our initial bracket,  narrowed down each time we compute
a new delta_lambda.  If our next step would take us outside that
bracket (or if too many steps have occurred, indicating slow
convergence),  we do a bisection step.

   Near the antipodes,  delta_lambda becomes ill-behaved,  and we
almost always go to the bisection code.  Elsewhere,  the newly-added
secant method makes for slightly faster convergence (though
convergence away from the antipodes was already pretty good).

   Also,  _really_ near the antipodes,  delta_lambda becomes
so poorly behaved that with mere 64-bit floats,  we can't get
delta < tolerance. That's why we break out of the bisection
code if it's clear that lambda isn't actually changing.  (Which
isn't a problem,  because we've determined lambda to sufficient
accuracy.)

   -- All this narrowed the "trouble zone" to the area where the
latitudes are of opposite sides,  and the longitudes are within
(flattening*180) degrees of being 180 degrees apart.  To evade this,
if we're within antipode_tol of this line segment,  the code just
bumps the second latitude by antipode_tol,  resulting in a possible
two-millimeter error.  Even much of this could be eliminated,  by
"adjusting" the result by the latitude correction multiplied by the
cosine of the azimuth... something for another day.

*/

static const double big_a_poly( const double usquared)
{
   return( 1. + usquared *
       (4096 + usquared * (-768 + usquared * (320 - 175 * usquared))) / 16384.);
}

static const double big_b_poly( const double usquared)
{
   return( usquared *
       (256. + usquared * (-128. + usquared * (74 - 47 * usquared))) / 1024.);
}

static const double cvt_lat( const double b, const double lat)
{
   return( atan2( b * sin( lat), cos( lat)));;
}

double vincenty_earth_dist( const double lat1, const double lon1,
                            const double lat2, const double lon2,
                            const double flat, double *azimuth)
{
   const double b = 1 - flat;       /* if a = 1 */
            /* u1, u2 are "reduced latitudes" */
   const double u1 = cvt_lat( b, lat1);
   const double u2 = cvt_lat( b, lat2);
   const double cos_u1 = cos( u1), sin_u1 = sin( u1);
   const double cos_u2 = cos( u2), sin_u2 = sin( u2);
   const double pi = 3.14159265358979323846264338327950288419716939937510582;
   double lambda, delta_lambda = 2.;
   double prev_lambda[2], prev_delta[2];
   double low_lambda = 0., high_lambda = pi;
   double dlon = fmod( lon2 - lon1, 2. * pi);
   const double antipode_tol = 3e-10;  /* corresponds to 2 mm */
   const double tolerance = 1e-12;     /* corresponds to .006 mm */
   double temp1, temp2, temp3;
   double sin_o, cos_o, o, cos2_a, cos_2om, big_c;
   int iter = 0, flipped_azimuth = 0;
   const int max_iter = 100;

   if( dlon > pi)          /* Things are a hair simpler if 0 < dlon < 180. */
      dlon -= 2. * pi;     /* We note if the azimuths have to be flipped.  */
   if( dlon < 0.)
      {
      dlon = -dlon;
      flipped_azimuth = 1;
      }
                       /* If we're too close to a line near the antipodes, */
                       /* we try again slightly offset from that line:     */
   if( fabs( u1 + u2) < .9 * antipode_tol && dlon > pi * b)
      {
      printf( "Using lat %.10lf\n", (antipode_tol - lat1) * 180. / pi);
      return( vincenty_earth_dist( lat1, lon1, antipode_tol - lat1,
                                 lon2, flat, azimuth));
      }
   if( lat1 == lat2 && !dlon)      /* points are identical */
      {
      if( *azimuth)
         *azimuth = 0.;
      return( 0.);
      }
   lambda = dlon;
   while( fabs( delta_lambda) > tolerance && iter < max_iter)
      {
      const double sin_lambda = sin( lambda);
      const double cos_lambda = cos( lambda);
      double sin_a;

      temp1 = cos_u2 * sin_lambda;
      temp2 = cos_u1 * sin_u2 - sin_u1 * cos_u2 * cos_lambda;
      cos_o = sin_u1 * sin_u2 + cos_u1 * cos_u2 * cos_lambda;
      sin_o = sqrt( temp1 * temp1 + temp2 * temp2);
      o = atan2( sin_o, cos_o);
      sin_a = cos_u1 * cos_u2 * sin_lambda / sin_o;
//    printf( "temp1 %.14lf temp2 %.14lf  sin_o %.14lf cos_o %.14lf sin_a %.14lf\n",
//             temp1, temp2, sin_o, cos_o, sin_a);
      cos2_a = 1 - sin_a * sin_a;  /* trig identity */
      if( !cos2_a)
         cos_2om = 0.;
      else
         cos_2om = cos_o - 2 * sin_u1 * sin_u2 / cos2_a;
      big_c = flat * cos2_a * (4. + flat * (4 - 3 * cos2_a)) / 16.;
      temp3 = 2. * cos_2om * cos_2om - 1.;
      delta_lambda = dlon - lambda + (1.-big_c) * flat * sin_a
                   * (o + big_c * sin_o * (cos_2om + big_c * cos_o * temp3));
      if( delta_lambda > 0.)     /* update our root brackets */
         low_lambda = lambda;
      else
         high_lambda = lambda;
      prev_lambda[1] = prev_lambda[0];
      prev_lambda[0] = lambda;       /* store previous lambda & delta_lambda */
      prev_delta[1] = prev_delta[0];   /* values for use in secant method */
      prev_delta[0] = delta_lambda;
//    printf( "Dist: %lf; lambda %.14lf; delta %.14lf\n", o * 180. / pi,
//              lambda * 180. / pi, delta_lambda * 180. / pi);
      if( !iter)
         lambda += delta_lambda;
      else
         {
         const double dx = prev_lambda[1] - prev_lambda[0];
         const double dy = prev_delta[1] - prev_delta[0];
         int auto_bisect;

//       printf( "delta_lambda = %.14lf; dx = %.14lf; dy = %.14lf\n",
//                               delta_lambda, dx, dy);
         if( !dy)          /* secant method would divide-by-zero */
            lambda += dy;
         else
            lambda -= delta_lambda * dx / dy;
                  /* following basically means "if the proposed lambda    */
                  /* lies outside the bracket,  or if we've done too many */
                  /* steps and obviously aren't converging briskly,  then*/
                  /* skip the secant steps and do a bisection step."    */
         auto_bisect = (iter > 10 && iter % 1);
         if( lambda <= low_lambda || lambda >= high_lambda || auto_bisect)
            {
            printf( "Gone outside brackets (%.14lf)!  Doing a bisection\n",
                           lambda * 180. / pi);
            printf( "%.14lf to %.14lf (%.7lg)\n",
                  low_lambda * 180. / pi, high_lambda * 180. / pi,
                  (high_lambda - low_lambda) * 180. / pi);
            lambda = low_lambda + (high_lambda - low_lambda) / 2.;
            if( lambda == prev_lambda[0])     /* if no change occurs  , */
               {                              /* we've  "maxed out" and */
               printf( "NO ACTUAL CHANGE\n"); /* should leave the loop  */
               delta_lambda = 0.;
               }
            }
         printf( "Iter %d: revised lambda %.14lf\n", iter, lambda * 180. / pi);
         }
      iter++;
      if( lambda > pi)
         lambda = pi;
      }

   const double usquared = cos2_a * (1.-b*b) / (b*b);
   const double big_b = big_b_poly( usquared);
   const double delta_o = big_b * sin_o *
              (cos_2om + (big_b / 4.) * (cos_o * temp3
            - (big_b / 6.) * cos_2om * (-3 + 4 * sin_o * sin_o)
            * (-3. + 4. * cos_o * cos_o)));
   const double rval = b * big_a_poly( usquared) * (o - delta_o);
   if( azimuth)
      {
      *azimuth = atan2( temp1, temp2);
      if( *azimuth < 0.)
         *azimuth += 2. * pi;
      if( flipped_azimuth)
         *azimuth = 2. * pi - *azimuth;
      }
   return( rval);
}

void vincenty_direct( const double lat1, const double lon1,
                            double *lat2, double *lon2,
                            const double flat, const double azimuth,
                            const double dist)
{
   const double pi = 3.14159265358979323846264338327950288419716939937510582;
   const double b = 1 - flat;       /* if a = 1 */
            /* u1, u2 are "reduced latitudes" */
   const double u1 = cvt_lat( b, lat1);
   const double cos_u1 = cos( u1), sin_u1 = sin( u1);
   const double sin_az = sin( azimuth), cos_az = cos( azimuth);
   const double sigma1 = atan2( sin_u1, cos_az * cos_u1);
   const double sin_alpha = cos_u1 * sin_az;
   const double cos2_alpha = 1. - sin_alpha * sin_alpha;
   const double u_squared = cos2_alpha * (1 - b * b) / (b * b);
   const double big_a = big_a_poly( u_squared);
   const double big_b = big_b_poly( u_squared);
   double sigma = dist / (b * big_a);
   double sigma_prime = 2. * pi;
   double cos_2sigma_m, sin_sigma, cos_sigma;
   const double tolerance = 1e-12;     /* corresponds to .006 mm */

   while( fabs( sigma - sigma_prime) > tolerance)
      {
      double delta_sigma, cos_2sigma_m2;

      cos_2sigma_m = cos( 2. * sigma1 + sigma);
      cos_2sigma_m2 = cos_2sigma_m * cos_2sigma_m;
      sin_sigma = sin( sigma);
      delta_sigma = big_b * sin_sigma *
               (cos_2sigma_m + (big_b / 4.) * (cos( sigma) *
               (2. * cos_2sigma_m2 - 1.) - (big_b / 6.)
               * cos_2sigma_m * (4. * sin_sigma * sin_sigma - 3.) *
               (4. * cos_2sigma_m2 - 3.)));

      sigma_prime = sigma;
      sigma = dist / (b * big_a) + delta_sigma;
      printf( "sigma = %.13lf\n", sigma);
      }

   sin_sigma = sin( sigma);
   cos_sigma = cos( sigma);
   cos_2sigma_m = cos( 2. * sigma1 + sigma);
   const double tval = sin_u1 * sin_sigma - cos_u1 * cos_sigma * cos_az;
   *lat2 = atan2( sin_u1 * cos_sigma + cos_u1 * sin_sigma * cos_az,
            b * sqrt( sin_alpha * sin_alpha + tval * tval));
   double lambda = atan2( sin_sigma * sin_az,
            cos_u1 * cos_sigma - sin_u1 * sin_sigma * cos_az);
// double big_c = flat / (16. * cos2_alpha * (4. + flat * (4. - 3 * cos2_alpha)));
   double big_c = (flat / 16) * cos2_alpha * (4. + flat * (4. - 3 * cos2_alpha));
   double big_l = lambda - (1. - big_c) * flat * sin_alpha *
         (sigma + big_c * sin_sigma * (cos_2sigma_m + big_c * cos_sigma *
         (2. * cos_2sigma_m * cos_2sigma_m - 1.)));

   *lon2 = lon1 + big_l;
}

/* For the 'spherical earth' formula,  we use a radius of 6371 km (close to
the average of the polar and equatorial radii.)  For the 'flattened earth',
the equatorial radius of 6378.140 km is used. */

int main( const int argc, const char **argv)
{
   const double pi = 3.14159265358979323846264338327950288419716939937510582;
   const double lat1 = atof( argv[1]) * pi / 180.;
   const double lon1 = atof( argv[2]) * pi / 180.;
   double lat2 = atof( argv[3]) * pi / 180.;
   double lon2 = atof( argv[4]) * pi / 180.;
   double azimuth, dist;
   const double flattening = 1. / 298.257223563;
   const double semimajor = 6378.137;
// const double flattening = 1. / 298.257;
// const double semimajor = 6378.14;

   double d1 = spherical_earth_dist( lat1, lon1, lat2, lon2);
   double d2 =           earth_dist( lat1, lon1, lat2, lon2, flattening);

   printf( "Distance from 'round earth': %lf\n", d1 * 6371.);
   printf( "Distance from H. Andoyer: %lf\n", d2 * semimajor);
   dist = vincenty_earth_dist( lat1, lon1, lat2, lon2, flattening, &azimuth);
   printf( "Distance from Vincenty: %.6lf\n", semimajor * dist);
   printf( "Azimuth %.10lf\n", azimuth * 180. / pi);

   lat2 = lon2 = -1.;
   vincenty_direct( lat1, lon1, &lat2, &lon2, flattening, azimuth, dist);
   printf( "Vincenty direct: %.10lf %.10lf\n",
            lat2 * 180. / pi, lon2 * 180. / pi);
   return( 0);
}
