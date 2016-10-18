/* ephem0.cpp: low-level funcs for ephemerides & pseudo-MPECs

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

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
#include "watdefs.h"
#include "afuncs.h"
#include "lunar.h"
#include "date.h"
#include "comets.h"
#include "mpc_obs.h"

#define J2000 2451545.0
#define EARTH_MAJOR_AXIS 6378140.
#define EARTH_MINOR_AXIS 6356755.
#define EARTH_AXIS_RATIO (EARTH_MINOR_AXIS / EARTH_MAJOR_AXIS)
#define EARTH_MAJOR_AXIS_IN_AU (EARTH_MAJOR_AXIS / (1000. * AU_IN_KM))
#define PI 3.1415926535897932384626433832795028841971693993751058209749445923
#define LOG_10  2.302585
#define LIGHT_YEAR_IN_KM    (365.25 * 86400. * 299792.456)

#define ROB_MATSON_TEST_CODE     1

char *fgets_trimmed( char *buff, size_t max_bytes, FILE *ifile); /*elem_out.c*/
int parallax_to_lat_alt( const double rho_cos_phi, const double rho_sin_phi,
       double *lat, double *ht_in_meters, const int planet_idx); /* ephem0.c */
double calc_obs_magnitude( const int is_comet, const double obj_sun,
          const double obj_earth, const double earth_sun, double *phase_ang);
int lat_alt_to_parallax( const double lat, const double ht_in_meters,
             double *rho_cos_phi, double *rho_sin_phi, const int planet_idx);
int write_residuals_to_file( const char *filename, const char *ast_filename,
          const int n_obs, const OBSERVE FAR *obs_data, const int format);
void put_observer_data_in_text( const char FAR *mpc_code, char *buff);
int compute_observer_loc( const double jde, const int planet_no,
               const double rho_cos_phi,                    /* mpc_obs.cpp */
               const double rho_sin_phi, const double lon, double FAR *offset);
int compute_observer_vel( const double jde, const int planet_no,
               const double rho_cos_phi,                    /* mpc_obs.cpp */
               const double rho_sin_phi, const double lon, double FAR *offset);
int make_pseudo_mpec( const char *mpec_filename, const char *obj_name);
                                              /* ephem0.cpp */
void remove_trailing_cr_lf( char *buff);      /* ephem0.cpp */
void create_obs_file( const OBSERVE FAR *obs, int n_obs, const int append);
const char *get_environment_ptr( const char *env_ptr);     /* mpc_obs.cpp */
void set_environment_ptr( const char *env_ptr, const char *new_value);
int text_search_and_replace( char FAR *str, const char *oldstr,
                                     const char *newstr);   /* ephem0.cpp */
void format_dist_in_buff( char *buff, const double dist_in_au); /* ephem0.c */
char int_to_mutant_hex_char( const int ival);              /* mpc_obs.c */
int debug_printf( const char *format, ...);                /* mpc_obs.c */
double planet_axis_ratio( const int planet_idx);            /* collide.cpp */
double planet_radius_in_meters( const int planet_idx);      /* collide.cpp */
double mag_band_shift( const char mag_band);                /* elem_out.c */

const char *residual_filename = "residual.txt";
const char *ephemeris_filename = "ephemeri.txt";
const char *elements_filename = "elements.txt";

int lat_alt_to_parallax( const double lat, const double ht_in_meters,
            double *rho_cos_phi, double *rho_sin_phi, const int planet_idx)
{
   const double axis_ratio = planet_axis_ratio( planet_idx);
   const double major_axis_in_meters = planet_radius_in_meters( planet_idx);
   const double u = atan( sin( lat) * axis_ratio / cos( lat));

   *rho_sin_phi = axis_ratio * sin( u) +
                            (ht_in_meters / major_axis_in_meters) * sin( lat);
   *rho_cos_phi = cos( u) + (ht_in_meters / major_axis_in_meters) * cos( lat);
   *rho_sin_phi *= major_axis_in_meters / AU_IN_METERS;
   *rho_cos_phi *= major_axis_in_meters / AU_IN_METERS;
   return( 0);
}

int parallax_to_lat_alt( const double rho_cos_phi, const double rho_sin_phi,
               double *lat, double *ht_in_meters, const int planet_idx)
{
   const double axis_ratio = planet_axis_ratio( planet_idx);
   const double tan_u = (rho_sin_phi / rho_cos_phi) / axis_ratio;
   const double radius = sqrt( rho_cos_phi * rho_cos_phi
                             + rho_sin_phi * rho_sin_phi);
   const double ang = atan( rho_sin_phi / rho_cos_phi);
   const double surface_radius = cos( ang) * cos( ang) +
                           sin( ang) * sin( ang) * axis_ratio;

   *ht_in_meters = (radius - surface_radius);
   *lat = atan( tan_u / axis_ratio);
   *ht_in_meters *= planet_radius_in_meters( planet_idx);
   return( 0);
}

/* format_dist_in_buff() formats the input distance (in AU) into a
seven-byte buffer.  It does this by choosing suitable units: kilometers
if the distance is less than a million km,  AU out to 10000 AU,  then
light-years.  In actual use,  light-years indicates some sort of error,
but I wanted the program to handle that problem without just crashing.

   28 Dec 2009:  To help puzzle out a strange bug,  I extended things so
that distances beyond 10000 light-years could be shown in kilo-light-years
(KLY),  mega,  giga,  etc. units,  making up prefixes after 'yocto'.  Thus,
one "ALY" = 10^78 light-years.  Since the universe is 13.7 billion years
old,  and nothing can be seen beyond that point,  such distances are
physically unreasonable by a factor of roughly 10^68.  But for debugging
purposes,  display of these distances was useful.  The bug is now fixed,
so with any luck,  nobody should see such distances again.

   2012 Feb 13:  there's some interest in determining orbits of _very_
close objects,  such as tennis balls.  To address this,  short distances
are now shown in millimeters,  centimeters,  meters,  or .1 km,
as appropriate.

NOTE:  It used to be that I followed MPC practice in showing all distances
between .01 and ten AU in the form d.dddd , but now,  that's only true for
one to ten AU.  For distances less than 1 AU,  I'm using .ddddd,  thereby
getting an extra digit displayed.   */

void format_dist_in_buff( char *buff, const double dist_in_au)
{
   if( dist_in_au < 0.)
      strcpy( buff, " <NEG!>");
   else
      {
      const double dist_in_km = dist_in_au * AU_IN_KM;
      const char *fmt;

                  /* for objects within a million km (about 2.5 times   */
                  /* the distance to the moon),  switch to km/m/cm/mm:  */
      if( dist_in_km < .0099)                 /* 0 to 9900 millimeters: */
         sprintf( buff, "%5.0lfmm", dist_in_km * 1e+6);    /* " NNNNmm" */
      else if( dist_in_km < .099)             /* 990 to 9900 centimeters: */
         sprintf( buff, "%5.0lfcm", dist_in_km * 1e+5);    /* " NNNNcm" */
      else if( dist_in_km < 99.)          /* 99 to 99000 meters: */
         sprintf( buff, "%6.0lfm", dist_in_km * 1e+3);     /* " NNNNNm" */
      else if( dist_in_km < 999.)         /* 99.0 to 999.9 kilometers: */
         sprintf( buff, "%6.1lfk", dist_in_km);            /* " NNN.Nk" */
      else if( dist_in_km < 999999.)      /* 999.9 to 999999 km: */
         sprintf( buff, "%7.0lf", dist_in_km);
      else if( dist_in_au > 9999.999)
         {
         double dist_in_light_years =
                               dist_in_au * AU_IN_KM / LIGHT_YEAR_IN_KM;

         if( dist_in_light_years > 9999.9)
            {
            const char *prefixes = "kMGTPEZYXWVUSRQONLJIHFDCBA";
            int idx = 0;

            dist_in_light_years /= 1000;
            for( idx = 0; prefixes[idx] && dist_in_light_years > 999.; idx++)
               dist_in_light_years /= 1000;
            if( !prefixes[idx])   /* can't represent this even in our */
               strcpy( buff, " <HUGE>");     /* largest made-up units */
            else
               {
               if( dist_in_light_years < 9.9)
                  sprintf( buff, "%4.1lfxLY", dist_in_light_years);
               else
                  sprintf( buff, "%4.0lfxLY", dist_in_light_years);
               buff[4] = prefixes[idx];
               }
            }
         else
            {
            if( dist_in_light_years > 99.999)  /* " 1234LY" */
               fmt = "%5.0lfLY";
            else if( dist_in_light_years > 9.999)   /* " 12.3LY" */
               fmt = "%5.1lfLY";
            else if( dist_in_light_years > .999)
               fmt = "%5.2lfLY";           /* " 1.23LY" */
            else
               fmt = "%5.3lfLY";           /* " .123LY" */
            sprintf( buff, fmt, dist_in_light_years);
            }
         }
      else
         {
         if( dist_in_au > 999.999)
            fmt = "%7.1lf";             /* " 1234.5" */
         else if( dist_in_au > 99.999)
            fmt = "%7.2lf";             /* " 123.45" */
         else if( dist_in_au > 9.999)
            fmt = "%7.3lf";             /* " 12.345" */
         else if( dist_in_au > .99)     /* used to be .01 */
            fmt = "%7.4lf";             /* " 1.2345" */
         else
            fmt = "%7.5lf";             /* " .12345" */
         sprintf( buff, fmt, dist_in_au);
         }
      *buff = ' ';                /* remove leading zero for small amts */
      }
}


static void format_velocity_in_buff( char *buff, const double vel)
{
   const char *format;

   if( vel < 9.999 && vel > -9.999)
      format = "%7.3lf";
   else if( vel < 99.999 && vel > -99.999)
      format = "%7.2lf";
   else if( vel < 999.9 && vel > -999.9)
      format = "%7.1lf";
   else
      format =  " !!!!!!";
   sprintf( buff, format, vel);
}

/* Rob Matson asked about having the program produce ECF (Earth-Centered
Fixed) coordinates,  in geometric lat/lon/altitude form.  I just hacked
it in,  then commented it out.  I'd have removed it,  but I'm thinking
this might be a useful option someday.  */

#ifdef ROB_MATSON_TEST_CODE
int find_lat_lon_alt( const double jd, const double *ivect,  /* collide.cpp */
                                       double *lat_lon_alt);
#endif

/* 'get_step_size' parses input text to get a step size in days,  so that */
/* '4h' becomes .16667 days,  '30m' becomes 1/48 day,  and '10s' becomes  */
/* 10/(24*60*60) days.  The units (days, hours, minutes,  or seconds) are */
/* returned in 'step_units' if the input pointer is non-NULL.  The number */
/* of digits necessary is returned in 'step_digits' if that pointer is    */
/* non-NULL.  Both are used to ensure correct time format in ephemeris    */
/* output;  that is,  if the step size is (say) .05d,  the output times   */
/* ought to be in a format such as '2009 Mar 8.34', two places in days.   */

double get_step_size( const char *stepsize, char *step_units, int *step_digits)
{
   double step = 0.;
   char units = 'd';

   if( sscanf( stepsize, "%lf%c", &step, &units) >= 1)
      if( step)
         {
         if( step_digits)
            {
            const char *tptr = strchr( stepsize, '.');

            *step_digits = 0;
            if( tptr)
               {
               tptr++;
               while( isdigit( *tptr++))
                  (*step_digits)++;
               }
#ifdef OBSOLETE_METHOD
            double tval = fabs( step);

            for( *step_digits = 0; tval < .999; (*step_digits)++)
               tval *= 10.;
#endif
            }
         units = tolower( units);
         if( step_units)
            *step_units = units;
         switch( units)
            {
            case 'd':
               break;
            case 'h':
               step /= 24.;
               break;
            case 'm':
               step /= 1440.;
               break;
            case 's':
               step /= 86400.;
               break;
            case 'w':
               step *= 7.;
               break;
            case 'y':
               step *= 365.25;
               break;
            }
         }
   return( step);
}

int find_precovery_plates( const char *filename, const double *orbit,
                           double epoch_jd)
{
   FILE *ofile, *ifile = fopen( "sky_cov.txt", "rb");
   double orbi[6];
   char buff[100];

   if( !ifile)
      return( -2);
   ofile = fopen( filename, "w");
   if( !ofile)
      {
      fclose( ifile);
      return( -1);
      }
   setvbuf( ofile, NULL, _IONBF, 0);
   memcpy( orbi, orbit, 6 * sizeof( double));
   while( fgets_trimmed( buff, sizeof( buff), ifile))
      {
      char tbuff[80];
      int line_no = 0;
      FILE *coverage_file = fopen( buff + 8, "rb");
      const long jd = dmy_to_day( 1, 1, atol( buff) / 1000, CALENDAR_GREGORIAN)
                                + atol( buff) % 1000;
      const double curr_jd = (double)jd + .5;
      double obs_posn[3], topo[3];
      double obj_ra, obj_dec;
      int i;

      integrate_orbit( orbi, epoch_jd, curr_jd);
      epoch_jd = curr_jd;
      compute_observer_loc( curr_jd, 3, 0., 0., 0., obs_posn);
      for( i = 0; i < 3; i++)
         topo[i] = orbi[i] - obs_posn[i];
      ecliptic_to_equatorial( topo);
      obj_ra = atan2( topo[1], topo[0]);
      obj_dec = asin( topo[2] / vector3_length( topo));

//    fprintf( ofile, "Looking at %s\n", buff);
//    fprintf( ofile, "JD %lf; %lf %lf\n", curr_jd, obj_ra * 180. / PI,
//             obj_dec * 180. / PI);
      if( coverage_file)
         {
         while( fgets( tbuff, sizeof( tbuff), coverage_file))
            {
            double ra_min, ra_max, dec_min, dec_max;

            line_no++;
            ra_min = ra_max = atof( tbuff);
            dec_min = dec_max = atof( tbuff + 10);
            for( i = 1; i < 4; i++)
               {
               double ra = atof( tbuff + i * 18 + 1);
               double dec = atof( tbuff + i * 18 + 10);

               while( ra - ra_min > PI)
                  ra -= PI + PI;
               while( ra - ra_max < -PI)
                  ra += PI + PI;
               if( ra_min > ra)
                  ra_min = ra;
               if( ra_max < ra)
                  ra_max = ra;
               if( dec_min > dec)
                  dec_min = dec;
               if( dec_max < dec)
                  dec_max = dec;
               }
            while( obj_ra - ra_min > PI)
               obj_ra -= PI + PI;
            while( obj_ra - ra_min < -PI)
               obj_ra += PI + PI;
            if (obj_ra > ra_min && obj_ra < ra_max && obj_dec > dec_min
                            && obj_dec < dec_max)
               fprintf( ofile, "%4d %s\n", line_no, buff + 8);
            }
         fclose( coverage_file);
         }
      }
   fclose( ifile);
   fclose( ofile);
   return( 0);
}

double vector_to_polar( double *lon, double *lat, const double *vector)
{
   const double r = vector3_length( vector);

   *lon = PI + atan2( -vector[1], -vector[0]);
   *lat =  asin( vector[2] / r);
   return( r);
}

static void format_motion( char *buff, const double motion)
{
   const char *motion_format;
   const double fabs_motion = fabs( motion);

   if( fabs_motion > 999999.)
      motion_format = "------";
   else if( fabs_motion > 999.)
      motion_format = "%6.0lf";
   else if( fabs_motion > 99.9)
      motion_format = "%6.1lf";
   else
      motion_format = "%6.2lf";
   sprintf( buff, motion_format, motion);
}

int ephemeris_in_a_file( const char *filename, const double *orbit,
         OBSERVE *obs, const int n_obs,
         const int planet_no,
         const double epoch_jd, const double jd_start, const char *stepsize,
         const double lon,
         const double rho_cos_phi, const double rho_sin_phi,
         const int n_steps, const char *note_text,
         const int options)
{
   double orbi[6], step;
   double prev_ephem_t = epoch_jd, prev_radial_vel = 0.;
   int i, j, date_format, hh_mm, n_step_digits;
   const int ephem_type = (options & 7);
   FILE *ofile;
   char step_units;
   const char *timescale = get_environment_ptr( "TT_EPHEMERIS");
   const char *override_date_format = get_environment_ptr( "DATE_FORMAT");
   const int show_alt_az = ((options & OPTION_ALT_AZ_OUTPUT)
              && rho_cos_phi && rho_sin_phi && ephem_type == OPTION_OBSERVABLES);
   double abs_mag = calc_absolute_magnitude( obs, n_obs);

   step = get_step_size( stepsize, &step_units, &n_step_digits);
   if( !step)
      return( -2);
   ofile = fopen( filename, "w");
   if( !ofile)
      return( -1);
   if( !abs_mag)
      abs_mag = atof( get_environment_ptr( "ABS_MAG"));
   memcpy( orbi, orbit, 6 * sizeof( double));
   setvbuf( ofile, NULL, _IONBF, 0);
   switch( step_units)
      {
      case 'd':
         hh_mm = 0;
         date_format = FULL_CTIME_DATE_ONLY;
         break;
      case 'h':
         hh_mm = 1;
         date_format = FULL_CTIME_FORMAT_HH;
         break;
      case 'm':
         hh_mm = 2;
         date_format = FULL_CTIME_FORMAT_HH_MM;
         break;
      case 's':
      default:
         hh_mm = 3;
         date_format = FULL_CTIME_FORMAT_SECONDS;
         break;
      }
   date_format |= FULL_CTIME_YEAR_FIRST | FULL_CTIME_MONTH_DAY
            | FULL_CTIME_MONTHS_AS_DIGITS | FULL_CTIME_LEADING_ZEROES
            | CALENDAR_JULIAN_GREGORIAN;
   date_format |= (n_step_digits << 4);
                     /* For ease of automated processing,  people may want */
   if( *override_date_format)     /* the time in some consistent format... */
      sscanf( override_date_format, "%x", &date_format);
   if( ephem_type == OPTION_STATE_VECTOR_OUTPUT ||
       ephem_type == OPTION_POSITION_OUTPUT ||
       ephem_type == OPTION_MPCORB_OUTPUT ||
       ephem_type == OPTION_8_LINE_OUTPUT)
      {
      timescale = "y";        /* force TT output */
      fprintf( ofile, "%.5lf %lf %d\n", jd_start, step, n_steps);
      }
   else if( ephem_type != OPTION_CLOSE_APPROACHES)
      {
      char hr_min_text[80];
      const char *pre_texts[4] = { "", " HH", " HH:MM", " HH:MM:SS" };

      strcpy( hr_min_text, pre_texts[hh_mm]);
      if( n_step_digits)
         {
         strcat( hr_min_text, ".");
         for( i = n_step_digits; i; i--)
            {
            char tbuff[2];

            tbuff[0] = step_units;
            tbuff[1] = '\0';
            strcat( hr_min_text, tbuff);
            }
         }
      if( note_text)
         fprintf( ofile, "#%s\n", note_text);
      fprintf( ofile, "Date (%s) %s   RA              ",
                     (*timescale ? "TT" : "UT"), hr_min_text);
      fprintf( ofile, "Dec         delta   r     elong ");
      if( options & OPTION_PHASE_ANGLE_OUTPUT)
         fprintf( ofile, " ph_ang  ");
      if( options & OPTION_PHASE_ANGLE_BISECTOR)
         fprintf( ofile, " ph_ang_bisector  ");
      if( options & OPTION_HELIO_ECLIPTIC)
         fprintf( ofile, " helio ecliptic   ");
      if( options & OPTION_TOPO_ECLIPTIC)
         fprintf( ofile, " topo ecliptic    ");
      if( abs_mag)
         fprintf( ofile, " mag");

      if( options & OPTION_MOTION_OUTPUT)
         fprintf( ofile, (options & OPTION_SEPARATE_MOTIONS) ? "  RA '/hr dec " : "  '/hr    PA  ");
      if( show_alt_az)
         fprintf( ofile, " alt  az");
      if( options & OPTION_RADIAL_VEL_OUTPUT)
         fprintf( ofile, "   rvel ");
      fprintf( ofile, "\n");

      for( i = 0; hr_min_text[i]; i++)
         if( hr_min_text[i] != ' ')
            hr_min_text[i] = '-';
      fprintf( ofile, "---- -- --%s  ------------   ",  hr_min_text);
      fprintf( ofile, "------------  ------ ------ ----- ");
      if( options & OPTION_PHASE_ANGLE_OUTPUT)
         fprintf( ofile, " ------  ");
      if( options & OPTION_PHASE_ANGLE_BISECTOR)
         fprintf( ofile, " ---------------  ");
      if( options & OPTION_HELIO_ECLIPTIC)
         fprintf( ofile, " ---------------  ");
      if( options & OPTION_TOPO_ECLIPTIC)
         fprintf( ofile, " ---------------  ");
      if( abs_mag)
         fprintf( ofile, " ---");
      if( options & OPTION_MOTION_OUTPUT)
         fprintf( ofile, " ------ ------");
      if( show_alt_az)
         fprintf( ofile, " --- ---");
      if( options & OPTION_RADIAL_VEL_OUTPUT)
         fprintf( ofile, " ------ ");
      fprintf( ofile, "\n");
      }

   for( i = 0; i < n_steps; i++)
      {
      double ephemeris_t, utc;
      double topo[3], r, obs_posn[3];
      double topo_ecliptic[3];
      double orbi_after_light_lag[3];
      double obs_vel[3], topo_vel[3];
#ifdef ROB_MATSON_TEST_CODE
      double lat_lon_alt[3];
#endif
      char buff[240];
      double curr_jd = jd_start + (double)i * step, delta_t;
      OBSERVE temp_obs;

      if( options & OPTION_ROUND_TO_NEAREST_STEP)
         curr_jd = floor( (curr_jd - .5) / step + .5) * step + .5;
      delta_t = td_minus_utc( curr_jd) / 86400.;
      if( *timescale)                     /* we want a TT ephemeris */
         {
         ephemeris_t = curr_jd;
         utc = curr_jd - delta_t;
         }
      else                                /* "standard" UTC ephemeris */
         {
         ephemeris_t = curr_jd + delta_t;
         utc = curr_jd;
         }
      integrate_orbit( orbi, prev_ephem_t, ephemeris_t);
      prev_ephem_t = ephemeris_t;
      compute_observer_loc( ephemeris_t, planet_no, rho_cos_phi, rho_sin_phi,
                              lon, obs_posn);
      compute_observer_vel( ephemeris_t, planet_no, rho_cos_phi, rho_sin_phi,
                              lon, obs_vel);

      for( j = 0; j < 3; j++)
         {
         topo[j] = orbi[j] - obs_posn[j];
         topo_vel[j] = orbi[j + 3] - obs_vel[j];
         }
      r = vector3_length( topo);
                 /* for "ordinary ephemeris" (not state vectors or */
                 /* orbital elements),  include light-time lag:    */
      if( ephem_type == OPTION_OBSERVABLES)
         for( j = 0; j < 3; j++)
            {
            const double diff = -orbi[j + 3] * r / AU_PER_DAY;

            orbi_after_light_lag[j] = orbi[j] + diff;
            topo[j] += diff;
            }

      memset( &temp_obs, 0, sizeof( OBSERVE));
      temp_obs.r = vector3_length( topo);
      for( j = 0; j < 3; j++)
         {
         temp_obs.vect[j] = topo[j] / r;
         temp_obs.obs_vel[j] = -topo_vel[j];
         }
      memcpy( temp_obs.obs_posn, obs_posn, 3 * sizeof( double));
      for( j = 0; j < 3; j++)
         temp_obs.obj_posn[j] = temp_obs.obs_posn[j] + topo[j];
//       temp_obs.obj_posn[j] = temp_obs.obs_posn[j] + r * temp_obs.vect[j];
                  /* rotate topo from ecliptic to equatorial */
      memcpy( topo_ecliptic, topo, 3 * sizeof( double));
      ecliptic_to_equatorial( topo);                           /* mpc_obs.cpp */
      ecliptic_to_equatorial( topo_vel);
      ecliptic_to_equatorial( obs_posn);

      if( ephem_type == OPTION_STATE_VECTOR_OUTPUT ||
            ephem_type == OPTION_POSITION_OUTPUT)
         {
         sprintf( buff, "%.5lf %15.10lf %15.10lf %15.10lf",
                           curr_jd, topo[0], topo[1], topo[2]);
         if( ephem_type == OPTION_STATE_VECTOR_OUTPUT)
            sprintf( buff + strlen( buff), "  %15.10lf %15.10lf %15.10lf",
                           topo_vel[0], topo_vel[1], topo_vel[2]);
         }
      else if( ephem_type == OPTION_8_LINE_OUTPUT
             || ephem_type == OPTION_MPCORB_OUTPUT)
         {
         FILE *ifile;

         write_out_elements_to_file( orbi, ephemeris_t, ephemeris_t,
                    obs, n_obs, "", 5, 0,
                    1 + (i == n_steps - 1 ? 0 : ELEM_OUT_NO_COMMENT_DATA));
                        /* only show comment data on last line */
         ifile = fopen( (ephem_type == OPTION_8_LINE_OUTPUT) ?
                  elements_filename : "mpc_fmt.txt", "rb");
         while( fgets( buff, sizeof( buff), ifile))
            fputs( buff, ofile);
         fclose( ifile);
         *buff = '\0';
         }
      else
         {
         DPT ra_dec;
         double ra, dec, sec, earth_r = 0.;
         double radial_vel = 0;
         int hr, min, deg, dec_sign = '+';
         char ra_buff[80], dec_buff[80], date_buff[80];
         char r_buff[10], solar_r_buff[10];
         double cos_elong, solar_r, v_dot_r = 0.;

         for( j = 0; j < 3; j++)
            v_dot_r += topo[j] * topo_vel[j];
         r = vector3_length( topo);
         radial_vel = v_dot_r / r;
         solar_r = vector3_length( orbi_after_light_lag);
         earth_r = vector3_length( obs_posn);
         cos_elong = r * r + earth_r * earth_r - solar_r * solar_r;
         cos_elong /= 2. * earth_r * r;

         ra_dec.x = atan2( topo[1], topo[0]);
         ra = ra_dec.x * 12. / PI;
         if( ra < 0.) ra += 24.;
         hr = (int)ra;
         min = (int)((ra - (double)hr) * 60.);
         sec =       (ra - (double)hr) * 3600. - 60. * (double)min;
         sprintf( ra_buff, "%02d %02d %6.3lf", hr, min, sec);
         if( ra_buff[6] == ' ')        /* leading zero */
            ra_buff[6] = '0';

         ra_dec.y = asin( topo[2] / r);
         dec = ra_dec.y * 180. / PI;
         if( dec < 0.)
            {
            dec = -dec;
            dec_sign = '-';
            }
         deg = (int)dec;
         min = (int)( (dec - (double)deg) * 60.);
         sec =        (dec - (double)deg) * 3600. - (double)min * 60.;
         sprintf( dec_buff, "%c%02d %02d %5.2lf", dec_sign, deg, min, sec);
         if( dec_buff[7] == ' ')        /* leading zero */
            dec_buff[7] = '0';

         full_ctime( date_buff, curr_jd, date_format);
         format_dist_in_buff( r_buff, r);
         format_dist_in_buff( solar_r_buff, solar_r);
         sprintf( buff, "%s  %s   %s %s%s %5.1lf",
               date_buff, ra_buff, dec_buff, r_buff, solar_r_buff,
               acose( cos_elong) * 180. / PI);
//       if( abs_mag)
            {
            char *endptr = buff + strlen( buff);
            double phase_ang;
            double curr_mag = abs_mag + calc_obs_magnitude( obs->flags & OBS_IS_COMET,
                       solar_r, r, earth_r, &phase_ang);  /* elem_out.cpp */

            if( curr_mag > 999.)       /* avoid overflow for objects     */
               curr_mag = 999.;        /* essentially at zero elongation */

            if( options & OPTION_PHASE_ANGLE_OUTPUT)
               {
               sprintf( endptr, " %8.4lf", phase_ang * 180. / PI);
               endptr += strlen( endptr);
               }

            if( options & OPTION_PHASE_ANGLE_BISECTOR)
               {
               double pab_vector[3], pab_lon, pab_lat;

               for( j = 0; j < 3; j++)
                  pab_vector[j] = topo_ecliptic[j] / r
                                       + orbi_after_light_lag[j] / solar_r;
#ifdef OBSOLETE_I_HOPE
               debug_printf( "%s: topo len %lf; %lf\n", date_buff,
                               r, vector3_length( topo_ecliptic));
               debug_printf( "        vect len %lf; %lf\n",
                              solar_r, vector3_length( orbi_after_light_lag));
               debug_printf( "   Ecliptic lat/lon from sun: %.5lf %.5lf\n",
                        180. + (180. / PI) * atan2( -orbi_after_light_lag[1], -orbi_after_light_lag[0]),
                        (180. / PI) * asin( orbi_after_light_lag[2] / solar_r));
               debug_printf( "   Ecliptic lat/lon from observer: %.5lf %.5lf\n",
                        180. + (180. / PI) * atan2( -topo_ecliptic[1], -topo_ecliptic[0]),
                        (180. / PI) * asin( topo_ecliptic[2] / r));
#endif
               vector_to_polar( &pab_lon, &pab_lat, pab_vector);
               sprintf( endptr, " %8.4lf %8.4lf", pab_lon * 180. / PI,
                                                  pab_lat * 180. / PI);
               endptr += strlen( endptr);
               }

            if( options & OPTION_HELIO_ECLIPTIC)
               {
               double eclip_lon, eclip_lat;

               vector_to_polar( &eclip_lon, &eclip_lat, orbi_after_light_lag);
               sprintf( endptr, " %8.4lf %8.4lf", eclip_lon * 180. / PI,
                                                  eclip_lat * 180. / PI);
               endptr += strlen( endptr);
               }

            if( options & OPTION_TOPO_ECLIPTIC)
               {
               double eclip_lon, eclip_lat;

               vector_to_polar( &eclip_lon, &eclip_lat, topo_ecliptic);
               sprintf( endptr, " %8.4lf %8.4lf", eclip_lon * 180. / PI,
                                                  eclip_lat * 180. / PI);
               endptr += strlen( endptr);
               }

            if( abs_mag)           /* don't show a mag if you dunno how bright */
               {                   /* the object really is! */
               if( curr_mag < 99)
                  sprintf( endptr, " %4.1lf", curr_mag + .05);
               else
                  sprintf( endptr, " %3d ", (int)( curr_mag + .5));
               }
            if( phase_ang > PI * 2. / 3.)    /* over 120 degrees */
               endptr[4] = '?';           /* signal doubtful magnitude */
            }
         if( options & OPTION_MOTION_OUTPUT)
            {
            MOTION_DETAILS m;
            char *end_ptr = buff + strlen( buff);

            compute_observation_motion_details( &temp_obs, &m);
            strcat( buff, " ");
            if( options & OPTION_SEPARATE_MOTIONS)
               {
               format_motion( end_ptr + 1, m.ra_motion);
               format_motion( end_ptr + 8, m.dec_motion);
               }
            else
               {
               format_motion( end_ptr + 1, m.total_motion);
               sprintf( end_ptr + 8, "%5.1lf ",
                               m.position_angle_of_motion);
               }
            end_ptr[7] = ' ';
            }
         if( show_alt_az)
            {
            DPT latlon, alt_az;
            double unused_ht_in_meters;

            ra_dec.x = -ra_dec.x;
            parallax_to_lat_alt( rho_cos_phi, rho_sin_phi, &latlon.y,
                                       &unused_ht_in_meters, planet_no);
            latlon.x = lon;
                        /* FIX someday:  this only works if planet_no == 3, */
                        /* i.e.,  topocentric ephemerides */
            full_ra_dec_to_alt_az( &ra_dec, &alt_az, NULL, &latlon, utc, NULL);
            alt_az.x += PI;
            while( alt_az.x < 0.)
               alt_az.x += PI + PI;
            while( alt_az.x > PI + PI)
               alt_az.x -= PI + PI;
            sprintf( buff + strlen( buff), " %c%02d %03d",
                                 (alt_az.y > 0. ? '+' : '-'),
                                 (int)( fabs( alt_az.y * 180. / PI) + .5),
                                 (int)( alt_az.x * 180. / PI + .5));
            }
         if( options & OPTION_RADIAL_VEL_OUTPUT)
            format_velocity_in_buff( buff + strlen( buff),
                                     radial_vel * AU_IN_KM / 86400.);
#ifdef ROB_MATSON_TEST_CODE
         if( options & OPTION_GROUND_TRACK)
            {
            find_lat_lon_alt( ephemeris_t, topo, lat_lon_alt);
            sprintf( buff + strlen( buff), "%9.4lf%9.4lf%10.3lf",
                  lat_lon_alt[0] * 180. / PI,
                  lat_lon_alt[1] * 180. / PI,
                  lat_lon_alt[2] * AU_IN_KM);
            }
#endif
         if( options & OPTION_CLOSE_APPROACHES)
            {
            if( (step > 0. && radial_vel >= 0. && prev_radial_vel < 0.)
               || (step < 0. && radial_vel <= 0. && prev_radial_vel > 0.))
               {
               double delta_t, v_squared = 0.;

               for( j = 0; j < 3; j++)
                  v_squared += topo_vel[j] * topo_vel[j];
               delta_t = -v_dot_r / v_squared;
//             full_ctime( date_buff, curr_jd, 0);
//             debug_printf( "At %s: delta_t = %lf, radial_vel = %lf, %lf\n",
//                         date_buff, delta_t, radial_vel, prev_radial_vel);
               full_ctime( date_buff, curr_jd + delta_t,
                        FULL_CTIME_FORMAT_HH_MM
                      | FULL_CTIME_YEAR_FIRST | FULL_CTIME_MONTH_DAY
                      | FULL_CTIME_MONTHS_AS_DIGITS
                      | FULL_CTIME_LEADING_ZEROES);
               sprintf( buff, "Close approach at %s: ", date_buff);
               for( j = 0; j < 3; j++)
                  topo[j] += delta_t * topo_vel[j];
               format_dist_in_buff( buff + strlen( buff),
                        vector3_length( topo));
               }
            else        /* suppress output */
               *buff = '\0';
            }
#ifdef NOT_IN_USE
         if( options & OPTION_SPACE_VEL_OUTPUT)
            {
                     /* get 'full' velocity; cvt AU/day to km/sec: */
            const double total_vel =
                       vector3_length( topo_vel) * AU_IN_KM / 86400.;

            format_velocity_in_buff( buff + strlen( buff), total_vel);

            }
#endif
         prev_radial_vel = radial_vel;
         }

      if( *buff)
         fprintf( ofile, "%s\n", buff);
      }
   fclose( ofile);
   return( 0);
}

static void output_angle_to_buff( char *obuff, const double angle,
                               const int precision)
{
   const int ihr = (int)angle;
   const double min = (angle - (double)ihr) * 60.;
   const int imin = (int)min;
   const double sec = (min - (double)imin) * 60.;

   if( precision >= 100)         /* decimal quantity */
      {
      char format_buff[10];
      int i;

      sprintf( format_buff, "%%%d.0%dlf\t", (precision - 100) + 3,
                                           (precision - 100));

      sprintf( obuff, format_buff, angle);
      if( *obuff == ' ')
         *obuff = '0';
      for( i = strlen( obuff); i < 12; i++)
         obuff[i] = ' ';
      obuff[12] = '\0';
      return;
      }

   sprintf( obuff, "%02d\t", ihr);
   switch( precision)
      {
      case -1:       /* hh mm,  integer minutes */
      case -2:       /* hh mm.m,  tenths of minutes */
      case -3:       /* hh mm.mm,  hundredths of minutes */
      case -4:       /* hh mm.mmm,  milliminutes */
         {
         static const char *format_text[4] = {
                              "%2.0lf\t      ",
                              "%4.1lf     ",
                              "%5.2lf    ",
                              "%6.3lf   " };

         sprintf( obuff + 3, format_text[ -1 - precision], min);
         }
         break;
      case 0:        /* hh mm ss,  integer seconds */
      case 1:        /* hh mm ss.s,  tenths of seconds */
      case 2:        /* hh mm ss.ss,  hundredths of seconds */
      case 3:        /* hh mm ss.sss,  thousands of seconds */
         {
         static const char *format_text[4] = {
                  "%02d\t%2.0lf    ",
                  "%02d\t%4.1lf  ",
                  "%02d\t%5.2lf ",
                  "%02d\t%6.3lf"   };

         sprintf( obuff + 3, format_text[precision], imin, sec);
         if( obuff[6] == ' ')
            obuff[6] = '0';
         }
         break;
      }
   if( obuff[3] == ' ')
      obuff[3] = '0';
}

/* 'put_residual_into_text( )' expresses a residual,  from 0 to 180 degrees, */
/* such that the text starts with a space,  followed by four characters,   */
/* and a sign:  six bytes in all.  This can be in forms such as:           */
/*                                                                         */
/*  179d-    (for values above 999 arcminutes = 16.6 degrees)              */
/*  314'+    (for values above 9999 arcsec but below 999 arcminutes)       */
/*  7821-    (for values below 9999 arcsec but above 99 arcsec)            */
/*  12.3+    (for values above one arcsec but below 99 arcsec)             */
/*  2.71-    (for values between one and ten arcsec,  if THREE_DIGITS=1)   */
/*   .87-    (for values under an arcsecond, usually)                      */
/*  .871-    (for values under an arcsecond, if THREE_DIGITS=1)            */

int precise_residuals = 0;

static void put_residual_into_text( char *text, const double resid)
{
   const double zval = fabs( resid);

   if( zval > 35999.9)             /* >999': show integer degrees */
      sprintf( text, "%4.0lfd", zval / 3600.);
   else if( zval > 9999.9)              /* 999' > x > 9999': show ###' arcmin */
      sprintf( text, "%4.0lf'", zval / 60.);
   else if( zval > 99.9)
      sprintf( text, "%5.0lf", zval);
   else if( zval > .99 && zval < 9.99 && precise_residuals)
      sprintf( text, "%5.2lf", zval);
   else if( zval > .99)
      sprintf( text, "%5.1lf", zval);
   else
      {
      sprintf( text, (precise_residuals? "%5.3lf" : "%5.2lf"), zval);
      text[1 - precise_residuals] = ' ';
      }
   if( !atof( text))
      text[5] = ' ';
   else
      text[5] = (resid > 0. ? '+' : '-');
   text[6] = '\0';
}

static void put_mag_resid( char *output_text, const double obs_mag,
                           const double computed_mag, const char mag_band)
{
   if( obs_mag && computed_mag)
      sprintf( output_text, "%6.2lf ",
               obs_mag - computed_mag);
//             obs_mag - computed_mag - mag_band_shift( mag_band);
   else
      strcpy( output_text, "------ ");
}

/* format_observation( ) takes an observation and produces text for it,
   suitable for display on a console (findorb) or in a Windoze scroll
   box (FIND_ORB),  or for writing to a file.  */

void format_observation( const OBSERVE FAR *obs, char *text,
                                        const int resid_format)
{
   double angle;
   char xresid[30], yresid[30];
   int i;
   const int base_format = (resid_format & 3);
   const int four_digit_years =
                    (resid_format & RESIDUAL_FORMAT_FOUR_DIGIT_YEARS);
   int month;
   long year;
   double day, utc;
   char *original_text_ptr = text;
   MOTION_DETAILS m;

   utc = obs->jd - td_minus_utc( obs->jd) / 86400.;
   day = decimal_day_to_dmy( utc, &year, &month, CALENDAR_JULIAN_GREGORIAN);

   if( base_format != RESIDUAL_FORMAT_SHORT)
      {
      const char *date_format_text[7] = { "%2.0lf       \t",
                                          "%4.1lf     \t",
                                          "%5.2lf    \t",
                                          "%6.3lf   \t",
                                          "%7.4lf  \t",
                                          "%8.5lf \t",
                                          "%9.6lf\t" };

      if( four_digit_years)
         sprintf( text, "%04ld\t%02d\t", year, month);
      else
         sprintf( text, "%02d\t%02d\t", abs( year % 100), month);
      text += strlen( text);
      if( resid_format & RESIDUAL_FORMAT_HMS)
         {
         const long seconds = (long)( day * 86400. + .001);

         sprintf( text, "%02ld %02ld:%02ld:%02ld \t", seconds / 86400,
                  (seconds / 3600) % 24, (seconds / 60) % 60, seconds % 60);
         }
      else
         sprintf( text, date_format_text[obs->time_precision], day);
      if( *text == ' ')       /* ensure a leading zero here: */
         *text = '0';
      sprintf( text + strlen( text), "%c\t%s\t",
                   (obs->is_included ? ' ' : 'X'), obs->mpc_code);
      angle = obs->ra * 12. / PI;
      angle = fmod( angle, 24.);
      if( angle < 0.) angle += 24.;
      output_angle_to_buff( text + strlen( text), angle, obs->ra_precision);
      strcat( text, (base_format == RESIDUAL_FORMAT_FULL_WITH_TABS) ?
                              "\t" : "\t ");
      }
   else        /* 'short' MPC format: */
      {
      if( four_digit_years)
         *text++ = int_to_mutant_hex_char( year / 100);
      sprintf( text, "%02d%02d%02d %s",
                       abs( year % 100L), month, (int)day, obs->mpc_code);
      }
   text += strlen( text);

   compute_observation_motion_details( obs, &m);        /* mpc_obs.cpp */
   if( resid_format & RESIDUAL_FORMAT_TIME_RESIDS)
      {
      const double abs_time_resid = fabs( m.time_residual);
      const char sign = (m.time_residual < 0. ? '-' : '+');

      if( abs_time_resid < .994)                 /* show as " +.31s " */
         sprintf( xresid, " %c.%02ds ", sign,
                     (int)( abs_time_resid * 100. + .5));
      else if( abs_time_resid < 9.9)             /* show as " -4.7s " */
         sprintf( xresid, " %+4.1lfs ", m.time_residual);
      else if( abs_time_resid < 999.)            /* show as " -217s " */
         sprintf( xresid, " %c%03ds ", sign,
                     (int)( abs_time_resid + .5));
      else if( abs_time_resid / 60. < 999.)      /* show as " +133m " */
         sprintf( xresid, " %c%03dm ", sign,
                     (int)( abs_time_resid / 60. + .5));
      else if( abs_time_resid / 3600. < 9999.)   /* show as " +027h " */
         sprintf( xresid, " %c%03dh ", sign,
                     (int)( abs_time_resid / 3600. + .5));
      else                                   /* Give up after 1000 hours; */
         strcpy( xresid, " !!!! ");          /* show "it's a long time"   */
      put_residual_into_text( yresid, m.cross_residual);
      }
   else
      {
      put_residual_into_text( xresid, m.xresid);
      put_residual_into_text( yresid, m.yresid);
      }
   if( base_format != RESIDUAL_FORMAT_SHORT)
      {
      const char *tab_separator =
             ((base_format == RESIDUAL_FORMAT_FULL_WITH_TABS) ? "\t" : "");

      angle = obs->dec * 180. / PI;
      if( angle < 0.)
         {
         angle = -angle;
         *text++ = '-';
         if( angle < -99.)       /* error prevention */
            angle = -99.;
         }
      else
         {
         *text++ = '+';
         if( angle > 99.)       /* error prevention */
            angle = 99.;
         }
      output_angle_to_buff( text, angle, obs->dec_precision);

      sprintf( text + strlen( text), "\t%s%s%s\t", xresid, tab_separator, yresid);
      format_dist_in_buff( xresid, obs->r);
      if( resid_format & RESIDUAL_FORMAT_MAG_RESIDS)
         put_mag_resid( yresid, obs->obs_mag, obs->computed_mag, obs->mag_band);
      else
         format_dist_in_buff( yresid, obs->solar_r);
      sprintf( text + strlen( text), "%s%s%s", xresid, tab_separator, yresid);
      }
   else        /* 'short' MPC format */
      {
      if( resid_format & RESIDUAL_FORMAT_MAG_RESIDS)
         {
         put_mag_resid( yresid, obs->obs_mag, obs->computed_mag, obs->mag_band);
         put_residual_into_text( xresid, sqrt( m.xresid * m.xresid
                                             + m.yresid * m.yresid));
         xresid[5] = ' ';        /* replace the '+' with a ' ' */
         }
      strcpy( text, xresid);
      strcpy( text + 6, yresid);
      text[0] = (obs->is_included ? ' ' : '(');
      text[12] = (obs->is_included ? ' ' : ')');
      text[13] = '\0';
      }
                       /* for all other formats, replace tabs w/spaces: */
   if( base_format != RESIDUAL_FORMAT_FULL_WITH_TABS)
      for( i = 0; original_text_ptr[i]; i++)
         if( original_text_ptr[i] == '\t')
            original_text_ptr[i] = ' ';
}

static const char *observe_filename = "observe.txt";

/* The MPC report format goes to only six decimal places in time,
a "microday".  If the reported time is more precise than that -- as can
happen with video observations -- a workaround is to make use of the object
motion data to adjust the position by up to half a microday.  We only do
this if the time given is more than a nanoday away from an integer
microday.  That simply avoids processing in the (overwhelmingly likely)
case that the data doesn't fall on a microday. */

static inline void set_obs_to_microday( OBSERVE FAR *obs)
{
   double delta_jd = obs->jd - floor( obs->jd);

   delta_jd = 1e+6 * delta_jd + .5;
   delta_jd = (delta_jd - floor( delta_jd) - .5) * 1e-6;
// if( delta_jd > 1e-9 || delta_jd < -1e-9)
      {
      MOTION_DETAILS m;
      const double cvt_motions_to_radians_per_day =
                  (PI / 180.) * 24. / 60.;

      compute_observation_motion_details( obs, &m);
      obs->jd -= delta_jd;
                  /* motions are given in '/hr, a.k.a. "/min: */
      obs->ra -= m.ra_motion * delta_jd * cvt_motions_to_radians_per_day
                        / cos( obs->dec);
      obs->dec -= m.dec_motion * delta_jd * cvt_motions_to_radians_per_day;
      }
}

void recreate_observation_line( char *obuff, const OBSERVE FAR *obs)
{
   char buff[100];
   int mag_digits_to_erase = 0;
   OBSERVE tobs = *obs;

   set_obs_to_microday( &tobs);
   format_observation( &tobs, buff, 4);
   memcpy( obuff, obs->packed_id, 12);
   obuff[12] = obs->discovery_asterisk;
   obuff[13] = obs->note1;
   obuff[14] = obs->note2;
   memcpy( obuff + 15, buff, 17);      /* date/time */
   memcpy( obuff + 32, buff + 24, 12);      /* RA */
   memcpy( obuff + 44, buff + 38, 13);      /* dec */
   sprintf( obuff + 57, "%13.2lf%c%c%s%s", obs->obs_mag,
              obs->mag_band, obs->mag_band2, obs->reference, obs->mpc_code);
   if( obs->obs_mag < 0.05)        /* no mag given;  clean out that value */
      mag_digits_to_erase = 5;
   else
      mag_digits_to_erase = 2 - obs->mag_precision;
   memset( obuff + 70 - mag_digits_to_erase, ' ', mag_digits_to_erase);
   if( !obs->is_included)
      obuff[64] = 'x';
}

#ifdef NOT_QUITE_READY_YET
void recreate_second_observation_line( char *buff, const OBSERVE FAR *obs)
{
   int i;
   double vect[3];

   buff[32] = '0' + obs->satellite_obs;
   for( i = 0; i < 3; i++)
      vect[i] = obs->obs_posn[j] - ?; (gotta get earths loc somewhere...)
   ecliptic_to_equatorial( vect);
   for( i = 0; i < 3; i++)
      sprintf( buff + 33 + i * 12, "%12.8lf", vect[i]);
   buff[69] = ' ';
}
#endif

void create_obs_file( const OBSERVE FAR *obs, int n_obs, const int append)
{
   FILE *ofile = fopen( observe_filename, append ? "ab" : "wb");
   double default_weight = 1.;

   while( n_obs--)
      {
      char obuff[81];

      if( obs->weight != default_weight)
         {
         fprintf( ofile, "#Weight %lg\n", obs->weight);
         default_weight = obs->weight;
         }
      recreate_observation_line( obuff, obs);
      fprintf( ofile, "%s\n", obuff);
      if( obs->second_line)
         fprintf( ofile, "%s\n", obs->second_line);
      obs++;
      }
   fclose( ofile);
}

static void add_final_period( char *buff)
{
   if( *buff && buff[strlen( buff) - 1] != '.')
      strcat( buff, ".");
}

static void tack_on_names( char *list, const char *names)
{
   while( *names)
      {
      int i, len, already_in_list = 0;

      while( *names == ' ')
         names++;
      for( len = 0; names[len] && names[len] != ','; len++)
         ;
      for( i = 0; list[i]; i++)
         if( !i || (i > 1 && list[i - 2] == ','))
            if( !memcmp( list + i, names, len))
               if( list[i + len] == ',' || !list[i + len])
                  already_in_list = 1;
      if( !already_in_list)
         {
         char *lptr;

         if( *list)
            strcat( list, ", ");
         lptr = list + strlen( list);
         memcpy( lptr, names, len);
         lptr[len] = '\0';
         }
      names += len;
      if( *names == ',')
         names++;
      }
}

static int get_observer_details( const char *observation_filename,
      const char *mpc_code, char *observers, char *measurers, char *scope)
{
   FILE *ifile = fopen( observation_filename, "rb");
   int rval = 0, no_codes_found = 1;

   *observers = *measurers = *scope = '\0';
   if( ifile)
      {
      char buff[90];

      while( fgets( buff, sizeof( buff), ifile))
         if( !memcmp( buff, "COD ", 4))
            {
            int new_code_found = 0;

            no_codes_found = 0;
            if( !memcmp( buff + 4, mpc_code, 3))
               while( fgets( buff, sizeof( buff), ifile) && !new_code_found)
                  {
                  char *add_to = NULL;

                  remove_trailing_cr_lf( buff);
                  if( !memcmp( buff, "OBS ", 4))
                     add_to = observers;
                  if( !memcmp( buff, "MEA ", 4))
                     add_to = measurers;
                  if( !memcmp( buff, "TEL ", 4))  /* allow for only one scope */
                     strcpy( scope, buff + 4);
                  if( add_to)
                     tack_on_names( add_to, buff + 4);
                  if( !memcmp( buff, "COD ", 4))
                     if( memcmp( buff + 4, mpc_code, 3))
                        new_code_found = 1;
                  }
            }
      fclose( ifile);
      add_final_period( observers);
      add_final_period( measurers);
      add_final_period( scope);
      if( !strcmp( observers, measurers))
         *measurers = '\0';
      }

   if( *observers)
      rval = 1;
   if( *measurers)
      rval |= 2;
   if( *scope)
      rval |= 4;
   if( no_codes_found)        /* we can just ignore this file completely, */
      rval = -1;              /* even for other observatory codes */
   return( rval);
}

#define REPLACEMENT_COLUMN 42

static void observer_link_substitutions( char *buff)
{
   FILE *ifile = fopen( "observer.txt", "rb");

   if( ifile)
      {
      char line[200], *loc;

      while( fgets( line, sizeof( line), ifile))
         if( *line != ';' && *line != '#')
            {
            line[REPLACEMENT_COLUMN - 1] = '\0';
            remove_trailing_cr_lf( line);
            loc = strstr( buff, line);
            if( loc)
               {
               int len = strlen( line), len2;

               if( loc[len] <= ' ' || loc[len] == '.' || loc[len] == ',')
                  {
                  remove_trailing_cr_lf( line + REPLACEMENT_COLUMN);
                  len2 = strlen( line + REPLACEMENT_COLUMN);
                  memmove( loc + len2, loc + len, strlen( loc + len) + 2);
                  memcpy( loc, line + REPLACEMENT_COLUMN, len2);
                  }
               }
            }
      fclose( ifile);
      }
}

static int write_observer_data_to_file( FILE *ofile, const char *ast_filename,
                 const int n_obs, const OBSERVE FAR *obs_data)
{
   int stations[500], n_stations = 0, i, j;
   int try_ast_file = 1, try_details_file = 1, try_scope_file = 1;

   for( i = 0; i < n_obs; i++)
      {
      int match_found = 0;

      for( j = 0; j < n_stations && !match_found; j++)
         match_found = !FSTRCMP( obs_data[i].mpc_code,
                                obs_data[stations[j]].mpc_code);
      if( !match_found)    /* new one:  add it to the array */
         stations[n_stations++] = i;
      }
               /* now do a simple bubblesort: */
   for( i = 0; i < n_stations; i++)
      for( j = 0; j < i; j++)
         if( FSTRCMP( obs_data[stations[i]].mpc_code,
                     obs_data[stations[j]].mpc_code) < 0)
            {
            int temp = stations[i];

            stations[i] = stations[j];
            stations[j] = temp;
            }

   for( i = 0; i < n_stations; i++)
      {
      char buff[200], tbuff[100];
      char details[3][300];
      int loc, j = 0, details_found = 0;

      FSTRCPY( tbuff, obs_data[stations[i]].mpc_code);
      put_observer_data_in_text( tbuff, buff);
      fprintf( ofile, "(%s) %s", tbuff, buff);

      if( try_ast_file)
         {
         details_found = get_observer_details( ast_filename, tbuff,
                                 details[0], details[1], details[2]);
         if( details_found == -1)         /* total washout; skip this file */
            {
            details_found = 0;
            try_ast_file = 0;
            }
         }

      if( !details_found && try_details_file)
         {
         details_found = get_observer_details( "details.txt", tbuff,
                                 details[0], details[1], details[2]);
         if( details_found == -1)         /* total washout; skip this file */
            {
            details_found = 0;
            try_details_file = 0;
            }
         }

      if( !details_found && try_scope_file)
         {
         details_found = get_observer_details( "scopes.txt", tbuff,
                                 details[0], details[1], details[2]);
         if( details_found == -1)         /* total washout; skip this file */
            {
            details_found = 0;
            try_scope_file = 0;
            }
         }

      fprintf( ofile, ".");
      loc = 7 + strlen( buff);
      for( j = 0; j < 3; j++)
         if( *details[j])
            {
            char inserted_text[15], *outtext = details[j];

            if( j == 2)
               strcpy( inserted_text, " ");
            else
               {
               strcpy( inserted_text, j ? " Measurer" : "  Observer");
               if( strchr( outtext, ','))
                  strcat( inserted_text, "s");
               strcat( inserted_text, " ");
               }
            memmove( outtext + strlen( inserted_text), outtext,
                              strlen( outtext) + 1);
            memcpy( outtext, inserted_text, strlen( inserted_text));
            while( *outtext)
               {
               int k, done;

               for( k = 0; outtext[k] && outtext[k] != ' '; k++)
                  ;
               done = !outtext[k];
               outtext[k] = '\0';
               if( loc + k > 78)    /* gotta go to a new line */
                  {
                  fprintf( ofile, "\n    %s", outtext);
                  loc = k + 4;
                  }
               else
                  {
                  fprintf( ofile, " %s", outtext);
                  loc += k + 1;
                  }
               outtext += k;
               if( !done)
                  outtext++;
               }
            }
      fprintf( ofile, "\n");
      }
   return( 0);
}

int write_residuals_to_file( const char *filename, const char *ast_filename,
       const int n_obs, const OBSERVE FAR *obs_data, const int resid_format)
{
   FILE *ofile = fopen( filename, "w");
   int rval = 0;

   if( ofile )
      {
      char buff[100];
      int number_lines = (n_obs + 2) / 3;
      int i;

      if( (resid_format & 3) == RESIDUAL_FORMAT_SHORT)
         for( i = 0; i < number_lines * 3; i++)
            {
            int num = (i % 3) * number_lines + i / 3;
            OBSERVE FAR *obs = ((OBSERVE FAR *)obs_data) + num;

            if( num < n_obs)
               format_observation( obs, buff, resid_format);
            else
               *buff = '\0';
            fprintf( ofile, "%s%s", buff, (i % 3 == 2) ? "\n" : "   ");
            }
      else
         for( i = 0; i < n_obs; i++)
            {
            format_observation( obs_data + i, buff, resid_format);
            fprintf( ofile, "%s\n", buff);
            }
      fprintf( ofile, "\nStation data:\n");
      write_observer_data_to_file( ofile, ast_filename, n_obs, obs_data);
      fclose( ofile);
      }
   else                    /* file not opened */
      rval = -1;
   return( rval);
}

#ifdef FUTURE_PROJECT_IN_WORKS

int create_residual_scattergram( const char *filename, const int n_obs,
                         const OBSERVE FAR *obs)
{
   const int tbl_height = 19, tbl_width = 71;
   const int xspacing = 12, yspacing = 5,
   short *remap_table = (short *)calloc( tbl_height * tbl_width, sizeof( short));
   int i, j;
   FILE *ofile = fopen( filename, "wb");

   for( i = 0; i < n_obs; i++, obs++)
      {
      const double yresid = 3600. * (180./pi) * (obs->dec - obs->computed_dec);
      const double xresid = 3600. * (180./pi) * (obs->ra - obs->computed_ra)
                                        * cos( obs->computed_dec);
      int xloc = (int)floor( xresid * (double)xspacing / .5)
      }

   fclose( ofile);
   free( remap_table);
}
#endif

void remove_trailing_cr_lf( char *buff)
{
   int i;

   for( i = 0; buff[i] && buff[i] != 13 && buff[i] != 10; i++)
      ;
   while( i && buff[i - 1] == ' ')
      i--;
   buff[i] = '\0';
}

int text_search_and_replace( char FAR *str, const char *oldstr,
                                     const char *newstr)
{
   int ilen = FSTRLEN( str), rval = 0;
   const int oldlen = strlen( oldstr);
   const int newlen = strlen( newstr);

   while( ilen >= oldlen)
      if( !FMEMCMP( str, oldstr, oldlen))
         {
         FMEMMOVE( str + newlen, str + oldlen, ilen - oldlen + 1);
         FMEMCPY( str, newstr, newlen);
         str += newlen;
         ilen -= oldlen;
         rval = 1;
         }
      else
         {
         str++;
         ilen--;
         }
   return( rval);
}

static long round_off( const double ival, const double prec)
{
   long rval = 0L, digit;
   int got_it = 0;
   double diff;

   for( digit = 10000000L; !got_it; digit /= 10)
      {
      rval = (((long)ival + digit / 2L) / digit) * digit;
      diff = fabs( (double)rval - ival);
      if( digit == 1 || diff < ival * prec)
         got_it = 1;
      }
   return( rval);
}

int make_pseudo_mpec( const char *mpec_filename, const char *obj_name)
{
   FILE *ofile = fopen( mpec_filename, "wb");
   FILE *ifile = fopen( "header.htm", "rb");
   FILE *elements_file;
   char buff[500], mpec_buff[7];
   int line_no = 0, rval = 0, total_lines = 0;
   int mpec_no = atoi( get_environment_ptr( "MPEC"));

   if( mpec_no)
      sprintf( mpec_buff, "_%02x", mpec_no % 256);
   else
      *mpec_buff = '\0';
   if( ifile)                 /* copy header data to pseudo-MPEC */
      {
      while( fgets( buff, sizeof( buff), ifile))
         if( *buff != '#')
            {
            char *tptr = strstr( buff, "_xx");

            if( tptr)
               {
               memmove( tptr + strlen( mpec_buff), tptr + 3, strlen( tptr));
               memcpy( tptr, mpec_buff, strlen( mpec_buff));
               }
            while( (tptr = strchr( buff, '$')) != NULL)
               {                       /* See comments in 'header.htm'.  */
               int got_it = 0, i;      /* code replaces text between $s  */
                                       /* in that file.                  */
               for( i = 1; tptr[i] && tptr[i] != '$'; i++)
                  ;        /* search for matching $ */
               if( i < 20 && tptr[i] == '$')
                  {
                  FILE *elem_file = fopen( elements_filename, "rb");
                  char search_str[80], replace_str[80];

                  memcpy( search_str, tptr, i);
                  search_str[i] = '\0';
                  if( elem_file)
                     {
                     char tbuff[300], *tptr2;

                     if( !strcmp( search_str, "$Tg"))
                        {
                        time_t t = time( NULL);

                        strcpy( replace_str, ctime( &t));
                        remove_trailing_cr_lf( replace_str);
                        got_it = 1;
                        }
                     else if( !strcmp( search_str, "$Name"))
                        {
                        strcpy( replace_str, obj_name);
                        got_it = 1;
                        }

                     while( !got_it &&
                             fgets_trimmed( tbuff, sizeof( tbuff), elem_file))
                        if( (tptr2 = strstr( tbuff, search_str)) != NULL
                                    && tptr2[i] == '=')
                           {
                           tptr2 += i + 1;
                           for( i = 0; tptr2[i] > ' '; i++)
                              replace_str[i] = tptr2[i];
                           replace_str[i] = '\0';
                           got_it = 1;
                           }
                     fclose( elem_file);
                     strcat( search_str, "$");
                     if( got_it)
                        text_search_and_replace( buff, search_str, replace_str);
                     }
                  }
               else        /* no matching '$' found */
                  *tptr = '!';
               }
            fputs( buff, ofile);
            }
      fclose( ifile);
      }

   if( mpec_no)
      {
      sprintf( buff, "%d", mpec_no % 255 + 1);
      set_environment_ptr( "MPEC", buff);
      }

   ifile = fopen( observe_filename, "rb");
   if( ifile)
      {
      int neocp_line_number = 0;

      while( fgets( buff, sizeof( buff), ifile))
         if( *buff != '#')       /* may be comments or 'weight' lines */
            {
            if( buff[14] == 's' || buff[14] == 'v' || buff[14] == 'r')
               fprintf( ofile, "%s", buff);
            else
               {
               char mpc_code[4];

               total_lines++;
               memcpy( mpc_code, buff + 77, 3);
               buff[12] = mpc_code[3] = buff[77] = '\0';
               fprintf( ofile, "<a name=\"o%s%03d\"></a><a href=\"#r%s%03d\">%s</a>",
                        mpec_buff, total_lines, mpec_buff, total_lines, buff);
               if( !memcmp( buff + 72, "NEOCP", 5))
                  {
                  const char *replacement_text;
                  int i;

                  strcpy( buff + 25, "<code class=\"neocp\">");
                  if( !neocp_line_number)
                     replacement_text = "~~~~</code>Astrometry<code class=\"neocp\">~~~~</code>redacted;<code class=\"neocp\">~~~~~~~~~~~~~~~~~~~~";
                  else if( neocp_line_number == 1)
                     replacement_text = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~</code>see<code class=\"neocp\">~~~~~~~</code><a href=\"http://www.cfa.harvard.edu/iau/NEO/ToConfirm.html\">NEOCP</a><code class=\"neocp\">";
                  else
                     replacement_text = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
                  strcat( buff + 25, replacement_text);
                  for( i = 25; buff[i]; i++)
                     if( buff[i] == '~')
                        {
                        const char *forbidden = "~ <>\"&";

                        while( strchr( forbidden, buff[i]))
                           buff[i] = (char)( ' ' + rand( ) % ('z' - ' '));
                        }
                  strcat( buff + 25, "NEOCP</code>");
                  neocp_line_number++;
                  }
               fprintf( ofile, " %s<a href=\"#stn_%s\">%s</a>\n",
                        buff + 13, mpc_code, mpc_code);
#ifdef NOW_OBSOLETE
               if( !strcmp( mpc_code, "247"))         /* roving observer */
                  {
                  extern double roving_lon, roving_lat, roving_ht_in_meters;

                  buff[12] = ' ';      /* repair line */
                  buff[14] = 'v';
                  sprintf( buff + 32, "%11.5lf%11.5lf%7d                247",
                           roving_lon * 180. / PI, roving_lat * 180. / PI,
                           (int)roving_ht_in_meters);
                  fprintf( ofile, "%s\n", buff);
                  }
#endif
               }
            }
      fclose( ifile);
      }
   else
      rval |= 1;

   ifile = fopen( residual_filename, "rb");
   if( ifile)
      {
      FILE *obslinks_file = fopen( "obslinks.htm", "rb");
      FILE *mpc_obslinks_file = fopen( "ObsCodesF.html", "rb");
      long obslinks_header_len = 0, mpc_obslinks_header_len = 0;
      char url[200];
              /* In making a pseudo-MPEC,  we attempt to include links     */
              /* to the Web sites of the observatory codes mentioned.  The */
              /* files 'obslinks.htm' (my effort at a complete list of Web */
              /* sites of observatory codes) and 'ObsLinksF.html' (the     */
              /* MPC's list) are both searched,  in that order.            */

      if( obslinks_file)
         {
         while( fgets( url, sizeof( url), obslinks_file)
                        && memcmp( url, "<a name=\"0\">", 12))
            ;
         obslinks_header_len = ftell( obslinks_file);
         }
      if( mpc_obslinks_file)
         {
         while( fgets( url, sizeof( url), mpc_obslinks_file)
                        && memcmp( url, "<pre>", 5))
            ;
         mpc_obslinks_header_len = ftell( mpc_obslinks_file);
         }
      while( fgets( buff, sizeof( buff), ifile) && memcmp( buff, "Station", 7))
         ;
      fprintf( ofile, "<a name=\"stations\"></a>\n");
      fprintf( ofile, "<b>%s</b>", buff);
      while( fgets( buff, sizeof( buff), ifile))
         if( *buff == ' ')
            {
            observer_link_substitutions( buff);
            fprintf( ofile, "%s", buff);
            }
         else
            {
            char tbuff[4];
            int compare = 1, i, url_index = 0;
            char *latlon, *remains;

            memcpy( tbuff, buff + 1, 3);
            tbuff[3] = '\0';
            remove_trailing_cr_lf( buff);
            fprintf( ofile, "<a name=\"stn_%s\"></a>", tbuff);

            for( i = 5; buff[i] && buff[i] != ')'; i++)
               ;

            if( buff[i])  /* found the closing paren of a lat/lon */
               {
               remains = buff + i;
               while( i > 2 && buff[i] != '(')
                  i--;
               buff[i - 2] = *remains = '\0';
               latlon = buff + i + 1;
               }
            else  /*  no lat/lon;  assume name ends with a . */
               {
               latlon = NULL;
               for( i = 0; buff[i] && buff[i] != '.'; i++)
                  ;
               if( buff[i] == '.' && buff[i + 1])        /* found the '.'! */
                  {
                  buff[i + 1] = '\0';
                  remains = buff + i + 1;
                  }
               else
                  remains = buff + i;
               }
            if( compare && obslinks_file)
               {
               url[19] = '\0';
               fseek( obslinks_file, obslinks_header_len, SEEK_SET);
               while( (compare = memcmp( url + 19, buff + 1, 3)) != 0 &&
                                 fgets( url, sizeof( url), obslinks_file))
                  ;
               url_index = 24;   /* if there is a link,  it starts in byte 24 */
               }
            if( compare && mpc_obslinks_file)
               {
               *url = '\0';
               fseek( mpc_obslinks_file, mpc_obslinks_header_len, SEEK_SET);
               while( (compare = memcmp( url, buff + 1, 3)) != 0 &&
                                 fgets( url, sizeof( url), mpc_obslinks_file))
                  ;
               if( !compare)     /* got the observatory code;  is there a link? */
                  if( url[33] != '<')   /* no;  don't bother with it */
                     compare = 1;
               url_index = 33;   /* if there is a link,  it starts in byte 33 */
               }

            if( !compare)   /* we got a link to an observatory code */
               {
               for( i = url_index; url[i] && url[i] != '>'; i++)
                  ;
               url[i + 1] = '\0';
               buff[5] = '\0';
               fprintf( ofile, "%s %s%s</a>",
                       buff, url + url_index, buff + 6);
               }
            else
               fprintf( ofile, "%s", buff);

            if( latlon)
               {
               double lat, lon;
               char lat_sign, lon_sign;

               sscanf( latlon, "%c%lf %c%lf", &lat_sign, &lat, &lon_sign, &lon);
               if( lat_sign == 'S')
                  lat = -lat;
               if( lon_sign == 'W')
                  lon = -lon;
               fprintf( ofile, " (<a title=\"Click for map\"");
               fprintf( ofile, " href=\"http://mappoint.msn.com/map.aspx?&amp;C=%.3lf,%.3lf&amp;A=1000\">",
                              lat, lon);
               fprintf( ofile, "%s</a>)", latlon);
               }
            observer_link_substitutions( remains + 1);
            fprintf( ofile, "%s\n", remains + 1);
            }
      if( obslinks_file)
         fclose( obslinks_file);
      if( mpc_obslinks_file)
         fclose( mpc_obslinks_file);
      }
   else
      rval |= 2;

   elements_file = fopen( elements_filename, "rb");
   if( elements_file)
      {
      fprintf( ofile, "<a name=\"elements%s\"></a>\n", mpec_buff);
      while( fgets( buff, sizeof( buff), elements_file))
         if( *buff != '#')
            {
            if( !memcmp( buff, "Orbital ele", 11))
               fprintf( ofile, "<b>%s</b>", buff);
            else if( *buff == 'P' && buff[19] == 'H')
               {
               const double abs_mag = atof( buff + 20);
                        /* H=4 indicates 420 to 940 km,  so: */
               double upper_size = 940. * exp( (4. - abs_mag) * LOG_10 / 5.);
               const char *units = "km";
               const char *size_url =
                   "href=\"http://www.cfa.harvard.edu/iau/lists/Sizes.html\">";
               char title[50];

               buff[19] = '\0';
               if( upper_size < .004)   /* under four meters,  use cm as units: */
                  {
                  upper_size *= 1000. * 100.;
                  units = "cm";
                  }
               else if( upper_size < 4.)  /* under four km,  use meters: */
                  {
                  upper_size *= 1000.;
                  units = "meters";
                  }
               sprintf( title, "\"Size is probably %ld to %ld %s\"\n",
                          round_off( upper_size / sqrt( 5.), .1),
                          round_off( upper_size, .1), units);
               fprintf( ofile, "%s<a title=%s%sH</a>%s",
                           buff, title, size_url, buff + 20);
               }
            else
               {
               text_search_and_replace( buff, "m^2", "m<sup>2</sup>");
               text_search_and_replace( buff, "   Find_Orb",
                                "   <a href=\"http://www.projectpluto.com/find_orb.htm\">Find_Orb</a>");
               fputs( buff, ofile);
               }
            }
      fclose( elements_file);
      }
   else
      rval |= 4;

               /* _now_ write out residuals: */
   if( ifile)
      {
      fseek( ifile, 0L, SEEK_SET);
      fprintf( ofile, "<a name=\"residuals%s\"></a>\n", mpec_buff);
      fprintf( ofile, "<b>Residuals in arcseconds:</b>\n");
      line_no = 0;
      while( fgets( buff, sizeof( buff), ifile) && *buff > ' ')
         {
         int i, column_off = (total_lines + 2) / 3, line;

         line_no++;
         for( i = 0; i < 3; i++)
            if( (line = line_no + column_off * i) <= total_lines)
               {
               char *tptr = buff + i * 26, tbuff[20];

               tptr[6] = '\0';         /* put out the YYMMDD... */
               fprintf( ofile, "<a name=\"r%s%03d\"></a><a href=\"#o%s%03d\">%s</a>",
                        mpec_buff, line, mpec_buff, line, tptr);

               memcpy( tbuff, tptr + 7, 3);     /* ...then the obs code.. */
               tbuff[3] = '\0';
               fprintf( ofile, " <a href=\"#stn_%s\">%s</a>", tbuff, tbuff);

               tptr[23] = '\0';        /* ...and finally,  the residuals */
               fprintf( ofile, "%s   ", tptr + 10);
               }
         fprintf( ofile, "\n");
         }
      fclose( ifile);
      }

               /* ...and now,  the ephemeris: */
   ifile = fopen( ephemeris_filename, "r");
   if( ifile && fgets( buff, sizeof( buff), ifile))
      {
      remove_trailing_cr_lf( buff);
      fprintf( ofile, "\n<a name=\"eph%s\"></a>", mpec_buff);
      if( !strcmp( buff, "#Geocentric"))
         fprintf( ofile, "<b>Ephemerides:</b>\n");
      else
         fprintf( ofile, "<b>Ephemerides (%s):</b>\n", buff + 1);

      while( fgets( buff, sizeof( buff), ifile))
         fputs( buff, ofile);
      fclose( ifile);
      }
   else
      rval |= 8;

   fprintf( ofile, "</pre></body></html>\n");
   fclose( ofile);
   return( 0);
}

