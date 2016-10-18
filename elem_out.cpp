/* elem_out.cpp: formatting elements into human-friendly form

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
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <assert.h>
#include "watdefs.h"
#include "comets.h"
#include "mpc_obs.h"
#include "date.h"
#include "afuncs.h"
#include "monte0.h"     /* for put_double_in_buff() proto */
#include "showelem.h"

#define J2000 2451545.
#define PI 3.1415926535897932384626433832795028841971693993751058209749445923
#define GAUSS_K .01720209895
#define SOLAR_GM (GAUSS_K * GAUSS_K)
#define JD_TO_YEAR(jd)  (2000. + ((jd)-J2000) / 365.25)

// void elements_in_tle_format( char *buff, const ELEMENTS *elem);
static int elements_in_mpcorb_format( char *buff, const char *packed_desig,
                const char *full_desig, const ELEMENTS *elem,
                const OBSERVE FAR *obs, const int n_obs);   /* orb_func.c */
int elements_in_guide_format( char *buff, const ELEMENTS *elem,
                     const char *obj_name);                 /* orb_func.c */
int find_worst_observation( const OBSERVE FAR *obs, const int n_obs);
double initial_orbit( OBSERVE FAR *obs, int n_obs, double *orbit);
int set_locs( const double *orbit, double t0, OBSERVE FAR *obs, int n_obs);
int text_search_and_replace( char FAR *str, const char *oldstr,
                                     const char *newstr);   /* ephem0.cpp */
double calc_obs_magnitude( const int is_comet, const double obj_sun,
          const double obj_earth, const double earth_sun, double *phase_ang);
int find_best_fit_planet( const double jd, const double *ivect,
                                 double *rel_vect);         /* runge.cpp */
const char *get_environment_ptr( const char *env_ptr);     /* mpc_obs.cpp */
void remove_trailing_cr_lf( char *buff);      /* ephem0.cpp */
int write_tle_from_vector( char *buff, const double *state_vect,
        const double epoch, const char *norad_desig, const char *intl_desig);
double find_moid( const ELEMENTS *elem1, const ELEMENTS *elem2); /* moid4.c */
void setup_planet_elem( ELEMENTS *elem, const int planet_idx,
                                          const double t_cen);   /* moid4.c */
void set_environment_ptr( const char *env_ptr, const char *new_value);
int store_defaults( const int ephemeris_output_options);    /* elem_out.c */
int get_defaults( int *ephemeris_output_options);           /* elem_out.c */
double find_collision_time( ELEMENTS *elem, double *latlon, const int is_impact);
char *fgets_trimmed( char *buff, size_t max_bytes, FILE *ifile); /*elem_out.c*/
int get_idx1_and_idx2( const int n_obs, const OBSERVE FAR *obs,
                                int *idx1, int *idx2);      /* elem_out.c */
int get_planet_posn_vel( const double jd, const int planet_no,
                     double *posn, double *vel);
char int_to_mutant_hex_char( const int ival);               /* mpc_obs.c */
double mag_band_shift( const char mag_band);                /* elem_out.c */
int get_jpl_ephemeris_info( int *de_version, double *jd_start, double *jd_end);

double asteroid_magnitude_slope_param = .15;
double comet_magnitude_slope_param = 10.;
char default_comet_magnitude_type = 'N';
int debug_printf( const char *format, ...);                /* runge.cpp */

char *fgets_trimmed( char *buff, size_t max_bytes, FILE *ifile)
{
   char *rval = fgets( buff, max_bytes, ifile);

   if( rval)
      {
      int i;

      for( i = 0; buff[i] && buff[i] != 10 && buff[i] != 13; i++)
         ;
      buff[i] = '\0';
      }
   return( rval);
}

static void get_first_and_last_included_obs( const OBSERVE *obs,
              const int n_obs, int *first, int *last)
{
   if( first)
      for( *first = 0; *first < n_obs - 1 && !obs[*first].is_included;
                                    (*first)++)
         ;
   if( last)
      for( *last = n_obs - 1; *last && !obs[*last].is_included; (*last)--)
         ;
}

void make_date_range_text( char *obuff, const double jd1, const double jd2)
{
   long year, year2;
   int month, month2;
   const int day1 = (int)decimal_day_to_dmy( jd1, &year,  &month,
                                    CALENDAR_JULIAN_GREGORIAN);
   const int day2 = (int)decimal_day_to_dmy( jd2, &year2, &month2,
                                    CALENDAR_JULIAN_GREGORIAN);
   static const char *month_names[] = { "Jan.", "Feb.", "Mar.", "Apr.", "May",
            "June", "July", "Aug.", "Sept.", "Oct.", "Nov.", "Dec." };

   if( year == year2)
      {
      sprintf( obuff, "%ld %s %d", year, month_names[month - 1], day1);
      obuff += strlen( obuff);
      if( month == month2 && day1 != day2)
         sprintf( obuff, "-%d", day2);
      else if( month != month2)
         sprintf( obuff, "-%s %d", month_names[month2 - 1], day2);
      }
   else              /* different years */
      sprintf( obuff, "%ld %s %d-%ld %s %d", year, month_names[month - 1],
                             day1, year2, month_names[month2 - 1], day2);

   obuff += strlen( obuff);
   if( jd2 - jd1 < 10. / 86400.)  /* less than 10 seconds: show to .01 sec */
      sprintf( obuff, " (%.2lf sec)", (jd2 - jd1) * 1440. * 60.);
   else if( jd2 - jd1 < 100. / 86400.) /* less than 100 seconds: show to .1 sec */
      sprintf( obuff, " (%.1lf sec)", (jd2 - jd1) * 1440. * 60.);
   else if( jd2 - jd1 < 100. / 1440.)     /* less than 100 minutes: show in min */
      sprintf( obuff, " (%.1lf min)", (jd2 - jd1) * 1440.);
   else if( jd2 - jd1 < 2.)
      sprintf( obuff, " (%.1lf hr)", (jd2 - jd1) * 24.);
}

/* observation_summary_data( ) produces the final line in an MPC report,
   such as 'From 20 observations 1997 Oct. 20-22;  mean residual 0".257.   '
   Note that the arcsecond mark comes before the decimal point;  this
   oddity is handled using the text_search_and_replace() function.
*/

static void observation_summary_data( char *obuff, const OBSERVE FAR *obs,
                              const int n_obs, const int options)
{
   int i, n_included, first_idx, last_idx;

   get_first_and_last_included_obs( obs, n_obs, &first_idx, &last_idx);

   for( i = n_included = 0; i < n_obs; i++)
      n_included += obs[i].is_included;
   if( (options & ELEM_OUT_ALTERNATIVE_FORMAT) && n_included != n_obs)
      sprintf( obuff, "%d of %d observations ", n_included, n_obs);
   else
      sprintf( obuff, "From %d observations ", n_included);
   if( n_included)
      {
      obuff += strlen( obuff);
      make_date_range_text( obuff, obs[first_idx].jd, obs[last_idx].jd);
      obuff += strlen( obuff);
      sprintf( obuff, "; mean residual %.3lf",
                           compute_rms( obs, n_obs));
      text_search_and_replace( obuff, ".", "\".");
      strcat( obuff, ".");
      }
}

static double centralize_ang( double ang)
{
   ang = fmod( ang, PI + PI);
   if( ang < 0.)
      ang += PI + PI;
   return( ang);
}

void convert_elements( const double epoch_from, const double epoch_to,
      double *incl, double *asc_node, double *arg_per);     /* conv_ele.cpp */

   /* Packed MPC designations have leading and/or trailing spaces.  This */
   /* function lets you get the designation minus those spaces.          */

static void packed_desig_minus_spaces( char *obuff, const char *ibuff)
{
   while( *ibuff && *ibuff == ' ')
      ibuff++;
   while( *ibuff && *ibuff != ' ')
      *obuff++ = *ibuff++;
   *obuff = '\0';
}

/* The 'standard' MPCORB format gives a mean anomaly and semimajor axis.
It can't handle parabolic or hyperbolic orbits,  and precision suffers
for high-eccentricity,  long-period objects.  If 'alt_mpcorb' != 0,
we switch to output wherein the time of perihelion is stored,  split into
an integer part,  where the mean anomaly used to go;  plus a fractional
part,  put where the mean motion used to be;  and the semimajor axis
is replaced with the perihelion distance.

   When loading such an orbit into another program,  we've the advantages
that (a) not much has changed,  so most of the code can be left the same,
and (b) it's easy to tell alternate from standard format;  the latter has
decimal points in columns 30 and 83,  whilst the alternate form has
digits there.

   We've the disadvantage that nearly everybody else uses just the standard
format.  So at least for the nonce,  this is mostly a theoretical improvement.

   Incidentally,  in the new format,  six bytes that used to contain
digits now contain spaces.  At some point,  they may contain an indicator
of a central object (so non-heliocentric orbits can be stored).  And perhaps
some other uses yet to be determined.        */

int alt_mpcorb;

static int elements_in_mpcorb_format( char *buff, const char *packed_desig,
                const char *full_desig, const ELEMENTS *elem,
                const OBSERVE FAR *obs, const int n_obs)   /* orb_func.c */
{
   extern int perturbers;
   int month, day, i, first_idx, last_idx, n_included_obs = 0;
   long year;
   const double jan_1970 = 2440587.5;
   const double rms_err = compute_rms( obs, n_obs);
   const unsigned hex_flags = 0;
            /* 'mpcorb' has four hexadecimal flags starting in column 162, */
            /* signifying if the object is in any of various classes such  */
            /* as Aten,  scattered-disk object,  PHA,  Jupiter Trojan,     */
            /*  etc.  None of those flags are set yet.                     */
   const int n_oppositions = 1;
            /* The above needs some work.  Problem is,  what constitutes   */
            /* an "opposition" for an NEO?  (It's more clearcut for MBOs.) */
            /* For the nonce,  we'll just say "one opposition".            */
   double arc_length;
   char packed_desig2[40];

   packed_desig_minus_spaces( packed_desig2, packed_desig);
   sprintf( buff, "%-8s%5.2lf  %4.2lf ", packed_desig2, elem->abs_mag,
                           asteroid_magnitude_slope_param);
   day = (int)( decimal_day_to_dmy( elem->epoch, &year,
                              &month, CALENDAR_JULIAN_GREGORIAN) + .0001);
   sprintf( buff + 20, "%c%02ld%X%c",
                  int_to_mutant_hex_char( year / 100),
                  year % 100L, month,
                  int_to_mutant_hex_char( day));
   sprintf( buff + 25, "%10.5lf%11.5lf%11.5lf%11.5lf%11.7lf",
           centralize_ang( elem->mean_anomaly) * 180. / PI,
           centralize_ang( elem->arg_per) * 180. / PI,
           centralize_ang( elem->asc_node) * 180. / PI,
           centralize_ang( elem->incl) * 180. / PI,
           elem->ecc);
   sprintf( buff + 79, "%12.8lf%12.7lf",
            (180 / PI) / elem->t0,        /* n */
            alt_mpcorb ? elem->q : elem->major_axis);
   if( alt_mpcorb)
      {
      const double fractional_day = elem->perih_time - floor( elem->perih_time);
      int microdays = (int)( fractional_day * 1000000. + .5);

      sprintf( buff + 26, "%7.0lf", floor( elem->perih_time));
      if( microdays >= 1000000)
         microdays = 999999;
      sprintf( buff + 80, " %06d", microdays);
      buff[33] = buff[34] = buff[87] = buff[88] = buff[89] = buff[90] = ' ';
      }
   for( i = 0; i < n_obs; i++)
      if( obs[i].is_included)
         n_included_obs++;
   day = (int)( decimal_day_to_dmy( jan_1970 + (double)time( NULL) / 86400.,
                         &year, &month, CALENDAR_JULIAN_GREGORIAN) + .0001);
   sprintf( buff + 103,
      "    FO %02d%02d%02d  %4d  %2d ****-**** ****         Find_Orb   %04x",
                  (int)( year % 100), month, (int)day,
                  n_included_obs, n_oppositions, hex_flags);
   get_first_and_last_included_obs( obs, n_obs, &first_idx, &last_idx);
   arc_length = obs[last_idx].jd - obs[first_idx].jd;
   if( arc_length < 99. / 86400.)
      sprintf( buff + 127, "%4.1lf sec ", arc_length * 86400.);
   else if( arc_length < 99. / 1440.)
      sprintf( buff + 127, "%4.1lf min ", arc_length * 1440.);
   else if( arc_length < 2.)
      sprintf( buff + 127, "%4.1lf hrs ", arc_length * 24.);
   else if( arc_length < 600.)
      sprintf( buff + 127, "%4d days", (int)arc_length + 1);
   else
      sprintf( buff + 127, "%4d-%4d",
                (int)JD_TO_YEAR( obs[first_idx].jd),
                (int)JD_TO_YEAR( obs[last_idx].jd));
   buff[136] = ' ';
   sprintf( buff + 165, " %-30s", full_desig);
   day = (int)( decimal_day_to_dmy( obs[last_idx].jd, &year,
                       &month, CALENDAR_JULIAN_GREGORIAN) + .0001);
   sprintf( buff + 194, "%04ld%02d%02d", year, month, day);
   if( rms_err < 9.9)
      sprintf( buff + 137, "%4.2lf", rms_err);
   else if( rms_err < 99.9)
      sprintf( buff + 137, "%4.1lf", rms_err);
   else if( rms_err < 9999.)
      sprintf( buff + 137, "%4.0lf", rms_err);
   buff[141] = ' ';
   if( (perturbers & 0x1fe) == 0x1fe)
      {    /* we have Mercury through Neptune,  at least */
      const char *coarse_perturb, *precise_perturb;

      if( perturbers & 0x700000)    /* asteroids included */
         {
         precise_perturb = (perturbers & 0x400 ? "3E" : "38");
         coarse_perturb = "M-v";
         }
      else        /* non-asteroid case */
         {
         precise_perturb = (perturbers & 0x400 ? "06" : "00");
         coarse_perturb = (perturbers & 0x200 ? "M-P" : "M-N");
         }
      memcpy( buff + 142, coarse_perturb, 3);
      memcpy( buff + 146, precise_perturb, 2);
      }
   return( 0);
}

int elements_in_guide_format( char *buff, const ELEMENTS *elem,
                     const char *obj_name)
{
   int month;
   double day;
   long year;

   day = decimal_day_to_dmy( elem->epoch, &year, &month,
                                              CALENDAR_JULIAN_GREGORIAN);
            /*      name day  mon yr MA      q      e */
   sprintf( buff, "%-43s%7.4lf%4d%5ld%10.5lf%14.7lf%12.7lf%11.6lf%12.6lf%12.6lf",
            obj_name, day, month, year,
            centralize_ang( elem->mean_anomaly) * 180. / PI,
            elem->q, elem->ecc,
            centralize_ang( elem->incl) * 180. / PI,
            centralize_ang( elem->arg_per) * 180. / PI,
            centralize_ang( elem->asc_node) * 180. / PI);
   if( elem->q < .01)
      {
      sprintf( buff + 71, "%12.10lf", elem->q);
      buff[71] = buff[83] = ' ';
      }
   sprintf( buff + strlen( buff), "  2000.0%7.1lf%6.2lf A",
            elem->abs_mag, elem->slope_param);
   if( elem->central_obj)
      sprintf( buff + strlen( buff), "  Center: %d", elem->central_obj);
   return( 0);
}

static int is_cometary( const char *constraints)
{
   const char *ecc = strstr( constraints, "e=1");

   return( ecc && atof( ecc + 2) == 1.);
}

#ifdef OBSOLETE_CODE
                   /* 'magic' number;  the CRC polynomial */
#define FM 0x4c11db7

static void fill_crc_table( uint32_t *tbl, const uint32_t polynomial)
{
   unsigned i, bit;
   uint32_t tblval;

   for( i = 0; i < 256; i++)
      {
      tblval = (uint32_t)i << 24;
      for( bit = 8; bit; bit--)
         if( tblval & 0x80000000)
            tblval = (tblval << 1) ^ polynomial;
         else
            tblval <<= 1;
      *tbl++ = tblval;
      }
}

static uint32_t crc_32( const unsigned char *buff, unsigned nbytes)
{
   static uint32_t *table = NULL;
   uint32_t rval;
   unsigned n;

   if( !table)    /* first pass through,  set up the CRC table */
      {           /* table can be made static,  no problem */
      const uint32_t crc_polynomial = 0x4c11db7;

      table = (uint32_t *)malloc( 256 * sizeof( uint32_t));
      fill_crc_table( table, crc_polynomial);
      }
   for( rval = 0xffffffff; nbytes; nbytes--)
      {
      n = ((unsigned)(rval >> 24) & 0xff) ^ *buff++;
      rval = (rval << 8) ^ table[n];
      }
   return( ~rval);
}

static inline uint32_t compute_orbit_crc( const double curr_epoch,
            const double epoch_shown,
            OBSERVE FAR *obs, const int n_obs,
            const int n_extra_params, const int perturbers)
{
   uint32_t rval;
   double buff[8];
   int i;

   buff[0] = curr_epoch;
   buff[1] = epoch_shown;
   rval = crc_32( (const unsigned char *)buff, 2 * sizeof( double));
   for( i = 0; i < n_obs; i++)
      if( obs[i].is_included)
         {
         buff[0] = obs[i].jd;
         buff[1] = obs[i].obs_posn[0];
         buff[2] = obs[i].obs_posn[1];
         buff[3] = obs[i].obs_posn[2];
         buff[4] = obs[i].weight;
         buff[5] = obs[i].obs_mag;
         rval ^= crc_32( (const unsigned char *)buff, 6 * sizeof( double));
         }
   rval ^= crc_32( (const unsigned char *)&n_extra_params, sizeof( int));
   rval ^= crc_32( (const unsigned char *)&perturbers, sizeof( int));
   return( rval);
}
#endif


int monte_carlo_object_count = 0;
int n_monte_carlo_impactors = 0;
int append_elements_to_element_file = 0;
int using_sr = 0;
char orbit_summary_text[80];
double max_monte_rms;

void set_statistical_ranging( const int new_using_sr)
{
   using_sr = new_using_sr;
   n_monte_carlo_impactors = monte_carlo_object_count = 0;
}


/* The results from write_out_elements_to_file() can be somewhat
varied.  The output for elliptical and parabolic/hyperbolic orbits are
very different.  Asteroids and comets differ in whether H and G are
given,  or M(N)/M(T) and K.  Or those fields can be blank if no
magnitude data is given.

   Constraints or Monte Carlo data can modify the perihelion line,  such as

   Perihelion 1998 Mar 16.601688 TT;  Constraint: a=2.5
   Perihelion 1998 Mar 16.601688 TT;  20.3% impact (203/1000)

   AMR data appears on the epoch line:

Epoch 1997 Oct 13.0 TT; AMR 0.034 m^2/kg

   MOIDs can carry on to next line,  wiping out P Q header and even
the (2000.0) text if there are enough extra MOIDs.

   Examples:

Orbital elements:
1996 XX1
   Perihelion 1998 Mar 8.969329 TT = 23:15:50 (JD 2450881.469329)
Epoch 1997 Oct 13.0 TT = JDT 2450734.5   Earth MOID: 0.0615   Ma: 0.0049
M 322.34260              (2000.0)            P               Q
n   0.25622619     Peri.   72.47318     -0.50277681     -0.86441429
a   2.45501241     Node    47.71084      0.79213697     -0.46159019
e   0.5741560      Incl.    0.14296      0.34602670     -0.19930484
P   3.85           H   15.8     U  8.4     q 1.04545221  Q 3.86457261
From 13 observations 1997 Oct. 12-22; mean residual 0".485.

Orbital elements:
1997 ZZ99
   Perilune 1997 Apr 22.543629 TT = 13:02:49 (JD 2450561.043629)
Epoch 1997 Apr 22.5 TT = JDT 2450561.0                  Find_Orb
q 24606.5028km           (2000.0)            P               Q
H   22.9  G 0.15   Peri.  111.16531      0.57204494     -0.81965341
                   Node   193.85513     -0.79189676     -0.54220560
e 659.1509995      Incl.    7.32769     -0.21369156     -0.18488202
From 35 observations 1997 Apr. 21-22; mean residual 1".128.

Possible new format,  allowing inclusion of uncertainties:

Orbital elements: 1996 XX1
   Perihelion 1998 Mar 8.977252 +/- 0.109 TT (JD 2450881.477252)
Epoch 1997 Oct 12.0 TT = JDT 2450733.5   Earth MOID: 0.0613   Ma: 0.0049
M     322.109164 +/-   0.0574     (More MOIDs or area/mass goes here)
Peri.  72.530699 +/-   0.0861     n   0.256058518 +/- 0.000205
Node   47.657952 +/-   0.0576     a   2.456084084 +/- 0.00131
Incl.   0.143216 +/-   0.000294   e   0.57438774 +/- 0.000178
H   16.6    G 0.15     U  7.1   P 3.85012 +/- .00021
q 1.045339477 +/- 0.000139      Q 3.866828691 +/- 0.0025
From 13 of 17 observations 1997 Oct. 12-22; mean residual 0".485.

   ...and the hyperbolic version:

Orbital elements: 1997 ZZ99
   Perilune 1997 Apr 22.543629 +/- 0.021 TT (JD 2450561.043629)
Epoch 1997 Apr 22.5 TT = JDT 2450561.0                  Find_Orb
Peri.  111.16531 +/- 0.0321       (More MOIDs or area/mass goes here)
Node   193.85513 +/- 0.145        z 934.123 +/- 3.14159
Incl.    7.32769 +/- 0.0056       e 659.1509995 +/- 3.456
H   22.9  G 0.15
q 24606.5028 +/- 3.345 km
From 35 of 35 observations 1997 Apr. 21-22; mean residual 1".128.

   ....wherein Q, a, P, n, U, and M are omitted as no longer meaningful,
but we show z=1/a,  which is sometimes helpful with near-parabolic cases.
Actually,  we should always show z if the uncertainty in a is comparable
to a itself,  for all types of orbits.

   We basically retain the first three lines of the existing format,
except that the second line needs a little revision to contain the
periapsis uncertainty.  And we retain the "from x to y" line.

   We should also "prune" quantities so the number of digits in the sigma
matches the number of digits in the quantity itself. */

int write_out_elements_to_file( const double *orbit,
            const double curr_epoch,
            const double epoch_shown,
            OBSERVE FAR *obs, const int n_obs, const char *constraints,
            const int precision, const int monte_carlo,
            const int options)
{
   static const double jan_1970 = 2440587.5;
   const char *file_permits = (append_elements_to_element_file ? "a" : "w+");
   extern const char *elements_filename;
   extern double planet_mass[];
   FILE *ofile = fopen( elements_filename, file_permits);
   double rel_orbit[6], orbit2[6];
   int planet_orbiting, n_lines, i, bad_elements;
   ELEMENTS elem, helio_elem;
   char *tptr, *tbuff = (char *)malloc( 80 * 9);
   char object_name[80], buff[100], more_moids[80];
   char impact_buff[80];
   int n_more_moids = 0;
   int output_format = (precision | SHOWELEM_PERIH_TIME_MASK);
   extern int n_extra_params, perturbers;
   double jd;
   double moids[9];
   const char *monte_carlo_permits;
   const int rms_ok = (compute_rms( obs, n_obs) < max_monte_rms);
   static int n_clones_accepted = 0;

   if( default_comet_magnitude_type == 'N')
      output_format |= SHOWELEM_COMET_MAGS_NUCLEAR;
   fprintf( ofile, "Orbital elements: %s",
         options & ELEM_OUT_ALTERNATIVE_FORMAT ? " " : "\n");
   memcpy( orbit2, orbit, 6 * sizeof( double));
   integrate_orbit( orbit2, curr_epoch, epoch_shown);
   planet_orbiting = find_best_fit_planet( epoch_shown, orbit2, rel_orbit);
   if( options & ELEM_OUT_HELIOCENTRIC_ONLY)
      {
      planet_orbiting = 0;       /* insist on heliocentric display */
      memcpy( rel_orbit, orbit2, 6 * sizeof( double));
      }
   if( planet_orbiting == 3)
      {
      const double obliq_2000 = 23.4392911 * PI / 180.;

                           /* Cvt ecliptic to equatorial: */
      rotate_vector( rel_orbit, obliq_2000, 0);
      rotate_vector( rel_orbit + 3, obliq_2000, 0);
      }

   calc_classical_elements( &elem, rel_orbit, epoch_shown, 1,
                                    SOLAR_GM * planet_mass[planet_orbiting]);
   if( elem.ecc < .9)
      sprintf( orbit_summary_text, "a=%.3lf, ", elem.major_axis);
   else
      sprintf( orbit_summary_text, "q=%.3lf, ", elem.q);
   sprintf( orbit_summary_text + strlen( orbit_summary_text),
            "e=%.3lf, i=%d", elem.ecc, (int)( elem.incl * 180. / PI + .5));
   elem.central_obj = planet_orbiting;
   elem.epoch = epoch_shown;
   elem.abs_mag = calc_absolute_magnitude( obs, n_obs);
   elem.is_asteroid = !(obs->flags & OBS_IS_COMET);
   elem.slope_param = (elem.is_asteroid ? asteroid_magnitude_slope_param
                                        : comet_magnitude_slope_param);
   helio_elem = elem;
   helio_elem.central_obj = 0;
   calc_classical_elements( &helio_elem, orbit2, epoch_shown, 1,
                                    SOLAR_GM * planet_mass[0]);
   get_object_name( object_name, obs->packed_id);
   n_lines = elements_in_mpc_format( tbuff, &elem, object_name,
               is_cometary( constraints) && fabs( elem.ecc - 1.) < 1.e-8,
               output_format);
   tptr = tbuff;
   *more_moids = '\0';
   for( i = 0; i < 9; i++)
      moids[i] = 0.;
   for( i = 0; *tptr && i < n_lines; i++)
      {
      char *tt_ptr;
      extern double solar_pressure[], uncertainty_parameter;
      const double SRP1AU = 2.3e-7;
             /* "Solar radiation pressure at 1 AU",  in             */
             /* kg*AU^3 / (m^2*d^2),  from a private communication  */
             /* from Steve Chesley; see orb_func.cpp for details    */

      strcpy( buff, tptr);
      tt_ptr = strstr( buff, "TT") + 2;
      if( !memcmp( buff, "   Peri", 7))
         {
         if( *constraints)
            sprintf( tt_ptr, ";  Constraint: %s", constraints);
         else if( n_extra_params == 2 || n_extra_params == 3)
            {
            char tbuff0[40], tbuff1[40];

            sprintf( tt_ptr, "; A1=%s, A2=%s",
                     put_double_in_buff( tbuff0, solar_pressure[0]),
                     put_double_in_buff( tbuff1, solar_pressure[1]));
            if( n_extra_params == 3)
               sprintf( tt_ptr + strlen( tt_ptr), ", A3=%s",
                  put_double_in_buff( tbuff1, solar_pressure[2]));
            }
         else if( n_monte_carlo_impactors && monte_carlo)
            sprintf( tt_ptr, ";  %.2lf%% impact (%d/%d)",
                100. * (double)n_monte_carlo_impactors /
                       (double)monte_carlo_object_count,
                       n_monte_carlo_impactors, monte_carlo_object_count);
         }

      else if( !memcmp( buff, "Epoch", 5))
         {
         int j;

         if( n_extra_params == 1)
            sprintf( tt_ptr, "; AMR %.5lg m^2/kg",
                                 solar_pressure[0] * SOLAR_GM / SRP1AU);
         else if( !planet_orbiting)
            for( j = 0; n_more_moids < 5 && j < 8; j++)
               {
               static const char moid_idx[] = { 3, 5, 2, 1, 4, 6, 7, 8 };
               double moid;
               ELEMENTS planet_elem;

               setup_planet_elem( &planet_elem, moid_idx[j],
                                (epoch_shown - J2000) / 36525.);
               moid = find_moid( &elem, &planet_elem);
               moids[(int)moid_idx[j]] = moid;
               if( (j < 2 && moid < 1.) || (j > 4 && moid < 1.)
                           || moid < .1)
                  {
                  char addendum[30];
                  static const char *moid_text[] = { "Earth MOID", "Ju",
                           "Ve", "Me", "Ma", "Sa", "Ur", "Ne" };

                  sprintf( addendum, "   %s: %.4lf", moid_text[j], moid);
                  if( strlen( addendum) + strlen( buff) < 79)
                     strcat( buff, addendum);
                  else
                     {
                     n_more_moids++;
                     strcat( more_moids, addendum);
                     }
                  if( !j && moid < .1)
                      sprintf( orbit_summary_text + strlen( orbit_summary_text),
                         " MOID %.3lf", moid);
                  }
            }
         if( strlen( buff) < 56)
            {
            const char *reference = get_environment_ptr( "REFERENCE");

            if( !*reference)
               {
               reference = "Find_Orb";
               set_environment_ptr( "REFERENCE", reference);
               }
            for( j = strlen( buff); j < 56; j++)
               buff[j] = ' ';
            strcpy( buff + 56, reference);
            }
         }
      else if( *more_moids)
         {
         int output_location = 38;

         if( n_more_moids == 4)
            output_location = 25;
         if( n_more_moids == 5)
            output_location = 12;
         strcpy( buff + output_location, more_moids);
         *more_moids = '\0';
         }
      if( !memcmp( buff + 31, " G ", 3) && uncertainty_parameter < 90.)
         {                       /* replace slope w/uncertainty */
         sprintf( buff + 32, "U%5.1lf  ", uncertainty_parameter);
         buff[40] = ' ';
         }
      fprintf( ofile, "%s\n", buff);
      tptr += strlen( tptr) + 1;
      }
   observation_summary_data( tbuff, obs, n_obs, options);
   fprintf( ofile, "%s\n", tbuff);
   if( elem.central_obj == 3 && elem.ecc < .99)
      {
      write_tle_from_vector( tbuff, rel_orbit, elem.epoch, NULL, NULL);
      tbuff[69] = tbuff[140] = '\0';
      fprintf( ofile, "%s\n%s\n", tbuff, tbuff + 71);
      }
   if( !(options & ELEM_OUT_NO_COMMENT_DATA))
      {
      fprintf( ofile, "# State vector (heliocentric ecliptic J2000):\n");
      fprintf( ofile, "# %17.12lf%17.12lf%17.12lf AU\n",
               orbit2[0], orbit2[1], orbit2[2]);
      fprintf( ofile, "# %17.12lf%17.12lf%17.12lf mAU/day\n",
               orbit2[3] * 1000., orbit2[4] * 1000., orbit2[5] * 1000.);
      if( !planet_orbiting)
         {
         fprintf( ofile, "# MOIDs: Me%8.4lf Ve%8.4lf Ea%8.4lf Ma%8.4lf\n",
                  moids[1], moids[2], moids[3], moids[4]);
         fprintf( ofile, "# MOIDs: Ju%8.4lf Sa%8.4lf Ur%8.4lf Ne%8.4lf\n",
                  moids[5], moids[6], moids[7], moids[8]);
         }
      }
   if( monte_carlo)
      {
      if( !monte_carlo_object_count)
         n_clones_accepted = 0;
      if( rms_ok)
         n_clones_accepted++;
      }

   monte_carlo_permits = (n_clones_accepted == 1 ? "wb" : "ab");
   if( monte_carlo && rms_ok)
      {
      FILE *ofile2 = fopen( "state.txt", monte_carlo_permits);

      if( ofile2)
         {
         fseek( ofile, 0L, SEEK_SET);
         while( fgets( buff, sizeof( buff), ofile))
            fwrite( buff, strlen( buff), 1, ofile2);
         fclose( ofile2);
         fseek( ofile, 0L, SEEK_END);
         }
      }

   jd = jan_1970 + time( NULL) / 86400.;
   full_ctime( buff, jd, CALENDAR_JULIAN_GREGORIAN);
   if( !(options & ELEM_OUT_NO_COMMENT_DATA))
      {
      fprintf( ofile, "# Elements written: %s (JD %lf)\n", buff, jd);
      make_date_range_text( buff, obs[0].jd, obs[n_obs - 1].jd);
      fprintf( ofile, "# Full range of obs: %s (%d observations)\n",
                              buff, n_obs);
      fprintf( ofile, "# Find_Orb ver: %s %s\n", __DATE__, __TIME__);
      fprintf( ofile, "# Perturbers: %08lx ", (unsigned long)perturbers);
      if( !perturbers)
         fprintf( ofile, "(unperturbed orbit)\n");
      else if( (perturbers & 0x3fe) == 0x3fe)
         fprintf( ofile, (perturbers & 0x400) ? "(Merc-Pluto plus Luna)\n" :
               "(Merc-Pluto, Earth & moon combined)\n");
      else
         fprintf( ofile, "\n");

      if( !elem.central_obj && elem.ecc != 1.)
         for( i = 0; i < 3; i++)           /* show Tisserand data for Jup & Nep */
            {
            const double semimajor_axes[3] = { 1., 5.2033, 30.069 };
            const char *names[3] = { "Earth", "Jupiter", "Neptune" };
            const double ratio =  semimajor_axes[i] / elem.major_axis;
            const double tisserand = ratio
               + 2. * sqrt( (1. - elem.ecc * elem.ecc) / ratio) * cos( elem.incl);

            fprintf( ofile, "# Tisserand relative to %s: %.5lf\n",
                     names[i], tisserand);
            }
      }

   *impact_buff = '\0';
   if( elem.central_obj < 15)
      {
      double latlon[2];
      const int is_an_impact = (obs->jd < elem.perih_time);
                                         /* basically means,  "if we */
                                         /* observed the object after */
                                         /* periapsis, must be a launch; */
                                         /* otherwise,  must be impact." */
//    const double saved_mean_anomaly = elem.mean_anomaly;
      const double t0 = find_collision_time( &elem, latlon, is_an_impact);

//    elem.mean_anomaly = saved_mean_anomaly;
      if( t0 < 1.)      /* t0 = 1 -> it was a miss after all */
         {
         char *end_ptr;
         const double lon = latlon[0] * 180. / PI;
         const double impact_time_td = elem.perih_time + t0;
         const double impact_time_utc = impact_time_td -
                        td_minus_utc( impact_time_td) / 86400.;

         full_ctime( buff, impact_time_utc,
                       FULL_CTIME_HUNDREDTH_SEC | CALENDAR_JULIAN_GREGORIAN);
         sprintf( impact_buff, " %s lat %+9.5lf lon ", buff,
               latlon[1] * 180. / PI);
         end_ptr = impact_buff + strlen( impact_buff);
                     /* 0 < longitude < 360;  for Earth,  show this in */
                     /* "conventional" East/West 0-180 degree format:  */
         if( elem.central_obj == 3)
            {
            sprintf( end_ptr, "%c%.5lf",
                  (lon < 180. ? 'E' : 'W'),
                  (lon < 180. ? lon : 360. - lon));
            fprintf( ofile, "%s at %s\n", (is_an_impact ? "IMPACT" : "LAUNCH"),
                                   impact_buff);
            }
                     /* Then show in 0-360 format,  for all other  */
                     /* planets, and for output to file:           */
         sprintf( end_ptr, "%9.5lf", lon);
         if( elem.central_obj != 3)
            fprintf( ofile, "%s at %s\n", (is_an_impact ? "IMPACT" : "LAUNCH"),
                             impact_buff);
         }
      }
   if( *get_environment_ptr( "PLANET_STATES"))
      {
      fprintf( ofile, "# Planet states as of JD %.2lf:\n", elem.epoch);
      for( i = 1; i < 11; i++)
         {
         double loc[3], vel[3];

         get_planet_posn_vel( elem.epoch, i, loc, vel);
         fprintf( ofile, "#%02d: %17.12lf%17.12lf%17.12lf AU\n",
                   i, loc[0], loc[1], loc[2]);
         fprintf( ofile, "#    %17.12lf%17.12lf%17.12lf mAU/day\n",
                   vel[0] * 1000., vel[1] * 1000., vel[2] * 1000.);
         }
      }
   if( !(options & ELEM_OUT_NO_COMMENT_DATA))
      {
      char time_buff[40];

      sprintf( tbuff, "#  $Name=%s", object_name);
      text_search_and_replace( tbuff + 2, " ", "%20");
               /* Epoch has to be given in YYYYMMDD.DDDDD format: */
      full_ctime( time_buff, helio_elem.perih_time,
               FULL_CTIME_YMD | FULL_CTIME_MONTHS_AS_DIGITS
               | FULL_CTIME_MICRODAYS | FULL_CTIME_LEADING_ZEROES);
      time_buff[4] = time_buff[7] = '\0';
      sprintf( tbuff + strlen( tbuff), "  $Ty=%s  $Tm=%s  $Td=%s",
               time_buff, time_buff + 5, time_buff + 8);
      sprintf( tbuff + strlen( tbuff), "  $MA=%.5lf",
                  centralize_ang( helio_elem.mean_anomaly) * 180. / PI);
      fprintf( ofile, "%s\n", tbuff);

      sprintf( tbuff, "#  $ecc=%.7lf  $Eqnx=2000.", helio_elem.ecc);
      fprintf( ofile, "%s\n", tbuff);

      sprintf( tbuff, "#  $a=%.7lf  $Peri=%.5lf  $Node=%.5lf",
                  helio_elem.major_axis,
                  centralize_ang( helio_elem.arg_per) * 180. / PI,
                  centralize_ang( helio_elem.asc_node) * 180. / PI);
      sprintf( tbuff + strlen( tbuff), "  $Incl=%.5lf",
                  helio_elem.incl * 180. / PI);
      fprintf( ofile, "%s\n", tbuff);

      sprintf( tbuff, "#  $EpJD=%.3lf  $q=%.6lf", helio_elem.epoch, helio_elem.q);
      sprintf( tbuff + strlen( tbuff), "  $T=%.6lf  $H=%.1lf",
               helio_elem.perih_time, helio_elem.abs_mag);
      fprintf( ofile, "%s\n", tbuff);
      }
   fclose( ofile);
         /* Return value indicates probable trouble if the eccentricity      */
         /* is greater than 1.2 for an heliocentric orbit.  If that happens, */
         /* the orbital elements ought to be shown in,  say,  flashing red.  */
   bad_elements = ( helio_elem.ecc < 1.2 || helio_elem.central_obj ? 0 : -1);
   if( helio_elem.q > 90.)
      bad_elements |= 2;

         /* Also,  write out elements in MPCORB-like format: */
   if( helio_elem.ecc < .999999 || alt_mpcorb)
      {
      elements_in_mpcorb_format( tbuff, obs->packed_id, object_name,
                              &helio_elem, obs, n_obs);
      ofile = fopen( "mpc_fmt.txt", file_permits);
      fprintf( ofile, "%s\n", tbuff);
      fclose( ofile);
      }

   if( monte_carlo)
      monte_carlo_object_count++;
   if( monte_carlo && rms_ok)
      {
      const char *element_filename = get_environment_ptr( "MONTE_CARLO");
      char name_buff[48], virtual_full_desig[40];

      if( !*element_filename)
         element_filename = "mpcorb.dat";
      if( *impact_buff)
         n_monte_carlo_impactors++;
      sprintf( name_buff, "%05d", monte_carlo_object_count);
      packed_desig_minus_spaces( virtual_full_desig, obs->packed_id);
      sprintf( virtual_full_desig + strlen( virtual_full_desig), " [%d]",
                                  monte_carlo_object_count);
      if( elem.central_obj || elem.ecc > .999999)
         {
         ofile = fopen( "virtual.txt", monte_carlo_permits);
         elements_in_guide_format( tbuff, &elem, virtual_full_desig);
         fprintf( ofile, "%s%s\n", tbuff, impact_buff);
         fclose( ofile);
         }

      if( helio_elem.ecc < .999999 || alt_mpcorb)
         {
         ofile = fopen( element_filename, monte_carlo_permits);
         if( !strcmp( monte_carlo_permits, "wb"))
            {        /* new file = write out a header for it */
            FILE *ifile = fopen( "mpcorb.hdr", "rb");
            time_t t0 = time( NULL);

            fprintf( ofile, "Monte Carlo orbits from Find_Orb\nComputed %s", ctime( &t0));
            fprintf( ofile, "Find_Orb version %s %s\n", __DATE__, __TIME__);
            fprintf( ofile, (using_sr ? "Statistical Ranging\n" : "Full Monte Carlo\n"));
            if( ifile)
               {
               while( fgets( tbuff, 80 * 9, ifile))
                  fputs( tbuff, ofile);
               fclose( ifile);
               }
            }
         elements_in_mpcorb_format( tbuff, name_buff, virtual_full_desig,
                              &helio_elem, obs, n_obs);
         fprintf( ofile, "%s%s\n", tbuff, impact_buff);
         fclose( ofile);
         }
//    if( helio_elem.ecc > 2. || helio_elem.q > 90.)
//       set_statistical_ranging( 1);
      }
   free( tbuff);
   return( bad_elements);
}

static const char *vector_data_file = "vectors.dat";

void set_solutions_found( OBJECT_INFO *ids, const int n_ids)
{
   FILE *ifile = fopen( vector_data_file, "rb");
   char buff[120];
   int i;

   for( i = 0; i < n_ids; i++)
      ids[i].solution_exists = 0;
   if( ifile)
      {
      while( fgets_trimmed( buff, sizeof( buff), ifile))
         {
         for( i = 0; i < n_ids; i++)
            if( !strcmp( buff + 11, ids[i].obj_name))
               ids[i].solution_exists = 1;
         fgets_trimmed( buff, sizeof( buff), ifile);
         fgets_trimmed( buff, sizeof( buff), ifile);
         }
      fclose( ifile);
      }
}

int fetch_previous_solution( OBSERVE *obs, const int n_obs, double *orbit,
               double *orbit_epoch, int *perturbers)
{
   const int ignore_prev_solns = atoi( get_environment_ptr( "NO_PREV_SOLNS"));
   FILE *ifile = (ignore_prev_solns ? NULL : fopen( vector_data_file, "rb"));
   int got_vectors = 0, i;

   if( ifile)
      {
      char buff[120], object_name[80];
      double jd1 = 0., jd2 = 0.;

      get_object_name( object_name, obs->packed_id);
      while( fgets_trimmed( buff, sizeof( buff), ifile))
         if( !FMEMCMP( object_name, buff + 11, FSTRLEN( object_name)))
            if( buff[ FSTRLEN( object_name) + 11] < ' ' && *buff == ' ')
               {
               extern int n_extra_params;
               double unused_step_size;
               extern double solar_pressure[];

               for( i = 0; i < 3; i++)
                  solar_pressure[i] = 0.;
               got_vectors = 1;
               *orbit_epoch = atof( buff);
               *perturbers = 0;
               fgets_trimmed( buff, sizeof( buff), ifile);
               sscanf( buff, "%lf%lf%lf%x%lf %lf %lf",
                             &orbit[0], &orbit[1], &orbit[2], perturbers,
                             solar_pressure,
                             solar_pressure + 1,
                             solar_pressure + 2);
               n_extra_params = 0;
               if( solar_pressure[2])
                  n_extra_params = 3;
               else if( solar_pressure[1])
                  n_extra_params = 2;
               else if( solar_pressure[0])
                  n_extra_params = 1;
               fgets_trimmed( buff, sizeof( buff), ifile);
               sscanf( buff, "%lf%lf%lf%lf%lf%lf",
                             orbit + 3, orbit + 4, orbit + 5,
                             &unused_step_size, &jd1, &jd2);
               for( i = 3; i < 6; i++)
                  orbit[i] *= .001;
               for( i = 0; i < n_obs; i++)
                  {
                  obs[i].computed_ra  = obs[i].ra;
                  obs[i].computed_dec = obs[i].dec;
                  }
               }
      if( got_vectors)
         {
         set_locs( orbit, *orbit_epoch, obs, n_obs);
         if( jd2)
            for( i = 0; i < n_obs; i++)
               if( obs[i].jd < jd1 - .00001 || obs[i].jd > jd2 + .00001)
                  obs[i].is_included = 0;
         }
      fclose( ifile);
      }
   if( !got_vectors)
      {
      perturbers = 0;
      *orbit_epoch = initial_orbit( obs, n_obs, orbit);
      }
   return( got_vectors);
}

int store_solution( const OBSERVE *obs, const int n_obs, const double *orbit,
       const double orbit_epoch, const int perturbers)
{
   FILE *ofile = fopen( vector_data_file, "ab");

   if( ofile)
      {
      char buff[80];
      int i, j;

      get_object_name( buff, obs->packed_id);
      fprintf( ofile, "%10.1lf %s\n", orbit_epoch, buff);
      for( i = 0; i < 3; i++)
         fprintf( ofile, "%17.12lf", orbit[i]);
      if( perturbers)
         {
         extern int n_extra_params;

         fprintf( ofile, " %04x", perturbers);
         if( n_extra_params)
            {
            extern double solar_pressure[];

            for( i = 0; i < n_extra_params; i++)
               fprintf( ofile, " %.5lg", solar_pressure[i]);
            }
         }
      get_first_and_last_included_obs( obs, n_obs, &i, &j);
      fprintf( ofile, "\n%17.12lf%17.12lf%17.12lf %.3lf %.5lf %.5lf\n",
              orbit[3] * 1000., orbit[4] * 1000., orbit[5] * 1000.,
              0., obs[i].jd, obs[j].jd);
      fclose( ofile);
      }
   return( ofile ? 0 : -1);
}


#define LOG_10 2.302585

double calc_obs_magnitude( const int is_comet, const double obj_sun,
            const double obj_earth, const double earth_sun, double *phase_ang)
{
   double rval;
   double ph_ang = obj_sun * obj_sun +
                  obj_earth * obj_earth - earth_sun * earth_sun;

   ph_ang /= 2. * obj_earth * obj_sun;
   ph_ang = acose( ph_ang);
   if( phase_ang)
      *phase_ang = ph_ang;

   if( is_comet)
      rval = comet_magnitude_slope_param * log( obj_sun);
   else
      {
      double phi1, phi2, log_tan_half_phase;

      log_tan_half_phase = log( sin( ph_ang / 2.) / cos( ph_ang / 2.));
      phi1 = exp( -3.33 * exp( log_tan_half_phase * 0.63));
      phi2 = exp( -1.87 * exp( log_tan_half_phase * 1.22));
      rval = 5. * log( obj_sun) - 2.5 *
                  log( (1. - asteroid_magnitude_slope_param) * phi1
                + asteroid_magnitude_slope_param * phi2);
      }
   rval += 5. * log( obj_earth);
   rval /= LOG_10;         /* allow for common logs,  not naturals */
   return( rval);
}

/* The following function,  as the comment indicates,  assumes that a */
/* "no band" case (obs->mag_band == ' ') must be an R mag.  That's     */
/* probably the best guess for most modern,  unfiltered CCD             */
/* observations.  MPC assumes a B (photographic) magnitude,  which is   */
/* probably the best guess for older observations.  I suppose one would */
/* ideally look at the second 'note'  which is C for CCD observations   */
/* and P for photographic observations.  The code could then assume a  */
/* default of R for CCD obs,  B for photographic,  and V for the      */
/* admittedly rare micrometer or encoder-based observations.         */
/* http://www.cfa.harvard.edu/iau/info/OpticalObs.html              */

/* Estimates for the mag band shifts in R, B, I, and U are from a post by */
/* Petr Pravec,  http://tech.groups.yahoo.com/group/mpml/message/24833,   */
/* in turn derived from data in Shevchenko and Lupishko, Solar System     */
/* Research 32, 220-232, 1998.                                            */

/* Note that it looks as if MPC uses V-w = .2,  just based on NEOCP  */
/* ephemeris V magnitudes based on PanSTARRS w data;  but,  as described */
/* in the MPML post referenced below,  this appears to be wrong.  */

double mag_band_shift( const char mag_band)
{
   double rval = 0.;

#ifdef OLD_VALUES
   if( mag_band == 'R' || mag_band == ' ')    /* V-R=+0.3 */
      rval = .3;           /* assume "no band" = "R" */
   if( mag_band == 'B')    /* B-V=-0.8 */
      rval = -.8;
#endif
   if( mag_band == 'R' || mag_band == ' ')    /* V-R=+0.3 */
      rval = .43;          /* assume "no band" = "R" */
   if( mag_band == 'B')    /* B-V=-0.8 */
      rval = -.77;
   if( mag_band == 'I')
      rval = .82;
   if( mag_band == 'U')
      rval = 1.16;
   if( mag_band == 'C')    /* similar for "C" = "clear"? */
      rval = .4;
   if( mag_band == 'w')    /* V-w=.5 to 1.0,  according to */
      rval = -.75; /* http://tech.groups.yahoo.com/group/mpml/message/24808 */
   return( rval);
}

double calc_absolute_magnitude( OBSERVE FAR *obs, int n_obs)
{
   int obs_no;
   double n_mags = 0.;
   double rval = 0.;

   for( obs_no = 0; obs_no < n_obs; obs_no++)
      {
      obs->computed_mag = 0.;
      if( obs->r && obs->solar_r)
         {
         const double earth_sun = vector3_length( obs->obs_posn);

         if( earth_sun)
            {
            double weight = 1.;

            if( obs->mag_precision == -1)    /* integer magnitude */
               weight = .1;
            if( obs->mag_precision == 2)  /* mag to .01;  assume it's good */
               weight = 5.;
            if( (obs->flags & OBS_IS_COMET) && obs->mag_band != default_comet_magnitude_type)
               weight = 0.;
            if( !obs->obs_mag || !obs->is_included)
               weight = 0.;
            obs->computed_mag = calc_obs_magnitude( (obs->flags & OBS_IS_COMET),
                     obs->solar_r, obs->r, earth_sun, NULL) - mag_band_shift( obs->mag_band);
            rval += weight * (obs->obs_mag - obs->computed_mag);
            n_mags += weight;
            }
         }
      obs++;
      }
   if( n_mags)
      rval /= n_mags;
   obs -= n_obs;
   for( obs_no = 0; obs_no < n_obs; obs_no++, obs++)
      if( n_mags)
         obs->computed_mag += rval;
      else
         obs->computed_mag = 0.;
   return( rval);
}

int find_worst_observation( const OBSERVE FAR *obs, const int n_obs)
{
   int i, rval = -1;
   double worst_rms = 0., rms;

   for( i = 0; i < n_obs; i++, obs++)
      if( obs->is_included)
         {
         rms = compute_rms( obs, 1);
         if( rms > worst_rms)
            {
            worst_rms = rms;
            rval = i;
            }
         }
   return( rval);
}

/* If you've got n_obs observations stored in the obs array,  the
   get_idx1_and_idx2( ) function will puzzle through them to find the first
   and last valid observation (those that haven't had their 'is_included'
   flags set to FALSE),  and will store the indices to them in *idx1 and
   *idx2.  These are shown near the top of the display,  and are used in
   the method of Herget.  Return value is the number of included obs.   */

int get_idx1_and_idx2( const int n_obs, const OBSERVE FAR *obs,
                                          int *idx1, int *idx2)
{
   int i, rval = 0;

   for( i = 0; i < n_obs && (!obs[i].is_included || !obs[i].r); i++)
      ;
   if( i == n_obs)
      *idx1 = *idx2 = 0;
   else
      {
      *idx1 = i;
      for( ; i < n_obs; i++)
         if( obs[i].is_included && obs[i].r)
            rval++;
      for( i = n_obs - 1; i && (!obs[i].is_included || !obs[i].r); i--)
         ;
      *idx2 = i;
      }
   return( rval);
}

int get_r1_and_r2( const int n_obs, const OBSERVE FAR *obs,
                           double *r1, double *r2)
{
   int idx1, idx2, rval = get_idx1_and_idx2( n_obs, obs, &idx1, &idx2);

   if( !rval)
      *r1 = *r2 = 0.;
   else
      {
      *r1 = obs[idx1].r;
      *r2 = obs[idx2].r;
      }
   return( rval);
}

static const char *defaults = "DEFAULTS";

int store_defaults( const int ephemeris_output_options)
{
   char buff[50];

   sprintf( buff, "%c%x", default_comet_magnitude_type,
                                        ephemeris_output_options);
   set_environment_ptr( defaults, buff);
   return( 0);
}

int get_defaults( int *ephemeris_output_options)
{
   const char *def_values = get_environment_ptr( defaults);

   *ephemeris_output_options = 0;
   sscanf( def_values, "%c%x", &default_comet_magnitude_type,
                               ephemeris_output_options);
   return( 0);
}

/* The following functions are used to "color" observations in the
console versions of Find_Orb.  The idea resembles that of the
four-color map problem,  except in this case, we'd like to show
observations in max_n_colors,  such that adjacent observations from
different MPC codes show up in different colors. You can't always do
this.  For example,  with observations from eight observatories,
mixed up so that each "pair" occurs,  you'd obviously need eight
different colors.  This code just tries a lot of possible colorings
and returns the one resulting in the fewest mismatches.

   To do this,  it uses an "annealing" sort of algorithm:  it first
sets the color for each MPC observatory at random (a value from zero
to max_n_colors - 1).  It then uses the improve_mpc_colors() to get
a better solution;  that function can make "obvious" improvements,
such as "if we change this MPC code to red,  there will be fewer
mismatches".  If the result has no mismatches,  we're home free and
stop looking for a "better" solution.  Otherwise,  we set a new set
of random colors and try to improve that... and repeat the procedure
for up to two seconds;  it's probably not worth spending much more time
on it than that.

   There are other ways to do this,  of course,  including a formal
pruned tree search among all possible color combinations.  But this
appears to work quite well,  and I thought it would result in simpler
code.  (I'm no longer so sure of that.  But I don't think I'll spend
the time to write a tree search version.)
*/

#define NO_MPC_COLOR_SET   -1

int find_mpc_color( const MPC_STATION *sdata, const char *mpc_code)
{
   int rval = NO_MPC_COLOR_SET;

   if( !mpc_code)       /* indicates 'just count colors */
      {
      rval = 0;
      while( sdata[rval].code[0])
         rval++;
      }
   else while( rval == NO_MPC_COLOR_SET && sdata->code[0])
      {
      if( mpc_code[0] == sdata->code[0] &&
             mpc_code[1] == sdata->code[1] &&
             mpc_code[2] == sdata->code[2])
         rval = sdata->color;
      sdata++;
      }
   return( rval);
}

static void set_mpc_colors_semirandomly( MPC_STATION *sdata,
               const int max_n_colors, unsigned long seed)
{
   srand( seed);
   while( sdata->code[0])
      {
      sdata->color = (char)( rand( ) % (unsigned long)max_n_colors);
      sdata++;
      }
}

/* After setting the colors at random,  we look for "simple" improvements:
for each MPC code,  we check the adjacent observations with different
MPC codes,  and see what colors they have.  We might see that (say)
there are four red neighbors,  three green,  and one blue.  In that case,
changing the color of the current MPC code to blue would result in only
one problem case,  instead of three or four.  (Ideally,  we'll find that
there are _no_ blue neighbors,  of course.)

   We keep trying this until no color changes are made.
*/

static void improve_mpc_colors( const int n_obs, const OBSERVE FAR *obs,
                   const int max_n_colors, MPC_STATION *sdata)
{
   int i, changes_made = 1, n_iterations = 0;

   while( changes_made)
      {
      changes_made = 0;
      for( i = 0; sdata[i].code[0]; i++)
         {
         int counts[20], j, color = sdata[i].color;

         assert( color >=0 && color < max_n_colors);
         for( j = 0; j < max_n_colors; j++)
            counts[j] = 0;
         for( j = 0; j < n_obs; j++)
            if( !strcmp( obs[j].mpc_code, sdata[i].code))
               {
               if( j > 0 && strcmp( obs[j - 1].mpc_code, sdata[i].code))
                  {
                  const int adjacent_color =
                          find_mpc_color( sdata, obs[j - 1].mpc_code);

                  assert( adjacent_color >=0 && adjacent_color < max_n_colors);
                  counts[adjacent_color]++;
                  }
               if( j < n_obs - 1 && strcmp( obs[j + 1].mpc_code, sdata[i].code))
                  {
                  const int adjacent_color =
                          find_mpc_color( sdata, obs[j + 1].mpc_code);

                  assert( adjacent_color >=0 && adjacent_color < max_n_colors);
                  counts[adjacent_color]++;
                  }
               }
         for( j = 0; j < max_n_colors; j++)
            if( counts[j] < counts[color])
               {
               color = j;
               sdata[i].color = (char)color;
               changes_made = 1;
               }
         assert( color >=0 && color < max_n_colors);
         }
      n_iterations++;
      assert( n_iterations < 10);
      }
}

extern int debug_level;
int debug_printf( const char *format, ...);                /* runge.cpp */

MPC_STATION *find_mpc_color_codes( const int n_obs, const OBSERVE FAR *obs,
                   const int max_n_colors)
{
   int n_codes = 0, i, j, n_alloced = 10;
   int best_score = 99999, best_seed = 0;
   clock_t t0;
   MPC_STATION *rval =
               (MPC_STATION *)calloc( n_alloced + 1, sizeof( MPC_STATION));

   for( i = 0; i < n_obs; i++)
      if( find_mpc_color( rval, obs[i].mpc_code) == NO_MPC_COLOR_SET)
         {
         int loc = 0;

         if( n_codes == n_alloced)
            {
            const int new_size = n_alloced * 2;
            MPC_STATION *new_array =
                   (MPC_STATION *)calloc( new_size + 1, sizeof( MPC_STATION));

            memcpy( new_array, rval, n_alloced * sizeof( MPC_STATION));
            n_alloced = new_size;
            free( rval);
            rval = new_array;
            }
         for( loc = 0; loc < n_codes
              && strcmp( rval[loc].code, obs[i].mpc_code) < 0; loc++)
            ;
                     /* move the rest of the array forward... */
         memmove( rval + loc + 1, rval + loc,
                          (n_codes - loc) * sizeof( MPC_STATION));
                     /* ...so we can copy in the new code: */
         strcpy( rval[loc].code, obs[i].mpc_code);
         n_codes++;
         }
   if( debug_level)
      {
      debug_printf( "%d obs, %d codes\n", n_obs, n_codes);
      for( i = 0; i < n_codes; i++)
         debug_printf( "%d: '%s'\n", i, rval[i].code);
      }
   t0 = clock( );
   for( i = 1; best_score && clock( ) < t0 + 2 * CLOCKS_PER_SEC; i++)
      {
      int score = 0;

      set_mpc_colors_semirandomly( rval, max_n_colors, (unsigned long)i);
      improve_mpc_colors( n_obs, obs, max_n_colors, rval);
      for( j = 0; j < n_obs - 1; j++)
         if( strcmp( obs[j].mpc_code, obs[j + 1].mpc_code))
            if( find_mpc_color( rval, obs[j].mpc_code) ==
               find_mpc_color( rval, obs[j + 1].mpc_code))
                  score++;
      if( score < best_score)   /* "lower" is "better" */
         {
         best_score = score;
         best_seed = i;
         }

      if( debug_level > 1 && (i < 10 || !(i % 100)))
          debug_printf( "Seed: %d, score %d, best %d\n", i, score, best_score);
      }
   if( debug_level)
      debug_printf( "Color setting: best score %d, i = %d\n", best_score, i);
   set_mpc_colors_semirandomly( rval, max_n_colors, (unsigned long)best_seed);
   improve_mpc_colors( n_obs, obs, max_n_colors, rval);
   return( rval);
}
