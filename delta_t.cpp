/* delta_t.cpp: computes Delta-T (UT-TD) for time system conversions

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
#include <stdint.h>
#include "watdefs.h"
#include "afuncs.h"

/* The following method for computing Delta-T works as follows.  For dates
between 1620 and 2012,  the value comes from the 'delta_t_table' (copied
in turn from Meeus' _Astronomical Algorithms_),  extended using data from
the USNO sites referenced below.  The table gives delta-t at two-year
intervals;  between those values,  a linear interpolation is used.  If
the year is before 1620,  one of two quadratic approximations is used
(see references below).   If the date is after the end of the table,
a linear extrapolation is used,  with a quadratic term added to that
which assumes an acceleration of 32.5 seconds/century^2.

   Updated 2 June 2001, 4 Jan 2005, 1 Nov 2007, 9 Jul 2008, 2010 Sep 19,
2012 May 3 to include new Delta-t data from:

ftp://maia.usno.navy.mil/ser7/deltat.data
ftp://maia.usno.navy.mil/ser7/deltat.preds
*/

#define DELTA_T_TABLE_SIZE 197

static const short delta_t_table[DELTA_T_TABLE_SIZE] =
 {12400, 11500, 10600, 9800, 9100, 8500, 7900,  /*  1620-1632 */
   7400, 7000, 6500, 6200, 5800, 5500, 5300,    /*  1634-1646 */
   5000, 4800, 4600, 4400, 4200, 4000, 3700,    /*  1648-1660 */
   3500, 3300, 3100, 2800, 2600, 2400, 2200,    /*  1662-1674 */
   2000, 1800, 1600, 1400, 1300, 1200, 1100,    /*  1676-1688 */
   1000,  900,  900,  900,  900,  900,  900,    /*  1690-1702 */
    900,  900, 1000, 1000, 1000, 1000, 1000,    /*  1704-1716 */
   1100, 1100, 1100, 1100, 1100, 1100, 1100,    /*  1718-1730 */
   1100, 1200, 1200, 1200, 1200, 1200, 1300,    /*  1732-1744 */
   1300, 1300, 1300, 1400, 1400, 1400, 1500,    /*  1746-1758 */
   1500, 1500, 1500, 1600, 1600, 1600, 1600,    /*  1760-1772 */
   1600, 1700, 1700, 1700, 1700, 1700, 1700,    /*  1774-1786 */
   1700, 1700, 1600, 1600, 1500, 1400, 1370,    /*  1788-1800 */
   1310, 1270, 1250, 1250, 1250, 1250, 1250,    /*  1802-1814 */
   1250, 1230, 1200, 1140, 1060,  960,  860,    /*  1816-1828 */
    750,  660,  600,  570,  560,  570,  590,    /*  1830-1842 */
    620,  650,  680,  710,  730,  750,  770,    /*  1844-1856 */
    780,  790,  750,  640,  540,  290,  160,    /*  1858-1870 */
   -100, -270, -360, -470, -540, -520, -550,    /*  1872-1884 */
   -560, -580, -590, -620, -640, -610, -470,    /*  1886-1898 */
   -270,    0,  260,  540,  770, 1050, 1340,    /*  1900-1912 */
   1600, 1820, 2020, 2120, 2240, 2350, 2390,    /*  1914-1926 */
   2430, 2400, 2390, 2390, 2370, 2400, 2430,    /*  1928-1940 */
   2530, 2620, 2730, 2820, 2910, 3000, 3070,    /*  1942-1954 */
   3140, 3220, 3310, 3400, 3500, 3650, 3830,    /*  1956-1968 */
   4020, 4220, 4450, 4650, 4850, 5050, 5220,    /*  1970-1982 */
   5380, 5490, 5580,                            /*  1984-1988 */
  5686,    /* 1990  1 1    56.8554     .3286        -24.6714 */
  5831,    /* 1992  1 1    58.3093    -.1253        -26.1253 */
  5998,    /* 1994  1 1    59.9847     .1993        -27.8007 */
  6163,    /* 1996  1 1    61.6287     .5553        -29.4447 */
  6297,    /* 1998  1 1    62.9659     .2181        -30.7819 */
  6383,    /* 2000  1 1    63.8285268  .3554732     -31.6445268 */
  6430,    /* 2002  1 1    64.2998152 -.1158152     -32.1158152 */
  6457,    /* 2004  1 1    64.5736400 -.3896400     -32.3896400 */
  6485,    /* 2006  1 1    64.8452                              */
  6546,    /* 2008  1 1:   65.4574                              */
  6607,    /* 2010  1 1:   66.0699                              */
  6660 };  /* 2012  1 1:   66.6030                              */

/* 8 Aug 2000:  Some people have expressed an interest in being able to
   insert their own formulae for Delta-T while running Guide.  I've
   enabled this by letting people add a string of coefficients in
   GUIDE.DAT;  gory details of how the string works are at

http://www.projectpluto.com/update7.htm#delta_t_def

   The following functions parse out the string of coefficients and
   let you set which string is used. */

static int evaluate_delta_t_string( const double year, double *delta_t,
                              const char *dt_string)
{
   int rval = -1;
   const char *tptr;

   for( tptr = dt_string; rval && *tptr; tptr++)
      if( tptr == dt_string || tptr[-1] == ';')
         {
         double year1, year2;
         int bytes_read;

         if( sscanf( tptr, "%lf,%lf:%n", &year1, &year2, &bytes_read) == 2
                  && year >= year1 && year < year2)
            {
            double offset = 2000., x, power = 1., coeff;

            tptr += bytes_read;
            if( *tptr == 'o')
               {
               tptr++;
               sscanf( tptr, "%lf,%n", &offset, &bytes_read);
               tptr += bytes_read;
               }
            x = (year - offset) / 100.;
            *delta_t = 0.;
            rval = 0;
            while( tptr[-1] != ';' && *tptr
                     && sscanf( tptr, "%lf%n", &coeff, &bytes_read))
               {
               *delta_t += power * coeff;
               power *= x;
               tptr += bytes_read;
               if( *tptr == ',')
                  tptr++;
               }
            tptr--;
            }
         }
   return( rval);
}

static const char *default_delta_t_string =
"-100000,-500:o1820,-20,0,32;\
-500,500:o0,10583.6,-1014.41,33.78311,-5.952053,-.1798452,.022174192,.0090316521;\
500,1600:o1000,1574.2,-556.01,71.23472,.319781,-.8503463,-.005050998,.0083572073;\
1600,1620:o1600,120,-98.08,-153.2,140.272";
// static const char *td_minus_dt_string = default_delta_t_string;
static const char *td_minus_dt_string = NULL;

/* If a user attempts to set a NULL or "blank" Delta-T definition,   */
/* we fall back on the above default_delta_t_string.                 */

void DLL_FUNC reset_td_minus_dt_string( const char *string)
{
   td_minus_dt_string = (string && *string ? string :
                        default_delta_t_string);
}


/* Niels Bohr noted that "prediction is hard,  especially about the future."
This is especially true for predicting Delta-T.  Morrison and Stephenson
arrived at a formula which is intended to match the long-term behavior of
Delta-T,  but it doesn't match current observations,  so it's pretty bad
for near-term results (say,  the next decade or two).  This code assumes
that the rate of change of Delta-T can be determined from the last two
entries in the above 'delta_t_table' array, with a quadratic term of 32.5
seconds/century^2 added to this.  This is likely to provide good
"near-term" results.  I don't know of anything likely to provide good
longer-term results.

   Prior to March 2012,  I used quadratic approximations for Delta-T for
years before 948 and years between 948 to 1620,  from F.R. Stephenson and
M. A. Houlden, "Atlas of Historical Eclipse Maps",  Cambridge University
Press,  England (1986), page x.  These can also be found on page 73 of
Meeus' _Astronomical Algorithms_,  1st edition.  They correspond to the
Delta-T string:

-10000,948:2715.6,573.36,46.5;948,1620:50.6,67.5,22.5

   Since then,  I've switched to quadratic and higher-order polynomials
from the _Five Millennium Catalogue of Solar Eclipses_ by Espenak and
Meeus,  as defined in the above default_delta_t_string.  */

double DLL_FUNC td_minus_ut( const double jd)
{
   double year, dt, rval;
   /* first convert t from JD to years */

   year = 2000. + (jd - 2451545.) / 365.25;
   if( td_minus_dt_string)
      if( !evaluate_delta_t_string( year, &rval, td_minus_dt_string))
         return( rval);
   if( !evaluate_delta_t_string( year, &rval, default_delta_t_string))
         return( rval);
   if( year < 0.)          /* fall back on simple quadratic */
      {
      dt = year - 1820.;
      return( .0032 * dt * dt - 20.);
      }
   dt = (year - 2000.) / 100.;
#ifdef NOW_OBSOLETE_MORRISON_STEPHENSON_FORMULA
   if( year < 948.)
      rval = 2715.6 + dt * (573.36 + 46.5 * dt);
   else if( year < 1620.)
      rval = 50.6 + dt * (67.5 + 22.5 * dt);
   else     /* from 1620 to +infinity */
#endif
      {
      double index_loc = (year - 1620.) / 2.;
      int index = (int)index_loc;
      const short *tptr;

      if( index > DELTA_T_TABLE_SIZE - 2)       /* running off end of table */
         index = DELTA_T_TABLE_SIZE - 2;
      dt = index_loc - (double)index;
      tptr = delta_t_table + index;
      rval = (double)tptr[0] + (double)(tptr[1] - tptr[0]) * dt;
      rval /= 100.;
      if( dt > 1.)            /* again, past end of table */
         {
         dt = (dt - 1.) / 50.;    /* cvt to centuries past end,  and add the */
         rval += 32.5 * dt * dt;  /* same 32.5 sec/cy^2 used by Stephenson   */
         }
      }
#ifdef TRY_OMITTING
   if( year < 1620.)       /* apply correction from _Astro Almanac_, */
      {                    /* 1991,  page K8: */
      const double n = -23.8946;   /* corrected lunar secular acceleration */

      rval -= 0.000091 * (n + 26.) * (year - 1955.) * (year - 1955.);
      }
#endif
   return( rval);
}

/* Some notes for other time systems that I may,  someday,  get around
to implementing:

TDT = TAI + 32.184 seconds

TAI & UTC differ by an integer number of seconds;  after 2012 Jul 1,
for example,  TAI-UTC = 35 seconds,  TDT-UTC = 67.184 seconds.
USNO & IERS bulletins often give TDT-UT1;  from this you can get
UT1-UTC (always between -1 and +1 second)

TCB - TDB = L * (jd - 2443144.5) * 86400,
   where L=1.550505e-08

TAI-GPS = 19 seconds ("fixed" as of 1980;  at that point,  GPS time was
identical to UTC.  But since then,  16 leap seconds have been added to
UTC (as of 2012),  so that GPS and UTC have slipped apart by that much.)

Note also that as of January 2012,  the last official leap second was
for July 2012.  Beyond December 2012,  UTC is indeterminate (until leap
seconds are announced,  or it's announced that they won't be inserted.)

Also note that leap seconds past 2038 April 22 = MJD 65535 will require
the 'leap_intervals' array to be switched from 16-bit integers to 32 bits.

Previously,  this function assumed no leap seconds except those which had
been officially announced.  The problem is that this causes UTC to deviate
from (predicted) UT.  Thus,  in computing TD-UTC for future dates (those
beyond the end of the leap second table),  we find out which in which
Gregorian half-year jd_utc falls,  then compute the current prediction for
Delta-T in the middle of that year.  Then

TD-UTC = td_minus_tai + floor( td_minus_ut + .5 - td_minus_tai)

   ...i.e.,  choose the value for TD-UTC for that half-year to match
Delta-T as closely as possible (within half a second at the center of
that half-year).

   You may not like this so-called "solution".  I don't,  either.  I'm
open to ideas,  but doubt there is a better "solution".

   Also,  note that before 1961 Jan 1 = MJD 37400, there was no "real" UTC,
 and we assume UTC=UT.  From then until 1972 Jan 1 = MJD 41317, the idea
was that UTC followed atomic clocks,  but at a slightly different rate to
keep UTC mostly in sync with the earth's rotation.  The end result is that
TD-UTC could be expressed as a linear relationship (or a series of 13
linear relationships given below,  each covering a given time span.)

   After 1972 Jan 1,  things switched over to the current system:
TAI-UTC must always be an integer,  with leap seconds inserted at
irregular intervals,  always at the end of June or December,  to keep
UTC within .9 seconds of UT.   */

#define N_LEAP_SECONDS 26

/* These macros determine the MJD of the given date in 'YEAR'.  */
/* They're valid for _non-negative_ years in the _Gregorian_ calendar. */

#define JAN_1( YEAR) (((YEAR) * 365 + ((YEAR) - 1) / 4 - ((YEAR) - 1) / 100 \
                         + ((YEAR) - 1) / 400) - 678940)
#define FEB_1( YEAR) (JAN_1( YEAR) + 31)
#define MAR_1( YEAR) (((YEAR)*365 + (YEAR)/4 - (YEAR)/100 + (YEAR)/400) - 678881)
#define APR_1( YEAR) (MAR_1( YEAR) + 31)
#define MAY_1( YEAR) (APR_1( YEAR) + 30)
#define JUN_1( YEAR) (MAY_1( YEAR) + 31)
#define JUL_1( YEAR) (JUN_1( YEAR) + 30)
#define AUG_1( YEAR) (JUL_1( YEAR) + 31)
#define SEP_1( YEAR) (AUG_1( YEAR) + 31)
#define OCT_1( YEAR) (SEP_1( YEAR) + 30)
#define NOV_1( YEAR) (OCT_1( YEAR) + 31)
#define DEC_1( YEAR) (NOV_1( YEAR) + 30)

double DLL_FUNC td_minus_utc( const double jd_utc)
{
   const double tdt_minus_tai = 32.184;
   const double mjd_utc = jd_utc - 2400000.5;
   int i;

   if( mjd_utc < (double)JAN_1( 1972))  /* between jan 1961 & dec 1972 */
      {
      static const unsigned short ranges[13] =  { JAN_1( 1961), AUG_1( 1961),
                      JAN_1( 1962), NOV_1( 1963), JAN_1( 1964), APR_1( 1964),
                      SEP_1( 1964), JAN_1( 1965), MAR_1( 1965), JUL_1( 1965),
                      SEP_1( 1965), JAN_1( 1966), FEB_1( 1968) };

      for( i = 12; i >= 0; i--)
         if( mjd_utc >= (double)ranges[i])
            {
            static const double offset[13] = {
                    1.4228180 - JAN_1( 1961) * 0.0012960,
                    1.3728180 - JAN_1( 1961) * 0.0012960,
                    1.8458580 - JAN_1( 1962) * 0.0011232,
                    1.9458580 - JAN_1( 1962) * 0.0011232,
                    3.2401300 - JAN_1( 1965) * 0.0012960,
                    3.3401300 - JAN_1( 1965) * 0.0012960,
                    3.4401300 - JAN_1( 1965) * 0.0012960,
                    3.5401300 - JAN_1( 1965) * 0.0012960,
                    3.6401300 - JAN_1( 1965) * 0.0012960,
                    3.7401300 - JAN_1( 1965) * 0.0012960,
                    3.8401300 - JAN_1( 1965) * 0.0012960,
                    4.3131700 - JAN_1( 1966) * 0.0025920,
                    4.2131700 - JAN_1( 1966) * 0.0025920  };
            static const short scale[13] =      { 12960, 12960, 11232,
                                    11232, 12960, 12960, 12960, 12960,
                                    12960, 12960, 12960, 25920, 25920 };
            const double tai_minus_utc = offset[i] +
                        mjd_utc * (double)scale[i] * 1.e-7;

            return( tdt_minus_tai + tai_minus_utc);
            }
      }
   else              /* integral leap seconds */
      {
      const int imjd_utc = (int)mjd_utc;
      static const unsigned short leap_intervals[N_LEAP_SECONDS] = {
                 JAN_1( 1972), JUL_1( 1972), JAN_1( 1973),
                 JAN_1( 1974), JAN_1( 1975), JAN_1( 1976),
                 JAN_1( 1977), JAN_1( 1978), JAN_1( 1979),
                 JAN_1( 1980), JUL_1( 1981), JUL_1( 1982),
                 JUL_1( 1983), JUL_1( 1985), JAN_1( 1988),
                 JAN_1( 1990), JAN_1( 1991), JUL_1( 1992),
                 JUL_1( 1993), JUL_1( 1994), JAN_1( 1996),
                 JUL_1( 1997), JAN_1( 1999), JAN_1( 2006),
                 JAN_1( 2009), JUL_1( 2012) };

      if( imjd_utc >= JUL_1( 2012))
         {
         int day = imjd_utc + 2400000 - 1721058;
         int year = (int)( (int64_t)day * (int64_t)400 / (int64_t)146097);
         int low, high, july_1;

         low = JAN_1( year);     /* The above value for 'year' is correct */
         if( imjd_utc < low)     /* more than 99% of the time.  But we    */
            {                    /* may find,  for 31 December,  that     */
            year--;              /* it's too high by one year.            */
            high = low;
            low = JAN_1( year);
            }
         else
            high = JAN_1( year + 1);
                      /*  jul  aug  sep  oct  nov  dec.. jul 1 is exactly 184 */
         july_1 = high - (31 + 31 + 30 + 31 + 30 + 31); /* days before jan 1  */
         if( imjd_utc < july_1) /* in first half of the year */
            high = july_1;
         else                /* in second half of the year */
            low = july_1;
         return( tdt_minus_tai + floor( td_minus_ut( 2400000.5 +
                        (double)( low + high) * .5) + .5 - tdt_minus_tai));
         }
      for( i = N_LEAP_SECONDS - 1; i >= 0; i--)
         if( imjd_utc >= (int)leap_intervals[i])
            return( (double)(i + 10) + tdt_minus_tai);
      }
                     /* still here?  Must be before jan 1961,  so UTC = UT1: */
   return( td_minus_ut( jd_utc));
}
