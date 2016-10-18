/* get_time.cpp: functions for parsing time/date text

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include "watdefs.h"
#include "date.h"

double DLL_FUNC split_time( double jd, long *year, int *month, int *day,
                                 int *hr, int *min, const int time_format)
{
   long int_jd;
   double seconds;

   jd += .5;
   int_jd = (long)floor( jd);
   jd = (jd - (double)int_jd) * 1440.;     /* t is now in minutes */
   *min = (int)jd;
   if( *min == 1440)       /* evade rounding errors: */
      {
      int_jd++;
      jd = 0.;
      *min = 0;
      }
   seconds = (jd - (double)*min) * 60.;
   day_to_dmy( int_jd, day, month, year, time_format & CALENDAR_MASK);
   *hr = *min / 60;
   *min %= 60;
   return( seconds);
}

#ifndef memicmp
#ifndef __WATCOMC__
static int memicmp( const char *s1, const char *s2, int n)
{
   int c1, c2;

   while( n--)
      {
      if( (c1 = tolower( *s1++)) != (c2 = tolower( *s2++)))
         return( c1 - c2);
      }
   return( 0);
}
#endif
#endif


/* Given a month name, month_name_to_index() will return its index from
1 to 12 (13 for thirteen-month-year calendars).  Three letters assure the
match,  at least for the Julian/Gregorian calendar,  but one can use a
partial month name such as "Fe" or even "O".  If the fragment isn't
unique,  you'll get whichever month comes first;  for example,  "M" = March,
"Ju" = June.   */

static int month_name_to_index( const char *str)
{
   int len = strlen( str);
   int rval = 0, i;
   const char *month_text;

   if( len > 3)      /* compare up to,  but no more than,  three bytes */
      len = 3;
   for( i = 1; !rval && i <= 13; i++)
      if( (month_text = set_month_name( i, NULL)) != NULL)
         if( !memicmp( month_text, str, len))
            rval = i;
   return( rval);
}

static int day_of_week_name_to_index( const char *str)
{
   int len = strlen( str);
   int rval = -1, i;
   const char *dow_text;

   if( len > 3)      /* compare up to,  but no more than,  three bytes */
      len = 3;
   for( i = 0; rval == -1 && i < 7; i++)
      if( (dow_text = set_day_of_week_name( i, NULL)) != NULL)
         if( !memicmp( dow_text, str, len))
            rval = i;
   return( rval);
}


/* A time string can have one or more time offsets at the end of it, such */
/* as '-10m' or '+3h'.  For example,  '13:14 -10m +3h' would mean 16:04.  */
/* The collect_time_offset() looks for the last such offset, inserts      */
/* a '\0' to remove it from the string,  and returns the offset in days.  */
/* The get_time_from_string() function calls it until 0 is returned.      */

static inline double collect_time_offset( char *istr)
{
   static const char *symbols = "smhdwlyc";
   static const double scales[8] = { 1. / 86400., 1. / 1440., 1. / 24.,
                     1., 7., 29.530588853, 365.25, 36525. };
   int  bytes_scanned;
   int len = strlen( istr);
   double rval = 0.;

   if( len > 1)
      {
      const char ending_char = istr[len - 1];

      if( strchr( symbols, ending_char) && istr[len - 2] == ' ')
         {
         char unused_char;
         int i;

         while( len >= 0 && istr[len] != '-' && istr[len] != '+')
            len--;
         if( len >= 0)         /* we found something... */
            if( sscanf( istr + len, "%lf %c%n", &rval, &unused_char,
                                    &bytes_scanned) == 2)
               if( !istr[len + bytes_scanned])
                  {
                  for( i = 0; symbols[i]; i++)
                     if( ending_char == symbols[i])
                        rval *= scales[i];
                  while( len && istr[len - 1] == ' ')
                     len--;
                  istr[len] = '\0';
                  }
         }
      else if( *istr == '-' || *istr == '+')  /* just add/subtract N days: */
         {
         if( sscanf( istr, "%lf%n", &rval, &bytes_scanned) == 1 &&
                  bytes_scanned == len)
            *istr = '\0';
         else
            rval = 0.;
         }
      }
   return( rval);
}

static size_t remove_trailing_spaces( char *istr)
{
   size_t len = strlen( istr);

   while( len && istr[len - 1] == ' ')
      len--;
   istr[len] = '\0';
   return( len);
}

static size_t remove_leading_and_trailing_spaces( char *istr)
{
   size_t i;

   for( i = 0; istr[i] == ' '; i++)
      ;
   if( i)
      memmove( istr, istr + i, strlen( istr + i) + 1);
   return( remove_trailing_spaces( istr));
}

static char *remove_substring( char *timestr, const char *substring)
{
   char *rval = strstr( timestr, substring);

   if( rval)
      {
      const size_t sublen = strlen( substring);

      memmove( rval, rval + sublen, strlen( rval + sublen) + 1);
      }
   return( rval);
}

/* You can enter times with an 'AD' or 'BC' anywhere in the string. In    */
/* the former case,  it's completely ignored and has no effect on the     */
/* result.  In the latter,  the year is negated and one added,  to bring  */
/* things in line with the astronomical convention:  for example,  the    */
/* year historians call BC 45 would be known,  to an astronomer,  as -44, */
/* because historians don't know anything about zeroes.                   */

static int check_for_bc( char *timestr)
{
   char *text_loc;

   remove_substring( timestr, "ad");     /* Just grab the 'ad' or 'a.d.';  */
   remove_substring( timestr, "a.d.");   /* they have no actual effect     */
   text_loc = remove_substring( timestr, "bc");
   if( !text_loc)
      text_loc = remove_substring( timestr, "b.c.");

   remove_leading_and_trailing_spaces( timestr);
   return( text_loc ? 1 : 0);
}

/* One can enter 'nm' followed by a number of days to get a particular */
/* lunar age;  e.g.,  'nm-3' would get the nearest time that is three  */
/* days before new moon,  'nm5.5' the nearest time 5.5 days after a    */
/* new moon,  and so on.  Or one can use 'fm-3' for three days before  */
/* full moon,  '1q-2' for two days before first quarter... or just use */
/* 'nm', '3q',  etc.  Uses formulae from Meeus, _Astronomical          */
/* Algorithms_,  chap 47, for a very approximate lunar age (I ignored  */
/* terms greater than about half an hour.)                             */

#define PHASE_IDX_UNDEFINED        -1
#define PHASE_IDX_NEW_MOON          0
#define PHASE_IDX_FIRST_QUARTER     1
#define PHASE_IDX_FULL_MOON         2
#define PHASE_IDX_THIRD_QUARTER     3

static int get_phase_idx( const char *istr)
{
   int i, rval = PHASE_IDX_UNDEFINED;

   for( i = 0; i < 4; i++)
      if( istr[0] == "n1f3"[i] && istr[1] == "mqmq"[i])
         rval = i;
   return( rval);
}

static double set_from_lunar( const int phase_idx, const double initial_jd)
{
   const double phase = (double)phase_idx * .25;
   const double pi = 3.1415926535897932384626433832795028841971693993751058209749445923;
   const double deg2rad = pi / 180.;
   const double jd0 = 2451550.09765;
   const double lunation = 29.530588853;
   const double k = floor((initial_jd - jd0) / lunation - phase + .5) + phase;
   const double amplit = ((phase_idx & 1) ? -.62801 : -.40720);
   const double jd = jd0 + k * lunation;

   return( jd
      + amplit * sin( 201.5643 * deg2rad + (385.81693528 * deg2rad) * k)
      + .17241 * sin(   2.5534 * deg2rad + (29.10535669 * deg2rad) * k));
}


/* get_time_from_string( ) first (*) checks for four simple types of input:
'J' or 'JD' followed by a Julian Day,  'y' followed by a decimal year,
'MJD' followed by a Modified Julian Day,  and a '+' or '-' followed
by a number and unit symbol.  That last case allows you to enter,  say,
'-3h' to change the time by three hours,  and is the reason why the
function takes the currently-set time as an input.

Next,  it looks for the common European usage d.m.y.  If it finds it,
both '.' characters are converted to '/' so that the punctuation isn't
misinterpreted as a decimal point.

Next, it looks for a time of day at the end of the string;  if found,
it's parsed and removed.  The code for interpreting the date becomes
much simpler if it doesn't have to worry about extraneous text.

Next,  we look for the usual date separators: space, '/', or '-'.
There may be two of these (separating three fields) or one (separating
two fields).  We check to see if any of the fields happens to match a
month name;  for example,  if we encounter 'oct',  that field is
unambiguously identified as a month,  October to be specific.

If any of the remaining fields is greater than 31 or less than zero,
it's probably a year.  A negative number is _definitely_ a year;
otherwise,  the largest field greater than 31 is assumed to be the
year.  (The code can handle "2008 132" to be the 132nd day of 2008;
that's why the logic has to be a bit convoluted.)

If there are only two fields remaining,  and one happens to be greater
than 12 (13 for the Jewish and Chinese calendars,  which can have thirteen
months in a year),  that larger field is assumed to be a day.

At the end of all this,  we may have identified both fields in a two-field
case,  or two out of three in a three-field case (in which case the third
is identified too,  by elimination).  That happens rather often;  you
can enter dates such as '17.08.2004' or '59 jun 3' with the day/month/year
in any order,  because it's clear what each field is.  (For modern-era
dates,  about 2/3 will be unambiguous in the three-number case,  and all
will be unambiguous if the month is spelt.)  If there is ambiguity,  the
remaining field(s) are puzzled out using the time_format input.  For
example,  with the FULL_CTIME_MONTH_DAY flag set,  '3 4' is interpreted
as the fourth day of March;  with it unset,  as the third day of April.

   (*) OK,  there are two caveats to that "first".  Before doing anything
else,  the function looks to see if there are trailing bits such as +3d,
-20m,  or +15.3s;  if so,  it removes them,  noting that three days should
be added to the final result,  or twenty minutes subtracted,  and so on.
Next, if the input string has 'nm', '1q', 'fm',  or '3q' at its end,  those
characters are removed,  the time for the _remainder_ of the string
computed,  and then the nearest lunar phase is found and returned.

   A final wrinkle:  if there's a default time zone,  "3:00" should be
assumed to be in that time zone.  But "JD 2451545.0" or an MJD should
be assumed to be UT.  (Or maybe UTC!)  To handle this,  look at the
value stored in *is_ut.  If it's zero,  the time should be assumed to
be local.  If it's non-zero,  the time should be assumed to be UT.
(All of which may be modified,  at some point,  to allow for input wherein
the zone is specified... say,  "3:00 UTC" or "4:57 EDT" or "10:00 TT".
In which case,  zero will continue to mean "time zone unspecified and
assumed local",  and non-zero values will indicate specific zones.)

   If you're working only with UTC (as many of my programs do),  you
can safetly pass a NULL for is_ut.
*/

#define AM_PM_UNSET        0
#define AM_PM_SET_TO_AM    1
#define AM_PM_SET_TO_PM    2

double DLL_FUNC get_time_from_string( double initial_jd,
         const char *time_str, const int time_format, int *is_ut)
{
   const double J2000 = 2451544.5;         /* 1.0 Jan 2000 = JD 2451544.5 */
   const int calendar = (time_format & CALENDAR_MASK);
            /* Certain solar-lunar calendars can have 13 months in a year: */
   const int max_month =
             ((calendar == CALENDAR_HEBREW || calendar == CALENDAR_CHINESE)
                    ? 13 : 12);
   int iday, month, hour, minute, n_bytes;
   int ival, i, colon_found = 0, is_bc = 0;
   int am_pm_indicator = AM_PM_UNSET;
   long year;
   double sec, dday;
   double rval = 0., offset = 0., tval;
   char buff[80];
   char symbol = 0;
   char *str = buff;

   if( is_ut)
      *is_ut = 0;
   while( *time_str == ' ')         /* skip leading spaces */
      time_str++;
   if( strlen( time_str) >= sizeof( buff) || !*time_str)
      return( initial_jd);       /* check/avoid possible buffer overflow */
             /* Ensure spaces between letters and digits.  For example,   */
             /* if time_str="11Nov 1918",  set str="11 Nov 1918".         */
             /* 'Q' is an exception to avoid '1q' becoming '1 q'.         */
   *str = tolower( *time_str++);
   i = 1;
   while( (size_t)i < sizeof( buff) - 1 && *time_str)
      {
      if( (isalpha( *time_str) && isdigit( time_str[-1]))
                         || (isdigit( *time_str) && isalpha( time_str[-1])))
         if( *time_str != 'q' && *time_str != 'Q')
            str[i++] = ' ';
      str[i++] = (char)tolower( *time_str++);
      }
   str[i] = '\0';
   remove_trailing_spaces( str);
   while( *str && (tval = collect_time_offset( str)) != 0.)
      {
      remove_trailing_spaces( str);
      offset += tval;
      }

   i = strlen( str);
   if( i > 1)
      {
      const int phase = get_phase_idx( str + i - 2);

      if( phase != PHASE_IDX_UNDEFINED)
         {
         str[i - 2] = '\0';
         rval = get_time_from_string( initial_jd, str, time_format, NULL);
         if( rval && is_ut)
            *is_ut = 1;
         return( set_from_lunar( phase, rval) + offset);
         }
      }

   is_bc = check_for_bc( str);

   if( *str == 'j')           /* decimal JD */
      {                       /* may begin with 'j' or 'jd' */
      if( str[1] == 'd')
         str++;
      rval = atof( str + 1);
      if( rval && is_ut)
         *is_ut = 1;
      }

   if( *str == 'y')      /* decimal years */
      rval = J2000 + (atof( str + 1) - 2000.) * 365.25;

   if( !memcmp( str, "mjd", 3))                 /* modified JD */
      {
      rval = atof( str + 3) + 2400000.5;
      if( is_ut)
         *is_ut = 1;
      }

   if( !memcmp( str, "now", 3))
      {
      static const double jan_1970 = 2440587.5;

      str += 3;
      while( *str == ' ')
         str++;
      initial_jd = jan_1970 + (double)time( NULL) / 86400.;
      }

   if( !*str)
      rval = initial_jd;

   if( rval)
      return( rval + offset);

               /* The common European format of separating day/month/year */
               /* with .s causes trouble,  because the code wants to see  */
               /* those as decimal numbers.  So if the input string starts */
               /* with three integers separated by dots,  we change both   */
               /* dots to '/' characters, then proceed normally:           */
   if( sscanf( str, "%d.%d.%d%n", &iday, &month, &hour, &n_bytes) == 3)
      for( i = 0; i < n_bytes; i++)
         if( str[i] == '.')
            str[i] = '/';

   sec = split_time( initial_jd, &year, &month, &iday, &hour,
                                         &minute, calendar);

               /* FITS times are always in the form YYYY-MM-DDTHH:MM:SS, */
               /* sometimes followed by .S.  This is handled separately, */
               /* in part to ensure that the month and day don't get     */
               /* swapped around:  they are _always_ in that order.      */
   i = strlen( str);
   if( i > 18 && str[10] == 't')
      if( str[4] == '-' && str[7] == '-' && str[13] == ':')
         {
         symbol = 'f';
         sscanf( str, "%ld-%d-%d", &year, &month, &iday);
         }

   if( i >= 4)
      {
      const char *search_text[4] = { " am", " a.m.", " pm", " p.m." };
      int j;

      for( j = 0; j < 4; j++)
         {
         const int slen = 3 + ((j & 1) << 1);

         if( !strcmp( str + i - slen, search_text[j]))
            {
            i -= slen;
            str[i] = '\0';
            am_pm_indicator = (j >= 2 ? AM_PM_SET_TO_PM :
                                        AM_PM_SET_TO_AM);
            j = 4;   /* break out of loop */
            }
         }
      }

            /* If the input text ends with something containing ':'s,     */
            /* assume there is a time to be extracted.  Back up along the */
            /* string,  looking for the start of the time string (which   */
            /* may be the beginning of the string,  or just after a space, */
            /* or (for FITS input) just after a 'T'... for simplicity,     */
            /* that last test just checks for any alphabetical char.)      */
   for( i = strlen( str); i && str[i - 1] != ' ' && !isalpha( str[i - 1]); i--)
      if( str[i - 1] == ':')
         colon_found = 1;

   if( strcmp( str + i, ":"))
      {
      const int saved_hour = hour;

      minute = hour = 0;
      sec = 0.;
      if( colon_found)
         {
         double dminute = 0.;

         if( str[i] != ':')
            {
            sscanf( str + i, "%d:%lf:%lf", &hour, &dminute, &sec);
            sec += dminute * 60.;
            if( am_pm_indicator == AM_PM_SET_TO_AM)
               if( hour == 12)
                  hour = 0;
            if( am_pm_indicator == AM_PM_SET_TO_PM)
               if( hour != 12)
                  hour += 12;
            }
         else      /* :MM:SS means "leave the hour unchanged" */
            {
            hour = saved_hour;
            sscanf( str + i + 1, "%lf:%lf", &dminute, &sec);
            sec += dminute * 60.;
            }
         }
      }
   if( colon_found)           /* lop the time off, leaving only the date: */
      str[i ? i - 1 : 0] = '\0';

   dday = (double)iday;
   i = 0;
   if( *str && symbol != 'f')
      {
      for( i = 1; str[i] && !strchr( "-:/ ", str[i]); i++)
         ;
      symbol = str[i];
      }
   switch( symbol)
      {
      case 'f':               /* FITS-format time: see above */
         break;
      case ':':               /* time of day */
         break;
      case '-':               /* dash-delimited such as '2009-01-20' */
      case ' ':               /* space-delimited format such as "25 dec 1980" */
      case '/':               /* common day/month/year dividing symbol */
         {
         int month_found = 0, n_fields_found = 2;
         int year_found = 0, day_found = 0;
         char tstr[80];
         double ivals[3];

         memcpy( tstr, str, i);
         tstr[i] = '\0';
         ival = month_name_to_index( tstr);
         if( ival)         /* month given first, such as 'jan 25' */
            {
            month_found = 1;
            ivals[0] = (double)ival;
            }
         else
            {
            ivals[0] = atof( tstr);
            if( strchr( tstr, '.'))   /* decimal day given */
               day_found = 1;
            }
         str += i + 1;
         for( i = 0; str[i] && str[i] != symbol && str[i] != ' '; i++)
            ;
         memcpy( tstr, str, i);
         tstr[i] = '\0';
         str += i;
         ival = month_name_to_index( tstr);
         if( ival)         /* month given second, such as '25-jan' */
            {
            month_found = 2;
            ivals[1] = (double)ival;
            }
         else
            {
            ivals[1] = atof( tstr);
            if( strchr( tstr, '.'))   /* decimal day given */
               day_found = 2;
            }

         if( *str == symbol)     /* maybe a third field was entered, but */
            {                       /* could be a time;  check for a ':' */
            str++;
            if( sscanf( str, "%s", tstr) == 1)
               if( (ival = month_name_to_index( tstr)) != 0)
                  {
                  month_found = 3;
                  n_fields_found = 3;
                  ivals[2] = (double)ival;
                  str += strlen( tstr);
                  }
            if( n_fields_found == 2)
               if( sscanf( str, "%lf%n", &ivals[2], &n_bytes) == 1)
                  if( str[n_bytes] != ':')
                     {
                     str += n_bytes;
                     if( strchr( tstr, '.'))
                        day_found = 3;
                     n_fields_found = 3;
                     }
            }
                     /* if one of the fields is negative, or if it's    */
                     /* greater than 32 and is the largest entry,  it   */
                     /* can be assumed to be the year:                  */
         for( i = 0; i < n_fields_found; i++)
            if( ivals[i] < 0.)
               {
               year_found = i + 1;   /* if we see a negative number, */
               i = n_fields_found;   /* we can stop looking further: */
               }
            else if( ivals[i] > 32.)
               if( !year_found || ivals[i] > ivals[year_found - 1])
                  year_found = i + 1;
         if( year_found || n_fields_found == 2)
            for( i = 0; i < n_fields_found; i++)
               if( ivals[i] > (double)max_month + .0001 && ivals[i] < 32.
                                      && i + 1 != year_found)
                  day_found = i + 1;

         if( n_fields_found == 2)
            {
            if( month_found)
               {
               double dval = ivals[2 - month_found];

               month = (int)ivals[month_found - 1];
               if( dval > .999 && dval < 32.)
                  dday = dval;
               else
                  year = (int)dval;
               }
            else if( year_found)         /* year/day of year format: */
               {
               year = (int)ivals[year_found - 1];
               month = 1;
               dday  = ivals[2 - year_found];
               }
            else if( day_found)     /* day/month, order is clear from */
               {                    /* the day being > 12 or having a decimal*/
               dday = ivals[day_found - 1];
               month = (int)ivals[2 - day_found];
               }
            else    /* can't tell what's day/month/year solely from input: */
               if( time_format & FULL_CTIME_MONTH_DAY)
                  {
                  month = (int)ivals[0];
                  dday = (int)ivals[1];
                  }
               else              /* day/month */
                  {
                  month = (int)ivals[1];
                  dday = (int)ivals[0];
                  }
            }
         else        /* three fields entered: */
            {
            const int year_first = (time_format & FULL_CTIME_YEAR_FIRST);

            if( !year_found)
               {
               if( !month_found)
                  {
                  if( !day_found || day_found == 2)
                     {
                     year_found = (year_first ? 1 : 3);
                     if( !day_found)   /* must rely solely on time format */
                        {              /* settings;  no fields autofound  */
                        day_found = (year_first ? 2 : 1);
                        if( time_format & FULL_CTIME_MONTH_DAY)
                           day_found++;      /* ymd or mdy case */
                        }
                     }
                  else     /* if day is 1st or last, year is last or 1st */
                     year_found = 4 - day_found;
                  }
               else if( !day_found)    /* only the month was found: */
                  {
                  if( month_found == 2)
                     year_found = (year_first ? 1 : 3);
                  else   /* if month is 1st or last, year must be last/1st */
                     year_found = 4 - month_found;
                  }
               }
            else        /* year_found... */
               if( !month_found && !day_found)  /* ...but nothing else */
                  {
                  if( time_format & FULL_CTIME_MONTH_DAY)
                     month_found = (year_found == 1 ? 2 : 1);
                  else
                     day_found = (year_found == 1 ? 2 : 1);
                  }
                  /* We now have the year nailed down.  If either the day */
                  /* or month is still not nailed down, we can find it    */
                  /* easily, since the 'found' values must sum up to 6:   */

            if( !day_found)
               day_found = 6 - year_found - month_found;
            else if( !month_found)
               month_found = 6 - year_found - day_found;
            year = (int)floor( ivals[year_found - 1] + .5);
            dday = ivals[day_found - 1];
            month = (int)( ivals[month_found - 1] + .5);
            }

         if( year > 0 && year < 100 && !is_bc)
            if( time_format & FULL_CTIME_TWO_DIGIT_YEAR)
               year += (year < 40 ? 2000 : 1900);
         }
         break;
      case '\0':       /* no dividing symbols found */
         if( *str)
            {
            ival = month_name_to_index( str);
            if( ival)
               month = ival;
            else if( (ival = day_of_week_name_to_index( str)) >= 0)
               {
               ival -= ((long)(initial_jd + 1.5)) % 7;
               if( ival < -3)
                  ival += 7;
               else if( ival > 3)
                  ival -= 7;
               dday += ival;
               }
            else
               {
               n_bytes = 0;
               if( sscanf( str, "%d%n", &ival, &n_bytes) == 1)
                  {
                  double tval = 0.;

                  str += n_bytes;
                  if( *str == '.')     /* also a fractional part to this: */
                     {
                     int n_bytes_in_fraction;

                     sscanf( str, "%lf%n", &tval, &n_bytes_in_fraction);
                     str += n_bytes_in_fraction;
                     }
                  tval += (double)ival;
                  switch( n_bytes)
                     {
                     case 1:
                     case 2:                 /* reset day */
                        dday = tval;
                        break;
                     case 3:                 /* reset day of year */
                        dday = tval;
                        month = 1;
                        break;
                     case 4:                 /* reset year, which may be */
                     case 5:                 /* four or five digits long */
                        if( (double)ival == tval)
                           {                    /* set 1 Jan of the year */
                           dday = 1.;
                           month = 1;
                           year = ival;
                           }
                        else        /* true decimal year */
                           return( J2000 + (tval - 2000.) * 365.25 + offset);
                        break;
                     case 7:               /* JD */
                        if( is_ut)
                           *is_ut = 1;
                        return( tval + offset);
                     case 6:     /* YYMMDD(.DD)   */
                     case 8:     /* YYYYMMDD(.DD) */
                        year = ival / 10000L;
                        if( n_bytes == 6)
                           year += (year < 40 ? 2000 : 1900);
                        month = (ival / 100) % 100L;
                        dday = (double)( ival % 100L);
                        dday += tval - (double)ival;
                        break;
                     }
                  }
               else        /* couldn't make sense of input text */
                  return( 0.);
               }
            }
         break;
      default:
         return( 0.);
//       break;
      }

   if( is_bc)
      year = 1 - year;
   rval = (double)dmy_to_day( 0, month, year, calendar) + dday - .5 +
                 (double)hour / 24. + (double)minute / 1440. + sec / 86400.;
   return( rval + offset);
}
