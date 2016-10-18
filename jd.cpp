/* jd.cpp: example of use of date/time functions & calendar conversions

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
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "watdefs.h"
#include "date.h"
#include "afuncs.h"

#ifdef OBSOLETE_FOR_REFERENCE_ONLY
/* At one point,  I had the months shown by name,  not number.  This is
no longer true,  but I saw no need to chop out the month names for
those who might be interested: */

extern const char *islamic_month_names[12] = {
      "Muharram", "Safar", "Rabi'a I",
      "Rabi'a II", "Jumada I", "Jumada II",
      "Rajab", "Sha'ban", "Ramadan",
      "Shawwal", "Dhu al-Q'adah", "Dhu al-Hijjah" };

extern const char *hebrew_month_names[13] = {
      "Tishri", "Heshvan", "Kislev", "Tevet", "Shevat", "Adar",
      "Adar II", "Nisan", "Iyar", "Sivan", "Tammuz", "Av", "Elul" };

extern const char *french_month_names[12] = {
  "Vend�miaire",           "Germinal",
  "Brumaire",              "Flor�al",
  "Frimaire",              "Prairial",
  "Niv�se",                "Messidor",
  "Pluvi�se",              "Thermidor",
  "Vent�se",               "Fructidor" };

extern const char *french_extra_day_names[6] = {
        "Jour de la vertu (Virtue Day)",
        "Jour du genie (Genius Day)",
        "Jour du travail (Labour Day)",
        "Jour de l'opinion (Reason Day)",
        "Jour des recompenses (Rewards Day)",
        "Jour de la revolution (Revolution Day)" };
#endif


/* The Chinese calendar involves added complications for two reasons.
First,  instead of being computed algorithmically (as all the other
calendars are),  it's computed using a pre-compiled data table.  So you
have to load the file into memory,  as shown in the following code;  the
data is stored using the set_chinese_calendar_data( ) function.

   Second,  as will be seen in the main( ) code,  you have to watch out
for the intercalary months.  */

static int load_chinese_calendar_data( const char *filename)
{
   FILE *ifile = fopen( filename, "rb");
   int rval = -1;

   if( ifile)
      {
      size_t filesize;
      static char *calendar_data = NULL;

      fseek( ifile, 0L, SEEK_END);
      filesize = (size_t)ftell( ifile);
      calendar_data = (char *)malloc( filesize);
      if( calendar_data)
         {
         fseek( ifile, 0L, SEEK_SET);
         if( !fread( calendar_data, filesize, 1, ifile))
            rval = -2;
         else
            {
            set_chinese_calendar_data( calendar_data);
            rval = 0;
            }
         }
      fclose( ifile);
      }
   return( rval);
}

static void error_exit( void)
{
   printf( "jd takes either a Julian Day or a year/month/day as command\n");
   printf( "line arguments.  It then expresses that date in several calendar\n");
   printf( "systems (Julian, Gregorian, Hebrew, Islamic, Persian, Chinese,\n");
   printf( "French Revolutionary).  Example usages:\n\n");
   printf( "jd 2450000.5   (to get calendar data for JD 2450000.5)\n");
   printf( "jd 2007 3 30   (to get calendar data for 2007 March 30)\n");
   printf( "jd 2007-mar-30   (same as above; 30-mar-2007 would also work)\n");
   exit( -1);
}

int main( int argc, char **argv)
{
   int calendar, err_code, i;
   static const double jan_1970 = 2440587.5;
               /* Set default initial time to "right now": */
   double jd = jan_1970 + time( NULL) / 86400.;
   char buff[90];

   if( argc < 2)
      error_exit( );
   strcpy( buff, argv[1]);
   for( i = 2; i < argc; i++)
      {
      strcat( buff, " ");
      strcat( buff, argv[i]);
      }

   jd = get_time_from_string( jd, buff,
        CALENDAR_JULIAN_GREGORIAN | FULL_CTIME_YMD | FULL_CTIME_TWO_DIGIT_YEAR, NULL);

   err_code = load_chinese_calendar_data( "chinese.dat");
   if( err_code)
      printf( "WARNING:  Chinese calendar data not loaded: %d\n", err_code);

   if( !jd)    /* no date found in command line;  show an error message: */
      error_exit( );
   else
      {
      full_ctime( buff, jd, CALENDAR_JULIAN_GREGORIAN | FULL_CTIME_YMD
                                    | FULL_CTIME_DAY_OF_WEEK_FIRST);
      printf( "%s = JD %.5lf\n", buff, jd);
      for( calendar = 0; calendar < 9; calendar++)
         {
         int day, month;
         long year;
         const long ljd = (long)floor( jd + .5);
         static const char *calendar_names[9] = {
                "Gregorian", "Julian", "Hebrew", "Islamic", "Revolutionary",
                "Persian (Jalali)", "Greg/Jul", "Chinese",  "Modern Persian" };
         char is_intercalary = ' ';

         day_to_dmy( ljd, &day, &month, &year, calendar);
         if( calendar == CALENDAR_CHINESE)
            {
            const int chinese_intercalary_month = get_chinese_intercalary_month( );

            if( chinese_intercalary_month)
               {
               if( month == chinese_intercalary_month)
                  is_intercalary = 'i';
               if( month >= chinese_intercalary_month)
                  month--;
               }
            }

         printf( "%-20s %3d%3d%c%4ld\n", calendar_names[calendar],
                        day, month, is_intercalary, year);
         }
      }
   printf( "Delta-T = TD - UT = %.4lf; TD - UTC = %.4lf; UTC - UT = %.4lf\n",
                            td_minus_ut( (double)jd),
                            td_minus_utc( (double)jd),
                            td_minus_ut( (double)jd) - td_minus_utc( (double)jd));
   return( 0);
}
