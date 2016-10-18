/* easter.cpp: functions for computing date of Easter, plus test code

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

/* See Meeus,  _Astronomical Algorithms_,  p 67.  He in turn states that
   "the following method has been given by Spencer Jones in his book
   _General Astronomy_ (p 73-74 of the 1922 edition).  It has been
   published again in the _Journal of the British Astronomical Association_,
   vol 88, page 91 (December 1977),  where it is said that it was devised
   in 1876 and appeared in Butcher's _Astronomical Calendar._

      Unlike the formula given by Gauss,  this method has no exception and
   is valid for all years in the Gregorian calendar,  hence from 1583 on."

   I modified the method to include negative years by taking advantage
of the fact that,  in the Gregorian calendar,  Easter recurs on a 5.7
million year cycle. Add 5.7 million to the year,  and all numbers are
positive (as long as year > -5700000,  of course!  I suppose one could
add,  say,  570 million instead,  extending the range to cover roughly
the Paleozoic era to a time 1.5 billion years hence...)

   The test program can be run with zero,  one,  or two command line
arguments,  to get three very different types of output.

   If run without command-line arguments,  the program computes the day
of Easter over an entire 5.7 million year cycle.  That lets it show the
frequency on which Easter occurs on a given day;  the result looks like

Mar 22: 0.483   Mar 29: 3.383  Apr  5: 3.383   Apr 12: 3.383  Apr 19: 3.867
Mar 23: 0.950   Mar 30: 3.325  Apr  6: 3.325   Apr 13: 3.325  Apr 20: 3.325
Mar 24: 1.425   Mar 31: 3.325  Apr  7: 3.325   Apr 14: 3.325  Apr 21: 2.850
Mar 25: 1.933   Apr  1: 3.383  Apr  8: 3.383   Apr 15: 3.383  Apr 22: 2.417
Mar 26: 2.333   Apr  2: 3.267  Apr  9: 3.267   Apr 16: 3.267  Apr 23: 1.867
Mar 27: 2.900   Apr  3: 3.383  Apr 10: 3.383   Apr 17: 3.383  Apr 24: 1.450
Mar 28: 3.267   Apr  4: 3.267  Apr 11: 3.267   Apr 18: 3.463  Apr 25: 0.737

   Run with a single command-line argument of a year,  you get the dates of
Easter for a 75-year time span.  For example,  'easter 2010' gives

2010 Apr  4    2025 Apr 20    2040 Apr  1    2055 Apr 18    2070 Mar 30
2011 Apr 24    2026 Apr  5    2041 Apr 21    2056 Apr  2    2071 Apr 19
2012 Apr  8    2027 Mar 28    2042 Apr  6    2057 Apr 22    2072 Apr 10
2013 Mar 31    2028 Apr 16    2043 Mar 29    2058 Apr 14    2073 Mar 26
2014 Apr 20    2029 Apr  1    2044 Apr 17    2059 Mar 30    2074 Apr 15
2015 Apr  5    2030 Apr 21    2045 Apr  9    2060 Apr 18    2075 Apr  7
2016 Mar 27    2031 Apr 13    2046 Mar 25    2061 Apr 10    2076 Apr 19
2017 Apr 16    2032 Mar 28    2047 Apr 14    2062 Mar 26    2077 Apr 11
2018 Apr  1    2033 Apr 17    2048 Apr  5    2063 Apr 15    2078 Apr  3
2019 Apr 21    2034 Apr  9    2049 Apr 18    2064 Apr  6    2079 Apr 23
2020 Apr 12    2035 Mar 25    2050 Apr 10    2065 Mar 29    2080 Apr  7
2021 Apr  4    2036 Apr 13    2051 Apr  2    2066 Apr 11    2081 Mar 30
2022 Apr 17    2037 Apr  5    2052 Apr 21    2067 Apr  3    2082 Apr 19
2023 Apr  9    2038 Apr 25    2053 Apr  6    2068 Apr 22    2083 Apr  4
2024 Mar 31    2039 Apr 10    2054 Mar 29    2069 Apr 14    2084 Mar 26

   Run with two command line arguments of a month (3=March or 4=April) and
a day,  one gets a list of years from 0 AD to 10000 AD when Easter would
fall on that day.  For example,  'easter 3 22' would give the following
list of years when Easter fell on 22 March,  the earliest day possible:

   15  387  482  539  607  854 1074 1131 1226 1503 1598 1693 1761 1818 2285
 2353 2437 2505 2972 3029 3401 3496 3564 3648 3716 4308 5299 5671 6043 6195
 6263 6415 6635 6703 6798 6882 6950 7322 7474 7542 7637 7789 7914 8161 8533
 8685 8753 8848 8905 9125 9220 9372 9440 9812 9964
55 found over 10000 years
*/

void easter_date( const long year, int *month, int *day)
{
   const long year2 = year + 5700000L;
   const long a = year2 % 19L, b = year2 / 100L, c = year2 % 100L;
   const long d = b / 4L, e = b % 4L, f = (b + 8L) / 25L;
   const long g = (b - f + 1L) / 3L, h = (19L * a + b - d - g + 15L) % 30L;
   const long i = c / 4L, k = c % 4L, l = (32L + e + e + i + i - h - k) % 7L;
   const long m = (a + 11L * h + 22L * l) / 451L, tval = h + l - 7L * m + 114L;

   *month = (int)( tval / 31L);
   *day = (int)( tval % 31L) + 1;
}

#ifdef TEST_CODE

#include <stdio.h>
#include <stdlib.h>

int main( int argc, char **argv)
{
   int month, day;
   long year;

   if( argc == 2)
      {
      int i;
      const int n_across = 5, n_down = 15;

      for( i = 0; i < n_across * n_down; i++)
         {
         year = atol( argv[1]) + (i % n_across) * n_down + i / n_across;
         easter_date( year, &month, &day);
         printf( "%ld %s %2d%s", year, (month == 3 ? "Mar" : "Apr"), day,
                  ((i + 1) % n_across) ? "    " : "\n");
         }
      }
   else if( argc == 3)
      {
      int n_found = 0;

      for( year = 0; year < 10000; year++)
         {
         easter_date( year, &month, &day);
         if( month == atoi( argv[1]) && day == atoi( argv[2]))
            {
            printf( "%5ld", year);
            if( n_found++ % 15 == 14)
               printf( "\n");
            }
         }
      printf( "\n%d found over 10000 years\n", n_found);
      }
   else
      {
      long march[32], april[32];

      for( day = 0; day < 32; day++)
         april[day] = march[day] = 0;
      for( year = 0; year < 5700000; year++)
         {
         easter_date( year, &month, &day);
         if( month == 3)
            march[day]++;
         else
            april[day]++;
         }
      for( day = 0; day < 32; day++)
         if( march[day])
            printf( "Mar %2d: %6.3lf\n", day, (double)march[day] / 57000.);
      for( day = 0; day < 32; day++)
         if( april[day])
            printf( "Apr %2d: %6.3lf\n", day, (double)april[day] / 57000.);

      printf( "Run 'easter' with a year on the command line to get the date\n");
      printf( "of Easter for that year.  For example,  'easter 2008' will\n");
      printf( "get the output 'March 23'.  Alternatively,  give a month and\n");
      printf( "day on the command line to get the years between 0 and 10000\n");
      printf( "when Easter will fall on that day.  For example,  'easter 3 23'\n");
      printf( "will produce a list of all years when Easter is on 23 March.\n");
      printf( "\nNote that Easter cannot occur before 22 March or after 25 April.\n");
      }
   return( 0);
}
#endif
