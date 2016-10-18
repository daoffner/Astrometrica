# Make file for console & non-interactive Find_Orb,  using PDCurses
# Note that this is very similar to 'linmake',  the make file for Linux,
# except that LINUX isn't defined and 'pdcurses.a' is linked to instead
# of the plain old 'curses' library.

all: find_orb.exe fo.exe

CC=g++

CFLAGS=-c -O3 -Wall

OBJS=b32_eph.o collide.o conv_ele.o eigen.o \
	elem2tle.o elem_out.o ephem0.o gauss.o get_pert.o \
	jpleph.o lsquare.o moid4.o monte0.o mpc_obs.o mt64.o \
	orb_func.o orb_fun2.o pl_cache.o roots.o runge.o \
	sm_vsop.o sr.o tle_out.o weight.o

LIBS=\lunar\posttest\lunar.a \pdcurses\win32a\pdcurses.a

find_orb.exe:      findorb.o $(OBJS) clipfunc.o $(LIBS) findorb.res
	$(CC) -o find_orb findorb.o $(OBJS) clipfunc.o $(LIBS) -l gdi32 -l user32 findorb.res

fo.exe:      fo.o $(OBJS) \lunar\posttest\lunar.a findorb.res
	$(CC) -o fo fo.o $(OBJS) \lunar\posttest\lunar.a

.cpp.o:
	$(CC) $(CFLAGS) $<

jpleph.o:   \jpl\pl\jpleph.cpp
	$(CC) $(CFLAGS) \jpl\pl\jpleph.cpp

lsquare.o:   \image\lsquare.cpp
	$(CC) $(CFLAGS) \image\lsquare.cpp

tle_out.o: \sattest\neoklis\tle_out.cpp
	$(CC) $(CFLAGS) \sattest\neoklis\tle_out.cpp

mt64.o: \prng\mt64.cpp
	$(CC) $(CFLAGS) \prng\mt64.cpp

findorb.res: res\find_orb.ico findorb.rc
	windres findorb.rc -O coff -o findorb.res
