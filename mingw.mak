all: astcheck.exe astephem.exe calendar.exe easter.exe get_test.exe \
 htc20b.exe jd.exe jsattest.exe lun_test.exe oblitest.exe persian.exe  \
 phases.exe ps_1996.exe testprec.exe relativi.exe tables.exe \
	ssattest.exe uranus1.exe

CC=g++

CFLAGS=-Wall -O3

.cpp.o:
	$(CC) $(CFLAGS) -c $<

OBJS= alt_az.o astfuncs.o big_vsop.o classel.o cospar.o date.o delta_t.o \
	de_plan.o dist_pa.o eart2000.o elp82dat.o getplane.o get_time.o \
	jsats.o lunar2.o miscell.o nutation.o obliquit.o pluto.o precess.o \
	showelem.o ssats.o triton.o vsopson.o

lunar.a: $(OBJS)
	del lunar.a
	ar rv lunar.a $(OBJS)

astcheck.exe:  astcheck.o mpcorb.o lunar.a
	$(CC) $(CFLAGS) -o astcheck astcheck.o mpcorb.o lunar.a

astephem.exe:  astephem.o mpcorb.o lunar.a
	$(CC) $(CFLAGS) -o astephem astephem.o mpcorb.o lunar.a

calendar.exe: calendar.o lunar.a
	$(CC) $(CFLAGS) -o calendar   calendar.o   lunar.a

easter.exe: easter.cpp lunar.a
	$(CC) $(CFLAGS) -o easter -DTEST_CODE easter.cpp lunar.a

get_test.exe: get_test.o lunar.a
	$(CC) $(CFLAGS) -o get_test get_test.o lunar.a

htc20b.exe: htc20b.cpp lunar.a
	$(CC) $(CFLAGS) -o htc20b -DTEST_MAIN htc20b.cpp lunar.a

jd.exe: jd.o lunar.a
	$(CC) $(CFLAGS) -o jd jd.o lunar.a

jsattest.exe: jsattest.o lunar.a
	$(CC) $(CFLAGS) -o jsattest jsattest.o lunar.a

lun_test.exe:                lun_test.o lun_tran.o riseset3.o lunar.a
	$(CC)	$(CFLAGS) -o lun_test lun_test.o lun_tran.o riseset3.o lunar.a

oblitest.exe: oblitest.o obliqui2.o spline.o lunar.a
	$(CC) $(CFLAGS) -o oblitest oblitest.o obliqui2.o spline.o lunar.a

persian.exe: persian.o solseqn.o lunar.a
	$(CC) $(CFLAGS) -o persian persian.o solseqn.o lunar.a

phases.exe: phases.o lunar.a
	$(CC) $(CFLAGS) -o phases   phases.o   lunar.a

ps_1996.exe: ps_1996.o lunar.a
	$(CC) $(CFLAGS) -o ps_1996   ps_1996.o   lunar.a

testprec.exe: testprec.o obliqui2.o spline.o lunar.a
	$(CC) $(CFLAGS) -o testprec testprec.o obliqui2.o spline.o lunar.a

relativi.exe: relativi.cpp lunar.a
	$(CC) $(CFLAGS) -o relativi -DTEST_CODE relativi.cpp lunar.a

ssattest.exe: ssattest.o lunar.a
	$(CC) $(CFLAGS) -o ssattest ssattest.o lunar.a

tables.exe: tables.o riseset3.o lunar.a
	$(CC) $(CFLAGS) -o tables tables.o riseset3.o lunar.a

uranus1.exe: uranus1.o gust86.o
	$(CC) $(CFLAGS) -o uranus1 uranus1.o gust86.o

