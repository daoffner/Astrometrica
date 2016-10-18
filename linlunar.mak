all: astcheck astephem calendar easter get_test \
 htc20b jd jsattest lun_test oblitest persian  \
 phases testprec tables \
	ssattest uranus1

CC=g++

CFLAGS=-Wall -O3

.cpp.o:
	$(CC) $(CFLAGS) -c $<

OBJS= alt_az.o astfuncs.o big_vsop.o classel.o cospar.o date.o delta_t.o \
	de_plan.o dist_pa.o eart2000.o elp82dat.o getplane.o get_time.o \
	jsats.o lunar2.o miscell.o nutation.o obliquit.o pluto.o precess.o \
	showelem.o ssats.o triton.o vsopson.o

lunar.a: $(OBJS)
	rm -f lunar.a
	ar rv lunar.a $(OBJS)

astcheck:  astcheck.o mpcorb.o lunar.a
	$(CC) $(CFLAGS) -o astcheck astcheck.o mpcorb.o lunar.a

astephem:  astephem.o mpcorb.o lunar.a
	$(CC) $(CFLAGS) -o astephem astephem.o mpcorb.o lunar.a

calendar: calendar.o lunar.a
	$(CC) $(CFLAGS) -o calendar   calendar.o   lunar.a

easter: easter.cpp lunar.a
	$(CC) $(CFLAGS) -o easter -DTEST_CODE easter.cpp lunar.a

get_test: get_test.o lunar.a
	$(CC) $(CFLAGS) -o get_test get_test.o lunar.a

htc20b: htc20b.cpp lunar.a
	$(CC) $(CFLAGS) -o htc20b -DTEST_MAIN htc20b.cpp lunar.a

jd: jd.o lunar.a
	$(CC) $(CFLAGS) -o jd jd.o lunar.a

jsattest: jsattest.o lunar.a
	$(CC) $(CFLAGS) -o jsattest jsattest.o lunar.a

lun_test:                lun_test.o lun_tran.o riseset3.o lunar.a
	$(CC)	$(CFLAGS) -o lun_test lun_test.o lun_tran.o riseset3.o lunar.a

oblitest: oblitest.o obliqui2.o spline.o lunar.a
	$(CC) $(CFLAGS) -o oblitest oblitest.o obliqui2.o spline.o lunar.a

persian: persian.o solseqn.o lunar.a
	$(CC) $(CFLAGS) -o persian persian.o solseqn.o lunar.a

phases: phases.o lunar.a
	$(CC) $(CFLAGS) -o phases   phases.o   lunar.a

testprec: testprec.o obliqui2.o spline.o lunar.a
	$(CC) $(CFLAGS) -o testprec testprec.o obliqui2.o spline.o lunar.a

relativi: relativi.cpp lunar.a
	$(CC) $(CFLAGS) -o relativi -DTEST_CODE relativi.cpp lunar.a

ssattest: ssattest.o lunar.a
	$(CC) $(CFLAGS) -o ssattest ssattest.o lunar.a

tables: tables.o riseset3.o lunar.a
	$(CC) $(CFLAGS) -o tables tables.o riseset3.o lunar.a

uranus1: uranus1.o gust86.o
	$(CC) $(CFLAGS) -o uranus1 uranus1.o gust86.o

