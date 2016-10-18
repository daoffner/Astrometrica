all: find_orb.exe fo.exe

OBJS=b32_eph.obj collide.obj conv_ele.obj eigen.obj \
  elem2tle.obj elem_out.obj ephem0.obj gauss.obj get_pert.obj \
  jpleph.obj lsquare.obj moid4.obj monte0.obj mpc_obs.obj mt64.obj \
  orb_func.obj orb_fun2.obj pl_cache.obj roots.obj runge.obj \
  sm_vsop.obj sr.obj tle_out.obj weight.obj

find_orb.exe:               findorb.obj $(OBJS) clipfunc.obj
     link /out:find_orb.exe findorb.obj $(OBJS) clipfunc.obj lunar.lib pdcurses.lib user32.lib

fo.exe:               fo.obj $(OBJS)
     link /out:fo.exe fo.obj $(OBJS) lunar.lib

CFLAGS=/c /Ot /W3 /nologo /MT /DCONSOLE

.cpp.obj:
   cl $(CFLAGS) $<

