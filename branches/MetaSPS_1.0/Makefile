ROOT := $(shell pwd)
XERCES_ROOT := ../src/xerces-c-3.1.1
#XERCES_ROOT := ../src/xerces-c-3.1.1
GNUPLOT_ROOT := ../src/gnuplot-4.2.5
BOOST_ROOT := ../src/boost_1_44_0

SPS_VERSION := $(shell date +%y%m%d)

CC := g++
CPPFLAGS := -DSPS_VERSION=$(SPS_VERSION) -I$(XERCES_ROOT)/src -Iinclude -lpthread

# 'make type=debug'
ifeq ($(type),debug)
  CPPFLAGS += -O0 -ggdb3
else
  CPPFLAGS += -O3 -ffunction-sections -fdata-sections
endif

# 'make type=32'
ifeq ($(type),32)
  LDFLAGS += -march=i686 -m32
  CPPFLAGS += -march=i686 -m32
  PREFIX := ../install/linux-g++

  # gcc 4
  ifeq ($(shell gcc -dumpversion | cut -d. -f 1), 4)
    LDFLAGS += -Wl,--gc-sections
  endif
endif

# 'make'
ifndef type
  # gcc 4
  ifeq ($(shell gcc -dumpversion | cut -d. -f 1), 4)
    LDFLAGS += -Wl,--gc-sections
    CPPFLAGS += -mtune=generic

  # gcc 3
  else
    CPPFLAGS += -march=$(shell uname -m)
  endif
endif

# Linux
ifneq (,$(findstring linux,$(shell gcc -dumpmachine)))
  ifneq (,$(findstring x86_64,$(shell gcc -dumpmachine)))
    PREFIX ?= ../install/linux-g++-64
  else
    PREFIX ?= ../install/linux-g++
  endif

  LDFLAGS += -L$(PREFIX)/lib

# Windows
else
  ifneq (,$(findstring x86_64,$(shell gcc -dumpmachine)))
    PREFIX ?= ../install/win32-g++-64
  else
    PREFIX ?= ../install/win32-g++
  endif

  LDFLAGS += -L$(PREFIX)/lib -L$(PREFIX)/bin
endif

CPPLIB_SRC = \
	abruijn.cpp \
	alignment_modmut.cpp \
	alignment_scoring.cpp \
	aminoacid.cpp \
	base64.cpp \
	batch.cpp \
	clusters.cpp \
	db_fasta.cpp \
	dekel_align.cpp \
	denovo.cpp \
	filters.cpp \
	graph.cpp \
	hash.cpp \
	inputParams.cpp \
	label.cpp \
	ms1.cpp \
	msn.cpp \
	PepnovoTags.cpp \
	setmerger.cpp \
	spectral_pairs.cpp \
	spectrum.cpp \
	spectrum_scoring.cpp \
	tags.cpp \
	ion.cpp \
	utils.cpp

CPPLIB_OBJ=$(CPPLIB_SRC:.cpp=.o)

EXECUTABLES = \
	filterpairs \
	starden \
	filterstarpairs \
	tagsearch \
	pathproj \
	protid \
	masab \
	matchma \
	mergefilter \
	homglue

SCRIPTS = \
	run_csps.sh \
	run_specnets.sh \
	grid.sh \
	grid_csps_start.sh \
	grid_sps_resume.sh \
	grid_sps_start.sh

SUBDIRS = \
	EulerOverlap

.PHONY: subdirs $(SUBDIRS) clean all

all: subdirs $(EXECUTABLES)

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@ LDFLAGS="$(LDFLAGS)" CPPFLAGS="$(CPPFLAGS)"

install: all
	mkdir -p $(PREFIX)/bin
	mkdir -p $(PREFIX)/share/gnuplot/4.2/PostScript
	cp -f $(EXECUTABLES) PepNovo_bin MsCluster_bin $(SCRIPTS) $(PREFIX)/bin
	cp -f $(GNUPLOT_ROOT)/term/PostScript/prologue.ps $(PREFIX)/share/gnuplot/4.2/PostScript
ifneq (,$(shell which mcc 2>/dev/null))
	cp -f csps specnets $(PREFIX)/bin
endif

clean:
	rm -f *.o $(EXECUTABLES); cd EulerOverlap; make clean; cd ..

distclean: clean
	$(MAKE) -C $(XERCES_ROOT) distclean
	$(MAKE) -C $(GNUPLOT_ROOT) distclean

.o.o: $(CPPLIB_SRC)
	$(CC) $(CPPFLAGS) $<

test_cs: test_cs.cpp $(CPPLIB_OBJ)
	$(CC) $(CPPFLAGS) $^ -o $@

test: test.cpp $(CPPLIB_OBJ)
	$(CC) $(CPPFLAGS) $^ -o $@

mergeConvert: mergeConvert.o $(CPPLIB_OBJ)
	$(CC) $(LDFLAGS) $^ -o $@

spectrum_statistics: spectrum_statistics.cpp inspect_parse.o spectrum_annot_statistics.o spectrum_annot_statistics_parse.o $(CPPLIB_OBJ)
	$(CC) $(CPPFLAGS) $^ -o $@
	
pklbin_to_mgf: pklbin_to_mgf.cpp $(CPPLIB_OBJ)
	$(CC) $(CPPFLAGS) $^ -o $@
	
spectrum_test_no_dll: spectrum_test.cpp $(CPPLIB_OBJ)
	$(CC) -I$(BOOST_ROOT) $(CPPFLAGS) $^ -o $@ -DNO_DLL

spectrum_test: spectrum_test.cpp $(CPPLIB_OBJ)
	$(CC) $(LDFLAGS) -I$(BOOST_ROOT) $(CPPFLAGS) $^ -o $@ -lboost_unit_test_framework

spectrum_annot_statistics_test_no_dll: spectrum_annot_statistics_test.cpp spectrum_annot_statistics.o $(CPPLIB_OBJ)
	$(CC) -I$(BOOST_ROOT) $(CPPFLAGS) $^ -o $@ -DNO_DLL

spectrum_annot_statistics_test: spectrum_annot_statistics_test.cpp spectrum_annot_statistics.o $(CPPLIB_OBJ)
	$(CC) $(LDFLAGS) -I$(BOOST_ROOT) $(CPPFLAGS) $^ -o $@ -lboost_unit_test_framework
	

# Linux
ifneq (,$(findstring linux,$(shell gcc -dumpmachine)))
convert: convert.o mzxml.o $(CPPLIB_OBJ)
	$(CC) $(LDFLAGS) $^ -o $@ -lxerces-c-3.0

# Windows
else
convert: convert.o mzxml.o $(CPPLIB_OBJ)
	$(CC) $(LDFLAGS) $^ -o $@ -lxerces-c-3-0
endif

filterpairs: filterpairs.o $(CPPLIB_OBJ)
	$(CC) $(LDFLAGS) $^ -o $@

#specalign: specalign.cpp $(CPPLIB_OBJ)
#	$(CC) $(CPPFLAGS) $^ -o $@

# starden is here temporarily as an interim replacement for specalign
starden: starden.o $(CPPLIB_OBJ)
	$(CC) $(LDFLAGS) $^ -o $@

filterstarpairs: filterstarpairs.o $(CPPLIB_OBJ)
	$(CC) $(LDFLAGS) $^ -o $@

tagsearch: tagsearch.o $(CPPLIB_OBJ)
	$(CC) $(LDFLAGS) $^ -o $@

pathproj: pathproj.o $(CPPLIB_OBJ)
	$(CC) $(LDFLAGS) $^ -o $@

findconnected: findconnected.o $(CPPLIB_OBJ)
	$(CC) $(LDFLAGS) $^ -o $@

msmod: msmod.o FloatMods.o FloatMods.h
	$(CC) $(LDFLAGS) $^ -o $@

# SPS modules
masab: masab.o $(CPPLIB_OBJ)
	$(CC) $(LDFLAGS) $^ -o $@

# CSPS modules
protid: protid.o $(CPPLIB_OBJ)
	$(CC) $(LDFLAGS) $^ -o $@

matchma: matchma.o $(CPPLIB_OBJ)
	$(CC) $(LDFLAGS) $^ -o $@

matchme: matchme.o $(CPPLIB_OBJ)
	$(CC) $(LDFLAGS) $^ -o $@

mergefilter: mergefilter.o $(CPPLIB_OBJ)
	$(CC) $(LDFLAGS) $^ -o $@

homglue: homglue.o $(CPPLIB_OBJ)
	$(CC) $(LDFLAGS) $^ -o $@

# NRP-sequencing modules
cycseq: cycseq.o $(CPPLIB_OBJ)
	$(CC) $(LDFLAGS) $^ -o $@

# Matlab
csps: csps.m
ifneq (,$(shell which mcc 2>/dev/null))
	mcc -R -nojvm -R -nojit -m $< -o $@
endif

specnets: specnets.m
ifneq (,$(shell which mcc 2>/dev/null))
	mcc -R -nojvm -R -nojit -m $< -o $@
endif

# specplot
%plot_l.o: %plot.cpp
	$(CC) $(CPPFLAGS) $< -o $@ -c -DSPECPLOT

$(GNUPLOT_ROOT)/Makefile:
	cd $(GNUPLOT_ROOT); ./configure --without-readline --without-lisp-files --without-tutorial --without-row-help --without-x --disable-wxwidgets; cd $(ROOT)

$(XERCES_ROOT)/Makefile:
	cd $(XERCES_ROOT); ./configure --disable-pretty-make --disable-network --disable-threads --without-curl --without-icu; cd $(ROOT)

$(PREFIX)/lib/libgnuplot.a: $(GNUPLOT_ROOT)/Makefile
	$(MAKE) -C $(GNUPLOT_ROOT)/src gnuplot CXXLINK="ar -rcs libgnuplot.a" gnuplot_LDADD="" LIBS="" CFLAGS="$(CPPFLAGS) -fPIC -Dmain=gnuplot_main"
	mv $(GNUPLOT_ROOT)/src/libgnuplot.a $(PREFIX)/lib

$(PREFIX)/lib/libxerces-c.a: $(XERCES_ROOT)/Makefile
	$(MAKE) -C $(XERCES_ROOT)/src LDFLAGS="$(LDFLAGS)" CPPFLAGS="$(CPPFLAGS)"
	mv $(XERCES_ROOT)/src/.libs/libxerces-c.a $(PREFIX)/lib

# Linux
ifneq (,$(findstring linux,$(shell gcc -dumpmachine)))
$(PREFIX)/lib/libgnuplot.so: $(GNUPLOT_ROOT)/Makefile
	$(MAKE) -C $(GNUPLOT_ROOT)/src gnuplot LDFLAGS="$(LDFLAGS) -shared -Wl,-o,libgnuplot.so" CFLAGS="$(CPPFLAGS) -fPIC -Dmain=gnuplot_main"
	mv $(GNUPLOT_ROOT)/src/libgnuplot.so $(PREFIX)/lib

$(PREFIX)/lib/libxerces-c-3.0.so: $(XERCES_ROOT)/Makefile
	$(MAKE) -C $(XERCES_ROOT)/src LDFLAGS="$(LDFLAGS)" CPPFLAGS="$(CPPFLAGS)"
	mv $(XERCES_ROOT)/src/.libs/libxerces-c-3.0.so $(PREFIX)/lib

ifneq (,$(shell locate libbfd.a))
%plot.o: %plot.cpp
	$(CC) $(CPPFLAGS) -ggdb $< -o $@ -c

spsplot: spsplot.o specplot.o contplot.o spsstats.o mzxml.o dbg.o $(CPPLIB_OBJ) #$(PREFIX)/lib/libgnuplot.so $(PREFIX)/lib/libxerces-c-3.0.so
	$(CC) $(LDFLAGS) $^ -o $@ -lgnuplot -lxerces-c-3.0 -lbfd -liberty
else
spsplot: spsplot.o specplot.o contplot.o spsstats.o mzxml.o $(CPPLIB_OBJ) #$(PREFIX)/lib/libgnuplot.so $(PREFIX)/lib/libxerces-c-3.0.so
	$(CC) $(LDFLAGS) $^ -o $@ -lgnuplot -lxerces-c-3.0
endif

specplot: spsplot_l.o specplot_l.o mzxml.o $(CPPLIB_OBJ) #$(PREFIX)/lib/libgnuplot.a $(PREFIX)/lib/libxerces-c.a
	$(CC) $(LDFLAGS) $^ -o $@ -Wl,-s -Wl,-Bstatic -lgnuplot -lxerces-c -Wl,-Bdynamic -lncurses -lgd

# Windows Cygwin
else
ifneq (,$(findstring cygwin,$(shell gcc -dumpmachine)))
$(PREFIX)/bin/gnuplot.dll: $(GNUPLOT_ROOT)/Makefile
	$(MAKE) -C $(GNUPLOT_ROOT)/src gnuplot.exe LDFLAGS="$(LDFLAGS) -shared -Wl,-o,gnuplot.dll" CFLAGS="$(CPPFLAGS) -fPIC -Dmain=gnuplot_main"
	mv $(GNUPLOT_ROOT)/src/gnuplot.dll $(PREFIX)/bin

$(PREFIX)/bin/cygxerces-c-3-0.dll: $(XERCES_ROOT)/Makefile
	$(MAKE) -C $(XERCES_ROOT)/src LDFLAGS="$(LDFLAGS)" CPPFLAGS="$(CPPFLAGS)"
	mv $(XERCES_ROOT)/src/.libs/cygxerces-c-3-0.dll $(PREFIX)/bin

spsplot: spsplot.o specplot.o contplot.o spsstats.o mzxml.o $(CPPLIB_OBJ) #$(PREFIX)/bin/gnuplot.dll $(PREFIX)/bin/cygxerces-c-3-0.dll
	$(CC) $(LDFLAGS) $^ -o $@ -lgnuplot -lxerces-c-3-0

specplot: spsplot_l.o specplot_l.o mzxml.o $(CPPLIB_OBJ) $(PREFIX)/lib/libgnuplot.a #$(PREFIX)/bin/cygxerces-c-3-0.dll
	$(CC) $(LDFLAGS) $^ -o $@
endif
endif

# dependencies
depends : $(CPPLIB_SRC) specplot.cpp contplot.cpp spsplot.cpp spsstats.cpp dbg.cpp mzxml.cpp
	@$(CC) -MM $^ > depends

# Linux
ifneq (,$(findstring linux,$(shell gcc -dumpmachine)))
-include depends
endif
