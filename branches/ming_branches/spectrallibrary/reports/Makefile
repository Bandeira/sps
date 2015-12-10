ROOT := $(shell pwd)
XERCES_ROOT := ../../src/xerces-c-3.1.1
GNUPLOT_ROOT := ../../src/gnuplot-4.2.5
SPS_VERSION := $(shell svnversion .)

CC := g++
# CC := i686-w64-mingw32-g++

CPPFLAGS := -fPIC -DSPS_VERSION=$(SPS_VERSION) -I$(XERCES_ROOT)/src -I. -I../include -I.. -I../ExecFramework -pthread
LDFLAGS := 

# 'make type=32'
ifeq ($(type),32)
  LDFLAGS += -march=i686 -m32
  CPPFLAGS += -march=i686 -m32
  PREFIX := ../../install/linux-g++

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
    PREFIX ?= ../../install/linux-g++-64
  else
    PREFIX ?= ../../install/linux-g++
  endif

  LDFLAGS += -L$(PREFIX)/lib

# Windows
else
  ifneq (,$(findstring x86_64,$(shell gcc -dumpmachine)))
    PREFIX ?= ../../install/win32-g++
  else
    PREFIX ?= ../../install/win32-g++
  endif

  LDFLAGS += -L$(PREFIX)/lib -L$(PREFIX)/bin/cygwin
endif




# type=debug
ifeq ($(type),debug)
  CPPFLAGS += -O0 -g
else
  CPPFLAGS += -O3 -ffunction-sections -fdata-sections -Wno-deprecated -Wno-write-strings 
endif




SUBDIRS = \
  spsReports.dir \
  contplot.dir \
  specplot.dir \
  uu64Encode.dir \
  clusterInfo.dir \
  spsReportServer.dir


.PHONY: subdirs $(SUBDIRS) clean all

all: subdirs
	cp spsReports.dir/spsReport spsReport
	cp spsReportServer.dir/spsReportServer spsReportServer
	cp specplot.dir/specplot specplot
	cp contplot.dir/contplot contplot
	cp uu64Encode.dir/uu64Encode uu64Encode
	cp clusterInfo.dir/clusterinfo clusterinfo


subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@
	
clean:
	cd spsReports.dir ; make clean ; cd ..
	cd spsReportServer.dir ; make clean ; cd ..
	cd specplot.dir ; make clean ; cd ..
	cd contplot.dir ; make clean ; cd ..
	cd uu64Encode.dir ; make clean ; cd ..
	cd clusterInfo.dir ; make clean ; cd ..
	rm common/*.o
	rm spsReport spsReportServer specplot contplot uu64Encode clusterinfo


