ROOT := $(shell pwd)
#XERCES_ROOT := ../../src/xerces-c-3.0.1
#XERCES_ROOT := ../../src/xerces-c-3.1.1
GNUPLOT_ROOT := ../../src/gnuplot-4.2.5
BOOST_ROOT := ../../src/boost_1_44_0
#SPS_VERSION := $(shell date +%y%m%d)
SPS_VERSION := $(shell cd ..; svnversion .; cd libCore)

CC := g++
# CC := i686-w64-mingw32-g++

CPPFLAGS := -fPIC -DSPS_VERSION=$(SPS_VERSION)

#type=debug
ifeq ($(type),debug)
  CPPFLAGS += -O0 -g
else
  CPPFLAGS += -O3 -ffunction-sections -fdata-sections -Wno-deprecated -Wno-write-strings
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
    PREFIX ?= ../install/win32-g++
  else
    PREFIX ?= ../install/win32-g++
  endif

  LDFLAGS += -L$(PREFIX)/lib -L$(PREFIX)/bin
endif

CPPLIB_SRC = \
	aminoacid.cpp \
	DeconvSpectrum.cpp \
	DelimitedTextReader.cpp \
	inputParams.cpp \
	IsoEnvelope.cpp \
	label.cpp \
	Logger.cpp \
	mzrange.cpp \
	PeptideSpectrumMatch.cpp \
	PeptideSpectrumMatchSet.cpp \
	SpecSet.cpp \
	spectrum_scoring.cpp \
	spectrum.cpp \
	utils.cpp

CPPLIB_OBJ=$(CPPLIB_SRC:.cpp=.o)

STATICLIB = libcore.a

DYNAMICLIB = libcore.so

all: $(STATICLIB)

clean:
	rm -f *.o *.a

$(STATICLIB): $(CPPLIB_OBJ)
	ar rcs $@ $^

$(DYNAMICLIB): $(CPPLIB_OBJ)
	$(CC) -shared  -Wl,-rpath,$(ROOT) -Wl,-soname,$@.1 -o $@.1.0 $^
	ln -s $@.1.0 $@.1
	ln -s $@.1 $@
