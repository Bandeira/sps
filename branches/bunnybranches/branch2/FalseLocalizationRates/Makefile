ROOT := $(shell pwd)
BOOST_ROOT := ../src/boost_1_44_0
LIBCORE_DIR := ../libCore
LIBCORE_SRC := $(LIBCORE_DIR)/libcore.a
#SPS_VERSION := $(shell date +%y%m%d)
#SPS_VERSION := $(shell svnversion .)

CC := g++
#CPPFLAGS := -fPIC -DSPS_VERSION=$(SPS_VERSION) -I.. -I$(LIBCORE_DIR)
CPPFLAGS := -fPIC -I.. -I$(LIBCORE_DIR)


type=debug
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
    PREFIX ?= ../install/win32-g++-64
  else
    PREFIX ?= ../install/win32-g++
  endif

  LDFLAGS += -L$(PREFIX)/lib -L$(PREFIX)/bin/cygwin
endif

CPPLIB_SRC = \
	FalseLocalizationRates.cpp \
	LinearEquation.cpp \
	Timer.cpp \
	../cosine_similarity.cpp 

CPPLIB_OBJ=$(CPPLIB_SRC:.cpp=.o)

STATICLIB = libflr.a

DYNAMICLIB = libflr.so


clean: 
	rm -f *.o $(STATICLIB)

.o.o: $(CPPLIB_SRC)
	$(CC) $(CPPFLAGS) $<

$(STATICLIB): $(CPPLIB_OBJ) $(LIBCORE_SRC)
	ar rcs $@ $^

#$(DYNAMICLIB): $(CPPLIB_OBJ) 
#	$(CC) -shared  -Wl,-rpath,$(ROOT) -Wl,-soname,$@.1 -o $@.1.0 $^
#	ln -s $@.1.0 $@.1
#	ln -s $@.1 $@