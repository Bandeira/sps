ROOT := $(shell pwd)
XERCES_ROOT := ../../../src/xerces-c-3.1.1
GNUPLOT_ROOT := ../../../src/gnuplot-4.2.5
SPS_VERSION := $(shell svnversion .)
INSTALL_DIR := ../../../install
REPORT_COMMON_FILES := ../common
COMMON_FILES := ../..


CC := g++
# CC := i686-w64-mingw32-g++

CPPFLAGS := -fPIC -DSPS_VERSION=$(SPS_VERSION) -I$(XERCES_ROOT)/src -I. -I$(REPORT_COMMON_FILES) -I$(COMMON_FILES)/include -I$(COMMON_FILES) -I$(COMMON_FILES)/ExecFramework -pthread 

LDFLAGS := 


# 'make type=32'
ifeq ($(type),32)
  LDFLAGS += -march=i686 -m32
  CPPFLAGS += -march=i686 -m32
  PREFIX := $(INSTALL_DIR)/linux-g++

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
    PREFIX ?= $(INSTALL_DIR)/linux-g++-64
  else
    PREFIX ?= $(INSTALL_DIR)/linux-g++
  endif

  LDFLAGS += -L$(PREFIX)/lib

# Windows
else
  ifneq (,$(findstring x86_64,$(shell gcc -dumpmachine)))
    PREFIX ?= $(INSTALL_DIR)/win32-g++
  else
    PREFIX ?= $(INSTALL_DIR)/win32-g++
  endif

  LDFLAGS += -L$(PREFIX)/lib -L$(PREFIX)/bin
endif




# type=debug
ifeq ($(type),debug)
  CPPFLAGS += -O0 -g
else
  CPPFLAGS += -O3 -ffunction-sections -fdata-sections -Wno-deprecated -Wno-write-strings 
endif


# type=debug
ifeq ($(type),debug)
  CPPFLAGS += -O0 -g
else
  CPPFLAGS += -O3 -ffunction-sections -fdata-sections -Wno-deprecated -Wno-write-strings 
endif



UU64_SRC = \
  uu64Encode.cpp \
  $(REPORT_COMMON_FILES)/base64.cpp


UU64_HDR = \
  uu64Encode.h \
  $(REPORT_COMMON_FILES)/base64.h


UU64_OBJ=$(UU64_SRC:.cpp=.o)


EXECUTABLES = \
	uu64Encode


all: $(EXECUTABLES)


clean:
	rm -f *.o $(EXECUTABLES) *.gch *~


# Linux

ifneq (,$(findstring linux,$(shell gcc -dumpmachine)))



uu64Encode: $(UU64_OBJ)
	$(CC) $(LDFLAGS) $^ -o $@ -Wall -static
	strip uu64Encode


# windows

else


uu64Encode: $(UU64_OBJ)
	$(CC) $(LDFLAGS) $^ -o $@ -Wall -static


# end
endif

