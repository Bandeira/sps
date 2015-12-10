ROOT := $(shell pwd)
XERCES_ROOT := ../../src/xerces-c-3.1.1
GNUPLOT_ROOT := ../../src/gnuplot-4.2.5
SPS_VERSION := $(shell svnversion .)
INSTALL_DIR := ../../../install
REPORT_COMMON_FILES := ../common
COMMON_FILES := ../..
LIBCORE_DIR := $(COMMON_FILES)/libCore
LIBCORE_SRC := $(LIBCORE_DIR)/libcore.a


CC := g++
# CC := i686-w64-mingw32-g++

CPPFLAGS := -fPIC -DSPS_VERSION=$(SPS_VERSION) -I$(XERCES_ROOT)/src -I. -I$(REPORT_COMMON_FILES) -I$(COMMON_FILES)/include -I$(COMMON_FILES) -I$(COMMON_FILES)/ExecFramework -I$(LIBCORE_DIR) -pthread 

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

  LDFLAGS += -L$(PREFIX)/lib -L$(PREFIX)/bin/cygwin
endif




# type=debug
ifeq ($(type),debug)
  CPPFLAGS += -O0 -g
else
  CPPFLAGS += -O3 -ffunction-sections -fdata-sections -Wno-deprecated -Wno-write-strings 
endif



CLUSTER_SRC = \
  ClusterInfo.cpp \
  ClusterInfoInterface.cpp \
  $(REPORT_COMMON_FILES)/ClusterData.cpp \
  $(COMMON_FILES)/ExecFramework/CommandLineParser.cpp \
  $(COMMON_FILES)/ExecFramework/ParameterList.cpp

CLUSTER_HDR = \
  ClusterInfo.h \
  ClusterInfoInterface.h \
  $(REPORT_COMMON_FILES)/ClusterData.h \
  $(COMMON_FILES)/ExecFramework/CommandLineParser.h \
  $(COMMON_FILES)/ExecFramework/ParameterList.h


CLUSTER_OBJ=$(CLUSTER_SRC:.cpp=.o)



EXECUTABLES = \
	clusterinfo



all: $(EXECUTABLES)
	strip  $(EXECUTABLES)


clean:
	rm -f *.o $(EXECUTABLES) *.gch *~



# Linux

ifneq (,$(findstring linux,$(shell gcc -dumpmachine)))

clusterinfo: $(CLUSTER_OBJ) $(LIBCORE_SRC)
	$(CC) $(LDFLAGS) $^ -o $@ -Wall -static


# windows

else

clusterinfo: $(CLUSTER_OBJ) $(LIBCORE_SRC)
	$(CC) $(LDFLAGS) $^ -o $@ -Wall -static


# end
endif

