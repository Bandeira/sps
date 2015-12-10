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



SPS_SERVER_SRC = \
	$(COMMON_FILES)/db_fasta.cpp \
  $(COMMON_FILES)/ExecFramework/CommandLineParser.cpp \
  $(COMMON_FILES)/ExecFramework/ParameterList.cpp  \
  $(COMMON_FILES)/aminoacid.cpp \
	$(COMMON_FILES)/spectrum.cpp \
  $(COMMON_FILES)/Logger.cpp \
	$(COMMON_FILES)/utils.cpp \
	$(COMMON_FILES)/mzrange.cpp \
	$(COMMON_FILES)/ms1.cpp \
	$(COMMON_FILES)/PeptideSpectrumMatch.cpp \
	$(COMMON_FILES)/PeptideSpectrumMatchSet.cpp \
	$(REPORT_COMMON_FILES)/ReportTableBase.cpp \
	$(REPORT_COMMON_FILES)/ReportTableHeader.cpp \
	$(REPORT_COMMON_FILES)/ReportTableProtein.cpp \
	$(REPORT_COMMON_FILES)/ReportTableProteinCoverage.cpp \
	$(REPORT_COMMON_FILES)/ReportTableContig.cpp \
	$(REPORT_COMMON_FILES)/ReportTableClusterConsensus.cpp \
	$(REPORT_COMMON_FILES)/ReportTableInputSpectra.cpp \
	ReportServerInterface.cpp \
	spsReportServer.cpp \
	$(COMMON_FILES)/spectrum.cpp \
	$(COMMON_FILES)/SpecSet.cpp

SPS_SERVER_HDR = \
	$(COMMON_FILES)/db_fasta.h \
  $(COMMON_FILES)/ExecFramework/CommandLineParser.h \
  $(COMMON_FILES)/ExecFramework/ParameterList.h  \
  $(COMMON_FILES)/aminoacid.h \
  $(COMMON_FILES)/alignment_scoring.h \
  $(COMMON_FILES)/DelimitedTextReader.h \
  $(COMMON_FILES)/inputParams.h \
  $(COMMON_FILES)/Logger.h \
	$(COMMON_FILES)/utils.h \
	$(COMMON_FILES)/mzrange.h \
	$(COMMON_FILES)/ms1.h \
	$(COMMON_FILES)/PeptideSpectrumMatch.h \
	$(COMMON_FILES)/PeptideSpectrumMatchSet.h \
	$(REPORT_COMMON_FILES)/ReportTableBase.h \
	$(REPORT_COMMON_FILES)/ReportTableHeader.h \
	$(REPORT_COMMON_FILES)/ReportTableProtein.h \
	$(REPORT_COMMON_FILES)/ReportTableProteinCoverage.h \
	$(REPORT_COMMON_FILES)/ReportTableContig.h \
	$(REPORT_COMMON_FILES)/ReportTableClusterConsensus.h \
	$(REPORT_COMMON_FILES)/ReportTableInputSpectra.h \
	ReportServerInterface.h \
	spsReportServer.h \
	$(COMMON_FILES)/spectrum.h \
	$(COMMON_FILES)/SpecSet.h




CPPLIB_SRC = \
  $(COMMON_FILES)/ExecFramework/CommandLineParser.cpp \
  $(COMMON_FILES)/ExecFramework/ParameterList.cpp  \
  $(COMMON_FILES)/alignment_scoring.cpp \
  $(COMMON_FILES)/aminoacid.cpp \
  $(COMMON_FILES)/base64.cpp \
  $(COMMON_FILES)/db_fasta.cpp \
  $(COMMON_FILES)/abruijn.cpp \
  $(COMMON_FILES)/inputParams.cpp \
  $(COMMON_FILES)/DelimitedTextReader.cpp \
	$(COMMON_FILES)/label.cpp \
  $(COMMON_FILES)/Logger.cpp \
  $(COMMON_FILES)/ms1.cpp \
  $(COMMON_FILES)/PeptideSpectrumMatch.cpp \
  $(COMMON_FILES)/PeptideSpectrumMatchSet.cpp \
	$(COMMON_FILES)/spectrum.cpp \
	$(COMMON_FILES)/SpecSet.cpp \
	$(COMMON_FILES)/SpectrumPair.cpp \
	$(COMMON_FILES)/SpectrumPairSet.cpp \
  $(COMMON_FILES)/spectrum_scoring.cpp \
	$(COMMON_FILES)/utils.cpp \
	$(COMMON_FILES)/mzxml.cpp \
	$(COMMON_FILES)/mzrange.cpp

CPPLIB_HDR = \
  $(COMMON_FILES)/ExecFramework/CommandLineParser.h \
  $(COMMON_FILES)/ExecFramework/ParameterList.h  \
  $(COMMON_FILES)/alignment_scoring.h \
  $(COMMON_FILES)/aminoacid.h \
  $(COMMON_FILES)/base64.h \
  $(COMMON_FILES)/db_fasta.h \
  $(COMMON_FILES)/abruijn.h \
  $(COMMON_FILES)/inputParams.h \
  $(COMMON_FILES)/DelimitedTextReader.h \
	$(COMMON_FILES)/label.h \
  $(COMMON_FILES)/Logger.h \
  $(COMMON_FILES)/ms1.h \
  $(COMMON_FILES)/PeptideSpectrumMatch.h \
  $(COMMON_FILES)/PeptideSpectrumMatchSet.h \
	$(COMMON_FILES)/spectrum.h \
	$(COMMON_FILES)/SpecSet.h \
	$(COMMON_FILES)/SpectrumPair.h \
	$(COMMON_FILES)/SpectrumPairSet.h \
  $(COMMON_FILES)/spectrum_scoring.h \
	$(COMMON_FILES)/utils.h \
	$(COMMON_FILES)/mzxml.h \
	$(COMMON_FILES)/mzrange.h



SPS_SERVER_OBJ=$(SPS_SERVER_SRC:.cpp=.o)

CPPLIB_OBJ=$(CPPLIB_SRC:.cpp=.o)



EXECUTABLES = \
	spsReportServer 



all: $(EXECUTABLES)
	strip  $(EXECUTABLES)


clean:
	rm -f *.o $(EXECUTABLES) *.gch *~


.cpp.o: $(CPPLIB_SRC)
	$(CC) $(CPPFLAGS) $< -c -o $@

.o.o: $(CPPLIB_SRC) $(CPPLIB_HDR)
	$(CC) $(CPPFLAGS) $<	


$(XERCES_ROOT)/Makefile:
	cd $(XERCES_ROOT); ./configure --disable-pretty-make --disable-network --disable-threads --without-curl --without-icu; cd $(ROOT)



# Linux

ifneq (,$(findstring linux,$(shell gcc -dumpmachine)))


spsReportServer: $(SPS_SERVER_OBJ)
	$(CC) $(CPPFLAGS) $(LDFLAGS) $^ -o $@ -static



# windows

else


spsReportServer: $(SPS_SERVER_OBJ)  $(CPPLIB_OBJ)
	$(CC) $(CPPFLAGS) $(LDFLAGS) $^ -o $@ -static -lxerces-c


# end
endif

