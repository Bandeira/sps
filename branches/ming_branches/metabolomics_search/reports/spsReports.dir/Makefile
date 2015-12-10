ROOT := $(shell pwd)
XERCES_ROOT := ../../../src/xerces-c-3.1.1
GNUPLOT_ROOT := ../../../src/gnuplot-4.2.5
SPS_VERSION := $(shell svnversion .)
INSTALL_DIR := ../../../install
REPORT_COMMON_FILES := ../common
COMMON_FILES := ../..

CONTPLOT_DIR := ../contplot.dir
SPECPLOT_DIR := ../specplot.dir


CC := g++
# CC := i686-w64-mingw32-g++

CPPFLAGS := -fPIC -DSPS_VERSION=$(SPS_VERSION) -I$(XERCES_ROOT)/src -I. -I$(REPORT_COMMON_FILES) -I$(COMMON_FILES)/include -I$(COMMON_FILES) -I$(COMMON_FILES)/ExecFramework -I$(CONTPLOT_DIR) -I$(SPECPLOT_DIR) -pthread 

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



SPS_SRC = \
	$(REPORT_COMMON_FILES)/ClusterData.cpp \
	$(REPORT_COMMON_FILES)/ReportTableGenerator.cpp \
	$(REPORT_COMMON_FILES)/ReportTableBase.cpp \
	$(REPORT_COMMON_FILES)/ReportTableHeader.cpp \
	$(REPORT_COMMON_FILES)/ReportTableProtein.cpp \
	$(REPORT_COMMON_FILES)/ReportTableProteinCoverage.cpp \
	$(REPORT_COMMON_FILES)/ReportTableContig.cpp \
	$(REPORT_COMMON_FILES)/ReportTableClusterConsensus.cpp \
	$(REPORT_COMMON_FILES)/ReportTableInputSpectra.cpp \
	ReportBase.cpp \
	ReportHeader.cpp \
	ReportProtein.cpp \
	ReportContig.cpp \
  ReportCluster.cpp \
	ReportInputSpectra.cpp \
	ReportRendererBase.cpp \
	ReportRendererHtml.cpp \
	ReportRendererHtmlDynamic.cpp \
	ReportModuleBase.cpp \
	ReportModuleFactory.cpp \
	ReportModuleSpecplot.cpp \
	ReportModuleContplot.cpp \
	$(REPORT_COMMON_FILES)/Tokenizer.cpp \
	$(COMMON_FILES)/abruijn.cpp \
	$(COMMON_FILES)/SpectrumAnnotStatistics.cpp \
	$(COMMON_FILES)/spectrum_scoring.cpp \
	ReportInterface.cpp \
	spsReport.cpp


SPS_HRD= \
	$(REPORT_COMMON_FILES)/ClusterData.h \
	$(REPORT_COMMON_FILES)/ReportTableGenerator.h \
	$(REPORT_COMMON_FILES)/ReportTableBase.h \
	$(REPORT_COMMON_FILES)/ReportTableHeader.h \
	$(REPORT_COMMON_FILES)/ReportTableProtein.h \
	$(REPORT_COMMON_FILES)/ReportTableProteinCoverage.h \
	$(REPORT_COMMON_FILES)/ReportTableContig.h \
	$(REPORT_COMMON_FILES)/ReportTableClusterConsensus.h \
	$(REPORT_COMMON_FILES)/ReportTableInputSpectra.h \
	ReportBase.h \
	ReportHeader.h \
	ReportProtein.h \
	ReportContig.h \
  ReportCluster.h \
	ReportInputSpectra.h \
	ReportRendererBase.h \
	ReportRendererHtml.h \
	ReportRendererHtmlDynamic.h \
	ReportModuleBase.h \
	ReportModuleFactory.h \
	ReportModuleSpecplot.h \
	ReportModuleContplot.h \
	$(REPORT_COMMON_FILES)/Tokenizer.h \
	$(COMMON_FILES)/abruijn.h \
	$(COMMON_FILES)/SpectrumAnnotStatistics.h \
	$(COMMON_FILES)/spectrum_scoring.h \
	ReportInterface.h \
	spsReport.h
	
	
CONT_SRC = \
	$(REPORT_COMMON_FILES)/Tokenizer.cpp \
  $(REPORT_COMMON_FILES)/base64.cpp \
	$(REPORT_COMMON_FILES)/RendererBase.cpp \
	$(REPORT_COMMON_FILES)/RendererGnu.cpp \
	$(REPORT_COMMON_FILES)/PlotBase.cpp \
	$(CONTPLOT_DIR)/ContigSequence.cpp \
	$(CONTPLOT_DIR)/ContigSpectrum.cpp \
	$(CONTPLOT_DIR)/PlotContig.cpp \
	$(CONTPLOT_DIR)/ContplotInterface.cpp

CONT_HDR = \
	$(REPORT_COMMON_FILES)/Tokenizer.h \
  $(REPORT_COMMON_FILES)/base64.h \
	$(REPORT_COMMON_FILES)/RendererBase.h \
	$(REPORT_COMMON_FILES)/RendererGnu.h \
	$(REPORT_COMMON_FILES)/PlotBase.h \
	$(CONTPLOT_DIR)/ContigSequence.h \
	$(CONTPLOT_DIR)/ContigSpectrum.h \
	$(CONTPLOT_DIR)/PlotContig.h \
	$(CONTPLOT_DIR)/contplot.h \
	$(CONTPLOT_DIR)/ContplotInterface.h
	
	

SPEC_SRC = \
  $(REPORT_COMMON_FILES)/base64.cpp \
	$(REPORT_COMMON_FILES)/RendererBase.cpp \
	$(REPORT_COMMON_FILES)/RendererGnu.cpp \
	$(REPORT_COMMON_FILES)/PlotBase.cpp \
	$(SPECPLOT_DIR)/PlotSpectrum.cpp \
	$(SPECPLOT_DIR)/SpecplotInterface.cpp

SPEC_HDR = \
  $(REPORT_COMMON_FILES)/base64.h \
	$(REPORT_COMMON_FILES)/RendererBase.h \
	$(REPORT_COMMON_FILES)/RendererGnu.h \
	$(REPORT_COMMON_FILES)/PlotBase.h \
	$(SPECPLOT_DIR)/PlotSpectrum.h \
	$(SPECPLOT_DIR)/SpecplotInterface.h



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



SPEC_OBJ=$(SPEC_SRC:.cpp=.o)

SPS_OBJ=$(SPS_SRC:.cpp=.o)

CONT_OBJ=$(CONT_SRC:.cpp=.o)

CPPLIB_OBJ=$(CPPLIB_SRC:.cpp=.o)



EXECUTABLES = \
	spsReport


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


$(PREFIX)/lib/libxerces-c.a: $(XERCES_ROOT)/Makefile
	$(MAKE) -C $(XERCES_ROOT)/src LDFLAGS="$(LDFLAGS)" CPPFLAGS="$(CPPFLAGS)"
	mv $(XERCES_ROOT)/src/.libs/libxerces-c.a $(PREFIX)/lib


$(GNUPLOT_ROOT)/Makefile:
	cd $(GNUPLOT_ROOT); ./configure --without-readline --without-lisp-files --without-tutorial --without-row-help --without-x --disable-wxwidgets; cd $(ROOT)




# Linux

ifneq (,$(findstring linux,$(shell gcc -dumpmachine)))



spsReport: $(SPS_OBJ) $(SPEC_OBJ) $(CPPLIB_OBJ) $(CONT_OBJ)
	$(CC) $(CPPFLAGS) $(LDFLAGS) $^ -o $@ -static -lxerces-c



# windows

else

$(PREFIX)/bin/gnuplot.dll: $(GNUPLOT_ROOT)/Makefile
	$(MAKE) -C $(GNUPLOT_ROOT)/src gnuplot.exe LDFLAGS="$(LDFLAGS) -shared -Wl,-o,gnuplot.dll" CFLAGS="$(CPPFLAGS) -fPIC -Dmain=gnuplot_main"
	mv $(GNUPLOT_ROOT)/src/gnuplot.dll $(PREFIX)/bin


spsReport: $(SPS_OBJ) $(SPEC_OBJ) $(CPPLIB_OBJ) $(CONT_OBJ) $(PREFIX)/bin/libxerces-c.a 
	$(CC) $(CPPFLAGS) $(LDFLAGS) $^ -o $@ -static -lxerces-c



# end
endif

