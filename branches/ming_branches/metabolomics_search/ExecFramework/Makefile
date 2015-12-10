ROOT := $(shell pwd)
XERCES_ROOT := ../../src/xerces-c-3.1.1
GNUPLOT_ROOT := ../../src/gnuplot-4.2.5
#SPS_VERSION := $(shell date +%y%m%d)
SPS_VERSION := $(shell svnversion .)
INSTALL_DIR := ../../install
REPORTS_DIR := ../reports
METASPS_DIR := ../MetaSPS
REPORT_COMMON_FILES := $(REPORTS_DIR)/common
CONTPLOT_DIR := $(REPORTS_DIR)/contplot.dir
SPECPLOT_DIR := $(REPORTS_DIR)/specplot.dir
REPORT_DIR := $(REPORTS_DIR)/spsReports.dir
COMMON_FILES := ..


CC := g++
# CC := i686-w64-mingw32-g++

CPPFLAGS := -DSPS_VERSION=$(SPS_VERSION) -pthread -I. -I$(COMMON_FILES) -I$(REPORT_COMMON_FILES) -I$(XERCES_ROOT)/src -I$(CONTPLOT_DIR) -I$(SPECPLOT_DIR) -I$(REPORT_DIR) -I$(METASPS_DIR)

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

  LDFLAGS += -L. -L.. -L$(PREFIX)/lib

# Windows
else
  ifneq (,$(findstring x86_64,$(shell gcc -dumpmachine)))
    PREFIX ?= $(INSTALL_DIR)/win32-g++
  else
    PREFIX ?= $(INSTALL_DIR)/win32-g++
  endif

  LDFLAGS += -L. -L.. -L$(PREFIX)/lib -L$(PREFIX)/bin
endif


LIBRARIES = -lpthread -lspecnets -lxerces-c 


CPPLIB_SRC = \
  AlignmentUtils.cpp \
  CommandLineParser.cpp \
  ExecAlignment.cpp \
  ExecAssembly.cpp \
  ExecBase.cpp \
  ExecGFNetwork.cpp \
  ExecSpecNetworkEval.cpp \
  ExecCreateSpectralLibrary.cpp \
  ExecDeconvoluteMS2.cpp \
  ExecFilterAligns.cpp \
  ExecFdrPeptide.cpp \
  ExecFilterContigPairs.cpp \
  ExecFilterPairs.cpp \
  ExecFilterStarPairs.cpp \
  ExecHomologyAssembly.cpp \
  ExecMainSpecnets.cpp \
  ExecMergeConvert.cpp \
  ExecMergePaired.cpp \
  ExecModuleFactory.cpp \
  ExecMsCluster.cpp \
  ExecReportProteinCoverage.cpp \
  ExecReportSpsplot.cpp \
  ExecReportSPSStats.cpp \
  ExecPrmScoring.cpp \
  ExecProjection.cpp \
  ExecProjectionStatistics.cpp \
  ExecSpectralLibrarySearch.cpp \
  ExecSpectralLibraryUnannotatedSearch.cpp \
  ExecSpectralLibraryDecoyCreation.cpp \
  ExecSpecNetsPropagation.cpp \
  ExecSpectraExtraction.cpp \
  ExecStatistics.cpp \
  ExecSvmStatistics.cpp \
  ExecTagSearch.cpp \
  ExecSpecProtAlign.cpp \
  ExecProtProtAlign.cpp \
  FileUtils.cpp \
  Filters.cpp \
  $(COMMON_FILES)/utils.cpp \
  $(COMMON_FILES)/Logger.cpp \
  $(COMMON_FILES)/mzxml.cpp \
  ParameterList.cpp \
  ParallelExecution.cpp \
  ParallelPbsExecution.cpp \
  ParallelSgeExecution.cpp \
  ParallelThreadedExecution.cpp \
	$(COMMON_FILES)/abruijn.cpp \
	$(COMMON_FILES)/SpectrumAnnotStatistics.cpp \
	$(COMMON_FILES)/spectrum_scoring.cpp \
	$(REPORT_COMMON_FILES)/ClusterData.cpp \
	$(REPORT_COMMON_FILES)/ReportTableGenerator.cpp \
	$(REPORT_COMMON_FILES)/ReportTableBase.cpp \
	$(REPORT_COMMON_FILES)/ReportTableHeader.cpp \
	$(REPORT_COMMON_FILES)/ReportTableProtein.cpp \
	$(REPORT_COMMON_FILES)/ReportTableProteinCoverage.cpp \
	$(REPORT_COMMON_FILES)/ReportTableContig.cpp \
	$(REPORT_COMMON_FILES)/ReportTableClusterConsensus.cpp \
	$(REPORT_COMMON_FILES)/ReportTableInputSpectra.cpp \
	$(REPORT_DIR)/ReportBase.cpp \
	$(REPORT_DIR)/ReportHeader.cpp \
	$(REPORT_DIR)/ReportProtein.cpp \
	$(REPORT_DIR)/ReportContig.cpp \
  $(REPORT_DIR)/ReportCluster.cpp \
	$(REPORT_DIR)/ReportInputSpectra.cpp \
	$(REPORT_DIR)/ReportRendererBase.cpp \
	$(REPORT_DIR)/ReportRendererHtml.cpp \
	$(REPORT_DIR)/ReportRendererHtmlDynamic.cpp \
	$(REPORT_DIR)/ReportModuleBase.cpp \
	$(REPORT_DIR)/ReportModuleFactory.cpp \
	$(REPORT_DIR)/ReportModuleSpecplot.cpp \
	$(REPORT_DIR)/ReportModuleContplot.cpp \
	$(REPORT_COMMON_FILES)/Tokenizer.cpp \
	$(REPORT_DIR)/ReportInterface.cpp \
  $(REPORT_COMMON_FILES)/base64.cpp \
	$(REPORT_COMMON_FILES)/RendererBase.cpp \
	$(REPORT_COMMON_FILES)/RendererGnu.cpp \
	$(SPECPLOT_DIR)/PlotSpectrum.cpp \
	$(SPECPLOT_DIR)/SpecplotInterface.cpp \
	$(REPORT_COMMON_FILES)/RendererBase.cpp \
	$(REPORT_COMMON_FILES)/PlotBase.cpp \
	$(CONTPLOT_DIR)/ContigSequence.cpp \
	$(CONTPLOT_DIR)/ContigSpectrum.cpp \
	$(CONTPLOT_DIR)/PlotContig.cpp \
	$(CONTPLOT_DIR)/ContplotInterface.cpp \
  SgeGridMonitor.cpp \
  ExecGenoMS.cpp
#  SetMerger.cpp \
#  SpectralPairs.cpp \
#  SpectrumPair.cpp \
#  SpectrumPairSet.cpp \

CPPLIB_HDR = \
  AlignmentUtils.h \
  CommandLineParser.h \
  ExecAlignment.h \
  ExecAssembly.h \
  ExecBase.h \
  ExecGFNetwork.h \
  ExecSpecNetworkEval.h \
  ExecCreateSpectralLibrary.h \
  ExecDeconvoluteMS2.h \
  ExecFilterAligns.h \
  ExecFdrPeptide.h \
  ExecFilterContigPairs.h \
  ExecFilterPairs.h \
  ExecFilterStarPairs.h \
  ExecHomologyAssembly.h \
  ExecMainSpecnets.h \
  ExecMergeConvert.h \
  ExecMergePaired.h \
  ExecModuleFactory.h \
  ExecMsCluster.h \
  ExecReportProteinCoverage.h \
  ExecReportSpsplot.h \
  ExecReportSPSStats.h \
  ExecProjection.h \
  ExecProjectionStatistics.h \
  ExecSpectralLibrarySearch.h \
  ExecSpectralLibraryUnannotatedSearch.h \
  ExecSpectralLibraryDecoyCreation.h \
  ExecSpectraExtraction.h \
  ExecPrmScoring.h \
  ExecTagSearch.h \
  ExecSpecProtAlign.h \
  ExecProtProtAlign.h \
  ExecSpecNetsPropagation.h \
  ExecStatistics.h \
  ExecSvmStatistics.h \
  FileUtils.h \
  Filters.h \
  $(COMMON_FILES)/utils.h \
  $(COMMON_FILES)/Logger.h \
	$(COMMON_FILES)/mzxml.h \
  ParameterList.h \
  ParallelExecution.h \
  ParallelPbsExecution.h \
  ParallelSgeExecution.h \
  ParallelThreadedExecution.h \
	$(COMMON_FILES)/abruijn.h \
	$(COMMON_FILES)/SpectrumAnnotStatistics.h \
	$(COMMON_FILES)/spectrum_scoring.h \
	$(REPORT_COMMON_FILES)/ClusterData.h \
	$(REPORT_COMMON_FILES)/ReportTableGenerator.h \
	$(REPORT_COMMON_FILES)/ReportTableBase.h \
	$(REPORT_COMMON_FILES)/ReportTableHeader.h \
	$(REPORT_COMMON_FILES)/ReportTableProtein.h \
	$(REPORT_COMMON_FILES)/ReportTableProteinCoverage.h \
	$(REPORT_COMMON_FILES)/ReportTableContig.h \
	$(REPORT_COMMON_FILES)/ReportTableClusterConsensus.h \
	$(REPORT_COMMON_FILES)/ReportTableInputSpectra.h \
	$(REPORT_DIR)/ReportBase.h \
	$(REPORT_DIR)/ReportHeader.h \
	$(REPORT_DIR)/ReportProtein.h \
	$(REPORT_DIR)/ReportContig.h \
  $(REPORT_DIR)/ReportCluster.h \
	$(REPORT_DIR)/ReportInputSpectra.h \
	$(REPORT_DIR)/ReportRendererBase.h \
	$(REPORT_DIR)/ReportRendererHtml.h \
	$(REPORT_DIR)/ReportRendererHtmlDynamic.h \
	$(REPORT_DIR)/ReportModuleBase.h \
	$(REPORT_DIR)/ReportModuleFactory.h \
	$(REPORT_DIR)/ReportModuleSpecplot.h \
	$(REPORT_DIR)/ReportModuleContplot.h \
	$(REPORT_COMMON_FILES)/Tokenizer.h \
	$(REPORT_DIR)/ReportInterface.h \
  $(REPORT_COMMON_FILES)/base64.h \
	$(REPORT_COMMON_FILES)/RendererBase.h \
	$(REPORT_COMMON_FILES)/RendererGnu.h \
	$(REPORT_COMMON_FILES)/PlotBase.h \
	$(SPECPLOT_DIR)/PlotSpectrum.h \
	$(SPECPLOT_DIR)/SpecplotInterface.h \
	$(CONTPLOT_DIR)/ContigSequence.h \
	$(CONTPLOT_DIR)/ContigSpectrum.h \
	$(CONTPLOT_DIR)/PlotContig.h \
	$(CONTPLOT_DIR)/contplot.h \
	$(CONTPLOT_DIR)/ContplotInterface.h \
  SgeGridMonitor.h \
  ExecGenoMS.h
#  SpectralPairs.h \
#  SetMerger.h \
#  SpectrumPair.h \
#  SpectrumPairSet.h \

CPPLIB_OBJ=$(CPPLIB_SRC:.cpp=.o)

EXECUTABLES = \
	main_diffspecs \
	main_specdump \
	main_specnets \
	main_execmodule

.PHONY: subdirs $(SUBDIRS) clean all

all: $(EXECUTABLES) main_specnets main_execmodule
	strip $(EXECUTABLES) main_specnets main_execmodule

clean:
	rm -f *.o $(EXECUTABLES) main_specnets main_execmodule 

.cpp.o: $(CPPLIB_SRC)
	$(CC) $(CPPFLAGS) $< -c -o $@

.o.o: $(CPPLIB_SRC) $(CPPLIB_HDR)
	$(CC) $(CPPFLAGS) $<

# Linux
ifneq (,$(findstring linux,$(shell gcc -dumpmachine)))
#convert: convert.o mzxml.o $(CPPLIB_OBJ)
#	$(CC) $(LDFLAGS) $^ -o $@ -lxerces-c-3.0

# Windows
else
#convert: convert.o mzxml.o $(CPPLIB_OBJ)
#	$(CC) $(LDFLAGS) $^ -o $@ -static -lxerces-c-3-0
endif

main_diffspecs: main_diffspecs.o ParameterList.o $(COMMON_FILES)/libspecnets.a
	$(CC) $(CPPFLAGS) $(LDFLAGS) -static $^ -o $@ $(LIBRARIES)

main_specdump: main_specdump.o  $(COMMON_FILES)/libspecnets.a
	$(CC) $(CPPFLAGS) $(LDFLAGS) -static $^ ParameterList.o -o $@ $(LIBRARIES)

main_specnets: main_specnets.o $(CPPLIB_OBJ) $(COMMON_FILES)/libspecnets.a
	$(CC) $(CPPFLAGS) $(LDFLAGS) -static $^ -o $@ $(LIBRARIES)

main_execmodule: main_execmodule.o $(CPPLIB_OBJ) $(COMMON_FILES)/libspecnets.a
	$(CC) $(CPPFLAGS) $(LDFLAGS) $^ -o $@ $(LIBRARIES)


# dependencies
depends : $(CPPLIB_SRC)  $(CPPLIB_HDR)
	@$(CC) $(CPPFLAGS) -MM $^ > depends

# Linux
ifneq (,$(findstring linux,$(shell gcc -dumpmachine)))
-include depends
endif

