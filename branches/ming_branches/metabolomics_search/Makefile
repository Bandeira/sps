ROOT := $(shell pwd)
#XERCES_ROOT := ../src/xerces-c-3.0.1
XERCES_ROOT := ../src/xerces-c-3.1.1
GNUPLOT_ROOT := ../src/gnuplot-4.2.5
BOOST_ROOT := ../src/boost_1_44_0
#SPS_VERSION := $(shell date +%y%m%d)
SPS_VERSION := $(shell svnversion .)

CC := g++
# CC := i686-w64-mingw32-g++

CPPFLAGS := -fPIC -DSPS_VERSION=$(SPS_VERSION) -I$(XERCES_ROOT)/src -IExecFramework -I$(ROOT)/MetaSPS -I$(ROOT)

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
	abruijn.cpp \
	alignment_modmut.cpp \
	alignment_scoring.cpp \
	aminoacid.cpp \
	base64.cpp \
	batch.cpp \
	clusters.cpp \
	db_fasta.cpp \
	DeconvSpectrum.cpp \
	DelimitedTextReader.cpp \
	denovo.cpp \
	FdrPeptide.cpp \
	filters.cpp \
	graph.cpp \
	inputParams.cpp \
	inspect_parse.cpp \
	label.cpp \
	MappedContig.cpp \
	MappedContigStatTable.cpp \
	MappedSpecnets.cpp \
	MappedSpectrum.cpp \
	MappedSPSStatTable.cpp \
	ms1.cpp \
	msn.cpp \
	mzrange.cpp \
	OutputTable.cpp \
	projectionutils.cpp \
	PepnovoTags.cpp \
	PeptideSpectrumMatch.cpp \
	PeptideSpectrumMatchNetwork.cpp \
	PeptideSpectrumMatchSet.cpp \
	SetMerger.cpp \
	SpecSet.cpp \
	SpectralPairs.cpp \
	spectrum.cpp \
	spectrum_window_filter.cpp \
	spectrum_scoring.cpp \
	SpectrumPair.cpp \
	SpectrumAnnotParameterList.cpp \
	SpectrumAnnotStatistics.cpp \
	SpectrumLibrary.cpp \
	SpectrumPairSet.cpp \
	SvmScaleParameterList.cpp \
	SpectralLibrary.cpp \
	SvmModel.cpp \
	tags.cpp \
	ion.cpp \
	Logger.cpp \
	mzxml.cpp \
	utils.cpp \
		MetaSPS/prm_alignment.cpp \
		MetaSPS/PRMSpecEdge.cpp \
		MetaSPS/PairedSpecSet.cpp \
		MetaSPS/Contig.cpp \
		MetaSPS/ContigNetwork.cpp
#	dekel_align.cpp \
#	hash.cpp \

CPPLIB_OBJ=$(CPPLIB_SRC:.cpp=.o)

STATICLIB = libspecnets.a

DYNAMICLIB = libspecnets.so

EXECUTABLES = \
	convert
#	MetaSPS/makeClustersFromAbinfo \
#	mergeConvert \
#	filterpairs \
#	starden \
#	filterstarpairs \
#	tagsearch \
#	protid \
#	matchma \
#	mergefilter \
#	homglue \
#	pathproj \
#	starden2 \

SCRIPTS = \
	run_csps.sh \
	run_specnets.sh \
	grid.sh \
	grid_csps_start.sh \
	grid_sps_resume.sh \
	grid_sps_start.sh

SUBDIRS = \
	pepnovo \
	mscluster \
	Clustal_w \
	ExecFramework \
	reports

.PHONY: subdirs $(SUBDIRS) clean all

all: subdirs $(STATICLIB) $(EXECUTABLES)
	cp pepnovo/PepNovo_bin PepNovo_bin
	cp mscluster/MsCluster/MsCluster_bin MsCluster_bin
	cp Clustal_w/clustalw clustalw

subdirs: $(STATICLIB) $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

install: all
	mkdir -p $(PREFIX)/bin
	mkdir -p $(PREFIX)/share/gnuplot/4.2/PostScript
	cp -f $(EXECUTABLES) PepNovo_bin MsCluster_bin $(SCRIPTS) $(PREFIX)/bin
	cp -f $(GNUPLOT_ROOT)/term/PostScript/prologue.ps $(PREFIX)/share/gnuplot/4.2/PostScript
ifneq (,$(shell which mcc 2>/dev/null))
	cp -f csps specnets $(PREFIX)/bin
endif

clean:
	rm -f *.o MetaSPS/*.o $(STATICLIB) $(EXECUTABLES) csps specnets PepNovo_bin MsCluster_bin pspectrum pcontig; cd pepnovo; make clean; cd ../mscluster; make clean; cd ../reports; make clean; cd ../ExecFramework;make clean;cd ..

distclean: clean
	$(MAKE) -C $(XERCES_ROOT) distclean
	$(MAKE) -C $(GNUPLOT_ROOT) distclean

.cpp.o: $(CPPLIB_SRC)
	$(CC) $(CPPFLAGS) $< -c -o $@

.o.o: $(CPPLIB_SRC)
	$(CC) $(CPPFLAGS) $<

$(STATICLIB): $(CPPLIB_OBJ)
	ar rcs $@ $^

$(DYNAMICLIB): $(CPPLIB_OBJ)
	$(CC) -shared  -Wl,-rpath,$(ROOT) -Wl,-soname,$@.1 -o $@.1.0 $^
	ln -s $@.1.0 $@.1
	ln -s $@.1 $@

test_cs: test_cs.cpp $(CPPLIB_OBJ)
	$(CC) $(CPPFLAGS) $^ -o $@

test: test.cpp $(CPPLIB_OBJ)
	$(CC) $(CPPFLAGS) $(LDFLAGS) $^ -o $@ -lxerces-c-3.1

# Linux
ifneq (,$(findstring linux,$(shell gcc -dumpmachine)))
convert: convert.o mzxml.o $(CPPLIB_OBJ)
	$(CC) $(LDFLAGS) $^ -o $@ -lxerces-c-3.1

MetaSPS/makeClustersFromAbinfo: MetaSPS/makeClustersFromAbinfo.o $(CPPLIB_OBJ)
	$(CC) $(LDFLAGS) $^ -o $@ -lxerces-c-3.1

mergeConvert: mergeConvert.o mzxml.o $(CPPLIB_OBJ)
	$(CC) $(LDFLAGS) $^ -o $@ -lxerces-c-3.1


# Windows
else
convert: convert.o mzxml.o $(CPPLIB_OBJ)
	$(CC) $(LDFLAGS) $^ -o $@ -lxerces-c-3-1

mergeConvert: mergeConvert.o mzxml.o $(CPPLIB_OBJ)
	$(CC) $(LDFLAGS) $^ -o $@ -lxerces-c-3-1
endif

filterpairs: filterpairs.o $(CPPLIB_OBJ)
	$(CC) $(LDFLAGS) $^ -o $@ -lpthread

starden: starden.o $(CPPLIB_OBJ)
	$(CC) $(LDFLAGS) $^ -o $@ -lpthread

filterstarpairs: filterstarpairs.o $(CPPLIB_OBJ)
	$(CC) $(LDFLAGS) $^ -o $@ -lpthread

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
	$(CC) $(LDFLAGS) $^ -o $@ -lpthread

homglue: homglue.o $(CPPLIB_OBJ)
	$(CC) $(LDFLAGS) $^ -o $@

# NRP-sequencing modules
cycseq: cycseq.o $(CPPLIB_OBJ)
	$(CC) $(LDFLAGS) $^ -o $@

#misc
pklbin_to_mgf: pklbin_to_mgf.cpp $(CPPLIB_OBJ)
	$(CC) $(CPPFLAGS) $^ -o $@

cosine_similarity: cosine_similarity_test.cpp cosine_similarity.o $(CPPLIB_OBJ)
	$(CC) $(LDFLAGS) -I$(BOOST_ROOT) $(CPPFLAGS) $^ -o $@


$(XERCES_ROOT)/Makefile:
	cd $(XERCES_ROOT); ./configure --disable-pretty-make --disable-network --disable-threads --without-curl --without-icu; cd $(ROOT)

$(PREFIX)/lib/libxerces-c.a: $(XERCES_ROOT)/Makefile
	$(MAKE) -C $(XERCES_ROOT)/src LDFLAGS="$(LDFLAGS)" CPPFLAGS="$(CPPFLAGS)"
	mv $(XERCES_ROOT)/src/.libs/libxerces-c.a $(PREFIX)/lib

# Linux
ifneq (,$(findstring linux,$(shell gcc -dumpmachine)))

$(PREFIX)/lib/libxerces-c-3.1.so: $(XERCES_ROOT)/Makefile
	$(MAKE) -C $(XERCES_ROOT)/src LDFLAGS="$(LDFLAGS)" CPPFLAGS="$(CPPFLAGS)"
	mv $(XERCES_ROOT)/src/.libs/libxerces-c-3.1.so $(PREFIX)/lib/libxerces-c-3.1.so


# Windows Cygwin
else
ifneq (,$(findstring cygwin,$(shell gcc -dumpmachine)))

$(PREFIX)/bin/cygxerces-c-3-1.dll: $(XERCES_ROOT)/Makefile
	$(MAKE) -C $(XERCES_ROOT)/src LDFLAGS="$(LDFLAGS)" CPPFLAGS="$(CPPFLAGS)"
	mv $(XERCES_ROOT)/src/.libs/cygxerces-c-3-1.dll $(PREFIX)/bin

endif
endif

# dependencies
depends : $(CPPLIB_SRC) spsstats.cpp dbg.cpp mzxml.cpp
	@$(CC) -MM $^ > depends

# Linux
ifneq (,$(findstring linux,$(shell gcc -dumpmachine)))
-include depends
endif
