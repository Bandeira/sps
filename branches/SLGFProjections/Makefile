################################################################################
# Definitions
################################################################################
BASE_DIR := ..
COMMON_FILES := $(BASE_DIR)/branches/SLGFProjections
#-------------------------------------------------------------------------------
include $(BASE_DIR)/Makefiles/Makefile.defs

################################################################################
# Conditions
################################################################################
include $(BASE_DIR)/Makefiles/Makefile.cond

################################################################################
# External Libraries
################################################################################
include $(BASE_DIR)/Makefiles/Makefile.libs

################################################################################
# Gnuplot
################################################################################
include $(BASE_DIR)/Makefiles/Makefile.gnuplot

################################################################################
# Compile rules - general
################################################################################
include $(BASE_DIR)/Makefiles/Makefile.general

################################################################################
# Source files - CPPLIB
################################################################################
CPPLIB_SRC = \
	abruijn.cpp \
	alignment_modmut.cpp \
	AlignmentPenaltyBased.cpp \
	base64.cpp \
	batch.cpp \
	clusters.cpp \
	cosine_similarity.cpp \
	denovo.cpp \
	FdrPeptide.cpp \
	filters.cpp \
	graph.cpp \
	inspect_parse.cpp \
	MappedContig.cpp \
	MappedContigStatTable.cpp \
	MappedSpecnets.cpp \
	MappedSpectrum.cpp \
	MappedSPSStatTable.cpp \
	ms1.cpp \
	msn.cpp \
	OutputTable.cpp \
	projectionutils.cpp \
	PenaltyMatrix.cpp \
	PepnovoTags.cpp \
	PeptideSpectrumMatchNetwork.cpp \
	SetMerger.cpp \
	SpectralPairs.cpp \
	spectrum_window_filter.cpp \
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
	mzxml.cpp \
#	dekel_align.cpp \
#	hash.cpp \

################################################################################
# Source files - META_SPS
################################################################################
META_SPS_SRC = MetaSPS/prm_alignment.cpp \
		MetaSPS/PRMSpecEdge.cpp \
		MetaSPS/PairedSpecSet.cpp \
		MetaSPS/CombineContigs.cpp \
		MetaSPS/Contig.cpp \
		MetaSPS/ContigNetwork.cpp

################################################################################
# Libraries - specnets
################################################################################
STATICLIB = libspecnets.a

DYNAMICLIB = libspecnets.so

################################################################################
# Executables & scripts
################################################################################
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

################################################################################
# Compile rules - install
################################################################################
install: all
	mkdir -p $(PREFIX)/bin
	mkdir -p $(PREFIX)/share/gnuplot/4.2/PostScript
	cp -f $(EXECUTABLES) $(SCRIPTS) $(PREFIX)/bin
	cp -f $(GNUPLOT_ROOT)/term/PostScript/prologue.ps $(PREFIX)/share/gnuplot/4.2/PostScript
ifneq (,$(shell which mcc 2>/dev/null))
	cp -f csps specnets $(PREFIX)/bin
endif

################################################################################
# Compile rules - cleanup
################################################################################
clean:
	rm -f *.o *.a MetaSPS/*.o PepNovo_bin GenoMS.jar MSGFDB.jar convert csps specnets MsCluster_bin pspectrum pcontig depends *.gch *~ *.bak; rm -rf DBs_GenoMS; cd libSPS; make clean; cd .. ; cd libCore; make clean; cd ../reports; make clean; cd ../ExecFramework;make clean;cd ..

distclean: clean
	$(MAKE) -C $(XERCES_ROOT) distclean
	$(MAKE) -C $(GNUPLOT_ROOT) distclean
	$(MAKE) -C $(PEPNOVO_DIR) clean
	$(MAKE) -C $(MSCLUSTER_DIR) clean
	$(MAKE) -C $(CLUSTAL_W_DIR) clean
	$(MAKE) -C $(PDFHARU_ROOT) distclean
	$(MAKE) -C $(LIBZ_ROOT) distclean
	$(MAKE) -C $(LIBPNG_ROOT) distclean
	rm -f $(PREFIX)/$(SUFFIX)/*.a

################################################################################
# Compile rules - specific
################################################################################
CPPLIB_OBJ=$(CPPLIB_SRC:.cpp=.o)

META_SPS_OBJ=$(META_SPS_SRC:.cpp=.o)

.PHONY: subdirs $(SUBDIRS) clean all

all: $(SUBDIRS) clustal_wTarget pepnovoTarget msclusterTarget ExternalLibs gnuplotExe libSPSTarget libCoreTarget $(STATICLIB) ReportsTarget ExecTarget  utilsTarget $(EXECUTABLES) GenoMS
	cp $(CLUSTAL_W_DIR)/clustalw .
	cp $(MSCLUSTER_DIR)/MsCluster/MsCluster_bin .
	cp $(PEPNOVO_DIR)/PepNovo_bin .
	cp $(GENOMS_DIR)/GenoMS.jar .
	cp $(GENOMS_DIR)/MSGFDB.jar .
	cp $(EXEC_DIR)/main_specnets .
	cp $(EXEC_DIR)/main_execmodule .

#subdirs: $(LIBCORE_SRC) $(STATICLIB) $(SUBDIRS)

SUBDIRS =

$(SUBDIRS):
	$(MAKE) -C $@

################################################################################
# Compile rules - libraries
################################################################################
ExecTarget: $(STATICLIB)
	$(MAKE) -C $(EXEC_DIR)

ReportsTarget: $(STATICLIB)
	$(MAKE) -C $(REPORTS_DIR)

libSPSTarget:
	$(MAKE) -C $(LIBSPS_DIR)

libCoreTarget:
	$(MAKE) -C $(LIBCORE_DIR)


$(STATICLIB): $(CPPLIB_OBJ) $(META_SPS_OBJ) $(LIBCORE_SRC)
	$(AR) rc $@ $^
	$(RANLIB) $@

$(LIBSPS_SRC):
	$(MAKE) -C $(LIBSPS_DIR)

$(LIBCORE_SRC):
	$(MAKE) -C $(LIBCORE_DIR)


$(DYNAMICLIB): $(CPPLIB_OBJ)
	$(CXX) -shared  -Wl,-rpath,$(ROOT) -Wl,-soname,$@.1 -o $@.1.0 $^
	ln -s $@.1.0 $@.1
	ln -s $@.1 $@

################################################################################
# Compile rules - external tools
################################################################################
pepnovoTarget:
	$(MAKE) -C $(PEPNOVO_DIR)

msclusterTarget:
	$(MAKE) -C $(MSCLUSTER_DIR)

clustal_wTarget:
	$(MAKE) -C $(CLUSTAL_W_DIR)

################################################################################
# Compile rules - utils directory
################################################################################
utilsTarget:
	$(MAKE) -C $(UTILS_DIR)

################################################################################
# Compile rules - test tools
################################################################################
test_cs: test_cs.cpp $(CPPLIB_OBJ) $(LIBCORE_SRC)
	$(CXX) $(CPPFLAGS) $^ -o $@ $(LIBS)

test: test.cpp $(CPPLIB_OBJ) $(LIBCORE_SRC)
	$(CXX) $(CPPFLAGS) $(LDFLAGS) $^ -o $@ -lxerces-c $(LIBS)

# Linux
ifneq (,$(findstring linux,$(shell gcc -dumpmachine)))
convert: convert.o mzxml.o $(CPPLIB_OBJ) $(LIBCORE_SRC) $(LIBSPS_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@ -lxerces-c $(LIBS)

MetaSPS/makeClustersFromAbinfo: MetaSPS/makeClustersFromAbinfo.o $(CPPLIB_OBJ) $(LIBCORE_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@ -lxerces-c $(LIBS)

mergeConvert: mergeConvert.o mzxml.o $(CPPLIB_OBJ) $(LIBCORE_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@ -lxerces-c $(LIBS)

# Windows
else
#convert: convert.o libCore/spectrum.o libCore/IsoEnvelope.o mzxml.o Base64.o alignment_scoring.o libCore/label.o libCore/Logger.o ms1.o libCore/utils.o libCore/aminoacid.o libCore/PeptideSpectrumMatch.o libCore/PeptideSpectrumMatchSet.o libCore/DelimitedTextReader.o libCore/SpecSet.o libCore/spectrum_scoring.o libCore/inputParams.o libCore/mzrange.o denovo.o clusters.o PepnovoTags.o abruijn.o graph.o SpectrumPair.o SetMerger.o SpectrumPairSet.o $(PREFIX)/bin/mingw/libxerces-c.a
#	$(CC) $(LDFLAGS) $^ -o $@ -l ws2_32 -lxerces-c
pklbin_to_mgf: pklbin_to_mgf.o $(CPPLIB_OBJ) $(LIBCORE_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@ -lxerces-c $(LIBS)

convert: convert.o mzxml.o $(CPPLIB_OBJ) $(LIBCORE_SRC)  $(LIBSPS_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@ -lxerces-c $(LIBS)

mergeConvert: mergeConvert.o mzxml.o $(CPPLIB_OBJ) $(LIBCORE_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@ -lxerces-c $(LIBS)
endif

filterpairs: filterpairs.o $(CPPLIB_OBJ) $(LIBCORE_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@ -lpthread $(LIBS)

starden: starden.o $(CPPLIB_OBJ) $(LIBCORE_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@ -lpthread $(LIBS)

filterstarpairs: filterstarpairs.o $(CPPLIB_OBJ) $(LIBCORE_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@ -lpthread $(LIBS)

tagsearch: tagsearch.o $(CPPLIB_OBJ) $(LIBCORE_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIBS)

pathproj: pathproj.o $(CPPLIB_OBJ) $(LIBCORE_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIBS)

findconnected: findconnected.o $(CPPLIB_OBJ) $(LIBCORE_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIBS)

msmod: msmod.o FloatMods.o FloatMods.h $(LIBCORE_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIBS)

#-------------------------------------------------------------------------------
# SPS modules
masab: masab.o $(CPPLIB_OBJ) $(LIBCORE_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIBS)

#-------------------------------------------------------------------------------
# CSPS modules
protid: protid.o $(CPPLIB_OBJ) $(LIBCORE_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIBS)
matchma: matchma.o $(CPPLIB_OBJ) $(LIBCORE_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIBS)

matchme: matchme.o $(CPPLIB_OBJ) $(LIBCORE_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIBS)

mergefilter: mergefilter.o $(CPPLIB_OBJ) $(LIBCORE_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@ -lpthread $(LIBS)

homglue: homglue.o $(CPPLIB_OBJ) $(LIBCORE_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIBS)

#-------------------------------------------------------------------------------
# NRP-sequencing modules
cycseq: cycseq.o $(CPPLIB_OBJ) $(LIBCORE_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIBS)

#-------------------------------------------------------------------------------
#misc
cosine_similarity: cosine_similarity_test.o cosine_similarity.o $(CPPLIB_OBJ) $(LIBCORE_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@ -lxerces-c $(LIBS)


################################################################################
# Compile rules - GenoMS
################################################################################
./DBs_GenoMS:
	mkdir ./DBs_GenoMS
	cp $(GENOMS_DIR)/DBs_GenoMS/* ./DBs_GenoMS/

GenoMS: ./DBs_GenoMS

################################################################################
# dependencies
################################################################################
depends : $(CPPLIB_SRC) spsstats.cpp dbg.cpp mzxml.cpp
	@$(CXX) $(CPPFLAGS) -MM $^ > depends

# Linux
ifneq (,$(findstring linux,$(shell gcc -dumpmachine)))
-include depends
endif

################################################################################
