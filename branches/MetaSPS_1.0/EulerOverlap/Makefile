
CC = g++

ifeq ($(type),debug)
  CPPFLAGS += -g -ffunction-sections -fdata-sections
else
  CPPFLAGS += -O3 -ffunction-sections -fdata-sections
endif

LINKER = -lpthread

CPPLIB_OBJ = ../abruijn.o\
	../alignment_modmut.o\
	../alignment_scoring.o\
	../aminoacid.o\
	../base64_.o\
	../batch.o\
	../clusters.o\
	../db_fasta.o\
	../dekel_align.o\
	../denovo.o\
	../filters.o\
	../graph.o\
	../hash.o\
	../inputParams.o\
	../inspect_parse.o\
	../label.o\
	../ms1.o\
	../msn.o\
	../mzrange.o\
	../PepnovoTags.o\
	../sequencing.o\
	../setmerger.o\
	../specnets_statistics.o\
	../spectral_pairs.o\
	../spectrum.o\
	../spectrum_scoring.o\
	../tags.o\
	../utils.o\
	../ion.o\
	OverlapOut.o\
	PruneEnvelopes.o\
	CombineContigs.o\
	parameter_names.o\
	SpsMergeHandler.o\
	SpsStatsHandler.o

CPPLIB_SRC = ../abruijn.cpp ../abruijn.h\
	../alignment_modmut.cpp ../alignment_modmut.h\
	../alignment_scoring.cpp ../alignment_scoring.h\
	../aminoacid.cpp ../aminoacid.h\
	../base64_.c ../base64_.h\
	../batch.cpp ../batch.h\
	../clusters.cpp ../clusters.h\
	../db_fasta.cpp ../db_fasta.h\
	../dekel_align.cpp ../dekel_align.h\
	../denovo.cpp ../denovo.h\
	../filters.cpp ../filters.h\
	../graph.cpp ../graph.h\
	../hash.cpp ../hash.h\
	../inputParams.cpp ../inputParams.h\
	../inspect_parse.cpp ../inspect_parse.h\
	../label.cpp ../label.h\
	../ms1.cpp ../ms1.h\
	../msn.cpp ../msn.h\
	../mzrange.cpp ../mzrange.h \
	../PepnovoTags.cpp ../PepnovoTags.h\
	../sequencing.cpp ../sequencing.h\
	../setmerger.cpp ../setmerger.h\
	../specnets_statistics.cpp ../specnets_statistics.h\
	../spectral_pairs.cpp ../spectral_pairs.h\
	../spectrum.cpp ../spectrum.h\
	../spectrum_scoring.cpp ../spectrum_scoring.h\
	../tags.cpp ../tags.h\
	../twovalues.h\
	../utils.cpp ../utils.h\
	../uncertain.h\
	../range.h\
	../ion.cpp ../ion.h\
	OverlapOut.cpp OverlapOut.h\
	PruneEnvelopes.cpp PruneEnvelopes.h\
	CombineContigs.cpp CombineContigs.h\
	parameter_names.cpp parameter_names.h\
	SpsMergeHandler.cpp SpsMergeHandler.h\
	SpsStatsHandler.cpp SpsStatsHandler.h

all: $(CPPLIB_SRC) spsoverlap spsmerge spsstats prunenv rankfilter

clean:
	rm -f *.o spsoverlap spsmerge spsstats prunenv getDBSearchStats filterMatchma mergeContigComponentSets filterOMSSA filterInspect rankfilter intensityFilter annotSpecs

.o.o: $(CPPLIB_SRC)
	$(CC) $(CPPFLAGS) $<

spsoverlap: spsoverlap.cpp $(CPPLIB_OBJ)
	$(CC) $(CPPFLAGS) $^ -o $@ $(LINKER)

prunenv: prunenv.cpp $(CPPLIB_OBJ)
	$(CC) $(CPPFLAGS) $^ -o $@ $(LINKER)

getDBSearchStats: getDBSearchStats.cpp $(CPPLIB_OBJ)
	$(CC) $(CPPFLAGS) $^ -o $@ $(LINKER)

filterMatchma: filterMatchma.cpp $(CPPLIB_OBJ)
	$(CC) $(CPPFLAGS) $^ -o $@ $(LINKER)

mergeContigComponentSets: mergeContigComponentSets.cpp $(CPPLIB_OBJ)
	$(CC) $(CPPFLAGS) $^ -o $@ $(LINKER)

filterOMSSA: filterOMSSA.cpp $(CPPLIB_OBJ)
	$(CC) $(CPPFLAGS) $^ -o $@ $(LINKER)

filterInspect: filterInspect.cpp $(CPPLIB_OBJ)
	$(CC) $(CPPFLAGS) $^ -o $@ $(LINKER)

spsmerge: spsmerge.cpp $(CPPLIB_OBJ)
	$(CC) $(CPPFLAGS) $^ -o $@ $(LINKER)

spsstats: spsstats.cpp $(CPPLIB_OBJ)
	$(CC) $(CPPFLAGS) $^ -o $@ $(LINKER)

mergeSpecs: mergeSpecs.cpp $(CPPLIB_OBJ)
	$(CC) $(CPPFLAGS) $^ -o $@ $(LINKER)

rankfilter: rankfilter.cpp $(CPPLIB_OBJ)
	$(CC) $(CPPFLAGS) $^ -o $@ $(LINKER)

intensityFilter: intensityFilter.cpp $(CPPLIB_OBJ)
	$(CC) $(CPPFLAGS) $^ -o $@ $(LINKER)

annotSpecs: annotSpecs.cpp $(CPPLIB_OBJ)
	$(CC) $(CPPFLAGS) $^ -o $@ $(LINKER)

