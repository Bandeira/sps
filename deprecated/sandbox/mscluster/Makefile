SHELL = /bin/sh
CC = g++
#CFLAGS = -O3 
LDFLAGS = -lm

OBJS = cluster_main.o AnnotatedSpecturm.o auxfun.o base64.o BasicDataStructs.o Config.o\
    ConversionTables.o DeNovoDp.o DeNovoSolutions.o DiscretePeakModel.o EdgeModel.o FileManagement.o\
    Fragmentation.o FragmentSelection.o FragProbTable.o Isotopes.o MassRangeMap.o ME_REG.o\
    ME_Regression_DataSet.o ME_Regression_Model.o Model.o MultiPath.o MZXML_parsing.o\
    PrmGraph.o PTM_stuff.o RegionalPepNovoModel.o RegularRankModel.o AminoAcidProbs.o\
    Spectrum.o PMCSQS.o PMCSQS_model.o QCBasicSpecReader.o QuickClusteringSpectra.o\
        QCoutputter.o PMC_rank.o RankBoost.o RankBoostSample.o AdvancedScoreModel.o AdvancedScoreModel_fill.o\
	AdvancedScoreModel_frag.o AdvancedScoreModel_regional.o PeptideComp.o PartitionModel.o PeakRankModel.o\
	PeakRankPartialDenovo.o PeakRank_combined.o PeakRank_advanced.o PeakRank_simple.o PeakRank_combined_dnv.o\
	DeNovoRankTrain.o DeNovoRankScore.o DeNovoPartModel.o ReScoreDB.o QCStatics.o QCDotProdSim.o QuickClustering.o QCDAT.o\
	CumulativeSeqProb.o

HDRS = AnnotatedSpectrum.h auxfun.h base64.h BasicDataStructs.h Config.h ConversionTables.h\
       DeNovoDp.h DeNovoSolutions.h DiscretePeakModel.h EdgeModel.h FileManagement.h Fragmentation.h\
       FragmentSelection.h FragProbTable.h includes.h Isotopes.h ME_REG.h Model.h\
       PrmGraph.h RegionalPepNovoModel.h RegularRankModel.h AminoAcidProbs.h CumulativeSeqProb.h\
       Spectrum.h TagModel.h PMCSQS.h RankBoost.h AdvancedScoreModel.h PeptideComp.h PeakRankModel.h DeNovoRankScore.h

EXE = MSCluster_bin

.SUFFIXES: .cpp 

all: $(EXE)

$(EXE): $(OBJS)
	$(CC) $(CPPFLAGS) $(LDFLAGS) -o $(EXE) $(OBJS)

.cpp.o: 
	$(CC) -c $(CPPFLAGS) -o $@ $<

clean-objs:
	@rm -f $(OBJS) 

clean-exe:
	@rm -f $(EXE) 

clean: clean-objs clean-exe
