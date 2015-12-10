################################################################################
# Definitions
################################################################################
BASE_DIR := ../../..
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
# Compile rules - general
################################################################################
include $(BASE_DIR)/Makefiles/Makefile.general

################################################################################
# Source files
################################################################################
DATACOMPARE_SRC = \
	DataComparer.cpp \
	DataCompareInterface.cpp

DATACOMPARE_HDR = \
	DataComparer.h \
	DataCompareInterface.h

CPPLIB_SRC = \
  $(COMMON_FILES)/ExecFramework/CommandLineParser.cpp \
  $(COMMON_FILES)/ExecFramework/ParameterList.cpp  \
  $(COMMON_FILES)/abruijn.cpp \
	$(COMMON_FILES)/SpectrumPair.cpp \
	$(COMMON_FILES)/SpectrumPairSet.cpp \

CPPLIB_HDR = \
  $(COMMON_FILES)/ExecFramework/CommandLineParser.h \
  $(COMMON_FILES)/ExecFramework/ParameterList.h  \
  $(COMMON_FILES)/abruijn.h \
  $(COMMON_FILES)/ms1.h \
	$(COMMON_FILES)/SpectrumPair.h \
	$(COMMON_FILES)/SpectrumPairSet.h \

################################################################################
# Compile rules - specific
################################################################################
DATACOMPARE_OBJ=$(DATACOMPARE_SRC:.cpp=.o)

CPPLIB_OBJ=$(CPPLIB_SRC:.cpp=.o)

dataCompare: $(PREFIX)/$(SUFFIX)/libxerces-c.a dataCompare.o $(DATACOMPARE_OBJ) $(CPPLIB_OBJ) $(LIBCORE_SRC)  $(LIBSPS_SRC) $(LIBSPSREPORTS_SRC) # $(PREFIX)/$(SUFFIX)/gnuplot.dll
	$(CXX) $(CPPFLAGS) $(LDFLAGS) dataCompare.o $(DATACOMPARE_OBJ) $(CPPLIB_OBJ) $(LIBSPS_SRC) $(LIBCORE_SRC) $(LIBSPSREPORTS_SRC) -o $@ -Wall -lxerces-c $(LIBS) $(PWIZ_LIBS) -lpthread

################################################################################
# Targets
################################################################################
EXECUTABLES = \
	dataCompare

all: $(EXECUTABLES)
	$(STRIP)  $(EXECUTABLES)

################################################################################
# Clean
################################################################################
clean:
	rm -f *.o $(EXECUTABLES) *.gch *~ *.bak

