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
# Compile rules - general
################################################################################
include $(BASE_DIR)/Makefiles/Makefile.general

################################################################################
# Source files - CLUSTER
################################################################################
CLUSTER_SRC = \
  ClusterInfo.cpp \
  ClusterInfoInterface.cpp \
  $(REPORT_COMMON_FILES)/ClusterData.cpp \
  $(COMMON_FILES)/ExecFramework/CommandLineParser.cpp \
  $(COMMON_FILES)/ExecFramework/ParameterList.cpp \
  $(COMMON_FILES)/Specific.cpp

CLUSTER_HDR = \
  ClusterInfo.h \
  ClusterInfoInterface.h \
  $(REPORT_COMMON_FILES)/ClusterData.h \
  $(COMMON_FILES)/ExecFramework/CommandLineParser.h \
  $(COMMON_FILES)/ExecFramework/ParameterList.h \
  $(COMMON_FILES)/Specific.h

################################################################################
# Targets
################################################################################
EXECUTABLES = \
	clusterinfo

all: $(EXECUTABLES)

################################################################################
# Compile rules - specific
################################################################################
CLUSTER_OBJ=$(CLUSTER_SRC:.cpp=.o)

clusterinfo: $(CLUSTER_OBJ) $(LIBCORE_SRC) $(LIBSPSREPORTS_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@ -Wall $(LIBS) $(PWIZ_LIBS) -lxerces-c -pthread

################################################################################
# Compile rules - cleanup
################################################################################
clean:
	rm -f *.o $(EXECUTABLES) depends *.gch *~ *.bak

################################################################################
