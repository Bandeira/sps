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
# Libraries
################################################################################
LIBRARIES = -lpthread -lspecnets -lflr -lxerces-c-3-1 -lcore

################################################################################
# Source files - CPPLIB
################################################################################
CPPLIB_SRC = \
  $(EXEC_DIR)/CommandLineParser.cpp \
  $(EXEC_DIR)/ExecBase.cpp \
  ExecModuleFactoryLP.cpp \
  ExecPairSpectra.cpp \
  ExecFlr.cpp \
  $(EXEC_DIR)/FileUtils.cpp \
  $(EXEC_DIR)/Filters.cpp \
  $(EXEC_DIR)/ParameterList.cpp \
  $(EXEC_DIR)/ParallelExecution.cpp \
  $(EXEC_DIR)/ParallelPbsExecution.cpp \
  $(EXEC_DIR)/ParallelSgeExecution.cpp \
  $(EXEC_DIR)/ParallelThreadedExecution.cpp \
  $(EXEC_DIR)/SgeGridMonitor.cpp
  
CPPLIB_HDR = \
  $(EXEC_DIR)/CommandLineParser.h \
  $(EXEC_DIR)/ExecBase.h \
  ExecModuleFactoryLP.h \
  ExecPairSpectra.h \
  ExecFlr.h \
  $(EXEC_DIR)/FileUtils.h \
  $(EXEC_DIR)/Filters.h \
  $(EXEC_DIR)/ParameterList.h \
  $(EXEC_DIR)/ParallelExecution.h \
  $(EXEC_DIR)/ParallelPbsExecution.h \
  $(EXEC_DIR)/ParallelSgeExecution.h \
  $(EXEC_DIR)/ParallelThreadedExecution.h \
  $(EXEC_DIR)/SgeGridMonitor.h

################################################################################
# Compile rules: Specific
################################################################################
CPPLIB_OBJ=$(CPPLIB_SRC:.cpp=.o)

main_execmodule: main_execmoduleLP.o $(CPPLIB_OBJ) $(CPPLIB_HDR)
	$(CXX) $(CPPFLAGS) $(LDFLAGS) $^ -o $@ $(LIBRARIES) $(LIBS)

main_flr: main_flr.o $(CPPLIB_OBJ) $(CPPLIB_HDR)
	$(CXX) $(CPPFLAGS) $(LDFLAGS) $^ -o $@ $(LIBRARIES) $(LIBS)

################################################################################
# Targets
################################################################################
EXECUTABLES = \
	main_execmodule \
	main_flr

.PHONY: subdirs $(SUBDIRS) clean all

all: $(EXECUTABLES) main_execmodule

################################################################################
# Compile rules: Cleanup
################################################################################
clean:
	rm -f *.o $(EXECUTABLES) 

################################################################################
# Compile rules: Dependencies
################################################################################
# dependencies
depends : $(CPPLIB_SRC) $(CPPLIB_HDR)
	@$(CXX) $(CPPFLAGS) -MM $^ > depends

# Linux
ifneq (,$(findstring linux,$(shell gcc -dumpmachine)))
-include depends
endif

################################################################################
