SHELL = /bin/sh
CXX ?= g++
CFLAGS += -ffunction-sections -fdata-sections -O3 -mtune=generic
LDFLAGS += -lm -static

OBJS := $(patsubst %.cpp,%.o,$(wildcard *.cpp))
ALL_MSCLUSTER_OBJS := $(patsubst %.cpp,%.o,$(wildcard ../MsCluster/*.cpp))
MSCLUSTER_OBJS := $(filter-out ../MsCluster/MsCluster_main.o,$(ALL_MSCLUSTER_OBJS))
ALL_MLLIB_OBJS := $(patsubst %.cpp,%.o,$(wildcard ../MLlib/*.cpp))
MLLIB_OBJS := $(filter-out ../MLlib/mllib.o,$(ALL_MLLIB_OBJS))
ALL_PEPNOVO_OBJS := $(patsubst %.cpp,%.o,$(wildcard ../PepNovo/*.cpp))
PEPNOVO_OBJS := $(filter-out ../PepNovo/PepNovo_main.o,$(ALL_PEPNOVO_OBJS))
COMMON_OBJS := $(patsubst %.cpp,%.o,$(wildcard ../Common/*.cpp))

EXE = MsLibrary_bin

.SUFFIXES: .cpp

all: $(EXE)

$(EXE): $(OBJS) $(COMMON_OBJS) $(MLLIB_OBJS) $(PEPNOVO_OBJS) $(MSCLUSTER_OBJS)
	$(CXX) $(CFLAGS) $(LDFLAGS) -o $(EXE) $(OBJS) $(COMMON_OBJS) $(MLLIB_OBJS) $(PEPNOVO_OBJS) $(MSCLUSTER_OBJS) $(LIBS)

.cpp.o:
	$(CXX) -c $(CFLAGS) -o $@ $<

clean-objs:
	@rm -f $(OBJS)

clean-exe:
	@rm -f $(EXE)

clean: clean-objs clean-exe

print_vars:
	echo "OBJS: " $(OBJS)
	echo "ALL_MLLIB_OBJS: " $(ALL_MLLIB_OBJS)
	echo "MLLIB_OBJS: " $(MLLIB_OBJS)
	echo "COMMIN_OBJS: " $(COMMON_OBJS)
	echo "PEPNOVO_OBJS: " $(PEPNOVO_OBJS)

