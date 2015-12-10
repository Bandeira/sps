SHELL = /bin/sh
CXX ?= g++
CFLAGS += -O3
LDFLAGS += -lm -static

OBJS := $(patsubst %.cpp,%.o,$(wildcard *.cpp))
COMMON_OBJS := $(patsubst %.cpp,%.o,$(wildcard ../Common/*.cpp))

EXE = mllib

.SUFFIXES: .cpp 

all: $(EXE)

$(EXE): $(OBJS) $(COMMON_OBJS)
	$(CXX) $(CFLAGS) $(LDFLAGS) -o $(EXE) $(OBJS) $(COMMON_OBJS) $(LIBS)

.cpp.o: 
	$(CXX) -c $(CFLAGS) -o $@ $<

clean-objs:
	@rm -f $(OBJS) 

clean-exe:
	@rm -f $(EXE) 

clean: clean-objs clean-exe
