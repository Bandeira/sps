SHELL = /bin/sh
CXX ?= g++
CFLAGS += -O3

OBJS := $(patsubst %.cpp,%.o,$(wildcard *.cpp))

.SUFFIXES: .cpp 

all:	$(OBJS)

.cpp.o:	$(OBJS) 
	$(CXX) -c $(CFLAGS) -o $@ $<

clean-objs:
	@rm -f $(OBJS) 

clean-exe:
	@rm -f $(EXE) 

clean: clean-objs clean-exe
