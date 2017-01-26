# ========================================
# Makefile for ReALLEN : Bamfilter
# (c) 2016 Ryo Kanno
# ========================================

ROOTDIR = $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
CXX	= g++
CFLAGS	= -O4 -Wall -I/usr/local/include -I$(ROOTDIR)/include/bamtools/include -I$(ROOTDIR)/include/cmdline
DEST	= /usr/local/bin
LDFLAGS	= -L/usr/local/lib -L$(ROOTDIR)/include/bamtools/lib
LIBS	= -lm -lz -lbamtools
SOURCE	= main.cc bamfilter-main.cc sequence-splitter.cc sequence-evaluator.cc
PROGRAM = bin/bamfilter

all:	$(PROGRAM)

$(PROGRAM):	$(SOURCE)
		$(CXX) $(CFLAGS) $(SOURCE) $(LDFLAGS) $(LIBS) -o $(PROGRAM)

clean:;	rm -f *.o *~ $(PROGRAM)

install:	$(PROGRAM)
		install -s $(PROGRAM) $(DEST)
