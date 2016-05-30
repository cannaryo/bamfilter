# ========================================
# Makefile for ReALLEN : Bamfilter
# (c) 2016 Ryo Kanno
# ========================================

CXX	= g++
CFLAGS	= -O4 -Wall -I/usr/local/include -I/home/clc/test/bamtools-master/include
DEST	= /usr/local/bin
LDFLAGS	= -L/usr/local/lib -L/home/clc/test/bamtools-master/lib
LIBS	= -lm -lz -lbamtools
SOURCE	= main.cc bamfilter-main.cc sequence-splitter.cc
PROGRAM = bin/bamfilter

all:	$(PROGRAM)

$(PROGRAM):	$(SOURCE)
		$(CXX) $(CFLAGS) $(SOURCE) $(LDFLAGS) $(LIBS) -o $(PROGRAM)

clean:;	rm -f *.o *~ $(PROGRAM)

install:	$(PROGRAM)
		install -s $(PROGRAM) $(DEST)
