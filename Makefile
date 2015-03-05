# Use the following on IRIX 6.5 and Linux 2.4

LIBS = $(shell root-config --glibs)
OPTCOMP = $(shell root-config --cflags)

#SPECIALFLAGS= --exceptions
#### TO TURN ON PDF WEIGHTS UNCOMMENT
#SPECIALFLAGS= --exceptions -D__USE_PDFS__
#SPECIALFLAGS= --exceptions -D__USE_PDFS_RESBOS__
SPECIALFLAGS= --exceptions -D__USE_PDFS_RESBOS__ -D__USE_PDFS__

CFLAGS = $(SPECIALFLAGS) -I- -I../ -I.
CFLAGS = $(SPECIALFLAGS) -iquote -I../ -I. -IAiUtil/
#LFLAGS = $(SPECIALFLAGS) -L../../lib/$(SRT_SUBDIR)/

RCXX=$(CFLAGS) $(ROOTCFLAGS) -ggdb
ROOTCINT=rootcint

#CC = KCC -n32 --exceptions --thread_safe -O $(OPTCOMP)
#CC = KCC +K0 --exceptions --thread_safe -O $(OPTCOMP)
CC = g++ $(RCXX) $(OPTCOMP) 

#all: tupleMaker tupleMaker2 tupleMaker3
all: tupleMaker3 tupleMaker_DYRES makeAiProfile


test: test.o Output.o
	$(CC) test.o Output.o $(LIBS) -o test -lEG

makeAiProfile: makeAiProfile.o AiMoments.o TLVUtils.o
	$(CC) $^ $(LIBS) -o $@ -lEG

tupleMaker_DYRES: tupleMaker_DYRES.o Output.o
	$(CC) tupleMaker_DYRES.o Output.o $(LIBS) -o tupleMaker_DYRES -lEG

tupleMaker3: tupleMaker3.o Output.o TLVUtils.o
	$(CC) $^ $(LIBS) -o tupleMaker3 -lEG

tupleMaker2: tupleMaker2.o Output.o
	$(CC) tupleMaker2.o Output.o $(LIBS) -o tupleMaker2 -lEG

tupleMaker: tupleMaker.o Output.o
	$(CC) tupleMaker.o Output.o $(LIBS) -o tupleMaker -lEG

Output.o: Output.cpp Output.hpp
	$(CC) -c Output.cpp -o Output.o

TLVUtils.o: AiUtil/TLVUtils.cxx AiUtil/TLVUtils.h
	$(CC) -c $< -o $@

AiMoments.o: AiUtil/AiMoments.cxx AiUtil/AiMoments.h AiUtil/TLVUtils.h
	$(CC) -c $< -o $@

tupleMaker.o: tupleMaker.cpp
	$(CC) -c tupleMaker.cpp -o tupleMaker.o

tupleMaker2.o: tupleMaker2.cpp
	$(CC) -c tupleMaker2.cpp -o tupleMaker2.o

tupleMaker3.o: tupleMaker3.cxx KinFile.C AiUtil/TLVUtils.h
	$(CC) -c tupleMaker3.cxx -o tupleMaker3.o

tupleMaker_DYRES.o: tupleMaker_DYRES.cxx
	$(CC) -c tupleMaker_DYRES.cxx -o tupleMaker_DYRES.o

makeAiProfile.o: makeAiProfile.cpp AiUtil/AiMoments.h ReadResbosROOT.C
	$(CC)  -c $< -o $@

test.o: test.cxx
	$(CC) -c test.cxx -o test.o

KinFileDict.cxx: KinFile.C
	$(ROOTCINT) -f $@ -c $^

clean:
	\rm -f *.o
	\rm -f *~
	\rm -fr test
	\rm -fr tupleMaker
	\rm -fr tupleMaker2
	\rm -fr tupleMaker3
	\rm -fr tupleMaker_DYRES
	\rm -fr makeAiProfile
