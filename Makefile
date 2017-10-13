#path var + dir commands                                                                                                         
path=$(PWD)
MKDIR_BIN = mkdir -p $(path)/bin
MKDIR_LIB = mkdir -p $(path)/lib
MKDIR_OUTTXT = mkdir -p $(path)/out

ROOT = `root-config --cflags --glibs`
CXX = g++
CXXFLAGS = -Wall -O2 -Wextra -Wno-unused-local-typedefs  -Werror -Wno-deprecated-declarations -std=c++11
ifeq "$(GCCVERSION)" "1"
  CXXFLAGS += -Wno-error=misleading-indentation
endif

FASTJET = `/afs/cern.ch/work/c/cmcginn/private/Generators/fastjet-install/bin/fastjet-config --cxxflags --libs --plugins`

#fortran var
FC := gfortran
FFLAGS := -g -static

#programs to make

all: mkdirBin mkdirLib mkdirOutTxt pythia-6.4.28.o main63.o main63.exe xsection63.o xsection63.exe parsePylist.exe clusterGenToJet.exe plotQGFraction.exe plotDijetFlavVsPthat.exe extractWeightAndErr.exe

mkdirBin:
	$(MKDIR_BIN)

mkdirLib:
	$(MKDIR_LIB)

mkdirOutTxt:
	$(MKDIR_OUTTXT)

pythia-6.4.28.o:
	$(FC) -c -o $(path)/lib/pythia-6.4.28.o $(path)/src/pythia-6.4.28.f

main63.o: pythia-6.4.28.o
	$(FC) -c -o $(path)/lib/main63.o $(path)/src/main63.f

main63.exe: $(path)/lib/pythia-6.4.28.o $(path)/lib/main63.o
	$(FC) -o $(path)/bin/$@ $^ 

xsection63.o: pythia-6.4.28.o
	$(FC) -c -o $(path)/lib/xsection63.o $(path)/src/xsection63.f

xsection63.exe: $(path)/lib/pythia-6.4.28.o $(path)/lib/xsection63.o
	$(FC) -o $(path)/bin/$@ $^ 

parsePylist.exe: src/parsePylist.C
	$(CXX) $(CXXFLAGS) $(ROOT) -o bin/parsePylist.exe src/parsePylist.C

clusterGenToJet.exe: src/clusterGenToJet.C
	$(CXX) $(CXXFLAGS) $(ROOT) $(FASTJET) -I$(PWD) -o bin/clusterGenToJet.exe src/clusterGenToJet.C

plotQGFraction.exe: src/plotQGFraction.C
	$(CXX) $(CXXFLAGS) $(ROOT) $(FASTJET) -I$(PWD) -o bin/plotQGFraction.exe src/plotQGFraction.C

plotDijetFlavVsPthat.exe: src/plotDijetFlavVsPthat.C
	$(CXX) $(CXXFLAGS) $(ROOT) $(FASTJET) -I$(PWD) -o bin/plotDijetFlavVsPthat.exe src/plotDijetFlavVsPthat.C

extractWeightAndErr.exe: src/extractWeightAndErr.C
	$(CXX) $(CXXFLAGS) $(ROOT) $(FASTJET) -I$(PWD) -o bin/extractWeightAndErr.exe src/extractWeightAndErr.C

clean:
	rm -f *~
	rm -f \#*.*#
	rm -f $(path)/src/#*.*#
	rm -f $(path)/src/*~
	rm -f $(path)/include/#*.*#
	rm -f $(path)/include/*~
	rm -f $(path)/lib/py*.o 
	rm -f $(path)/lib/main*.o 
	rm -f $(path)/bin/*.exe

	rmdir bin
	rmdir lib

.PHONY: all

