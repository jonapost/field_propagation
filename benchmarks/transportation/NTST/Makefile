CURRENTDIR:=`pwd`
CUSTOMFLAGS := -O2 -ftree-vectorize -ftree-vectorizer-verbose=2 -I $(G4INCLUDE) -I $(CURRENTDIR)/include
SOURCES := $(wildcard src/*.cc)
SOURCES += testNTST.cc
OBJECTS := $(SOURCES:.cc=.o)
GEANT4LIBFLAGS := `geant4-config --libs`

%.o: %.cc
	$(CXX) $(CUSTOMFLAGS) -c -o $@ $< 

testNTST: $(OBJECTS)
	$(CXX) $(OBJECTS) -o $@ $(GEANT4LIBFLAGS)


clean:
	rm $(OBJECTS)
	rm testNTST