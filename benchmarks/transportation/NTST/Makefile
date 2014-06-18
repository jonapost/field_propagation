CURRENTDIR:=`pwd`
CUSTOMFLAGS := -O2 -ftree-vectorize -ftree-vectorizer-verbose=2 -I $(G4INCLUDE) -I $(CURRENTDIR)/include
SOURCES := $(wildcard src/*.cc)
SOURCES += testNTST.cc
OBJECTS := $(SOURCES:.cc=.o)
GEANT4LIBFLAGS := `geant4-config --libs`

%.o: %.cc
	$(CXX) $(CUSTOMFLAGS) -c -o $@ $< 

testNTST: $(OBJECTS)
	$(CC) $(OBJECTS) -o NTSTbench $(GEANT4LIBFLAGS)


clean:
	rm $(OBJECTS)
	rm testNTST