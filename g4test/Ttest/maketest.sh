g++ -DINLINEDUMBSTEPPERS -o testPropagateMagField -g -O2 -ftree-vectorize -ftree-vectorizer-verbose=0  -I ${G4INCLUDE} testPropagateMagField.cc `geant4-config --libs` 
