#ifndef ChordFinder_HH
#define ChordFinder_HH

#include "G4FieldTrack.hh"
#include <iomanip>
#include "G4SystemOfUnits.hh"

template <class T_Field, class T_Equation, class T_Driver>
class ChordFinder{

public:
    ChordFinder():
            fDeltaChord(0.25 * mm  ),    //   Internal parameters
            fIntgrDriver(new T_Driver(new T_Equation))             // Dependent objects

    {
    }
private:
    G4double fDeltaChord;
    T_Driver* fIntgrDriver;

public:
    G4double AdvanceChordLimited(G4FieldTrack& yCurrent,
                                  G4double      stepMax,
                                  G4double      epsStep,
                                  const G4ThreeVector latestSafetyOrigin,
                                  G4double       latestSafetyRadius );



     virtual G4double FindNextChord( const  G4FieldTrack& yStart,
                                    G4double     stepMax,
                                    G4FieldTrack&   yEnd, // Endpoint
                                    G4double&   dyErrPos, // Error of endpoint
                                    G4double    epsStep,
                                    G4double*  pStepForAccuracy,
                                    const  G4ThreeVector, //  latestSafetyOrigin,
                                    G4double);       //  latestSafetyRadius


};

#endif
