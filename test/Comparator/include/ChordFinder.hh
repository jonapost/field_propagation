#ifndef ChordFinder_HH
#define ChordFinder_HH

#include "G4FieldTrack.hh"
#include <iomanip>
#include "G4SystemOfUnits.hh"

#include "G4EquationOfMotion.hh"

template <class T_Field, class T_Driver>
class ChordFinder{

public:
    ChordFinder(G4EquationOfMotion* equation);


    G4double AdvanceChordLimited(G4FieldTrack& yCurrent,
                                  G4double      hstep,
                                  G4double      epsStep);



     virtual G4double FindNextChord(const  G4FieldTrack& trackIn,
                                    G4double     hstep,
                                    G4FieldTrack&   trackOut, // Endpoint
                                    G4double&   yErr, // Error of endpoint
                                    G4double    epsStep);       //  latestSafetyRadius

private:
    G4double fDeltaChord;
    T_Driver* fIntgrDriver;


};

template <class T_Field, class T_Driver>
ChordFinder<T_Field, T_Driver>::ChordFinder(G4EquationOfMotion* equation):
            fDeltaChord(0.25 * mm),    //   Internal parameters
            fIntgrDriver(new T_Driver(equation))             // Dependent objects

    {
    }


template <class T_Field, class T_Driver>
G4double ChordFinder<T_Field, T_Driver>::AdvanceChordLimited(G4FieldTrack& yCurrent,
                                                                         G4double hstep,
                                                                         G4double epsStep){


     G4double stepPossible;
     G4double yErr;
     G4FieldTrack yEnd(yCurrent);

     stepPossible = FindNextChord(yCurrent, hstep, yEnd, yErr, epsStep);
     G4cout<<"AdvanceChordLimited: stepPossible "<<stepPossible<<G4endl;

     if (yErr < epsStep * stepPossible){
     // Accept this accuracy.
         yCurrent = yEnd;
     }

     else{
             fIntgrDriver->AccurateAdvance(yCurrent, stepPossible, epsStep);
    }

    return stepPossible;
 }

template <class T_Field, class T_Driver>
G4double ChordFinder<T_Field, T_Driver>::FindNextChord(const  G4FieldTrack& trackIn,
                                G4double     hstep,
                                G4FieldTrack&   trackOut,
                                G4double&   yErr, // Error of endpoint
                                G4double    epsStep){


     // Returns Length of Step taken
     G4FieldTrack yCurrent =  trackIn;
     //  1.)  Try to "leap" to end of interval
     //  2.)  Evaluate if resulting chord gives d_chord that is good enough.
     // 2a.)  If d_chord is not good enough, find one that is.

     G4bool validEndPoint = false;
     G4double dChordStep;

     do{
         yCurrent = trackIn;    // Always start from initial point

         yErr = fIntgrDriver->QuickAdvance(yCurrent, hstep, dChordStep, epsStep);

         //  We check whether the criterion is met here.
         validEndPoint = dChordStep < fDeltaChord;
         G4cout<<"dChordStep "<<dChordStep<<"  fDeltaChord "<<fDeltaChord<<" hstep "<<hstep<<G4endl;

         if( ! validEndPoint ) hstep *= 0.1;

     }while(!validEndPoint);   // End of do-while  RKD

     trackOut =  yCurrent;

     return hstep;
}

#endif
