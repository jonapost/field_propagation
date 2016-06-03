#ifndef ChordFinder_HH
#define ChordFinder_HH

#include "G4FieldTrack.hh"
#include <iomanip>
#include "G4SystemOfUnits.hh"

#include "G4EquationOfMotion.hh"

template <class T_Driver>
class ChordFinder{

public:
    ChordFinder(G4EquationOfMotion* equation, G4double deltaChord = 0.25*mm);


    G4double AdvanceChordLimited(G4FieldTrack& yCurrent,
                                 G4double      hstep,
                                 G4double      epsStep);



    G4double FindNextChord(G4FieldTrack& track,
                                   G4double hstep, // Error of endpoint
                                   G4double epsStep);

private:
    G4double fDeltaChord;
    T_Driver* fIntgrDriver;


};

template <class T_Driver>
ChordFinder<T_Driver>::ChordFinder(G4EquationOfMotion* equation, G4double deltaChord):
                                   fDeltaChord(deltaChord),
                                   fIntgrDriver(new T_Driver(equation))
    {
    }


template <class T_Driver>
G4double ChordFinder<T_Driver>::AdvanceChordLimited(G4FieldTrack& track,
                                                    G4double hstep,
                                                    G4double epsStep){

    //G4cout<<"AdvanceChordLimited: hstep "<<hstep<<G4endl;
    G4FieldTrack tmpTrack(track);
    G4double dostep = FindNextChord(tmpTrack, hstep, epsStep);
    //G4cout<<"do step "<<dostep<<G4endl;

    fIntgrDriver->AccurateAdvance(track, dostep, epsStep);

    return dostep;
 }

template <class T_Driver>
G4double ChordFinder<T_Driver>::FindNextChord(G4FieldTrack& track,
                                              G4double  hstep,
                                              G4double epsStep){


    //  1)  Try to "leap" to end of interval
    //  2)  Evaluate if resulting chord gives d_chord that is good enough.
    // 2a)  If d_chord is not good enough, find one that is.

     G4FieldTrack yCurrent =  track;

     G4bool validEndPoint = false;
     G4double chordStep;

     do{
         yCurrent = track;    // Always start from initial point

         hstep = fIntgrDriver->QuickAdvance(yCurrent, hstep, chordStep, epsStep);

         // Check whether the chord is small enough.
         validEndPoint = chordStep < fDeltaChord;
         //G4cout<<"chordStep "<<chordStep<<"  fDeltaChord "<<fDeltaChord<<" hstep "<<hstep<<G4endl;

         if(!validEndPoint ) hstep *= 0.1;

     }while(!validEndPoint);

     track =  yCurrent;

     return hstep;
}

#endif
