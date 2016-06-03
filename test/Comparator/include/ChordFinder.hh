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

    bool AcceptableMissDist(G4double dChordStep);

    G4double NewStep(G4double  stepTrialOld,
                     G4double  dChordStep, // Curr. dchord achieved
                     G4double& stepEstimate_Unconstrained );


private:
    G4double fDeltaChord;
    T_Driver* fIntgrDriver;

    G4double fFractionNextEstimate;
    G4double fFractionLast;
    G4double fLastStepEstimate_Unconstrained;
    G4double fFirstFraction;
};

template <class T_Driver>
ChordFinder<T_Driver>::ChordFinder(G4EquationOfMotion* equation, G4double deltaChord):
                                   fDeltaChord(deltaChord),
                                   fIntgrDriver(new T_Driver(equation))
    {
        fFirstFraction = 0.999;
        fFractionLast = 1.00;
        fFractionNextEstimate = 0.98;
        fLastStepEstimate_Unconstrained = DBL_MAX;
    }


template <class T_Driver>
G4double ChordFinder<T_Driver>::AdvanceChordLimited(G4FieldTrack& track,
                                                    G4double hstep,
                                                    G4double epsStep){

    G4cout<<"AdvanceChordLimited: hstep "<<hstep<<G4endl;
    G4FieldTrack tmpTrack(track);
    G4double dostep = FindNextChord(tmpTrack, hstep, epsStep);
    G4cout<<"do step "<<dostep<<G4endl;

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

     const G4double safetyFactor= fFirstFraction; //  0.975 or 0.99 ? was 0.999
     G4double stepTrial = std::min(hstep, safetyFactor*fLastStepEstimate_Unconstrained);

     G4bool validEndPoint = false;
     G4double dChordStep;
     G4double stepForChord;
     G4double newStepEst_Uncons= 0.0;


     do{
         yCurrent = track;    // Always start from initial point


         stepTrial = fIntgrDriver->QuickAdvance(yCurrent, stepTrial, dChordStep, epsStep);

         // Check whether the chord is small enough.
         validEndPoint = AcceptableMissDist(dChordStep);
         G4cout<<"stepTrial "<<stepTrial<<" valid "<<validEndPoint<<G4endl;

         // This method estimates to step size for a good chord.
         stepForChord = NewStep(stepTrial, dChordStep, newStepEst_Uncons );
         G4cout<<"stepForChord "<<stepForChord<<G4endl;
         if(!validEndPoint)
         {
            if(stepTrial<=0.0 )
            {
              stepTrial = stepForChord;
            }
            else if (stepForChord <= stepTrial)
            {
              // Reduce by a fraction, possibly up to 20%
              stepTrial = std::min( stepForChord, fFractionLast * stepTrial);
              G4cout<<"if "<<stepTrial<<G4endl;
            }
            else
            {
              stepTrial *= 0.1;
            }
         }
         //G4cout<<"  "<<stepTrial<<G4endl;

     }while(!validEndPoint);

     fLastStepEstimate_Unconstrained= newStepEst_Uncons;

     track =  yCurrent;

     return stepTrial;
}

template <class T_Driver>
bool ChordFinder<T_Driver>::AcceptableMissDist(G4double dChordStep){
    return dChordStep < fDeltaChord;
}

template <class T_Driver>
G4double ChordFinder<T_Driver>::NewStep(G4double  stepTrialOld,
                                G4double  dChordStep, // Curr. dchord achieved
                                G4double& stepEstimate_Unconstrained )
{
  // Is called to estimate the next step size, even for successful steps,
  // in order to predict an accurate 'chord-sensitive' first step
  // which is likely to assist in more performant 'stepping'.

  G4double stepTrial;

#if 1

  if (dChordStep > 0.0)
  {
    stepEstimate_Unconstrained =
                 stepTrialOld*std::sqrt( fDeltaChord / dChordStep );
    stepTrial =  fFractionNextEstimate * stepEstimate_Unconstrained;
  }
  else
  {
    // Should not update the Unconstrained Step estimate: incorrect!
    stepTrial =  stepTrialOld * 2.;
  }

  if( stepTrial <= 0.001 * stepTrialOld)
  {
     if ( dChordStep > 1000.0 * fDeltaChord )
     {
        stepTrial= stepTrialOld * 0.03;
     }
     else
     {
        if ( dChordStep > 100. * fDeltaChord )
        {
          stepTrial= stepTrialOld * 0.1;
        }
        else   // Try halving the length until dChordStep OK
        {
          stepTrial= stepTrialOld * 0.5;
        }
     }
  }
  else if (stepTrial > 1000.0 * stepTrialOld)
  {
     stepTrial= 1000.0 * stepTrialOld;
  }

  if( stepTrial == 0.0 )
  {
     stepTrial= 0.000001;
  }

#else

  if ( dChordStep > 1000. * fDeltaChord )
  {
        stepTrial= stepTrialOld * 0.03;
  }
  else
  {
     if ( dChordStep > 100. * fDeltaChord )
     {
        stepTrial= stepTrialOld * 0.1;
     }
     else  // Keep halving the length until dChordStep OK
     {
        stepTrial= stepTrialOld * 0.5;
     }
  }

#endif

  // A more sophisticated chord-finder could figure out a better
  // stepTrial, from dChordStep and the required d_geometry
  //   e.g.
  //      Calculate R, r_helix (eg at orig point)
  //      if( stepTrial < 2 pi  R )
  //          stepTrial = R arc_cos( 1 - fDeltaChord / r_helix )
  //      else
  //          ??

  return stepTrial;
}


#endif
