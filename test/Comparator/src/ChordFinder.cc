#include "ChordFinder.hh"
template <class T_Field, class T_Equation, class T_Driver>
G4double ChordFinder<T_Field, T_Equation, T_Driver>::AdvanceChordLimited(G4FieldTrack& yCurrent,
                              G4double      stepMax,
                              G4double      epsStep,
                              const G4ThreeVector latestSafetyOrigin,
                              G4double       latestSafetyRadius ){


     G4double stepPossible;
     G4double dyErr;
     G4FieldTrack yEnd( yCurrent);
     G4double  startCurveLen= yCurrent.GetCurveLength();
     G4double nextStep;


     stepPossible= FindNextChord(yCurrent, stepMax, yEnd, dyErr, epsStep,
                                 &nextStep, latestSafetyOrigin, latestSafetyRadius);


     G4bool good_advance;


     if (dyErr < epsStep * stepPossible ){
                // Accept this accuracy.

         yCurrent = yEnd;
         good_advance = true;}

     else{
                // Advance more accurately to "end of chord"
                //                           ***************

             good_advance = fIntgrDriver->AccurateAdvance(yCurrent, stepPossible,
                                                          epsStep, nextStep);

             if ( ! good_advance ){
                    // In this case the driver could not do the full distance
                    stepPossible= yCurrent.GetCurveLength()-startCurveLen;
                }
    }

    return stepPossible;

 }

/*
G4double ChordFinder::FindNextChord( const  G4FieldTrack& yStart,
                                G4double     stepMax,
                                G4FieldTrack&   yEnd, // Endpoint
                                G4double&   dyErrPos, // Error of endpoint
                                G4double    epsStep,
                                G4double*  pStepForAccuracy,
                                const  G4ThreeVector, //  latestSafetyOrigin,
                                G4double){       //  latestSafetyRadius


     // Returns Length of Step taken
     G4FieldTrack yCurrent=  yStart;
     G4double    stepTrial, stepForAccuracy;
     G4double    dydx[G4FieldTrack::ncompSVEC];

     //  1.)  Try to "leap" to end of interval
     //  2.)  Evaluate if resulting chord gives d_chord that is good enough.
     // 2a.)  If d_chord is not good enough, find one that is.

     G4bool validEndPoint= false;
     G4double dChordStep, lastStepLength; //  stepOfLastGoodChord;
     fIntgrDriver-> GetDerivatives( yCurrent, dydx );
     G4int noTrials=0;

     const G4double safetyFactor= fFirstFraction; //  0.975 or 0.99 ? was 0.999

     stepTrial = std::min( stepMax, safetyFactor*fLastStepEstimate_Unconstrained );
     G4double newStepEst_Uncons= 0.0;

     do{
         G4double stepForChord;
         yCurrent = yStart;    // Always start from initial point

         fIntgrDriver->QuickAdvance( yCurrent, dydx, stepTrial, dChordStep, dyErrPos);

         //  We check whether the criterion is met here.
         validEndPoint = AcceptableMissDist(dChordStep);

         lastStepLength = stepTrial;

         // This method estimates to step size for a good chord.
         stepForChord = NewStep(stepTrial, dChordStep, newStepEst_Uncons );

         if( ! validEndPoint ){

             if( stepTrial<=0.0 ){
                 stepTrial = stepForChord;
             }

             else if (stepForChord <= stepTrial){

                 // Reduce by a fraction, possibly up to 20%
                 stepTrial = std::min( stepForChord, fFractionLast * stepTrial);
             }

             else{
                 stepTrial *= 0.1;
             }
         }
         noTrials++;
     }while( ! validEndPoint );   // End of do-while  RKD

     yEnd=  yCurrent;
     return stepTrial;
}*/
