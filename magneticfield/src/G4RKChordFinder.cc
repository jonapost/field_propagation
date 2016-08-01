#include "G4RKChordFinder.hh"

#include "G4MagIntegratorDriver.hh"
#include "G4ClassicalRK4.hh"
#include "G4Mag_UsualEqRhs.hh"

//#define G4DEBUG_FIELD 1
//#define DEBUG_PRINTS


G4RKChordFinder::G4RKChordFinder(G4VIntegrationDriver *pIntegrationDriver,
                                 G4int staticsticsVerbosity):
    G4VRevisedChordFinder(pIntegrationDriver, staticsticsVerbosity),
    fEquation(nullptr),
    fStepper(nullptr),
    fAllocatedStepper(false),
    fAllocatedEquation(false)
{
}

G4RKChordFinder::G4RKChordFinder(G4MagneticField *magField,
                                 G4double stepMinimum,
                                 G4MagIntegratorStepper *pItsStepper,
                                 G4int statisticsVerbosity):
    G4VRevisedChordFinder(nullptr, statisticsVerbosity)
{
    fEquation = new G4Mag_UsualEqRhs(magField);
    fAllocatedEquation = true;
    fStepper = pItsStepper;

    if(fStepper == nullptr)
    {
       fStepper =  new G4ClassicalRK4(fEquation);   // The old default
       fAllocatedStepper = true;
    }
    else
    {
       fAllocatedStepper = false;
    }
    SetIntegrationDriver(new G4MagInt_Driver(stepMinimum, fStepper, fStepper->GetNumberOfVariables()));
}


G4RKChordFinder::~G4RKChordFinder()
{
    delete fEquation; // fIntgrDriver->pIntStepper->theEquation_Rhs;
    if (fAllocatedStepper) delete fStepper;
}

G4RKChordFinder::G4RKChordFinder(const G4RKChordFinder& other):
    G4VRevisedChordFinder(nullptr, other.GetVerbose())
{
    G4Exception("G4RKChordFinder::G4RKChordFinder()", "GeomField0003",
                JustWarning, "Copy constructor is not implemented!");
}

const G4RKChordFinder& G4RKChordFinder::operator =(const G4RKChordFinder& /*other*/){
    G4Exception("G4RKChordFinder::operator =()", "GeomField0003",
                JustWarning, "Assignment operator is not implemented!");
    return *this;
}

void G4RKChordFinder::reset(){
}

G4bool G4RKChordFinder::DoStepForIntersection(G4FieldTrack &track, G4double stepLen, G4double eps){
    return GetIntegrationDriver()->AccurateAdvance(track, stepLen, eps);
}


G4double G4RKChordFinder::AdvanceChordLimited( G4FieldTrack& yCurrent, G4double  stepMax, G4double epsStep)
{
  G4double stepPossible;
  G4double dyErr;
  G4FieldTrack yEnd( yCurrent);
  G4double  startCurveLen= yCurrent.GetCurveLength();
  G4double nextStep;
  //            *************
  stepPossible = FindNextChord(yCurrent, stepMax, yEnd, dyErr, epsStep,&nextStep);
  //            *************

  G4bool good_advance;

  if ( dyErr < epsStep * stepPossible )
  {
     // Accept this accuracy.

     yCurrent = yEnd;
     good_advance = true;
  }
  else
  {
     // Advance more accurately to "end of chord"
     //                           ***************
     good_advance = GetIntegrationDriver()->AccurateAdvance(yCurrent, stepPossible,
                                                  epsStep, nextStep);
     if ( ! good_advance )
     {
       // In this case the driver could not do the full distance
       stepPossible= yCurrent.GetCurveLength()-startCurveLen;
     }
  }
  return stepPossible;
}


// ............................................................................

G4double
G4RKChordFinder::FindNextChord( const  G4FieldTrack& yStart,
                                     G4double     stepMax,
                                     G4FieldTrack&   yEnd, // Endpoint
                                     G4double&   dyErrPos, // Error of endpoint
                                     G4double    epsStep,
                                     G4double*  pStepForAccuracy)
{
  // Returns Length of Step taken

  G4FieldTrack yCurrent=  yStart;
  G4double    stepTrial, stepForAccuracy;
  G4double    dydx[G4FieldTrack::ncompSVEC];

  //  1.)  Try to "leap" to end of interval
  //  2.)  Evaluate if resulting chord gives d_chord that is good enough.
  // 2a.)  If d_chord is not good enough, find one that is.

  G4bool    validEndPoint= false;
  G4double  dChordStep, lastStepLength; //  stepOfLastGoodChord;

  GetIntegrationDriver()-> GetDerivatives( yCurrent, dydx );

  unsigned int        noTrials=0;
  const unsigned int  maxTrials= 300; // Avoid endless loop for bad convergence

  const G4double safetyFactor= GetFirstFraction(); //  0.975 or 0.99 ? was 0.999

  stepTrial = std::min( stepMax, safetyFactor* GetLastStepEstimateUnc());

  G4double newStepEst_Uncons= 0.0;
  G4double stepForChord;
  do
  {
     yCurrent = yStart;    // Always start from initial point

     //            ************
     GetIntegrationDriver()->QuickAdvance( yCurrent, dydx, stepTrial,
                                 dChordStep, dyErrPos);
     //            ************

     //  We check whether the criterion is met here.
     validEndPoint = AcceptableMissDist(dChordStep);

     lastStepLength = stepTrial;

     // This method estimates to step size for a good chord.
     stepForChord = NewStep(stepTrial, dChordStep, newStepEst_Uncons );

     if( ! validEndPoint )
     {
        if( stepTrial<=0.0 )
        {
          stepTrial = stepForChord;
        }
        else if (stepForChord <= stepTrial)
        {
          // Reduce by a fraction, possibly up to 20%
          stepTrial = std::min( stepForChord, GetFractionLast() * stepTrial);
        }
        else
        {
          stepTrial *= 0.1;
        }
     }
     noTrials++;
  }
  while( (! validEndPoint) && (noTrials < maxTrials) );   // End of do-while  RKD

  if( noTrials >= maxTrials )
  {
      std::ostringstream message;
      message << "Exceeded maximum number of trials= " << maxTrials << G4endl
              << "Current sagita dist= " << dChordStep << G4endl
              << "Step sizes (actual and proposed): " << G4endl
              << "Last trial =         " << lastStepLength  << G4endl
              << "Next trial =         " << stepTrial  << G4endl
              << "Proposed for chord = " << stepForChord  << G4endl
              ;
      G4Exception("G4ChordFinder::FindNextChord()", "GeomField0003",
                  JustWarning, message);
  }

  if( newStepEst_Uncons > 0.0  )
  {
      SetLastStepEstimateUnc(newStepEst_Uncons);
  }

  AccumulateStatistics( noTrials );

  if( pStepForAccuracy )
  {
     // Calculate the step size required for accuracy, if it is needed
     //
     G4double dyErr_relative = dyErrPos/(epsStep*lastStepLength);
     if( dyErr_relative > 1.0 )
     {
        stepForAccuracy =
                GetIntegrationDriver()->ComputeNewStepSize( dyErr_relative, lastStepLength );
     }
     else
     {
        stepForAccuracy = 0.0;   // Convention to show step was ok
     }
     *pStepForAccuracy = stepForAccuracy;
  }

#ifdef  TEST_CHORD_PRINT
  static int dbg=0;
  if( dbg )
  {
    G4cout << "ChordF/FindNextChord:  NoTrials= " << noTrials
           << " StepForGoodChord=" << std::setw(10) << stepTrial << G4endl;
  }
#endif
  yEnd =  yCurrent;
  return stepTrial;
}


