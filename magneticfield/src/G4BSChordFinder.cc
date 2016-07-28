#include "G4BSChordFinder.hh"

#include "G4SystemOfUnits.hh"
#include "G4LineSection.hh"

#include <cstdio>
#include <iomanip>

#include "G4MagIntegratorDriver.hh"
#include "G4ClassicalRK4.hh"

#include "BulirschStoerDriver.hh"

//#define G4DEBUG_FIELD 1
//#define DEBUG_PRINTS

#define ncomp G4FieldTrack::ncompSVEC

G4BSChordFinder::G4BSChordFinder(G4double hminimum,
                             G4EquationOfMotion* equation,
                             G4int numberOfComponents,
                             G4int  statisticsVerbosity):
    G4VChordFinder(new BulirschStoerDriver(hminimum, equation), statisticsVerbosity),
    fIntgrDriver(equation, numberOfComponents, 0, true),
    fMinimumStep(hminimum),
    fverb(statisticsVerbosity),
    tBegin(0),tEnd(0),
    nextStepSize(DBL_MAX),
    lastStepSize(0),
    eps_prev(0)

{
}

G4BSChordFinder::~G4BSChordFinder()
{
}

void G4BSChordFinder::reset()
{
    tBegin = tEnd = 0;
}

G4double G4BSChordFinder::OneGoodStep(G4FieldTrack& track, G4double stepLen, G4double eps)
{
    fIntgrDriver.set_max_relative_error(eps);
    track.DumpToArray(yIn);
    G4double time = track.GetCurveLength();
    tBegin = time;
    GetEquationOfMotion()->RightHandSide(yIn,dydx);
    eps_prev = eps;

    G4int result = step_result::fail;
    G4int ncalls = 0, max_calls = 10000;

    nextStepSize = stepLen;

    while(result == step_result::fail)
    {
        nextStepSize = std::min(stepLen,nextStepSize);
        //nextStepSize = std::max(nextStepSize,fMinimumStep);
        result = fIntgrDriver.try_step(yIn,dydx,time,yOut,nextStepSize);
        ++ncalls;
        if (ncalls > max_calls)
        {
            char buff[256];
            sprintf(buff, "ncalls exeeded max_calls(%i)! \n "
                          "stepLen: %g, nextStepSize: %g, fminStep: %g",max_calls, stepLen, nextStepSize, fMinimumStep);
            G4Exception("G4BSChordFinder::OneGoodStep()", "GeomField0001",
                        FatalException, buff);
        }
    }

    //G4cout<<"ncalls: "<<ncalls<<G4endl;
    //G4cout<<fIntgrDriver.GetCurrentOrder()<<G4endl;

    memcpy(yOut + GetNumberOfVariables(),
           yIn + GetNumberOfVariables(),
           sizeof(G4double)*(ncomp-GetNumberOfVariables()));

    track.LoadFromArray(yOut,ncomp);
    track.SetCurveLength(time);
    tEnd = time;
    //G4cout<<"hdid: "<<tEnd - tBegin<<G4endl;
    return tEnd - tBegin;
}

void G4BSChordFinder::DoInterpolation(G4FieldTrack &track, G4double hstep, G4double eps)
{
    G4double time = track.GetCurveLength();
    G4double twant = time + hstep;

    //little upperflow, allow.
    G4double upperflow = 0;
    if (twant > tEnd)
        upperflow = (twant - tEnd)/twant;
    if (upperflow > 0 && upperflow < perMillion){
        twant = tEnd;
#ifdef DEBUG_PRINTS
        G4cout<<"twant: "<<twant<<" tEnd: "<<tEnd<<G4endl;
#endif
    }

    if (time >= tBegin && twant <= tEnd)
    {
        track.DumpToArray(yOut);
        fIntgrDriver.do_interpolation(twant, yOut);
        track.LoadFromArray(yOut, ncomp);
        track.SetCurveLength(twant);
        if (eps != 0 && eps != eps_prev)
        {
            char buff[256];
            sprintf(buff,"Accuracy changed. eps: %g, eps_prev: %g Interpolation is not accurate!",eps,eps_prev);
            G4Exception("G4BSChordFinder::DoInterpolation()", "GeomField0001",
                        FatalException, buff);
        }
    }
    else
    {
        char buff[256];
        sprintf(buff,"time = %g is out of interpolation interval (%g,%g)!",twant,tBegin,tEnd);
        G4Exception("G4BSChordFinder::DoInterpolation()", "GeomField0001", FatalException, buff);
    }
}


void G4BSChordFinder::DoStepForChord(G4FieldTrack &track, G4double stepLen)
{
    DoInterpolation(track, stepLen);
}

G4bool G4BSChordFinder::DoStepForIntersection(G4FieldTrack &track, G4double stepLen, G4double eps)
{
    //DoInterpolation(track, stepLen, eps);
    //return true;

    return GetIntegrationDriver()->AccurateAdvance(track, stepLen, eps);
}


G4double G4BSChordFinder::AdvanceChordLimited(G4FieldTrack& trackCurrent,
                                            G4double stepMax, G4double eps)
{
    const G4double time = trackCurrent.GetCurveLength();
    G4double twant = time + stepMax;
    G4double stepPossible = DBL_MAX;

    // initialize interpolation. Do a step.
    G4FieldTrack tmpTrack(trackCurrent);
    if (time < tBegin || (time + fMinimumStep) >= tEnd)
    {
        //G4FieldTrack tmpTrack(trackCurrent);
        stepPossible = std::min(stepMax, nextStepSize);
        stepPossible = std::max(stepPossible, fMinimumStep);
        lastStepSize = OneGoodStep(tmpTrack, stepPossible, eps);
    }


    if (twant > tEnd) twant = tEnd;
    stepPossible = twant - time;

    tmpTrack = trackCurrent;

    //find step size using interpolation
    G4double hdid =  FindNextChord(trackCurrent, stepPossible);

    //G4cout<<"lastStepSize: "<<lastStepSize<<" hdid: "<<hdid<<" eps: "<<eps<<G4endl;

    return hdid;
}


// Returns Length of Step taken
G4double G4BSChordFinder::FindNextChord(G4FieldTrack& trackCurrent, G4double stepMax)
{
    //  1.)  Try to "leap" to end of interval
    //  2.)  Evaluate if resulting chord gives d_chord that is good enough.
    // 2a.)  If d_chord is not good enough, find one that is.

    G4double tCurrent = trackCurrent.GetCurveLength();
    G4FieldTrack trackEnd =  trackCurrent, trackMid = trackCurrent;
    G4double    stepTrial;

    G4bool    validEndPoint = false;
    G4double  dChordStep = DBL_MAX, lastStepLength; //  stepOfLastGoodChord;


  unsigned int        noTrials = 0;
  const unsigned int  maxTrials = 300; // Avoid endless loop for bad convergence

  const G4double safetyFactor = GetFirstFraction(); //  0.975 or 0.99 ? was 0.999

  stepTrial = std::min( stepMax, safetyFactor*GetLastStepEstimateUnc());

  G4double newStepEst_Uncons = 0.0;
  G4double stepForChord;
  do
  {
     // Always start from initial point
     trackEnd.SetCurveLength(tCurrent);
     trackMid.SetCurveLength(tCurrent);


     DoStepForChord(trackEnd, stepTrial);
     DoStepForChord(trackMid, stepTrial/2.);


     //  We check whether the criterion is met here.
     dChordStep = G4LineSection::Distline(trackMid.GetPosition(),
                                          trackCurrent.GetPosition(),
                                          trackEnd.GetPosition());
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
      G4Exception("G4BSChordFinder::FindNextChord()", "GeomField0003",
                  JustWarning, message);
  }

  if( newStepEst_Uncons > 0.0  )
  {
      SetLastStepEstimateUnc(newStepEst_Uncons);
  }

  AccumulateStatistics( noTrials );


  //change track
  trackCurrent = trackEnd;

  return stepTrial;
}
