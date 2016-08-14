#include <iomanip>
#include "G4RevisedChordFinder.hh"
#include "G4SystemOfUnits.hh"
#include "G4LineSection.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4ClassicalRK4.hh"
#include "G4Mag_UsualEqRhs.hh"

//#define G4DEBUG_FIELD 1
//#define ncomp G4FieldTrack::ncompSVEC

G4RevisedChordFinder::G4RevisedChordFinder(G4VIntegrationDriver *pIntDriver, G4int  VerboseLevel):
    fpIntDriver(pIntDriver),
    fVerboseLevel(VerboseLevel),
    fDefaultDeltaChord( 0.25 * mm ),      // Parameters
    fDeltaChord( fDefaultDeltaChord ),    //   Internal parameters
    fFirstFraction(0.999), fFractionLast(1.00),  fFractionNextEstimate(0.98),
    fMultipleRadius(15.0),
    fStatsVerbose(0),
    fLastStepEstimate_Unconstrained(DBL_MAX),          // Should move q, p to
    fTotalNoTrials_FNC(0), fNoCalls_FNC(0), fmaxTrials_FNC(0),
    fpEquation(nullptr),fpStepper(nullptr)

{
    // check the values and set the other parameters
    SetFractions_Last_Next( fFractionLast, fFractionNextEstimate);
}

G4RevisedChordFinder::G4RevisedChordFinder(G4MagneticField* magField, G4double stepMinimum,
                                           G4MagIntegratorStepper* pItsStepper,G4int VerboseLevel):
    fpIntDriver(nullptr),
    fVerboseLevel(VerboseLevel),
    fDefaultDeltaChord( 0.25 * mm ),      // Parameters
    fDeltaChord( fDefaultDeltaChord ),    //   Internal parameters
    fFirstFraction(0.999), fFractionLast(1.00),  fFractionNextEstimate(0.98),
    fMultipleRadius(15.0),
    fStatsVerbose(0),
    fLastStepEstimate_Unconstrained(DBL_MAX),          // Should move q, p to
    fTotalNoTrials_FNC(0), fNoCalls_FNC(0), fmaxTrials_FNC(0),
    fpEquation(nullptr),fpStepper(nullptr)

{
       fpEquation = new G4Mag_UsualEqRhs(magField);
       fAllocatedEquation = true;
       fpStepper = pItsStepper;
       if(fpStepper == nullptr)
       {
          fpStepper =  new G4ClassicalRK4(fpEquation);   // The old default
          fAllocatedStepper = true;
       }
       else
       {
          fAllocatedStepper = false;
       }
       SetIntegrationDriver(new G4MagInt_Driver(stepMinimum, fpStepper, fpStepper->GetNumberOfVariables()));

}

G4RevisedChordFinder::~G4RevisedChordFinder()
{
    delete fpEquation;
    if (fAllocatedStepper) delete fpStepper;

    delete fpIntDriver;
    if( fStatsVerbose ) { PrintStatistics(); }
}


void G4RevisedChordFinder::SetFractions_Last_Next( G4double fractLast, G4double fractNext )
{
  // Use -1.0 as request for Default.
  if( fractLast == -1.0 )   fractLast = 1.0;   // 0.9;
  if( fractNext == -1.0 )   fractNext = 0.98;  // 0.9;

  // fFirstFraction  = 0.999; // Orig 0.999 A safe value, range: ~ 0.95 - 0.999
  // fMultipleRadius = 15.0;  // For later use, range: ~  2 - 20

  if( fStatsVerbose )
  {
    G4cout << " ChordFnd> Trying to set fractions: "
           << " first " << fFirstFraction
           << " last " <<  fractLast
           << " next " <<  fractNext
           << " and multiple " << fMultipleRadius
           << G4endl;
  }

  if( (fractLast > 0.0) && (fractLast <=1.0) )
  {
    fFractionLast= fractLast;
  }
  else
  {
    G4cerr << "G4RevisedChordFinder::SetFractions_Last_Next: Invalid "
           << " fraction Last = " << fractLast
           << " must be  0 <  fractionLast <= 1 " << G4endl;
  }
  if( (fractNext > 0.0) && (fractNext <1.0) )
  {
    fFractionNextEstimate = fractNext;
  }
  else
  {
    G4cerr << "G4RevisedChordFinder:: SetFractions_Last_Next: Invalid "
           << " fraction Next = " << fractNext
           << " must be  0 <  fractionNext < 1 " << G4endl;
  }
}


G4double G4RevisedChordFinder::NewStep(G4double  stepTrialOld,
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



G4FieldTrack G4RevisedChordFinder::ApproxCurvePointS( const G4FieldTrack&  CurveA_PointVelocity,
                                  const G4FieldTrack&  CurveB_PointVelocity,
                                  const G4FieldTrack&  ApproxCurveV,
                                  const G4ThreeVector& CurrentE_Point,
                                  const G4ThreeVector& CurrentF_Point,
                                  const G4ThreeVector& PointG,
                                       G4bool first, G4double eps_step)
{
  // ApproxCurvePointS is 2nd implementation of ApproxCurvePoint.
  // Use Brent Algorithm (or InvParabolic) when possible.
  // Given a starting curve point A (CurveA_PointVelocity), curve point B
  // (CurveB_PointVelocity), a point E which is (generally) not on the curve
  // and  a point F which is on the curve (first approximation), find new
  // point S on the curve closer to point E.
  // While advancing towards S utilise 'eps_step' as a measure of the
  // relative accuracy of each Step.

  G4FieldTrack EndPoint(CurveA_PointVelocity);
  if(!first){EndPoint= ApproxCurveV;}

  G4ThreeVector Point_A,Point_B;
  Point_A=CurveA_PointVelocity.GetPosition();
  Point_B=CurveB_PointVelocity.GetPosition();

  G4double xa,xb,xc,ya,yb,yc;

  // InverseParabolic. AF Intersects (First Part of Curve)

  if(first)
  {
    xa=0.;
    ya=(PointG-Point_A).mag();
    xb=(Point_A-CurrentF_Point).mag();
    yb=-(PointG-CurrentF_Point).mag();
    xc=(Point_A-Point_B).mag();
    yc=-(CurrentE_Point-Point_B).mag();
  }
  else
  {
    xa=0.;
    ya=(Point_A-CurrentE_Point).mag();
    xb=(Point_A-CurrentF_Point).mag();
    yb=(PointG-CurrentF_Point).mag();
    xc=(Point_A-Point_B).mag();
    yc=-(Point_B-PointG).mag();
    if(xb==0.)
    {
      EndPoint=
      ApproxCurvePointV(CurveA_PointVelocity, CurveB_PointVelocity,
                        CurrentE_Point, eps_step);
      return EndPoint;
    }
  }

  const G4double tolerance = 1.e-12;
  if(std::abs(ya)<=tolerance||std::abs(yc)<=tolerance)
  {
    ; // What to do for the moment: return the same point as at start
      // then PropagatorInField will take care
  }
  else
  {
    G4double test_step = InvParabolic(xa,ya,xb,yb,xc,yc);
    G4double curve;
    if(first)
    {
      curve=std::abs(EndPoint.GetCurveLength()
                    -ApproxCurveV.GetCurveLength());
    }
    else
    {
      test_step=(test_step-xb);
      curve=std::abs(EndPoint.GetCurveLength()
                    -CurveB_PointVelocity.GetCurveLength());
      xb=(CurrentF_Point-Point_B).mag();
    }

    if(test_step<=0)    { test_step=0.1*xb; }
    if(test_step>=xb)   { test_step=0.5*xb; }
    if(test_step>=curve){ test_step=0.5*curve; }

    if(curve*(1.+eps_step)<xb) // Similar to ReEstimate Step from
    {                          // G4VIntersectionLocator
      test_step=0.5*curve;
    }

    if (fpIntDriver->isDense())
    {
        fpIntDriver->DoInterpolation(EndPoint,test_step, eps_step);
    }
    else
    {
        fpIntDriver->AccurateAdvance(EndPoint, test_step, eps_step);
    }

#ifdef G4DEBUG_FIELD
    // Printing Brent and Linear Approximation
    //
    G4cout << "G4RevisedChordFinder::ApproxCurvePointS() - test-step ShF = "
           << test_step << "  EndPoint = " << EndPoint << G4endl;

    //  Test Track
    //
    G4FieldTrack TestTrack( CurveA_PointVelocity);
    TestTrack = ApproxCurvePointV( CurveA_PointVelocity,
                                   CurveB_PointVelocity,
                                   CurrentE_Point, eps_step );
    G4cout.precision(14);
    G4cout << "G4RevisedChordFinder::BrentApprox = " << EndPoint  << G4endl;
    G4cout << "G4RevisedChordFinder::LinearApprox= " << TestTrack << G4endl;
#endif
  }
  return EndPoint;
}



G4FieldTrack G4RevisedChordFinder::ApproxCurvePointV(const G4FieldTrack& trackPointA,
                                              const G4FieldTrack& trackPointB,
                                              const G4ThreeVector& pointE,
                                              G4double eps_step)
{
  // If r=|AE|/|AB|, and s=true path lenght (AB)
  // return the point that is r*s along the curve!

  G4FieldTrack   trackPointE = trackPointA;

  G4ThreeVector  pointA = trackPointA.GetPosition();
  G4ThreeVector  pointB = trackPointB.GetPosition();

  G4ThreeVector  vectorAB = pointB - pointA;
  G4ThreeVector  vectorAE = pointE - pointA;

  G4double  ABdist = vectorAB.mag();

  //curve length  of AB
  G4double  curve_length = trackPointB.GetCurveLength() - trackPointA.GetCurveLength();

  //G4cout<<"cl(A): "<<trackPointA.GetCurveLength()<<" cl(B): "<<trackPointB.GetCurveLength()
  //      <<" cl(B - A): "<<curve_length<<" ABdist: "<<ABdist<<G4endl;



  G4double  integrationInaccuracyLimit = std::max(perMillion, 0.5*eps_step );

  if( curve_length < ABdist * (1. - integrationInaccuracyLimit) )
  {
#ifdef G4DEBUG_FIELD
    G4cerr << " Warning in G4RevisedChordFinder::ApproxCurvePointV(): \n"
           << " The two points are further apart than the curve length \n"
           << " Dist = "         << ABdist
           << " curve length = " << curve_length
           << " relativeDiff = " << (curve_length-ABdist)/ABdist
           << G4endl;
    if( curve_length < ABdist * (1. - 10*eps_step) )
    {
        char message[] = "Unphysical curve length.\n "
                         "The size of the above difference exceeds allowed limits. \n "
                         "Aborting.";
        G4Exception("G4RevisedChordFinder::ApproxCurvePointV()", "GeomField0003",
                  FatalException, message);
    }
#endif
    // Take default corrective action: adjust the maximum curve length.
    // NOTE: this case only happens for relatively straight paths.
    // curve_length = ABdist;
  }

  G4double  AE_fraction = 0.5;   // Guess .. ?;
  if ( ABdist > 0.0 )
  {
     AE_fraction = vectorAE.mag() / ABdist;
  }
  else
  {
#ifdef G4DEBUG_FIELD
     G4cout << "Warning in G4RevisedChordFinder::ApproxCurvePointV():"
               " A and B are the same point! \n"
               " Chord AB length = " << vectorAE.mag() << "\n\n";
#endif
  }

  if(AE_fraction > 1.0 + perMillion)
  {
#ifdef G4DEBUG_FIELD
    G4cerr << " G4RevisedChordFinder::ApproxCurvePointV() - Warning:"
           << " Anomalous condition:AE > AB or AE/AB <= 0 " << G4endl
           << "   AE_fraction = " <<  AE_fraction << G4endl
           << "   Chord AE length = " << vectorAE.mag() << G4endl
           << "   Chord AB length = " << ABdist << G4endl << G4endl;
    G4cerr << " OK if this condition occurs after a recalculation of 'B'"
           << G4endl << " Otherwise it is an error. " << G4endl ;
#endif
     // This course can now result if B has been re-evaluated,
     // without E being recomputed (1 July 99).
     // In this case this is not a "real error" - but it is undesired
     // and we cope with it by a default corrective action ...
     //
     AE_fraction = 0.5;                         // Default value
  }


  G4double stepLength = AE_fraction * curve_length;
  /*G4cout<<"tBegin: "<<trackPointE.GetCurveLength()<<
          " tWant: "<<trackPointE.GetCurveLength() + stepLength<<
          " AE_fraction: "<<AE_fraction<<G4endl;
*/
  /*
  if (AE_fraction < G4ThreeVector::getTolerance() || 1 - AE_fraction < G4ThreeVector::getTolerance())
  {
      char buff[256];
      sprintf(buff, "AE_fraction = %g",AE_fraction);
      G4Exception("G4RevisedChordFinder::ApproxCurvePointV()","123",
                  JustWarning, buff);
  }
*/
  if ( AE_fraction > 0.0 )
  {
      //calculate approximation of the intersection point
      if (fpIntDriver->isDense())
      {
        fpIntDriver->DoInterpolation(trackPointE, stepLength, eps_step);
      }
      else
      {
        fpIntDriver->AccurateAdvance(trackPointE, stepLength, eps_step);
      }
  }


  return trackPointE;
}



void G4RevisedChordFinder::PrintStatistics()
{
  // Print Statistics

  G4cout << "G4RevisedChordFinder statistics report: " << G4endl;
  G4cout
    << "  No trials: " << fTotalNoTrials_FNC
    << "  No Calls: "  << fNoCalls_FNC
    << "  Max-trial: " <<  fmaxTrials_FNC
    << G4endl;
  G4cout
    << "  Parameters: "
    << "  fFirstFraction "  << fFirstFraction
    << "  fFractionLast "   << fFractionLast
    << "  fFractionNextEstimate " << fFractionNextEstimate
    << G4endl;
}



void G4RevisedChordFinder::TestChordPrint( G4int    noTrials,
                                    G4int    lastStepTrial,
                                    G4double dChordStep,
                                    G4double nextStepTrial )
{
     G4int oldprec= G4cout.precision(5);
     G4cout << " ChF/fnc: notrial " << std::setw( 3) << noTrials
            << " this_step= "       << std::setw(10) << lastStepTrial;
     if( std::fabs( (dChordStep / fDeltaChord) - 1.0 ) < 0.001 )
     {
       G4cout.precision(8);
     }
     else
     {
       G4cout.precision(6);
     }
     G4cout << " dChordStep=  " << std::setw(12) << dChordStep;
     if( dChordStep > fDeltaChord ) { G4cout << " d+"; }
     else                           { G4cout << " d-"; }
     G4cout.precision(5);
     G4cout <<  " new_step= "       << std::setw(10)
            << fLastStepEstimate_Unconstrained
            << " new_step_constr= " << std::setw(10)
            << lastStepTrial << G4endl;
     G4cout << " nextStepTrial = " << std::setw(10) << nextStepTrial << G4endl;
     G4cout.precision(oldprec);
}


G4double G4RevisedChordFinder::AdvanceChordLimited(G4FieldTrack& trackCurrent,
                                                   G4double  stepMax,
                                                   G4double epsStep)
{
    if (!fpIntDriver->isDense())
    {
      G4double dyErr;
      G4FieldTrack yEnd = trackCurrent;
      G4double startCurveLen = trackCurrent.GetCurveLength();
      G4double nextStep;

      G4double stepPossible = FindNextChord(trackCurrent, stepMax, yEnd, dyErr, epsStep,&nextStep);
      //G4cout<<"G4RevisedChordFinder::AdvanceChordLimited \n stepPossible: "<<stepPossible<<G4endl;

      G4bool good_advance = false;

      if ( dyErr < epsStep * stepPossible )
      {
         // Accept this accuracy.

         trackCurrent = yEnd;
         good_advance = true;
      }
      else
      {
         // Advance more accurately to "end of chord"
         //                           ***************
         good_advance = GetIntegrationDriver()->AccurateAdvance(trackCurrent, stepPossible,
                                                                epsStep, nextStep);
         //G4cout<<"AccurateAdvance: "<<good_advance<<G4endl;
         if ( ! good_advance )
         {
           // In this case the driver could not do the full distance
           stepPossible = trackCurrent.GetCurveLength()-startCurveLen;
         }
      }
      return stepPossible;
    }
    else
    {
        const G4double curveLength = trackCurrent.GetCurveLength();
        G4double clWant = curveLength + stepMax;
        G4double stepPossible = DBL_MAX;

        // initialize interpolation. Do a big step.
        G4FieldTrack tmpTrack(trackCurrent);
        interpolationInterval& interval = fpIntDriver->GetInterpolationInterval();
        if (curveLength < interval.first || (curveLength + GetIntegrationDriver()->GetMinimumStep()) > interval.second)
        {
            stepPossible = std::max(stepMax, GetIntegrationDriver()->GetMinimumStep());
            fpIntDriver->DoStep(tmpTrack, stepPossible, epsStep);
            //const G4double lineDist = (tmpTrack.GetPosition() - trackCurrent.GetPosition()).mag();
            //G4cout<<"lineDist: "<<lineDist<<" curveDist: "<<lastStepSize<<"\n\n\n";
            //CheckInterpol(trackCurrent, 1000);
            //++i;
            //G4cout<<"interval: "<<interval.first<<","<<interval.second<<"\n\n";
        }


         if (clWant > interval.second) clWant = interval.second;
         stepPossible = clWant - curveLength;

         //tmpTrack = trackCurrent;

         //find step size using interpolation, moves trackCurrent
         tmpTrack = trackCurrent;
         G4double dyErr;
         G4double hdid =  FindNextChord(trackCurrent, stepPossible, tmpTrack, dyErr, epsStep, nullptr);
         trackCurrent = tmpTrack;
         //G4cout<<"hdid: "<<hdid<<G4endl;

         //G4cout<<"lastStepSize: "<<lastStepSize<<" hdid: "<<hdid<<" eps: "<<eps<<G4endl;

         return hdid;
    }
}


// Returns Length of Step taken
G4double G4RevisedChordFinder::FindNextChord(const  G4FieldTrack& trackStart, G4double stepMax,
                                             G4FieldTrack& trackEnd, G4double& dyErrPos,
                                             G4double epsStep, G4double* pStepForAccuracy)
{
    //  1.)  Try to "leap" to end of interval
    //  2.)  Evaluate if resulting chord gives d_chord that is good enough.
    // 2a.)  If d_chord is not good enough, find one that is.

    G4FieldTrack trackMid = trackStart;
                 trackEnd = trackStart;

    G4double clStart = trackStart.GetCurveLength();

    G4double stepForAccuracy;

    G4bool    validEndPoint = false;
    G4double  dChordStep, lastStepLength; //  stepOfLastGoodChord;

    G4double dydx[G4FieldTrack::ncompSVEC];
    if (!fpIntDriver->isDense())
    {
        GetIntegrationDriver()-> GetDerivatives(trackEnd, dydx);
    }

    unsigned int        noTrials = 0;
    const unsigned int  maxTrials = 300; // Avoid endless loop for bad convergence

    const G4double safetyFactor = GetFirstFraction(); //  0.975 or 0.99 ? was 0.999


    G4double stepTrial = std::min(stepMax, safetyFactor* GetLastStepEstimateUnc());

    G4double newStepEst_Uncons =  0.0;
    G4double stepForChord;

    do
    {
        // Always start from initial point
        if (fpIntDriver->isDense())
        {
            trackEnd.SetCurveLength(clStart);
            trackMid.SetCurveLength(clStart);

            fpIntDriver->DoInterpolation(trackEnd, stepTrial);
            fpIntDriver->DoInterpolation(trackMid, stepTrial/2.);



            dChordStep = G4LineSection::Distline(trackMid.GetPosition(),
                                                 trackStart.GetPosition(),
                                                 trackEnd.GetPosition());
/*
            G4cout<<"start: "<<trackStart.GetPosition()<<" ("<<trackStart.GetCurveLength()<<") \n"
                    " mid: "<<trackMid.GetPosition()<<" ("<<trackMid.GetCurveLength()<<") \n"
                    " end: "<<trackEnd.GetPosition()<<" ("<<trackEnd.GetCurveLength()<<") \n"
                    " dist: "<<dChordStep<<G4endl<<G4endl;
*/
        }
        else
        {
            trackEnd = trackStart;
            GetIntegrationDriver()->QuickAdvance(trackEnd, dydx, stepTrial, dChordStep, dyErrPos);
        }



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
                << "Proposed for chord = " << stepForChord  << G4endl;
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
  return stepTrial;
}
