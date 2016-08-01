#include <iomanip>
#include "G4VRevisedChordFinder.hh"
#include "G4SystemOfUnits.hh"
#include "G4LineSection.hh"

//#define G4DEBUG_FIELD 1

#define ncomp G4FieldTrack::ncompSVEC

G4VRevisedChordFinder::G4VRevisedChordFinder(G4VIntegrationDriver *pIntDriver, G4int  statisticsVerbosity):
    fIntDriver(pIntDriver),
    fverb(statisticsVerbosity),
    fDefaultDeltaChord( 0.25 * mm ),      // Parameters
    fDeltaChord( fDefaultDeltaChord ),    //   Internal parameters
    fFirstFraction(0.999), fFractionLast(1.00),  fFractionNextEstimate(0.98),
    fMultipleRadius(15.0),
    fStatsVerbose(0),
    fLastStepEstimate_Unconstrained(DBL_MAX),          // Should move q, p to
    fTotalNoTrials_FNC(0), fNoCalls_FNC(0), fmaxTrials_FNC(0)

{
    // check the values and set the other parameters
    SetFractions_Last_Next( fFractionLast, fFractionNextEstimate);
}

G4VRevisedChordFinder::~G4VRevisedChordFinder()
{
    delete fIntDriver;
    if( fStatsVerbose ) { PrintStatistics(); }
}


void G4VRevisedChordFinder::SetFractions_Last_Next( G4double fractLast, G4double fractNext )
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
    G4cerr << "G4VRevisedChordFinder::SetFractions_Last_Next: Invalid "
           << " fraction Last = " << fractLast
           << " must be  0 <  fractionLast <= 1 " << G4endl;
  }
  if( (fractNext > 0.0) && (fractNext <1.0) )
  {
    fFractionNextEstimate = fractNext;
  }
  else
  {
    G4cerr << "G4VRevisedChordFinder:: SetFractions_Last_Next: Invalid "
           << " fraction Next = " << fractNext
           << " must be  0 <  fractionNext < 1 " << G4endl;
  }
}


G4double G4VRevisedChordFinder::NewStep(G4double  stepTrialOld,
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



G4FieldTrack G4VRevisedChordFinder::ApproxCurvePointS( const G4FieldTrack&  CurveA_PointVelocity,
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

    DoStepForIntersection(EndPoint,test_step, eps_step);

#ifdef G4DEBUG_FIELD
    // Printing Brent and Linear Approximation
    //
    G4cout << "G4VRevisedChordFinder::ApproxCurvePointS() - test-step ShF = "
           << test_step << "  EndPoint = " << EndPoint << G4endl;

    //  Test Track
    //
    G4FieldTrack TestTrack( CurveA_PointVelocity);
    TestTrack = ApproxCurvePointV( CurveA_PointVelocity,
                                   CurveB_PointVelocity,
                                   CurrentE_Point, eps_step );
    G4cout.precision(14);
    G4cout << "G4VRevisedChordFinder::BrentApprox = " << EndPoint  << G4endl;
    G4cout << "G4VRevisedChordFinder::LinearApprox= " << TestTrack << G4endl;
#endif
  }
  return EndPoint;
}



G4FieldTrack G4VRevisedChordFinder::ApproxCurvePointV(const G4FieldTrack& trackPointA,
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

  G4double  integrationInaccuracyLimit = std::max(perMillion, 0.5*eps_step );
  //G4cout<<"curve_length: "<<curve_length<<" ABdist: "<<ABdist<<G4endl;
  if( curve_length < ABdist * (1. - integrationInaccuracyLimit) )
  {
#ifdef G4DEBUG_FIELD
    G4cerr << " Warning in G4VRevisedChordFinder::ApproxCurvePointV(): \n"
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
        G4Exception("G4VRevisedChordFinder::ApproxCurvePointV()", "GeomField0003",
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
     G4cout << "Warning in G4VRevisedChordFinder::ApproxCurvePointV():"
               " A and B are the same point! \n"
               " Chord AB length = " << vectorAE.mag() << "\n\n";
#endif
  }

  if(AE_fraction > 1.0 + perMillion)
  {
#ifdef G4DEBUG_FIELD
    G4cerr << " G4VRevisedChordFinder::ApproxCurvePointV() - Warning:"
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
  if ( AE_fraction > 0.0 )
  {
      //calculate approximation of the intersection point
      DoStepForIntersection(trackPointE, stepLength, eps_step);
  }

  return trackPointE;
}



void G4VRevisedChordFinder::PrintStatistics()
{
  // Print Statistics

  G4cout << "G4VRevisedChordFinder statistics report: " << G4endl;
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



void G4VRevisedChordFinder::TestChordPrint( G4int    noTrials,
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
