#ifndef TChordFinder_HH
#define TChordFinder_HH

#include "G4FieldTrack.hh"
#include <iomanip>
#include "G4SystemOfUnits.hh"

template
< class Driver>
class TChordFinder : public G4ChordFinder
{
      typedef Driver T_Driver;
      typedef typename T_Driver::T_Stepper T_Stepper;
      typedef typename T_Stepper::T_Equation T_Equation;
      typedef typename T_Equation::T_Field T_Field;

    public:  // with description

        // ..........................................................................

        TChordFinder(T_Driver* pIntegrationDriver)
            : fDefaultDeltaChord( 0.25 * mm ),      // Parameters
            fDeltaChord( fDefaultDeltaChord ),    //   Internal parameters
            fFirstFraction(0.999), fFractionLast(1.00),  fFractionNextEstimate(0.98), 
            fMultipleRadius(15.0), 
            fStatsVerbose(0),
            fDriversStepper(0),                    // Dependent objects 
            fAllocatedStepper(false),
            fEquation(0),      
            fTotalNoTrials_FNC(0), fNoCalls_FNC(0), fmaxTrials_FNC(0), 
            G4ChordFinder(pIntegrationDriver)
    {
        // Simple constructor -- it does not create equation
        fIntgrDriver= pIntegrationDriver;
        fAllocatedStepper= false;

        fLastStepEstimate_Unconstrained = DBL_MAX;          // Should move q, p to

        SetFractions_Last_Next( fFractionLast, fFractionNextEstimate);  
        // check the values and set the other parameters
    }


        // ..........................................................................

        TChordFinder( T_Field*        theMagField,
                G4double                stepMinimum, 
                T_Stepper* pItsStepper )
            : fDefaultDeltaChord( 0.25 * mm ),     // Constants 
            fDeltaChord( fDefaultDeltaChord ),   // Parameters
            fFirstFraction(0.999), fFractionLast(1.00),  fFractionNextEstimate(0.98), 
            fMultipleRadius(15.0), 
            fStatsVerbose(0),
            fDriversStepper(0),                  //  Dependent objects
            fAllocatedStepper(false),
            fEquation(0), 
            fTotalNoTrials_FNC(0), fNoCalls_FNC(0), fmaxTrials_FNC(0),  // State - stats
            G4ChordFinder(theMagField, stepMinimum, pItsStepper)
    {
        //  Construct the Chord Finder
        //  by creating in inverse order the  Driver, the Stepper and EqRhs ...

        T_Equation *pEquation = new T_Equation(theMagField);
        fEquation = pEquation;                            
        fLastStepEstimate_Unconstrained = DBL_MAX;          // Should move q, p to
        //    G4FieldTrack ??

        SetFractions_Last_Next( fFractionLast, fFractionNextEstimate);  
        // check the values and set the other parameters

        // --->>  Charge    Q = 0 
        // --->>  Momentum  P = 1       NOMINAL VALUES !!!!!!!!!!!!!!!!!!

        if( pItsStepper == 0 )
        { 
            pItsStepper = fDriversStepper = new T_Stepper(pEquation);
            fAllocatedStepper= true;
        }
        else
        {
            fAllocatedStepper= false; 
        }
        fIntgrDriver = new T_Driver(stepMinimum, pItsStepper, 
                pItsStepper->GetNumberOfVariables() );
    }


        // ......................................................................

        virtual ~TChordFinder()
        {
            delete   fEquation; // fIntgrDriver->pIntStepper->theEquation_Rhs;
            if( fAllocatedStepper)
            { 
                delete fDriversStepper; 
            }
            delete   fIntgrDriver; 

            if( fStatsVerbose ) { PrintStatistics(); }
        }


        // ......................................................................

        G4double 
            AdvanceChordLimited( G4FieldTrack& yCurrent,
                    G4double      stepMax,
                    G4double      epsStep,
                    const G4ThreeVector latestSafetyOrigin,
                    G4double       latestSafetyRadius )
            {
                G4double stepPossible;
                G4double dyErr;
                G4FieldTrack yEnd( yCurrent);
                G4double  startCurveLen= yCurrent.GetCurveLength();
                G4double nextStep;
                //            *************
                stepPossible= FindNextChord(yCurrent, stepMax, yEnd, dyErr, epsStep,
                        &nextStep, latestSafetyOrigin, latestSafetyRadius
                        );
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
                    good_advance = fIntgrDriver->AccurateAdvance(yCurrent, stepPossible,
                            epsStep, nextStep);
                    if ( ! good_advance )
                    { 
                        // In this case the driver could not do the full distance
                        stepPossible= yCurrent.GetCurveLength()-startCurveLen;
                    }
                }
                return stepPossible;
            }

        // ...........................................................................

        G4FieldTrack
            ApproxCurvePointS( const G4FieldTrack&  CurveA_PointVelocity, 
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

                const G4double tolerance= 1.e-12;
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

                    fIntgrDriver->AccurateAdvance(EndPoint,test_step, eps_step);

#ifdef G4DEBUG_FIELD
                    // Printing Brent and Linear Approximation
                    //
                    G4cout << "ApproxCurvePointS() - test-step ShF = "
                        << test_step << "  EndPoint = " << EndPoint << G4endl;

                    //  Test Track
                    //
                    G4FieldTrack TestTrack( CurveA_PointVelocity);
                    TestTrack = ApproxCurvePointV( CurveA_PointVelocity, 
                            CurveB_PointVelocity, 
                            CurrentE_Point, eps_step );
                    G4cout.precision(14);
                    G4cout << "BrentApprox = " << EndPoint  << G4endl;
                    G4cout << "LinearApprox= " << TestTrack << G4endl; 
#endif
                }
                return EndPoint;
            }


        // ...........................................................................

        G4FieldTrack 
            ApproxCurvePointV( const G4FieldTrack& CurveA_PointVelocity, 
                    const G4FieldTrack& CurveB_PointVelocity, 
                    const G4ThreeVector& CurrentE_Point,
                    G4double eps_step)
            {
                // If r=|AE|/|AB|, and s=true path lenght (AB)
                // return the point that is r*s along the curve!

                G4FieldTrack   Current_PointVelocity = CurveA_PointVelocity; 

                G4ThreeVector  CurveA_Point= CurveA_PointVelocity.GetPosition();
                G4ThreeVector  CurveB_Point= CurveB_PointVelocity.GetPosition();

                G4ThreeVector  ChordAB_Vector= CurveB_Point   - CurveA_Point;
                G4ThreeVector  ChordAE_Vector= CurrentE_Point - CurveA_Point;

                G4double       ABdist= ChordAB_Vector.mag();
                G4double  curve_length;  //  A curve length  of AB
                G4double  AE_fraction; 

                curve_length= CurveB_PointVelocity.GetCurveLength()
                    - CurveA_PointVelocity.GetCurveLength();  

                G4double  integrationInaccuracyLimit= std::max( perMillion, 0.5*eps_step ); 
                if( curve_length < ABdist * (1. - integrationInaccuracyLimit) )
                { 
#ifdef G4DEBUG_FIELD
                    G4cerr << " Warning in ApproxCurvePoint: "
                        << G4endl
                        << " The two points are further apart than the curve length "
                        << G4endl
                        << " Dist = "         << ABdist
                        << " curve length = " << curve_length 
                        << " relativeDiff = " << (curve_length-ABdist)/ABdist 
                        << G4endl;
                    if( curve_length < ABdist * (1. - 10*eps_step) )
                    {
                        std::ostringstream message;
                        message << "Unphysical curve length." << G4endl
                            << "The size of the above difference exceeds allowed limits."
                            << G4endl
                            << "Aborting.";
                        G4Exception("ApproxCurvePointV()", "GeomField0003",
                                FatalException, message);
                    }
#endif
                    // Take default corrective action: adjust the maximum curve length. 
                    // NOTE: this case only happens for relatively straight paths.
                    // curve_length = ABdist; 
                }

                G4double  new_st_length; 

                if ( ABdist > 0.0 )
                {
                    AE_fraction = ChordAE_Vector.mag() / ABdist;
                }
                else
                {
                    AE_fraction = 0.5;                         // Guess .. ?; 
#ifdef G4DEBUG_FIELD
                    G4cout << "Warning in ApproxCurvePointV():"
                        << " A and B are the same point!" << G4endl
                        << " Chord AB length = " << ChordAE_Vector.mag() << G4endl
                        << G4endl;
#endif
                }

                if( (AE_fraction> 1.0 + perMillion) || (AE_fraction< 0.) )
                {
#ifdef G4DEBUG_FIELD
                    G4cerr << " ApproxCurvePointV() - Warning:"
                        << " Anomalous condition:AE > AB or AE/AB <= 0 " << G4endl
                        << "   AE_fraction = " <<  AE_fraction << G4endl
                        << "   Chord AE length = " << ChordAE_Vector.mag() << G4endl
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

                new_st_length= AE_fraction * curve_length; 

                if ( AE_fraction > 0.0 )
                { 
                    fIntgrDriver->AccurateAdvance(Current_PointVelocity, 
                            new_st_length, eps_step );
                    //
                    // In this case it does not matter if it cannot advance the full distance
                }

                // If there was a memory of the step_length actually required at the start 
                // of the integration Step, this could be re-used ...

                G4cout.precision(14);

                return Current_PointVelocity;
            }


        // ......................................................................

        inline 
            void  SetIntegrationDriver(T_Driver* IntegrationDriver)
            {
                fIntgrDriver=IntegrationDriver;
            }

        inline
            T_Driver*  GetIntegrationDriver()
            {
                return fIntgrDriver;
            }


        inline
            G4double  GetDeltaChord() const
            {
                return fDeltaChord;
            }

        inline
            void  SetDeltaChord(G4double newval)
            {
                fDeltaChord=newval;
            }

        // ......................................................................


        inline
            void  ResetStepEstimate()
            {
                fLastStepEstimate_Unconstrained = DBL_MAX;    
            }

        // ......................................................................
        inline G4int  GetNoCalls()     { return fNoCalls_FNC; }
        inline G4int  GetNoTrials()    { return fTotalNoTrials_FNC; }
        inline G4int  GetNoMaxTrials() { return fmaxTrials_FNC; } 

        inline G4double  GetFirstFraction() { return fFirstFraction; } 
        inline G4double  GetFractionLast()  { return fFractionLast; } 
        inline G4double  GetFractionNextEstimate() 
        { return fFractionNextEstimate; } 
        inline G4double  GetMultipleRadius()
        { return fMultipleRadius; } 

        inline void  SetFirstFraction(G4double val){ fFirstFraction=val; }

        inline G4int  SetVerbose( G4int newvalue )
        { 
            G4int oldval= fStatsVerbose; 
            fStatsVerbose = newvalue;
            return oldval; 
        }


        // A  member that calculates the inverse parabolic through
        // the three points (x,y) and returns the value x that, for the
        // inverse parabolic, corresponds to y=0.
        //
        inline G4double  InvParabolic ( const G4double xa, const G4double ya,
                const G4double xb, const G4double yb,
                const G4double xc, const G4double yc )
        {   	const G4double R = yb/yc,
            S = yb/ya,
            T = ya/yc;
            const G4double Q = (T-1)*(R-1)*(S-1);
            if (std::fabs(Q) <DBL_MIN ) return  DBL_MAX;

            const G4double P = S*(T*(R-T)*(xc-xb) - (1-R)*(xb-xa));
            return xb + P/Q;
        }

        virtual void
            PrintStatistics()
            {
                // Print Statistics

                G4cout << "TChordFinder statistics report: " << G4endl;
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


        // ...........................................................................


        void   
            SetFractions_Last_Next( G4double fractLast=0.90, G4double fractNext=0.95 )
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
                    G4cerr << "SetFractions_Last_Next: Invalid "
                        << " fraction Last = " << fractLast
                        << " must be  0 <  fractionLast <= 1 " << G4endl;
                }
                if( (fractNext > 0.0) && (fractNext <1.0) )
                {
                    fFractionNextEstimate = fractNext;
                }
                else
                {
                    G4cerr << " SetFractions_Last_Next: Invalid "
                        << " fraction Next = " << fractNext
                        << " must be  0 <  fractionNext < 1 " << G4endl;
                }
            }




    public:  // without description

        void TestChordPrint( G4int    noTrials, 
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

    protected:   // .........................................................

        inline  
            void     AccumulateStatistics( G4int noTrials ) 
            {
                // Statistics 
                fTotalNoTrials_FNC += noTrials; 
                fNoCalls_FNC++; 
                // if( noTrials >= fmaxTrials_FNC ){
                if (noTrials > fmaxTrials_FNC ) { 
                    fmaxTrials_FNC=noTrials; 
                    // fnoTimesMaxTrFNC=0; 
                } else { 
                    // fnoTimesMaxTrFNC++; 
                } 
                // } 
            }

        inline
            G4bool  AcceptableMissDist(G4double dChordStep) const
            { 
                return (dChordStep <= fDeltaChord) ;
            }

        // ...........................................................................

        G4double NewStep(G4double  stepTrialOld, 
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

        // ............................................................................
        virtual 
            G4double
            FindNextChord( const  G4FieldTrack& yStart,
                    G4double     stepMax,
                    G4FieldTrack&   yEnd, // Endpoint
                    G4double&   dyErrPos, // Error of endpoint
                    G4double    epsStep,
                    G4double*  pStepForAccuracy, 
                    const  G4ThreeVector, //  latestSafetyOrigin,
                    G4double       //  latestSafetyRadius 
                    )
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
                fIntgrDriver-> GetDerivatives( yCurrent, dydx );
                G4int     noTrials=0;
                const G4double safetyFactor= fFirstFraction; //  0.975 or 0.99 ? was 0.999

                stepTrial = std::min( stepMax, safetyFactor*fLastStepEstimate_Unconstrained );
                G4double newStepEst_Uncons= 0.0; 
                do
                { 
                    G4double stepForChord;  
                    yCurrent = yStart;    // Always start from initial point

                    //            ************
                    fIntgrDriver->QuickAdvance( yCurrent, dydx, stepTrial, 
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
                            stepTrial = std::min( stepForChord, fFractionLast * stepTrial);
                        }
                        else
                        {
                            stepTrial *= 0.1;
                        }
                    }
                    noTrials++; 
                }
                while( ! validEndPoint );   // End of do-while  RKD 

                if( newStepEst_Uncons > 0.0  )
                {
                    fLastStepEstimate_Unconstrained= newStepEst_Uncons;
                }

                AccumulateStatistics( noTrials );

                if( pStepForAccuracy )
                { 
                    // Calculate the step size required for accuracy, if it is needed
                    //
                    G4double dyErr_relative = dyErrPos/(epsStep*lastStepLength);
                    if( dyErr_relative > 1.0 )
                    {
                        stepForAccuracy = fIntgrDriver->ComputeNewStepSize( dyErr_relative,
                                lastStepLength );
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
                yEnd=  yCurrent;  
                return stepTrial; 
            }


        // ...........................................................................


        void     PrintDchordTrial(G4int     noTrials, 
                G4double  stepTrial, 
                G4double  oldStepTrial, 
                G4double  dChordStep);

        inline 
            G4double  GetLastStepEstimateUnc()
            {
                return fLastStepEstimate_Unconstrained;   
            } 
        inline 
            void  SetLastStepEstimateUnc( G4double stepEst )
            {
                fLastStepEstimate_Unconstrained = stepEst;    
            }

    private:  // ............................................................

        TChordFinder(const TChordFinder&);
        TChordFinder& operator=(const TChordFinder&);
        // Private copy constructor and assignment operator.

    private:  // ............................................................
        // G4int    nOK, nBAD;

        // Constants
        const G4double fDefaultDeltaChord;  // SET in TChordFinder.cc = 0.25 mm

        //  PARAMETERS 
        //  ---------------------
        G4double  fDeltaChord;               //  Maximum miss distance 
        //    Internal parameters
        G4double  fFirstFraction, fFractionLast, fFractionNextEstimate;
        G4double  fMultipleRadius; 
        G4int     fStatsVerbose;  // if > 0, print Statistics in destructor

        //  DEPENDENT Objects
        //  ---------------------
        T_Driver*        fIntgrDriver;
        T_Stepper* fDriversStepper; 
        G4bool                  fAllocatedStepper;  // Bookkeeping of dependent object
        T_Equation*     fEquation; 

        //  STATE information
        //  --------------------
        G4double    fLastStepEstimate_Unconstrained;
        //  State information for efficiency

        // For Statistics
        // -- G4int   fNoTrials, fNoCalls;
        G4int   fTotalNoTrials_FNC,  fNoCalls_FNC, fmaxTrials_FNC; // fnoTimesMaxTrFNC; 
};

// Inline function implementation:

#endif  // TChordFinder_HH
