//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: FSALMagIntegratorDriver.cc
//
// 
//
// Implementation for class FSALMagInt_Driver
// Tracking in space dependent magnetic field
//
// History of major changes: /To be filled/

#include "G4FSALIntegrationDriver.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4GeometryTolerance.hh"
#include "G4FieldTrack.hh"
#include "Utils.hh"

#include <iomanip>

using namespace magneticfield;

namespace {
const G4int NCOMP = G4FieldTrack::ncompSVEC;

} // namespace

//  Stepsize can increase by no more than 5.0
//           and decrease by no more than 1/10. = 0.1
//
const G4double G4FSALIntegrationDriver::max_stepping_increase = 5.0;
const G4double G4FSALIntegrationDriver::max_stepping_decrease = 0.1;

//  The (default) maximum number of steps is Base
//  divided by the order of Stepper
//
const G4int  G4FSALIntegrationDriver::fMaxStepBase = 100000;  // Was 5000, was 250

#ifndef G4NO_FIELD_STATISTICS
#define G4FLD_STATS  1
#endif

#ifndef G4DEBUG_FIELD
#define G4DEBUG_FIELD 1
#endif


G4FSALIntegrationDriver::G4FSALIntegrationDriver(
    G4double hminimum,
    G4VFSALIntegrationStepper* pStepper,
    G4int numComponents,
    G4int statisticsVerbose):
    fSmallestFraction( 1.0e-12 ),
    fNoIntegrationVariables(numComponents),
    fMinNoVars(12),
    fNoVars(std::max(fNoIntegrationVariables, fMinNoVars)),
    fStatisticsVerboseLevel(statisticsVerbose),
    fNoTotalSteps(0),
    fNoBadSteps(0),
    fNoSmallSteps(0),
    fNoInitialSmallSteps(0),
    fDyerr_max(0.0), fDyerr_mx2(0.0),
    fDyerrPos_smTot(0.0),
    fDyerrPos_lgTot(0.0),
    fDyerrVel_lgTot(0.0),
    fSumH_sm(0.0),
    fSumH_lg(0.0),
    fVerboseLevel(0),
    TotalNoStepperCalls(0),
    fObserver(fVerboseLevel, hminimum, numComponents)
{
    // In order to accomodate "Laboratory Time", which is [7], fMinNoVars=8
    // is required. For proper time of flight and spin,  fMinNoVars must be 12
    
    RenewStepperAndAdjust(pStepper);
    fMinimumStep = hminimum;
    fMaxNoSteps = fMaxStepBase / pIntStepper->IntegratorOrder();
#ifdef G4DEBUG_FIELD
    fVerboseLevel = 2;
#endif
    
    if (fVerboseLevel > 0 || fStatisticsVerboseLevel > 1) {
        G4cout << "MagIntDriver version: Accur-Adv: "
        << "invE_nS, QuickAdv-2sqrt with Statistics "
#ifdef G4FLD_STATS
        << " enabled "
#else
        << " disabled "
#endif
        << G4endl;
    }
}

G4FSALIntegrationDriver::~G4FSALIntegrationDriver()
{
    if (fStatisticsVerboseLevel > 1) {
        //fObserver.PrintStatisticsReport();
        PrintStatisticsReport();
    }
}

// To add much printing for debugging purposes, uncomment the following
// and set verbose level to 1 or higher value !
// #define  G4DEBUG_FIELD 1


// Runge-Kutta driver with adaptive stepsize control. Integrate starting
// values at y_current over hstep x2 with accuracy eps.
// On output ystart is replaced by values at the end of the integration
// interval. RightHandSide is the right-hand side of ODE system.
// The source is similar to odeint routine from NRC p.721-722 .
G4bool G4FSALIntegrationDriver::AccurateAdvance(
    G4FieldTrack& y_current,
    G4double hstep,
    G4double eps,
    G4double hinitial)
{
    G4int nstp, i, no_warnings=0;
    G4double x, hnext, hdid, h;
    
#ifdef G4DEBUG_FIELD
    static G4int dbg = 10;
    static G4int nStpPr = 50;   // For debug printing of long integrations
    G4double ySubStepStart[NCOMP];
    G4FieldTrack  yFldTrkStart(y_current);
#endif
    
    G4double y[NCOMP], dydx[NCOMP];
    G4double ystart[NCOMP], yEnd[NCOMP];
    G4double  x1, x2;
    G4bool succeeded = true, lastStepSucceeded;
    
    G4double startCurveLength;
    
    G4int  noFullIntegr=0, noSmallIntegr = 0 ;
    static G4ThreadLocal G4int  noGoodSteps =0 ;  // Bad = chord > curve-len
    const  G4int  nvar= fNoVars;
    
    G4FieldTrack yStartFT(y_current);
    
    //  Ensure that hstep > 0
    //
    if( hstep <= 0.0 )
    {
        if(hstep == 0.0)
        {
            std::ostringstream message;
            message << "Proposed step is zero; hstep = " << hstep << " !";
            G4Exception("G4FSALIntegrationDriver::AccurateAdvance()",
                        "GeomField1001", JustWarning, message);
            return succeeded;
        }
        else
        {
            std::ostringstream message;
            message << "Invalid run condition." << G4endl
            << "Proposed step is negative; hstep = " << hstep << "." << G4endl
            << "Requested step cannot be negative! Aborting event.";
            G4Exception("G4FSALIntegrationDriver::AccurateAdvance()",
                        "GeomField0003", EventMustBeAborted, message);
            return false;
        }
    }
    
    y_current.DumpToArray( ystart );
    
    startCurveLength= y_current.GetCurveLength();
    x1= startCurveLength;
    x2= x1 + hstep;
    
    if ( (hinitial > 0.0) && (hinitial < hstep)
        && (hinitial > perMillion * hstep) )
    {
        h = hinitial;
    }
    else  //  Initial Step size "h" defaults to the full interval
    {
        h = hstep;
    }
    
    x = x1;
    
    for (i=0;i<nvar;i++) {
        y[i] = ystart[i];
    }
    
    G4bool lastStep = false;
    nstp = 1;
    
    pIntStepper->ComputeRightHandSide(y, dydx);

    
    do
    {
        G4ThreeVector StartPos( y[0], y[1], y[2] );
        
#ifdef G4DEBUG_FIELD
        G4double xSubStepStart = x;
        for (i = 0; i < nvar; i++)  { ySubStepStart[i] = y[i]; }
        yFldTrkStart.LoadFromArray(y, fNoIntegrationVariables);
        yFldTrkStart.SetCurveLength(x);
#endif
        
        
        
        
        fNoTotalSteps++;
        
        // Perform the Integration
        //
        if( h > fMinimumStep )
        {
            OneGoodStep(y,dydx,x,h,eps,hdid,hnext) ;
            //--------------------------------------
            lastStepSucceeded= (hdid == h);
#ifdef G4DEBUG_FIELD
            if (dbg>2) {
                fObserver.PrintStatus( ySubStepStart, xSubStepStart, y, x, h,  nstp); // Only
            }
#endif
        }
        else
        {
            G4FieldTrack yFldTrk( G4ThreeVector(0,0,0),
                                 G4ThreeVector(0,0,0), 0., 0., 0., 0. );
            G4double dchord_step, dyerr, dyerr_len;   // What to do with these ?
            yFldTrk.LoadFromArray(y, fNoIntegrationVariables);
            yFldTrk.SetCurveLength( x );
            
            QuickAdvance( yFldTrk, dydx, h, dchord_step, dyerr_len );
            //-----------------------------------------------------
            
            yFldTrk.DumpToArray(y);

#ifdef G4FLD_STATS
            fNoSmallSteps++;
            if ( dyerr_len > fDyerr_max)  { fDyerr_max= dyerr_len; }
            fDyerrPos_smTot += dyerr_len;
            fSumH_sm += h;  // Length total for 'small' steps
            if (nstp<=1)  { fNoInitialSmallSteps++; }
#endif
#ifdef G4DEBUG_FIELD
            if (dbg > 1) {
                if (fNoSmallSteps < 2) {
                    fObserver.PrintStatus(ySubStepStart, x1, y, x, h, -nstp);
                }
                G4cout << "Another sub-min step, no " << fNoSmallSteps
                << " of " << fNoTotalSteps << " this time " << nstp << G4endl;
                fObserver.PrintStatus(ySubStepStart, x1, y, x, h,  nstp);   // Only this
                G4cout << " dyerr= " << dyerr_len << " relative = " << dyerr_len / h
                << " epsilon= " << eps << " hstep= " << hstep
                << " h= " << h << " hmin= " << fMinimumStep << G4endl;
            }
#endif
            if( h == 0.0 )
            {
                G4Exception("G4FSALIntegrationDriver::AccurateAdvance()",
                            "GeomField0003", FatalException,
                            "Integration Step became Zero!");
            }
            dyerr = dyerr_len / h;
            hdid= h;
            x += hdid;
            
            // Compute suggested new step
            hnext= ComputeNewStepSize( dyerr/eps, h);
            
            // .. hnext= ComputeNewStepSize_WithinLimits( dyerr/eps, h);
            lastStepSucceeded= (dyerr<= eps);
        }
        
        if (lastStepSucceeded)  { noFullIntegr++; }
        else                    { noSmallIntegr++; }
        
        G4ThreeVector EndPos( y[0], y[1], y[2] );
        
#ifdef  G4DEBUG_FIELD
        if( (dbg>0) && (dbg<=2) && (nstp>nStpPr))
        {
            if( nstp==nStpPr )  { G4cout << "***** Many steps ****" << G4endl; }
            G4cout << "MagIntDrv: " ;
            G4cout << "hdid="  << std::setw(12) << hdid  << " "
            << "hnext=" << std::setw(12) << hnext << " "
            << "hstep=" << std::setw(12) << hstep << " (requested) "
            << G4endl;
            fObserver.PrintStatus( ystart, x1, y, x, h, (nstp==nStpPr) ? -nstp: nstp);
        }
#endif
        
        // Check the endpoint
        G4double endPointDist= (EndPos-StartPos).mag();
        if ( endPointDist >= hdid*(1.+perMillion) )
        {
            fNoBadSteps++;
            
            // Issue a warning only for gross differences -
            // we understand how small difference occur.
            if ( endPointDist >= hdid*(1.+perThousand) )
            {
#ifdef G4DEBUG_FIELD
                if (dbg)
                {
                    fObserver.WarnEndPointTooFar ( endPointDist, hdid, eps, dbg );
                    G4cerr << "  Total steps:  bad " << fNoBadSteps
                    << " good " << noGoodSteps << " current h= " << hdid
                    << G4endl;
                    fObserver.PrintStatus( ystart, x1, y, x, hstep, no_warnings?nstp:-nstp);
                }
#endif
                no_warnings++;
            }
        }
        else
        {
            noGoodSteps ++;
        }
        // #endif
        
        //  Avoid numerous small last steps
        if( (h < eps * hstep) || (h < fSmallestFraction * startCurveLength) )
        {
            // No more integration -- the next step will not happen
            lastStep = true;
        }
        else
        {
            // Check the proposed next stepsize
            if(std::fabs(hnext) <= Hmin())
            {
#ifdef  G4DEBUG_FIELD
                // If simply a very small interval is being integrated, do not warn
                if( (x < x2 * (1-eps) ) &&        // The last step can be small: OK
                   (std::fabs(hstep) > Hmin()) ) // and if we are asked, it's OK
                {
                    if(dbg>0)
                    {
                        fObserver.WarnSmallStepSize( hnext, hstep, h, x-x1, nstp );
                        fObserver.PrintStatus( ystart, x1, y, x, hstep, no_warnings?nstp:-nstp);
                    }
                    no_warnings++;
                }
#endif
                // Make sure that the next step is at least Hmin.
                h = Hmin();
            }
            else
            {
                h = hnext;
            }
            
            //  Ensure that the next step does not overshoot
            if ( x+h > x2 )
            {                // When stepsize overshoots, decrease it!
                h = x2 - x ;   // Must cope with difficult rounding-error
            }                // issues if hstep << x2
            
            if ( h == 0.0 )
            {
                // Cannot progress - accept this as last step - by default
                lastStep = true;
#ifdef G4DEBUG_FIELD
                if (dbg>2)
                {
                    int prec= G4cout.precision(12);
                    G4cout << "Warning: FSALMagIntegratorDriver::AccurateAdvance"
                    << G4endl
                    << "  Integration step 'h' became "
                    << h << " due to roundoff. " << G4endl
                    << " Calculated as difference of x2= "<< x2 << " and x=" << x
                    << "  Forcing termination of advance." << G4endl;
                    G4cout.precision(prec);
                }
#endif
            }
        }
    } while ( ((nstp++)<=fMaxNoSteps) && (x < x2) && (!lastStep) );
    // Loop checking, 07.10.2016, J. Apostolakis

    // Have we reached the end ?
    // --> a better test might be x-x2 > an_epsilon
    
    succeeded=  (x>=x2);  // If it was a "forced" last step
    
    for (i=0;i<nvar;i++)  { yEnd[i] = y[i]; }
    
    // Put back the values.
    y_current.LoadFromArray( yEnd, fNoIntegrationVariables );
    y_current.SetCurveLength( x );
    
    if(nstp > fMaxNoSteps)
    {
        no_warnings++;
        succeeded = false;
#ifdef G4DEBUG_FIELD
        if (dbg)
        {
            fObserver.WarnTooManyStep( x1, x2, x );  //  Issue WARNING
            fObserver.PrintStatus( yEnd, x1, y, x, hstep, -nstp);
        }
#endif
    }
    
#ifdef G4DEBUG_FIELD
    if( dbg && no_warnings )
    {
        G4cerr << "FSALMagIntegratorDriver exit status: no-steps " << nstp <<G4endl;
        fObserver.PrintStatus( yEnd, x1, y, x, hstep, nstp);
    }
#endif
    
    return succeeded;
}  // end of AccurateAdvance ...........................

G4double G4FSALIntegrationDriver::ShrinkStep(G4double error, G4double hstep)
{
    G4double htemp = GetSafety() * hstep * std::pow(error, GetPshrnk());
    return std::max(htemp, 0.1 * hstep);
}

G4double G4FSALIntegrationDriver::GrowStep(G4double error, G4double hstep)
{
    G4double hnext;
    if (error > errcon) {
        hnext = GetSafety() * hstep * std::pow(error, GetPgrow());
    } else {
        hnext = max_stepping_increase * hstep;
    }

    return hnext;
}

void G4FSALIntegrationDriver::OneGoodStep(G4double y[],             // InOut
                                  G4double dydx[],
                                  G4double& trackLength,   // InOut
                                  G4double htry,
                                  G4double eps_rel_max,
                                  G4double& hdid,          // Out
                                  G4double& hnext )        // Out
{
    G4double yerr[NCOMP], ytemp[NCOMP], dydxOut[NCOMP];

    static G4ThreadLocal G4int tot_no_trials = 0;
    const G4int maxTrials = 100;
    // Set stepsize to the initial trial value
    G4double hstep = htry;
    G4double error;
    for (G4int iter = 0; iter < maxTrials; ++iter) {
        ++tot_no_trials;
        hstep = std::max(hstep, fMinimumStep);
        pIntStepper->Stepper(y, dydx, hstep, ytemp, yerr, dydxOut);

        error = relativeError(y, yerr, hstep, eps_rel_max, fObserver);

         // Step succeeded
        if (error <= 1.0) {
            break;
        }

        // Step failed, compute the size of retrial step
        hstep = ShrinkStep(error, hstep);
    }

    // Compute size of next step
    hnext = GrowStep(error, hstep);
    trackLength += (hdid = hstep);

    for(G4int k = 0; k < fNoIntegrationVariables; ++k) {
        y[k] = ytemp[k];
        dydx[k] = dydxOut[k];
    }
}


// QuickAdvance just tries one Step - it does not ensure accuracy
//
G4bool  G4FSALIntegrationDriver::QuickAdvance(
                                      G4FieldTrack& y_posvel,         // INOUT
                                      G4double     dydx[],
                                      G4double     hstep,       // In
                                      G4double&    dchord_step,
                                      G4double&    dyerr_pos_sq,
                                      G4double&    dyerr_mom_rel_sq )
{
    G4Exception("G4FSALIntegrationDriver::QuickAdvance()", "GeomField0001",
                FatalException, "Not yet implemented.");
    
    // Use the parameters of this method, to please compiler
    dchord_step = dyerr_pos_sq = hstep * hstep * dydx[0];
    dyerr_mom_rel_sq = y_posvel.GetPosition().mag2();
    return true;
}

G4bool  G4FSALIntegrationDriver::QuickAdvance(G4FieldTrack& track,
                                              const G4double dydx[],
                                              G4double hstep,
                                              G4double& dchord_step,
                                              G4double& dyerr)
{
    G4double yError[NCOMP], yIn[NCOMP], yOut[NCOMP], dydxOut[NCOMP];
    
    static G4ThreadLocal G4int no_call = 0;
    ++no_call;
    
    track.DumpToArray(yIn);
    G4double curveLengthIn = track.GetCurveLength();
    
    pIntStepper->Stepper(yIn, dydx, hstep, yOut, yError, dydxOut);
    dchord_step = pIntStepper->DistChord();
    
    track.LoadFromArray(yOut, fNoIntegrationVariables);
    track.SetCurveLength(curveLengthIn + hstep);
    
#ifdef  G4DEBUG_FIELD
    if (fVerboseLevel > 2) {
        G4cout << "G4MagIntDrv: Quick Advance" << G4endl;
        fObserver.PrintStatus(yIn, curveLengthIn, yOut, curveLengthIn + hstep, hstep,  1);
    }
#endif
    

    // A single measure of the error
    //      TO-DO :  account for  energy,  spin, ... ?
    G4double positionError2 = extractValue2(yError, Value3D::Position);

    G4double momentumError2 = extractValue2(yError, Value3D::Momentum);
    G4double momentum2 = extractValue2(yOut, Value3D::Momentum);
    G4double relativeMomentumError2 = momentumError2 / momentum2;

    // Calculate also the change in the momentum squared also ???
    // G4double veloc_square = track.GetVelocity().mag2();
    // ...
    
#ifdef RETURN_A_NEW_STEP_LENGTH
    // The following step cannot be done here because "eps" is not known.
    dyerr_len = std::sqrt( dyerr_len_sq );
    dyerr_len_sq /= eps ;
    
    // Look at the velocity deviation ?
    //  sqr(yError[3])+sqr(yError[4])+sqr(yError[5]));
    
    // Set suggested new step
    hstep = ComputeNewStepSize( dyerr_len, hstep);
#endif
    
    dyerr = std::max(positionError2, relativeMomentumError2 * sqr(hstep));
    dyerr = std::sqrt(dyerr);
    
    return true;
}

void G4FSALIntegrationDriver::GetDerivatives(const G4FieldTrack& track,
                                             G4double dydx[]) const
{
    G4double  y[NCOMP];
    track.DumpToArray(y);
    pIntStepper->RightHandSide(y, dydx);
}

void G4FSALIntegrationDriver::SetVerboseLevel(G4int level)
{
    fVerboseLevel = level;
}

G4int G4FSALIntegrationDriver::GetVerboseLevel() const
{
    return fVerboseLevel;
}

void G4FSALIntegrationDriver::SetEquationOfMotion(G4EquationOfMotion *equation)
{
    pIntStepper->SetEquationOfMotion(equation);
}

G4EquationOfMotion* G4FSALIntegrationDriver::GetEquationOfMotion()
{
    return pIntStepper->GetEquationOfMotion();
}

// --------------------------------------------------------------------------

#ifdef QUICK_ADV_ARRAY_IN_AND_OUT
G4bool  G4FSALIntegrationDriver::QuickAdvance(
                                      G4double     yarrin[],    // In
                                      const G4double     dydx[],
                                      G4double     hstep,       // In
                                      G4double     yarrout[],
                                      G4double&    dchord_step,
                                      G4double&    dyerr )      // In length
{
    G4Exception("G4FSALIntegrationDriver::QuickAdvance()", "GeomField0001",
                FatalException, "Not yet implemented.");
    dyerr = dchord_step = hstep * yarrin[0] * dydx[0];
    yarrout[0]= yarrin[0];
}
#endif

// --------------------------------------------------------------------------

//  This method computes new step sizes - but does not limit changes to
//   within  certain factors
G4double G4FSALIntegrationDriver::ComputeNewStepSize(G4double  errMaxNorm,    // max error  (normalised)
                                                     G4double  hstepCurrent)  // current step size
{
    if (errMaxNorm > 1) { // Step failed; compute the size of retrial Step.
        return GetSafety() * hstepCurrent * std::pow(errMaxNorm, GetPshrnk());
    } else if (errMaxNorm > 0) { // Compute size of next Step for a successful step
        return GetSafety() * hstepCurrent * std::pow(errMaxNorm, GetPgrow());
    }

    // if error estimate is zero (possible) or negative (dubious)
    return max_stepping_increase * hstepCurrent;
}

// ---------------------------------------------------------------------------

// This method computes new step sizes limiting changes within certain factors
//
// It shares its logic with AccurateAdvance.
// They are kept separate currently for optimisation.
//
G4double G4FSALIntegrationDriver::ComputeNewStepSize_WithinLimits(
    G4double  errMaxNorm,    // max error  (normalised)
    G4double  hstepCurrent)  // current step size
{
    G4double hnew;
    
    // Compute size of next Step for a failed step
    if (errMaxNorm > 1.0) {
        // Step failed; compute the size of retrial Step.
        hnew = GetSafety() * hstepCurrent * std::pow(errMaxNorm,GetPshrnk()) ;
        
        if (hnew < max_stepping_decrease * hstepCurrent) {
            hnew = max_stepping_decrease * hstepCurrent;
            // reduce stepsize, but no more
            // than this factor (value= 1/10)
        }
    } else {
        // Compute size of next Step for a successful step
        if (errMaxNorm > errcon) {
            hnew = GetSafety() * hstepCurrent * std::pow(errMaxNorm, GetPgrow());
        } else { // No more than a factor of 5 increase
            hnew = max_stepping_increase * hstepCurrent;
        }
    }

    return hnew;
}

void G4FSALIntegrationDriver::PrintStatisticsReport()
{
    G4int noPrecBig= 6;
    G4int oldPrec= G4cout.precision(noPrecBig);
    
    G4cout << "G4FSALIntegrationDriver Statistics of steps undertaken. " << G4endl;
    G4cout << "G4FSALIntegrationDriver: Number of Steps: "
    << " Total= " <<  fNoTotalSteps
    << " Bad= "   <<  fNoBadSteps 
    << " Small= " <<  fNoSmallSteps 
    << " Non-initial small= " << (fNoSmallSteps - fNoInitialSmallSteps)
    << G4endl;
    
#ifdef G4FLD_STATS
    G4cout << "MID dyerr: " 
    << " maximum= " << fDyerr_max 
    << " Sum small= " << fDyerrPos_smTot 
    << " std::sqrt(Sum large^2): pos= " << std::sqrt(fDyerrPos_lgTot)
    << " vel= " << std::sqrt( fDyerrVel_lgTot )
    << " Total h-distance: small= " << fSumH_sm 
    << " large= " << fSumH_lg
    << G4endl;
    
#if 0
    G4int noPrecSmall=4; 
    // Single line precis of statistics ... optional
    G4cout.precision(noPrecSmall);
    G4cout << "MIDnums: " << fMinimumStep
    << "   " << fNoTotalSteps 
    << "  "  <<  fNoSmallSteps
    << "  "  << fNoSmallSteps-fNoInitialSmallSteps
    << "  "  << fNoBadSteps         
    << "   " << fDyerr_max
    << "   " << fDyerr_mx2 
    << "   " << fDyerrPos_smTot 
    << "   " << fSumH_sm
    << "   " << fDyerrPos_lgTot
    << "   " << fDyerrVel_lgTot
    << "   " << fSumH_lg
    << G4endl;
#endif 
#endif 
    
    G4cout.precision(oldPrec);
}

// ---------------------------------------------------------------------------

void G4FSALIntegrationDriver::SetSmallestFraction(G4double newFraction)
{
    if (newFraction > 1.e-16 && newFraction < 1e-8) {
        fSmallestFraction = newFraction;
    } else {
        G4cerr << "Warning: SmallestFraction not changed. " << G4endl
        << "  Proposed value was " << newFraction << G4endl
        << "  Value must be between 1.e-8 and 1.e-16" << G4endl;
    }
}

