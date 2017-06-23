#include "G4IntegrationObserver.hh"

#include "globals.hh"
#include "G4GeometryTolerance.hh"

#include <iomanip>

G4IntegrationObserver::G4IntegrationObserver(
    G4int verboseLevel,
    G4double minimumStepLength,
    G4double numberOfIntegrationVariables):
    fVerboseLevel(verboseLevel),
    fMinimumStep(minimumStepLength),
    fNoIntegrationVariables(numberOfIntegrationVariables),
    fNoTotalSteps(0),
    fNoBadSteps(0),
    fNoSmallSteps(0)
{
}

G4IntegrationObserver::~G4IntegrationObserver()
{
}

void G4IntegrationObserver::WarnSmallStepSize(
    G4double hnext,
    G4double hstep,
    G4double h,
    G4double xDone,
    G4int stepNumber)
{
    static G4ThreadLocal G4int noWarningsIssued = 0;
    const  G4int maxNoWarnings =  10;   // Number of verbose warnings
    std::ostringstream message;
    if ( noWarningsIssued < maxNoWarnings || fVerboseLevel > 10) {
        message << "The stepsize for the next iteration, " << hnext
        << ", is too small - in Step number " << stepNumber << "." << G4endl
        << "The minimum for the driver is " << fMinimumStep  << G4endl
        << "Requested integr. length was " << hstep << " ." << G4endl
        << "The size of this sub-step was " << h     << " ." << G4endl
        << "The integrations has already gone " << xDone;
    } else {
        message << "Too small 'next' step " << hnext
        << ", step-no: " << stepNumber << G4endl
        << ", this sub-step: " << h
        << ",  req_tot_len: " << hstep
        << ", done: " << xDone << ", min: " << fMinimumStep;
    }
    G4Exception("G4FSALIntegrationDriver::WarnSmallStepSize()", "GeomField1001",
                JustWarning, message);
    ++noWarningsIssued;
}

void G4IntegrationObserver::WarnTooManyStep(
    G4double x1start,
    G4double x2end,
    G4double xCurrent)
{
    std::ostringstream message;
    message << "The number of steps used in the Integration driver"
    << " (Runge-Kutta) is too many." << G4endl
    << "Integration of the interval was not completed !" << G4endl
    << "Only a " << (xCurrent-x1start) * 100 / (x2end - x1start)
    << " % fraction of it was done.";
    G4Exception("G4FSALIntegrationDriver::WarnTooManyStep()", "GeomField1001",
                JustWarning, message);
}


void G4IntegrationObserver::WarnEndPointTooFar (
    G4double endPointDist,
    G4double h,
    G4double eps,
    G4int dbg)
{
    static G4ThreadLocal G4double maxRelError=0.0;
    G4bool isNewMax, prNewMax;

    isNewMax = endPointDist > (1.0 + maxRelError) * h;
    prNewMax = endPointDist > (1.0 + 1.05 * maxRelError) * h;
    if(isNewMax) {
        maxRelError = endPointDist / h - 1.0;
    }

    if (dbg && (h > G4GeometryTolerance::GetInstance()->GetSurfaceTolerance())
       && ( (dbg>1) || prNewMax || (endPointDist >= h*(1.+eps) ) ) )
    {
        static G4ThreadLocal G4int noWarnings = 0;
        std::ostringstream message;
        if( (noWarnings ++ < 10) || (dbg>2) )
        {
            message << "The integration produced an end-point which " << G4endl
            << "is further from the start-point than the curve length."
            << G4endl;
        }
        message << "  Distance of endpoints = " << endPointDist
        << ", curve length = " << h << G4endl
        << "  Difference (curveLen-endpDist)= " << (h - endPointDist)
        << ", relative = " << (h-endPointDist) / h
        << ", epsilon =  " << eps;
        G4Exception("G4FSALIntegrationDriver::WarnEndPointTooFar()", "GeomField1001",
                    JustWarning, message);
    }
}

void G4IntegrationObserver::PrintStatus(
    const G4double* StartArr,
    G4double xstart,
    const G4double* CurrentArr,
    G4double xcurrent,
    G4double requestStep,
    G4int subStepNo)
// Potentially add as arguments:
//                                 <dydx>           - as Initial Force
//                                 stepTaken(hdid)  - last step taken
//                                 nextStep (hnext) - proposal for size
{
    G4FieldTrack StartFT(G4ThreeVector(0,0,0),
                          G4ThreeVector(0,0,0), 0., 0., 0., 0. );
    G4FieldTrack  CurrentFT (StartFT);

    StartFT.LoadFromArray( StartArr, fNoIntegrationVariables);
    StartFT.SetCurveLength( xstart);
    CurrentFT.LoadFromArray( CurrentArr, fNoIntegrationVariables);
    CurrentFT.SetCurveLength( xcurrent );

    PrintStatus(StartFT, CurrentFT, requestStep, subStepNo );
}

void G4IntegrationObserver::PrintStatus(
    const G4FieldTrack&  StartFT,
    const G4FieldTrack&  CurrentFT,
    G4double requestStep,
    G4int subStepNo)
{
    G4int verboseLevel= fVerboseLevel;
    static G4ThreadLocal G4int noPrecision= 5;
    G4int oldPrec= G4cout.precision(noPrecision);
    // G4cout.setf(ios_base::fixed,ios_base::floatfield);

    const G4ThreeVector StartPosition=       StartFT.GetPosition();
    const G4ThreeVector StartUnitVelocity=   StartFT.GetMomentumDir();
    const G4ThreeVector CurrentPosition=     CurrentFT.GetPosition();
    const G4ThreeVector CurrentUnitVelocity= CurrentFT.GetMomentumDir();

    G4double  DotStartCurrentVeloc= StartUnitVelocity.dot(CurrentUnitVelocity);

    G4double step_len= CurrentFT.GetCurveLength() - StartFT.GetCurveLength();
    G4double subStepSize = step_len;

    if( (subStepNo <= 1) || (verboseLevel > 3) )
    {
        subStepNo = - subStepNo;        // To allow printing banner

        G4cout << std::setw( 6)  << " " << std::setw( 25)
        << " G4FSALIntegrationDriver: Current Position  and  Direction" << " "
        << G4endl;
        G4cout << std::setw( 5) << "Step#" << " "
        << std::setw( 7) << "s-curve" << " "
        << std::setw( 9) << "X(mm)" << " "
        << std::setw( 9) << "Y(mm)" << " "
        << std::setw( 9) << "Z(mm)" << " "
        << std::setw( 8) << " N_x " << " "
        << std::setw( 8) << " N_y " << " "
        << std::setw( 8) << " N_z " << " "
        << std::setw( 8) << " N^2-1 " << " "
        << std::setw(10) << " N(0).N " << " "
        << std::setw( 7) << "KinEner " << " "
        << std::setw(12) << "Track-l" << " "   // Add the Sub-step ??
        << std::setw(12) << "Step-len" << " "
        << std::setw(12) << "Step-len" << " "
        << std::setw( 9) << "ReqStep" << " "
        << G4endl;
    }

    if( (subStepNo <= 0) )
    {
        PrintStat_Aux( StartFT,  requestStep, 0.,
                      0,        0.0,         1.0);
        //*************
    }

    if( verboseLevel <= 3 )
    {
        G4cout.precision(noPrecision);
        PrintStat_Aux( CurrentFT, requestStep, step_len,
                      subStepNo, subStepSize, DotStartCurrentVeloc );
        //*************
    }

    else // if( verboseLevel > 3 )
    {
        //  Multi-line output

        // G4cout << "Current  Position is " << CurrentPosition << G4endl
        //    << " and UnitVelocity is " << CurrentUnitVelocity << G4endl;
        // G4cout << "Step taken was " << step_len
        //    << " out of PhysicalStep= " <<  requestStep << G4endl;
        // G4cout << "Final safety is: " << safety << G4endl;
        // G4cout << "Chord length = " << (CurrentPosition-StartPosition).mag()
        //        << G4endl << G4endl;
    }
    G4cout.precision(oldPrec);
}

void G4IntegrationObserver::PrintStat_Aux(
    const G4FieldTrack&  aFieldTrack,
    G4double requestStep,
    G4double step_len,
    G4int subStepNo,
    G4double subStepSize,
    G4double dotVeloc_StartCurr)
{
    const G4ThreeVector Position = aFieldTrack.GetPosition();
    const G4ThreeVector UnitVelocity = aFieldTrack.GetMomentumDir();

    if( subStepNo >= 0)
    {
        G4cout << std::setw( 5) << subStepNo << " ";
    }
    else
    {
        G4cout << std::setw( 5) << "Start" << " ";
    }
    G4double curveLen= aFieldTrack.GetCurveLength();
    G4cout << std::setw( 7) << curveLen;
    G4cout << std::setw( 9) << Position.x() << " "
    << std::setw( 9) << Position.y() << " "
    << std::setw( 9) << Position.z() << " "
    << std::setw( 8) << UnitVelocity.x() << " "
    << std::setw( 8) << UnitVelocity.y() << " "
    << std::setw( 8) << UnitVelocity.z() << " ";
    G4int oldprec= G4cout.precision(3);
    G4cout << std::setw( 8) << UnitVelocity.mag2()-1.0 << " ";
    G4cout.precision(6);
    G4cout << std::setw(10) << dotVeloc_StartCurr << " ";
    G4cout.precision(oldprec);
    G4cout << std::setw( 7) << aFieldTrack.GetKineticEnergy();
    G4cout << std::setw(12) << step_len << " ";

    static G4ThreadLocal G4double oldCurveLength= 0.0;
    static G4ThreadLocal G4double oldSubStepLength= 0.0;
    static G4ThreadLocal G4int oldSubStepNo = -1;

    G4double subStep_len = 0.0;
    if( curveLen > oldCurveLength )
    {
        subStep_len= curveLen - oldCurveLength;
    }
    else if (subStepNo == oldSubStepNo)
    {
        subStep_len= oldSubStepLength;
    }
    oldCurveLength= curveLen;
    oldSubStepLength= subStep_len;

    G4cout << std::setw(12) << subStep_len << " ";
    G4cout << std::setw(12) << subStepSize << " ";
    if( requestStep != -1.0 )
    {
        G4cout << std::setw( 9) << requestStep << " ";
    }
    else
    {
        G4cout << std::setw( 9) << " InitialStep " << " ";
    }
    G4cout << G4endl;
}
/*
void G4IntegrationObserver::PrintStatisticsReport()
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
*/

