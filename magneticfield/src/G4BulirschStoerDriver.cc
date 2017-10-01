#include "G4BulirschStoerDriver.hh"

#include "G4LineSection.hh"
#include "G4SystemOfUnits.hh"
#include "G4FieldUtils.hh"

#define ncomp G4FieldTrack::ncompSVEC

G4BulirschStoerDriver::G4BulirschStoerDriver(
    G4double hminimum,
    G4EquationOfMotion* equation,
    G4int integratedComponents,
    G4int verb)
        : fEquation(equation)
        , fMinimumStep(hminimum)
        , fNumberOfComponents(integratedComponents)
        , fVerbosity(verb)
        , fMidpointMethod(equation, integratedComponents)
        , bulirschStoer(equation, integratedComponents, 0, DBL_MAX)
        , interval_sequence{2,4}
{
}

G4bool  G4BulirschStoerDriver::AccurateAdvance(
    G4FieldTrack& track,
    G4double hstep,
    G4double eps,
    G4double hinitial)
{
    G4int fNoTotalSteps = 0;
    G4int fMaxNoSteps = 10000;
    G4double fNoBadSteps = 0;
    G4double fSmallestFraction = 1.0e-12;

    // Driver with adaptive stepsize control. Integrate starting
    // values at y_current over hstep x2 with accuracy eps.
    // On output ystart is replaced by values at the end of the integration
    // interval. RightHandSide is the right-hand side of ODE system.
    // The source is similar to odeint routine from NRC p.721-722 .


    //  Ensure that hstep > 0
    if(hstep == 0)
    {
      std::ostringstream message;
      message << "Proposed step is zero; hstep = " << hstep << " !";
      G4Exception("G4MagInt_Driver::AccurateAdvance()",
                  "GeomField1001", JustWarning, message);

      return true;
    }
    if(hstep < 0)
    {
        std::ostringstream message;
        message << "Invalid run condition." << G4endl
                << "Proposed step is negative; hstep = " << hstep << "." << G4endl
                << "Requested step cannot be negative! Aborting event.";
        G4Exception("G4MagInt_Driver::AccurateAdvance()",
                    "GeomField0003", EventMustBeAborted, message);

        return false;
    }



    //init first step size
    G4double h;
    if ( (hinitial > 0) && (hinitial < hstep)
      && (hinitial > perMillion * hstep) )
    {
       h = hinitial;
    }
    else  //  Initial Step size "h" defaults to the full interval
    {
       h = hstep;
    }

    //integration variables
    track.DumpToArray(yCurrent);
    //copy non-integration variables to out array
    memcpy(yOut+fNumberOfComponents,
           yCurrent+fNumberOfComponents,
           sizeof(G4double)*(ncomp-fNumberOfComponents));

    G4double startCurveLength = track.GetCurveLength();
    G4double curveLength = startCurveLength;
    G4double endCurveLength = startCurveLength + hstep;

    //loop variables
    G4int nstp = 1, no_warnings = 0;
    G4double hnext, hdid;


    G4bool succeeded = true, lastStepSucceeded;

    G4int  noFullIntegr = 0, noSmallIntegr = 0 ;
    static G4ThreadLocal G4int  noGoodSteps = 0 ;  // Bad = chord > curve-len

    G4bool lastStep = false;

    //BulirschStoer.reset();

    G4FieldTrack yFldTrk(track);

    do
    {
        G4ThreeVector StartPos(yCurrent[0], yCurrent[1], yCurrent[2]);
        GetEquationOfMotion()->RightHandSide(yCurrent, dydxCurrent);
        fNoTotalSteps++;

        // Perform the Integration
        if(h == 0){
            G4Exception("G4MagInt_Driver::AccurateAdvance()",
                        "GeomField0003", FatalException,
                        "Integration Step became Zero!");
        }
        else if(h > fMinimumStep){
            //step size if Ok
            OneGoodStep(yCurrent,dydxCurrent,curveLength,h,eps,hdid,hnext);
            lastStepSucceeded = (hdid == h);
        }
        else{
            // for small steps call QuickAdvance for speed

            G4double dchord_step, dyerr, dyerr_len;   // What to do with these ?
            yFldTrk.LoadFromArray(yCurrent, ncomp);
            yFldTrk.SetCurveLength(curveLength);

            QuickAdvance(yFldTrk, dydxCurrent, h, dchord_step, dyerr_len);

            yFldTrk.DumpToArray(yCurrent);


            dyerr = dyerr_len / h;
            hdid = h;
            curveLength += hdid;

            // Compute suggested new step
            //hnext = ComputeNewStepSize(dyerr/eps, h);
            hnext = h;

            //hnext= ComputeNewStepSize_WithinLimits( dyerr/eps, h);
            lastStepSucceeded = (dyerr <= eps);
        }


        lastStepSucceeded ? ++noFullIntegr : ++noSmallIntegr;


        G4ThreeVector EndPos(yCurrent[0], yCurrent[1], yCurrent[2]);

        // Check the endpoint
        G4double endPointDist = (EndPos - StartPos).mag();
        if (endPointDist >= hdid*(1. + perMillion))
        {
            ++fNoBadSteps;

            // Issue a warning only for gross differences -
            // we understand how small difference occur.
            if (endPointDist >= hdid*(1.+perThousand))
            {
                ++no_warnings;
            }
        }
        else
        {
            ++noGoodSteps;
        }

        //  Avoid numerous small last steps
        if((h < eps * hstep) || (h < fSmallestFraction * startCurveLength))
        {
            // No more integration -- the next step will not happen
            lastStep = true;
        }
        else
        {
            // Check the proposed next stepsize
            if(std::fabs(hnext) < fMinimumStep)
            {
              // Make sure that the next step is at least Hmin.
              h = fMinimumStep;
            }
            else
            {
              h = hnext;
            }

            //  Ensure that the next step does not overshoot
            if (curveLength + h > endCurveLength)
            {
              h = endCurveLength - curveLength;
            }

            if (h == 0)
            {
              // Cannot progress - accept this as last step - by default
              lastStep = true;
            }
        }
    } while (((nstp++) <= fMaxNoSteps) && (curveLength < endCurveLength) && (!lastStep));
       // Have we reached the end ?
       // --> a better test might be x-x2 > an_epsilon

    succeeded = (curveLength >= endCurveLength);  // If it was a "forced" last step

    //copy integrated vars to output array
    memcpy(yOut,yCurrent,sizeof(G4double)*fNumberOfComponents);

    // upload new state
    track.LoadFromArray(yOut, ncomp);
    track.SetCurveLength(curveLength);

    if(nstp > fMaxNoSteps)
    {
      ++no_warnings;
      succeeded = false;
    }

    return succeeded;
}

G4bool G4BulirschStoerDriver::QuickAdvance(
    G4FieldTrack& track,
    const G4double dydx[],
    G4double hstep,
    G4double& missDist,
    G4double& dyerr)
{
    track.DumpToArray(yIn);
    const G4double curveLength = track.GetCurveLength();

    fMidpointMethod.SetSteps(interval_sequence[0]);
    fMidpointMethod.DoStep(yIn, dydx, yOut, hstep, yMid, derivs[0]);

    fMidpointMethod.SetSteps(interval_sequence[1]);
    fMidpointMethod.DoStep(yIn, dydx, yOut2, hstep, yMid2, derivs[1]);

    //extrapolation
    static const G4double coeff =
        1. / (sqr(static_cast<G4double>(interval_sequence[1]) /
                     static_cast<G4double>(interval_sequence[0])) - 1.);
    for (G4int i = 0; i < fNumberOfComponents; ++i){
        yOut[i] = yOut2[i] + (yOut2[i] - yOut[i]) * coeff;
        yMid[i] = yMid2[i] + (yMid2[i] - yMid[i]) * coeff;
    }

    //calc chord lenght
    const auto mid = field_utils::makeVector(yMid, field_utils::Value3D::Position);
    const auto in = field_utils::makeVector(yIn, field_utils::Value3D::Position);
    const auto out = field_utils::makeVector(yOut, field_utils::Value3D::Position);

    missDist = G4LineSection::Distline(mid, in, out);

    //calc error
    for (G4int i = 0; i < fNumberOfComponents; ++i){
        yError[i] = yOut[i] - yOut2[i];
    }

    G4double errPos2 = sqr(yError[0]) + sqr(yError[1]) + sqr(yError[2]);
    G4double errMom2 = sqr(yError[3]) + sqr(yError[4]) + sqr(yError[5]);
    G4double Mom2 = sqr(yOut[3]) + sqr(yOut[4]) + sqr(yOut[5]);
    errMom2 /= Mom2;

    dyerr = std::sqrt(std::max(errPos2, errMom2*sqr(hstep)));


    //copy non-integrated variables to output array
    memcpy(yOut + fNumberOfComponents,
           yIn + fNumberOfComponents,
           sizeof(G4double) * (ncomp - fNumberOfComponents));

    //set new state
    track.LoadFromArray(yOut, ncomp);
    track.SetCurveLength(curveLength + hstep);

    return true;
}

void  G4BulirschStoerDriver::OneGoodStep(
    G4double y[],
    const G4double dydx[],
    G4double& curveLength,
    G4double htry,
    G4double eps,
    G4double& hdid,
    G4double& hnext)
{
    hnext = htry;
    G4double curveLengthBegin = curveLength;

    // set maximum allowed error
    bulirschStoer.set_max_relative_error(eps);

    step_result res = step_result::fail;
    while (res == step_result::fail)
    {
        res = bulirschStoer.try_step(y, dydx, curveLength, yOut, hnext);
    }

    memcpy(y, yOut, sizeof(G4double) * fNumberOfComponents);
    hdid = curveLength - curveLengthBegin;
}

void G4BulirschStoerDriver::GetDerivatives(
    const G4FieldTrack& track,
    G4double dydx[]) const
{
    G4double y[G4FieldTrack::ncompSVEC];
    track.DumpToArray(y);
    fEquation->RightHandSide(y, dydx);
}

void G4BulirschStoerDriver::SetVerboseLevel(G4int level)
{
    fVerbosity = level;
}

G4int G4BulirschStoerDriver::GetVerboseLevel() const
{
    return fVerbosity;
}

G4double G4BulirschStoerDriver::ComputeNewStepSize(
    G4double /* errMaxNorm*/,
    G4double  hstepCurrent)
{
    return hstepCurrent;
}

G4EquationOfMotion* G4BulirschStoerDriver::GetEquationOfMotion()
{
    return fEquation;
}

void G4BulirschStoerDriver::SetEquationOfMotion(
    G4EquationOfMotion* equation)
{
    fEquation = equation;
}
