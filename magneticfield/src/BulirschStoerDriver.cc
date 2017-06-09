#include "BulirschStoerDriver.hh"
#include "G4LineSection.hh"
#include "G4SystemOfUnits.hh"
#include "ModifiedMidpoint.hh"
#include "ModifiedMidpointDenseOut.hh"

#include "G4CashKarpRKF45.hh"

#include <iostream>

//#include <boost/numeric/odeint.hpp>

#define ncomp G4FieldTrack::ncompSVEC
//typedef std::array<G4double,ncomp> state_type;

//using namespace boost::numeric::odeint;

BulirschStoerDriver::BulirschStoerDriver(G4double hminimum,
                                         G4EquationOfMotion* equation,
                                         G4int integratedComponents,
                                         G4int verb):
    G4VIntegrationDriver(hminimum,equation,integratedComponents,verb),
    dummyStepper(new BSStepper(equation)),
    modifiedMidpoint(equation,integratedComponents),
    denseMidpoint(equation,integratedComponents),
    bulirschStoer(equation,integratedComponents,0,DBL_MAX),
    tBegin(0),tEnd(0),
    interval_sequence{2,4},
    fcoeff(1./(sqr(G4double(interval_sequence[1])/G4double(interval_sequence[0]))-1.))

{
//    system = [this](const state_type& y, state_type& dydx, double /*t*/){
//        fequation->RightHandSide(y.data(),dydx.data());
//    };
}

BulirschStoerDriver::~BulirschStoerDriver()
{
    delete dummyStepper;
}


G4bool  BulirschStoerDriver::AccurateAdvance(G4FieldTrack&  track,
                         G4double hstep,
                         G4double eps,
                         G4double hinitial)
{
    //reset QuickAdvance
    tBegin = tEnd = 0;

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
    if( hstep < 0)
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
    memcpy(yOut+GetNumberOfVariables(),
           yCurrent+GetNumberOfVariables(),
           sizeof(G4double)*(ncomp-GetNumberOfVariables()));

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
        else if(h > GetMinStep()){
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
            hnext = ComputeNewStepSize(dyerr/eps, h);

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
            if(std::fabs(hnext) < GetMinStep())
            {
              // Make sure that the next step is at least Hmin.
              h = GetMinStep();
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
    memcpy(yOut,yCurrent,sizeof(G4double)*GetNumberOfVariables());

    // upload new state
    track.LoadFromArray(yOut, ncomp);
    track.SetCurveLength(curveLength);

    if(nstp > fMaxNoSteps)
    {
      ++no_warnings;
      succeeded = false;
    }

    return succeeded;

/*
    if (hstep > fMinimumStep){
        //G4cout<<"BulirschStoerDriver::AccurateAdvance \n";
        G4double dydx[ncomp];
        G4double y[ncomp];
        track.DumpToArray(y);
        G4double hdid = 0;
        G4double hnext = beginStep != 0 ? beginStep : hstep;
        //G4double hnext = beginStep;
        G4double hrest = hstep;
        G4double curveLength = track.GetCurveLength();
        do{
            fequation->RightHandSide(y,dydx);
            OneGoodStep(y,dydx,curveLength,hnext,eps,hdid,hnext);
            hrest -= hdid;
            hnext = std::min(hrest, hnext);
            //G4cout<<"hdid "<<hdid<<" hrest "<<hrest<<" hnext "<<hnext<<G4endl;
        }while(hrest > eps*hstep);
        //G4cout<<"driver: did step "<<stepInitial<<G4endl;
        track.SetCurveLength(curveLength);
        track.LoadFromArray(y,ncomp);
        //G4cout<<y[0]<<"  "<<y[1]<<"  "<<y[2]<<G4endl;
    }

    return  true;*/
}

/*
 * //version with G4ChashKarp for quick advance
G4bool  BulirschStoerDriver::QuickAdvance(G4FieldTrack& track,
                                          const G4double dydx[],
                                          G4double hstep,
                                          G4double& missDist,
                                          G4double& dyerr)
{

    track.DumpToArray(yIn);
    const G4double curveLength = track.GetCurveLength();

    G4CashKarpRKF45 cashKarp(GetEquationOfMotion());

    cashKarp.Stepper(yIn,dydx,hstep,yOut,yError);
    missDist = cashKarp.DistChord();


    G4double errPos2 = sqr(yError[0]) + sqr(yError[1]) + sqr(yError[2]);
    G4double errMom2 = sqr(yError[3]) + sqr(yError[4]) + sqr(yError[5]);
    G4double Mom2 = sqr(yOut[3]) + sqr(yOut[4]) + sqr(yOut[5]);
    errMom2 /= Mom2;

    dyerr = std::sqrt(std::max(errPos2, errMom2*sqr(hstep)));


    //copy non-integrated variables to output array
    memcpy(yOut+GetNumberOfVariables(),
           yIn+GetNumberOfVariables(),
           sizeof(G4double)*(ncomp-GetNumberOfVariables()));

    //set new state
    track.LoadFromArray(yOut,ncomp);
    track.SetCurveLength(curveLength + hstep);


    return true;
}
*/

//#define dense
#ifndef dense
G4bool  BulirschStoerDriver::QuickAdvance(G4FieldTrack& track,
                                          const G4double dydx[],
                                          G4double hstep,
                                          G4double& missDist,
                                          G4double& dyerr)
{

    track.DumpToArray(yIn);
    const G4double curveLength = track.GetCurveLength();

    /*
    G4int interval1 = 2, interval2 = 4;
    G4double coeff = 1./(sqr(static_cast<G4double>(interval2)/static_cast<G4double>(interval1))-1.);
    */

    denseMidpoint.set_steps(interval_sequence[0]);
    denseMidpoint.do_step(yIn, dydx, yOut, hstep, yMid, derivs[0]);

    denseMidpoint.set_steps(interval_sequence[1]);
    denseMidpoint.do_step(yIn, dydx, yOut2, hstep, yMid2, derivs[1]);

    //extrapolation
    for (G4int i = 0; i < GetNumberOfVariables(); ++i){
        yOut[i] = yOut2[i] + (yOut2[i] - yOut[i])*fcoeff;
        yMid[i] = yMid2[i] + (yMid2[i] - yMid[i])*fcoeff;
    }



    //calc chord lenght
    missDist = G4LineSection::Distline(G4ThreeVector(yMid[0],yMid[1],yMid[2]),
                                       G4ThreeVector(yIn[0],yIn[1],yIn[2]),
                                       G4ThreeVector(yOut[0],yOut[1],yOut[2]));

    //calc error
    for (G4int i = 0; i < GetNumberOfVariables(); ++i){
        yError[i] = yOut[i] - yOut2[i];
    }

    G4double errPos2 = sqr(yError[0]) + sqr(yError[1]) + sqr(yError[2]);
    G4double errMom2 = sqr(yError[3]) + sqr(yError[4]) + sqr(yError[5]);
    G4double Mom2 = sqr(yOut[3]) + sqr(yOut[4]) + sqr(yOut[5]);
    errMom2 /= Mom2;

    dyerr = std::sqrt(std::max(errPos2, errMom2*sqr(hstep)));


    //copy non-integrated variables to output array
    memcpy(yOut+GetNumberOfVariables(),
           yIn+GetNumberOfVariables(),
           sizeof(G4double)*(ncomp-GetNumberOfVariables()));

    //set new state
    track.LoadFromArray(yOut,ncomp);
    track.SetCurveLength(curveLength + hstep);


    return true;
}
#endif

#ifdef dense
//Dense version
G4bool  BulirschStoerDriver::QuickAdvance(G4FieldTrack& track,
                                          const G4double dydx[],
                                          G4double hstep,
                                          G4double& missDist,
                                          G4double& dyerr)
{

    track.DumpToArray(yIn);
    const G4double curveLength = track.GetCurveLength();


    if ((curveLength >= tBegin) && ((curveLength+hstep) <= tEnd))
    {
        //do interpolation
        interpolate(curveLength+hstep,yOut);
        interpolate(curveLength+0.5*hstep,yMid);

        //error is estimated as highest coeff
        memcpy(yError,diffs[3][0],sizeof(G4double)*ncomp);
        //G4cout<<"interpolate! \n";

/*

        G4int interval1 = 2, interval2 = 6;
        G4double coeff = 1./(sqr(static_cast<G4double>(interval2)/static_cast<G4double>(interval1))-1.);

        G4double tmpOut[12],tmpOut2[12],tmpMid[12],tmpMid2[12];

        denseMidpoint.set_steps(interval1);
        denseMidpoint.do_step(yIn, dydx, tmpOut, hstep, tmpMid, derivs[0]);

        denseMidpoint.set_steps(interval2);
        denseMidpoint.do_step(yIn, dydx, tmpOut2, hstep, tmpMid2, derivs[1]);


        G4double tmpError[12];
        //extrapolation
        for (G4int i = 0; i < GetNumberOfVariables(); ++i)
        {
            // extrapolate solution (4th order)
            tmpOut[i] = tmpOut2[i] + (tmpOut2[i] - tmpOut[i])*coeff;
            tmpMid[i] = tmpMid2[i] + (tmpMid2[i] - tmpMid[i])*coeff;
            tmpError[i] = tmpOut[i] - yOut[i];
        }


        G4double errPos2 = sqr(tmpError[0]) + sqr(tmpError[1]) + sqr(tmpError[2]);
        G4double errMom2 = sqr(tmpError[3]) + sqr(tmpError[4]) + sqr(tmpError[5]);
        G4double Mom2 = sqr(tmpOut[3]) + sqr(tmpOut[4]) + sqr(tmpOut[5]);
        errMom2 /= Mom2;

        G4double err = std::max(sqrt(errPos2), sqrt(errMom2)*hstep);
        //G4cout<<"interp err: "<<err<<" hstep "<<hstep<<G4endl;

        if (err > 1)
        {
            memcpy(yOut,tmpOut,sizeof(G4double)*GetNumberOfVariables());
            memcpy(yMid,tmpMid,sizeof(G4double)*GetNumberOfVariables());
            G4cout<<"interpol is bad, err: "<<err<<G4endl;
        }
        else
        {
            //G4cout<<"interpol is good \n";
        }
*/
    }
    else
    {
        //G4cout<<"calculate \n";
        //set curve length interval in which interpolation is avaliable
        tBegin = curveLength;
        tEnd = tBegin + hstep;

        G4int interval1 = 2, interval2 = 6;
        G4double coeff = 1./(sqr(static_cast<G4double>(interval2)/static_cast<G4double>(interval1))-1.);


        denseMidpoint.set_steps(interval1);
        denseMidpoint.do_step(yIn, dydx, yOut, hstep, yMid, derivs[0]);

        denseMidpoint.set_steps(interval2);
        denseMidpoint.do_step(yIn, dydx, yOut2, hstep, yMid2, derivs[1]);

        memcpy(diffs[0][0],derivs[0][0],sizeof(G4double)*GetNumberOfVariables());
        memcpy(diffs[0][1],derivs[1][2],sizeof(G4double)*GetNumberOfVariables());

        for (G4int i = 0; i < GetNumberOfVariables(); ++i)
        {
            diffs[1][0][i] = (derivs[0][1][i] - dydx[i])*1./hstep;
            diffs[1][1][i] = (derivs[1][3][i] - derivs[1][1][i])*3./hstep;
            diffs[2][0][i] = (derivs[1][4][i] - 2*derivs[1][2][i] + derivs[1][0][i])*9./(hstep*hstep);
            diffs[3][0][i] = (derivs[1][5][i] - 3*derivs[1][3][i] + 3*derivs[1][1][i] - dydx[i])*
                             27./(hstep*hstep*hstep);
        }


        //extrapolation
        for (G4int i = 0; i < GetNumberOfVariables(); ++i)
        {
            // extrapolate solution (4th order)
            yOut[i] = yOut2[i] + (yOut2[i] - yOut[i])*coeff;
            yMid[i] = yMid2[i] + (yMid2[i] - yMid[i])*coeff;

            // extrapolate finite differences (4th order)
            diffs[0][0][i] = diffs[0][1][i] + (diffs[0][1][i] - diffs[0][0][i])*coeff;
            diffs[1][0][i] = diffs[1][1][i] + (diffs[1][1][i] - diffs[1][0][i])*coeff;
        }


        G4double d = hstep/2.;
        for(G4int i = 0; i < 4; ++i)
        {
            // we need these terms for dense output interpolation
            for (G4int j = 0; j < GetNumberOfVariables(); ++j)
            {
                diffs[i][0][j] *= d;
            }
            d *= hstep/(2.*(i+2.));
        }


        //calc error
        for (G4int i = 0; i < GetNumberOfVariables(); ++i){
            yError[i] = yOut[i] - yOut2[i];
        }

    }


    //calc chord lenght
    missDist = G4LineSection::Distline(G4ThreeVector(yMid[0],yMid[1],yMid[2]),
                                       G4ThreeVector(yIn[0],yIn[1],yIn[2]),
                                       G4ThreeVector(yOut[0],yOut[1],yOut[2]));



    G4double errPos2 = sqr(yError[0]) + sqr(yError[1]) + sqr(yError[2]);
    G4double errMom2 = sqr(yError[3]) + sqr(yError[4]) + sqr(yError[5]);
    G4double Mom2 = sqr(yOut[3]) + sqr(yOut[4]) + sqr(yOut[5]);
    errMom2 /= Mom2;

    dyerr = std::max(sqrt(errPos2), sqrt(errMom2)*hstep);


    //copy non-integrated variables to output array
    memcpy(yOut+GetNumberOfVariables(),
           yIn+GetNumberOfVariables(),
           sizeof(G4double)*(ncomp-GetNumberOfVariables()));

    //set new state
    track.LoadFromArray(yOut,ncomp);
    track.SetCurveLength(curveLength + hstep);


    return true;
}
#endif


void BulirschStoerDriver::interpolate(G4double time, G4double out[])
{
    G4double theta = 2 * (time - tBegin)/(tEnd - tBegin) - 1;
    memcpy(out,yMid,sizeof(G4double)*ncomp);
    G4double theta_pow = theta;

    //4th order interpolation
    for (G4int i = 0; i < 4; ++i)
    {
        for (G4int j = 0; j < GetNumberOfVariables(); ++j)
        {
            out[j] += diffs[i][0][j]*theta_pow;
        }
        theta_pow *= theta;
    }
}


void BulirschStoerDriver::GetDerivatives(const G4FieldTrack& track, G4double dydx[] )
{
    track.DumpToArray(yIn);
    GetEquationOfMotion()->RightHandSide(yIn,dydx);
}

#define use_boost 0

void  BulirschStoerDriver::OneGoodStep(G4double  y[],
                                       const G4double  dydx[],
                                       G4double& curveLength,
                                       G4double htry,
                                       G4double eps,
                                       G4double& hdid,
                                       G4double& hnext)
{
    hnext = htry;
    hdid = 0;
#if use_boost
    static bulirsch_stoer<state_type> stepper(0,eps,1,0);
    stepper.set_max_relative_error(eps);
    state_type yStateIn,dydxState,yStateOut;
    memcpy(yStateIn.data(),y,sizeof(G4double)*GetNumberOfVariables());
    memcpy(dydxState.data(),dydx,sizeof(G4double)*GetNumberOfVariables());
    memset(dydxState.data()+GetNumberOfVariables(),0,sizeof(G4double)*(ncomp-GetNumberOfVariables()));
    stepper.try_step(
            [this](const state_type ySt, state_type& dydxSt, G4double /*t*/){
                this->GetEquationOfMotion()->RightHandSide(ySt.data(),dydxSt.data());},
            yStateIn, dydxState, hdid, yStateOut, hnext);
    memcpy(y,yStateOut.data(),sizeof(G4double)*GetNumberOfVariables());

#else
    //bulirschStoer.reset();
    //bulirschStoer.set_max_dt(htry);

    // set maximum allowed error
    bulirschStoer.set_max_relative_error(eps);

    // try step, changes hdid & hnext
    step_result res = bulirschStoer.try_step(y,dydx,hdid,yOut,hnext);

    // if step succeed copy integrated vars to output array
    if (res == step_result::success)
        memcpy(y,yOut,sizeof(G4double)*GetNumberOfVariables());
#endif


    curveLength += hdid;

}





G4double BulirschStoerDriver::ComputeNewStepSize(double /*error*/,
                                                 double h)
{
    /*
    G4double STEPFAC1( 0.65 ),
            STEPFAC2( 0.94 ),
            STEPFAC3( 0.02 ),
            STEPFAC4( 4.0 );

    G4int k = 1;
    G4double expo( 1.0/(2*k+1) );
    G4double facmin = std::pow( STEPFAC3 , expo );

    G4double fac;
    if (error == 0.0)
        fac=1.0/facmin;
    else
    {
        fac = STEPFAC2 / std::pow( error / STEPFAC1 , expo );
        fac = std::max(facmin/STEPFAC4, std::min(1.0/facmin, fac));
    }

    return h*fac;
    */
    return h;
}

G4MagIntegratorStepper* BulirschStoerDriver::GetStepper()
{
    return dummyStepper;
}
