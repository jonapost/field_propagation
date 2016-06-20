#include "BulirschStoerDriver.hh"
#include "G4LineSection.hh"
#include "G4SystemOfUnits.hh"
#include "ModifiedMidpoint.hh"

#include <boost/numeric/odeint.hpp>

#define ncomp G4FieldTrack::ncompSVEC

using namespace boost::numeric::odeint;

BulirschStoerDriver::BulirschStoerDriver(G4double hminimum,
                                         G4EquationOfMotion* equation,
                                         G4int integratedComponents,
                                         G4int verb):
    G4VIntegrationDriver(hminimum,equation,integratedComponents,verb),
    dummyStepper(new BSStepper(equation)),
    modifiedMidpoint(equation,integratedComponents),
    bulirschStoer(equation,integratedComponents,0,0)

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
    G4ThreeVector StartPos, EndPos;

    G4FieldTrack yFldTrk(track);

    do
    {
        StartPos = {yCurrent[0], yCurrent[1], yCurrent[2]};

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


        EndPos = {yCurrent[0], yCurrent[1], yCurrent[2]};

        // Check the endpoint
        G4double endPointDist = (EndPos - StartPos).mag();
        if (endPointDist >= hdid*(1. + perMillion))
        {
            fNoBadSteps++;

            // Issue a warning only for gross differences -
            // we understand how small difference occur.
            if (endPointDist >= hdid*(1.+perThousand))
            {
                ++no_warnings;
            }
        }
        else
        {
            noGoodSteps ++;
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
      no_warnings++;
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


G4bool  BulirschStoerDriver::QuickAdvance(G4FieldTrack& track,
                                          const G4double dydx[],
                                          G4double hstep,
                                          G4double& missDist,
                                          G4double& dyerr){



    track.DumpToArray(yIn);
    const G4double curveLength = track.GetCurveLength();
/*
    // boost version is not used
    // yIn,dydxIn... should be state_type not G4double*

    modified_midpoint<state_type> mp;
    mp.do_step(system,yIn,dydxIn,curveLength,yOut2,hstep);

    mp.do_step(system,yIn,dydxIn,curveLength,yOutMid,hstep/2.);
    system(yOutMid,dydxMid,curveLength+hstep/2.);
    mp.do_step(system,yOutMid,dydxMid,curveLength+hstep/2.,yOut,hstep/2.);
*/

    //do full step by modified midpoint
    modifiedMidpoint.do_step(yIn,dydx,yOut2,hstep);

    //do two steps by h/2
    modifiedMidpoint.do_step(yIn,dydx,yMid,0.5*hstep);
    GetEquationOfMotion()->RightHandSide(yMid,dydxMid);
    modifiedMidpoint.do_step(yMid,dydxMid,yOut,0.5*hstep);

    //calc chord lenght
    missDist = G4LineSection::Distline(G4ThreeVector(yMid[0],yMid[1],yMid[2]),
                                       G4ThreeVector(yIn[0],yIn[1],yIn[2]),
                                       G4ThreeVector(yOut[0],yOut[1],yOut[2]));

    //calc error
    for (int i = 0; i < 6; ++i){
        yError[i] = yOut[i] - yOut2[i];
    }
    G4double errPos2 = sqr(yError[0]) + sqr(yError[1]) + sqr(yError[2]);
    G4double errMom2 = sqr(yError[3]) + sqr(yError[4]) + sqr(yError[5]);
    G4double Mom2 = sqr(yOut[3]) + sqr(yOut[4]) + sqr(yOut[5]);
    errMom2 /= Mom2;
    errPos2 /= (hstep*hstep);

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

void BulirschStoerDriver::GetDerivatives(const G4FieldTrack& track, G4double dydx[] )
{
    track.DumpToArray(yIn);
    GetEquationOfMotion()->RightHandSide(yIn,dydx);
}


void  BulirschStoerDriver::OneGoodStep(G4double  y[],
                                       const G4double  dydx[],
                                       G4double& curveLength,
                                       G4double htry,
                                       G4double eps,
                                       G4double& hdid,
                                       G4double& hnext)
{
    //G4cout<<"OneGoodStep "<<htry<<G4endl;
    //bulirsch_stoer<state_type> stepper(0,eps,1,0,htry);

    //bulirschStoer.reset();
    bulirschStoer.set_max_dt(htry);
    bulirschStoer.set_max_relative_error(eps);

    hnext = htry;
    hdid = 0;
    bulirschStoer.try_step(y,dydx,hdid,yOut,hnext);

    memcpy(y,yOut,sizeof(G4double)*GetNumberOfVariables());
    curveLength += hdid;
/*
    //copy input to std::array and set non-integrated variables 0
    memcpy(yInOut.data(),y,sizeof(G4double)*fnvar);
    memset(yInOut.data()+fnvar,0,sizeof(G4double)*(ncomp - fnvar));

    memcpy(dydxIn.data(),dydx,sizeof(G4double)*fnvar);
    memset(dydxIn.data()+fnvar,0,sizeof(G4double)*(ncomp - fnvar));

    //do Bulisrsch-Stoer step
    bulirschStoer.try_step(system, yInOut, dydxIn, hdid, hnext);

    //copy integrated variables to output array
    memcpy(y,yInOut.data(),sizeof(G4double)*fnvar);

    curveLength += hdid;
    */
}

G4double BulirschStoerDriver::ComputeNewStepSize(double /*dyErr_relative*/,
                                                 double lastStepLength )
{
    return lastStepLength;
}

G4MagIntegratorStepper* BulirschStoerDriver::GetStepper()
{
    return dummyStepper;
}
