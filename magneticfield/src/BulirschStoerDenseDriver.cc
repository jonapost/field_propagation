#include "BulirschStoerDenseDriver.hh"
#include "G4LineSection.hh"

#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <functional>

#define ncomp G4FieldTrack::ncompSVEC


void print(const state_type& state){
    for (int i = 0; i < 6; ++i){
        G4cout<<state[i]<<"  ";
    }
    G4cout<<G4endl;
}

BulirschStoerDenseDriver::BulirschStoerDenseDriver(G4EquationOfMotion* pequation,
                                                   G4int numberOfComponents,
                                                   G4int statisticsVerbosity):
    BaseDriver(pequation,numberOfComponents,statisticsVerbosity),
    quickEps(1e50),
    tBegin(DBL_MAX),
    tEnd(DBL_MIN),
    theStepper(0,0,1,0,0,true)

{
}

BulirschStoerDenseDriver::~BulirschStoerDenseDriver(){

}

G4bool  BulirschStoerDenseDriver::QuickAdvance(G4FieldTrack& track,
                                               const G4double /*dydx*/[],
                                                 G4double hstep,
                                                 G4double& missDist,
                                                 G4double& dyerr){


    G4double eps = std::sqrt(sqr(quickEps) + sqr(theStepper.m_error_checker.m_eps_rel)); //geometric mean
    dyerr = eps*hstep;
    theStepper.m_error_checker.m_eps_rel = eps;
    theStepper.m_max_dt = hstep;
    state_type yIn, yMid, yOut;
    G4ThreeVector inVec,midVec,outVec;
    track.DumpToArray(yIn.data());
    G4double time = track.GetCurveLength();
    if (time >= tBegin && (time+hstep) <= tEnd){
        theStepper.calc_state(time + hstep/2.,yMid);

        theStepper.calc_state(time+hstep,yOut);

        for (int i = 0; i < 3; ++i){
            inVec[i] = yIn[i];
            midVec[i] = yMid[i];
            outVec[i] = yOut[i];
        }

        missDist = G4LineSection::Distline(midVec,inVec,outVec);

        track.SetCurveLength(time+hstep);
        track.LoadFromArray(yOut.data(),ncomp);

    }
    else{

        theStepper.initialize(yIn,time,hstep);
        auto tInterval = theStepper.do_step(system);
        tBegin = tInterval.first;
        tEnd = tInterval.second;
        G4double tMid = 0.5*(tBegin + tEnd);

        //interpolate to h/2
        theStepper.calc_state(tMid, yMid);

        yOut = theStepper.current_state();

        for (int i = 0; i < 3; ++i){
            inVec[i] = yIn[i];
            midVec[i] = yMid[i];
            outVec[i] = yOut[i];
        }

        missDist = G4LineSection::Distline(midVec,inVec,outVec);


        track.SetCurveLength(tEnd);
        track.LoadFromArray(yOut.data(),ncomp);
    }

    return true;
}

G4bool  BulirschStoerDenseDriver::AccurateAdvance(G4FieldTrack&  track,
                        G4double hstep,
                        G4double eps,
                        G4double /*beginStep = 0*/){  // Suggested 1st interval

    //G4cout<<"BulirschStoerDenseDriver::AccurateAdvance \n";

    theStepper.m_error_checker.m_eps_rel = eps;
    theStepper.m_max_dt = hstep;

    state_type y;
    track.DumpToArray(y.data());
    G4double curveLength = track.GetCurveLength();
    G4double clEnd = curveLength + hstep;
    G4double hrest = hstep;

    //G4cout<<"in: ";
    //print(y);

    theStepper.initialize(y,curveLength,hrest);
    while(hrest > eps){
        theStepper.do_step(system);
        curveLength = theStepper.current_time();
        hrest = clEnd - curveLength;
        theStepper.m_max_dt = hrest;
        theStepper.m_dt = hrest;
    }
    y = theStepper.current_state();
    //G4cout<<" out: ";
    //print(y);

    track.SetCurveLength(curveLength);
    track.LoadFromArray(y.data(),ncomp);

    return true;
}

void BulirschStoerDenseDriver::GetDerivatives(const G4FieldTrack &/*track*/, G4double /*dydx*/[]){
    //no need to calculate dydx
    //it is calculated internaly
}
