#include "BulirschStoerDriver.hh"
#include "G4LineSection.hh"

#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>

#define ncomp G4FieldTrack::ncompSVEC
#define nvar 6

using namespace boost::numeric::odeint;
typedef boost::array<G4double,6> state_type;

BulirschStoerDriver::BulirschStoerDriver(G4EquationOfMotion* _equation,G4double  _hmin,G4int _integComp, G4int  _verb):
    equation(_equation),hmin(_hmin),integComp(_integComp),verb(_verb)
{

}

BulirschStoerDriver::~BulirschStoerDriver(){

}


G4bool  BulirschStoerDriver::AccurateAdvance(G4FieldTrack&  track,
                         G4double hstep,
                         G4double eps,
                         G4double beginStep){
    G4double dydx[ncomp];
    G4double y[ncomp];
    track.DumpToArray(y);
    G4double hdid = 0;
    G4double hnext = beginStep;
    G4double hrest = hstep;
    G4double curveLength = track.GetCurveLength();
    do{
        equation->RightHandSide(y,dydx);
        OneGoodStep(y,dydx,curveLength,hnext,eps,hdid,hnext);
        hrest -= hdid;
        hnext = std::min(hrest, hnext);
        //G4cout<<"hdid "<<hdid<<" hrest "<<hrest<<" hnext "<<hnext<<G4endl;
    }while(hrest > eps*hstep);
    //G4cout<<"driver: did step "<<stepInitial<<G4endl;
    track.SetCurveLength(curveLength);
    track.LoadFromArray(y,ncomp);

    return  true;
}

G4bool  BulirschStoerDriver::QuickAdvance(G4FieldTrack& track,
                                          const G4double dydx[],
                                          G4double hstep,
                                          G4double& missDist,
                                          G4double& dyerr){

    state_type yIn, yOutMid, yOut1, yOut2;
    state_type dydxMid, dydxIn;
    track.DumpToArray(yIn.data());
    for (int i = 0; i < ncomp; ++i){
        dydxIn[i] = dydx[i];
    }
    G4double curveLength1 = track.GetCurveLength();
    G4double curveLength2 = curveLength1;
    const auto system = [this](const state_type & _y, state_type & _dydx, const double /*t*/){
        // there is access out of range!
        //solution: copy _y and _dydx to array of 8 elements
        equation->RightHandSide(_y.data(),_dydx.data());
    };

    modified_midpoint<state_type> mp;
    mp.do_step(system,yIn,dydxIn,curveLength2,yOut2,hstep);

    mp.do_step(system,yIn,dydxIn,curveLength1,yOutMid,hstep/2.);
    system(yOutMid,dydxMid,1);
    mp.do_step(system,yOutMid,dydxMid,curveLength1,yOut1,hstep/2.);

    missDist = G4LineSection::Distline(G4ThreeVector(yOutMid[0],yOutMid[1],yOutMid[2]),
                                       G4ThreeVector(yIn[0],yIn[1],yIn[2]),
                                       G4ThreeVector(yOut1[0],yOut1[1],yOut1[2]));


    G4double yErr[ncomp];
    for (int i = 0; i < nvar; ++i){
        yErr[i] = yOut1[i] - yOut2[i];
    }

    G4double errPos2 = sqr(yErr[0]) + sqr(yErr[1]) + sqr(yErr[2]);
    G4double errMom2 = sqr(yErr[3]) + sqr(yErr[4]) + sqr(yErr[5]);
    G4double Mom2 = sqr(yOut1[3]) + sqr(yOut1[4]) + sqr(yOut1[5]);

    dyerr = std::max(std::sqrt(errPos2)/hstep, std::sqrt(errMom2/Mom2));

    return true;

}


void BulirschStoerDriver::GetDerivatives(const G4FieldTrack& track, G4double dydx[] ){
    G4double y[ncomp];
    track.DumpToArray(y);
    equation->RightHandSide(y,dydx);
}

void  BulirschStoerDriver::OneGoodStep(G4double  y[],
                                       const G4double  dydx[],
                                       G4double& curveLength,
                                       G4double htry,
                                       G4double eps,
                                       G4double& hdid,
                                       G4double& hnext){
    bulirsch_stoer<state_type> stepper(0,eps,1,0,htry); //(epsStep_Relative, epsStep_Relative);

    hnext = htry;
    hdid = 0;

    const auto system = [this](const state_type & _y, state_type & _dydx, const double /*t*/){
        // there is access out of range!
        //solution: copy _y and _dydx to array of 8 elements
        equation->RightHandSide(_y.data(),_dydx.data());
    };
    state_type yState, dydxState;
    for (int i = 0; i < 6; ++i){
        yState[i] = y[i];
        dydxState[i] = dydx[i];
    }
    // System system , const StateIn &in , time_type &t , StateOut &out , time_type &dt
    stepper.try_step(system, yState, dydxState, hdid, hnext);
    for (int i = 0; i < 6; ++i){
        y[i] = yState[i];
    }
    curveLength += hdid;
}
