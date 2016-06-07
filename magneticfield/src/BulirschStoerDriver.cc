#include "BulirschStoerDriver.hh"
#include "G4LineSection.hh"

#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>

#define ncomp G4FieldTrack::ncompSVEC

using namespace boost::numeric::odeint;

BulirschStoerDriver::BulirschStoerDriver(G4EquationOfMotion* equation,
                                         G4int integratedComponents,
                                         G4int  verb):
    BaseDriver(equation,integratedComponents,verb)

{
}

BulirschStoerDriver::~BulirschStoerDriver(){

}


G4bool  BulirschStoerDriver::AccurateAdvance(G4FieldTrack&  track,
                         G4double hstep,
                         G4double eps,
                         G4double beginStep){
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

    return  true;
}

G4bool  BulirschStoerDriver::QuickAdvance(G4FieldTrack& track,
                                          const G4double dydx[],
                                          G4double hstep,
                                          G4double& missDist,
                                          G4double& dyerr){

    //G4cout<<"BulirschStoerDriver::QuickAdvance \n";
    state_type yIn, yOutMid, yOut1, yOut2;
    state_type dydxMid, dydxIn;
    track.DumpToArray(yIn.data());
    memcpy(dydxIn.data(),dydx,sizeof(G4double)*ncomp);
    G4double curveLength1 = track.GetCurveLength();
    G4double curveLength2 = curveLength1;

    modified_midpoint<state_type> mp;
    mp.do_step(system,yIn,dydxIn,curveLength2,yOut2,hstep);

    mp.do_step(system,yIn,dydxIn,curveLength1,yOutMid,hstep/2.);
    system(yOutMid,dydxMid,1);
    mp.do_step(system,yOutMid,dydxMid,curveLength1,yOut1,hstep/2.);

    missDist = G4LineSection::Distline(G4ThreeVector(yOutMid[0],yOutMid[1],yOutMid[2]),
                                       G4ThreeVector(yIn[0],yIn[1],yIn[2]),
                                       G4ThreeVector(yOut1[0],yOut1[1],yOut1[2]));


    G4double yErr[ncomp];
    for (int i = 0; i < fnvar; ++i){
        yErr[i] = yOut1[i] - yOut2[i];
    }

    G4double errPos2 = sqr(yErr[0]) + sqr(yErr[1]) + sqr(yErr[2]);
    G4double errMom2 = sqr(yErr[3]) + sqr(yErr[4]) + sqr(yErr[5]);
    G4double Mom2 = sqr(yOut1[3]) + sqr(yOut1[4]) + sqr(yOut1[5]);
    errMom2 /= Mom2;
    errPos2 /= (hstep*hstep);

    track.LoadFromArray(yOut1.data(),ncomp);

    dyerr = std::max(sqrt(errPos2), sqrt(errMom2)*hstep);
    //G4cout<<"dyerr "<<dyerr<<G4endl;

    //G4cout<<yOut2[0]<<"  "<<yOut2[1]<<"  "<<yOut2[2]<<G4endl;

    return true;
}

void BulirschStoerDriver::GetDerivatives(const G4FieldTrack& track, G4double dydx[] ){
    G4double y[G4FieldTrack::ncompSVEC];
    track.DumpToArray(y);
    fequation->RightHandSide(y,dydx);
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

    state_type yState, dydxState;
    memcpy(yState.data(),y,sizeof(G4double)*ncomp);
    memcpy(dydxState.data(),dydx,sizeof(G4double)*ncomp);

    // System system , const StateIn &in , time_type &t , StateOut &out , time_type &dt
    stepper.try_step(system, yState, dydxState, hdid, hnext);
    memcpy(y,yState.data(),sizeof(G4double)*ncomp);

    curveLength += hdid;
}
