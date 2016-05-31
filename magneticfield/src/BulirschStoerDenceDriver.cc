#include "BulirschStoerDenceDriver.hh"
#include "G4LineSection.hh"

#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <functional>

#define ncomp G4FieldTrack::ncompSVEC
#define nvar 6

using namespace boost::numeric::odeint;
typedef boost::array<G4double,6> state_type;

void print(const state_type& state){
    for (int i = 0; i < 6; ++i){
        G4cout<<state[i]<<"  ";
    }
    G4cout<<G4endl;
}

BulirschStoerDenceDriver::BulirschStoerDenceDriver(G4EquationOfMotion* _equation, G4int _integComp, G4int  _verb):
    equation(_equation),integComp(_integComp),verb(_verb)
{
    system = [_equation](const state_type& y, state_type& dydx, double /*t*/){
        G4double arrayY[ncomp],arrayDydx[ncomp];
        memcpy(arrayY,y.data(),sizeof(G4double)*nvar);
        memcpy(arrayDydx,dydx.data(),sizeof(G4double)*nvar);
        _equation->RightHandSide(arrayY,arrayDydx);
        memcpy(dydx.data(),arrayDydx,sizeof(G4double)*nvar);
        //print(dydx);
    };
}

BulirschStoerDenceDriver::~BulirschStoerDenceDriver(){

}

G4double  BulirschStoerDenceDriver::QuickAdvance(G4FieldTrack& track,
                                              G4double hstep,
                                              G4double& missDist,
                                              G4double epsStep){

    G4double eps = std::sqrt(sqr(quickEps) + sqr(epsStep));
    bulirsch_stoer_dense_out<state_type> stepper(0,eps,1,0);

    G4double arrayY[nvar];
    state_type yIn, yMid, yOut;
    track.DumpToArray(arrayY);
    memcpy(yIn.data(),arrayY,sizeof(G4double)*6);
    G4double curveLength = track.GetCurveLength();
    //print(yIn);

    //G4cout<<"QuickAdvance: "<<"curverLength "<<curveLength<<" hstep "<<hstep<<G4endl;

    stepper.initialize(yIn,curveLength,hstep);

    stepper.do_step(system);
    //interpolate to h/2
    stepper.calc_state(curveLength + hstep/2., yMid);

    yOut = stepper.current_state();

    G4ThreeVector mid(yMid[0],yMid[1],yMid[2]);
    G4ThreeVector in(yIn[0],yIn[1],yIn[2]);
    G4ThreeVector out(yOut[0],yOut[1],yOut[2]);

    //G4cout<<"in "<<in<<" mid "<<mid<<" out "<<out<<G4endl;

    missDist = G4LineSection::Distline(mid,in,out);

    curveLength = stepper.current_time();

    track.SetCurveLength(curveLength);
    memcpy(arrayY,yOut.data(),sizeof(G4double)*6);
    track.LoadFromArray(arrayY,ncomp);

    return eps;
}


G4bool  BulirschStoerDenceDriver::AccurateAdvance(G4FieldTrack&  track,
                         G4double hstep,
                         G4double eps){


    bulirsch_stoer_dense_out<state_type> stepper(0,eps,1,0);
    G4double arrayY[ncomp];
    state_type yIn, yOut;
    track.DumpToArray(arrayY);
    memcpy(yIn.data(),arrayY,sizeof(G4double));
    G4double curveLength = track.GetCurveLength();

    stepper.initialize(yIn,curveLength,hstep);

    stepper.do_step(system);

    yOut = stepper.current_state();

    curveLength = stepper.current_time();

    track.SetCurveLength(curveLength);
    memcpy(arrayY,yOut.data(),sizeof(G4double)*6);
    track.LoadFromArray(yOut.data(),ncomp);

    return  true;
}


void BulirschStoerDenceDriver::GetDerivatives(const G4FieldTrack& track, G4double dydx[] ){
    G4double y[ncomp];
    track.DumpToArray(y);
    equation->RightHandSide(y,dydx);
}

