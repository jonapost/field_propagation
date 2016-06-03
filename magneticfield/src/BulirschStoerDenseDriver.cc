#include "BulirschStoerDenseDriver.hh"
#include "G4LineSection.hh"

#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <functional>

#define ncomp G4FieldTrack::ncompSVEC

using namespace boost::numeric::odeint;
typedef boost::array<G4double,6> state_type;

void print(const state_type& state){
    for (int i = 0; i < 6; ++i){
        G4cout<<state[i]<<"  ";
    }
    G4cout<<G4endl;
}

BulirschStoerDenseDriver::BulirschStoerDenseDriver(G4EquationOfMotion* pequation,
                                                   G4int numberOfComponents,
                                                   G4int statisticsVerbosity):
    equation(pequation),
    nvar(numberOfComponents),
    verb(statisticsVerbosity),
    quickEps(0.01)
{
    system = [this](const state_type& y, state_type& dydx, double /*t*/){
        G4double arrayY[ncomp],arrayDydx[ncomp];
        memcpy(arrayY,y.data(),sizeof(G4double)*nvar);
        memcpy(arrayDydx,dydx.data(),sizeof(G4double)*nvar);
        equation->RightHandSide(arrayY,arrayDydx);
        memcpy(dydx.data(),arrayDydx,sizeof(G4double)*nvar);
        //print(dydx);
    };
}

BulirschStoerDenseDriver::~BulirschStoerDenseDriver(){

}

G4double  BulirschStoerDenseDriver::QuickAdvance(G4FieldTrack& track,
                                                 G4double hstep,
                                                 G4double& chord,
                                                 G4double eps){

    eps = std::sqrt(sqr(quickEps) + sqr(eps)); //geometric mean
    bulirsch_stoer_dense_out<state_type> stepper(0,eps,1,0);

    G4double arrayY[ncomp];
    state_type yIn, yMid, yOut;
    track.DumpToArray(arrayY);
    memcpy(yIn.data(),arrayY,sizeof(G4double)*nvar);
    G4double curveLengthBegin = track.GetCurveLength();
    //print(yIn);

    //G4cout<<"QuickAdvance: "<<"curverLength "<<curveLength<<" hstep "<<hstep<<G4endl;

    stepper.initialize(yIn,curveLengthBegin,hstep);
    stepper.do_step(system);
    G4double curveLengthEnd = stepper.current_time();
    G4double curverLengthMid = 0.5*(curveLengthBegin + curveLengthEnd);

    //interpolate to h/2
    stepper.calc_state(curverLengthMid, yMid);

    yOut = stepper.current_state();

    G4ThreeVector mid(yMid[0],yMid[1],yMid[2]);
    G4ThreeVector in(yIn[0],yIn[1],yIn[2]);
    G4ThreeVector out(yOut[0],yOut[1],yOut[2]);

    //G4cout<<"in "<<in<<" mid "<<mid<<" out "<<out<<G4endl;

    chord = G4LineSection::Distline(mid,in,out);


    track.SetCurveLength(curveLengthEnd);
    memcpy(arrayY,yOut.data(),sizeof(G4double)*nvar);
    track.LoadFromArray(arrayY,ncomp);

    return stepper.current_time_step();
}


void BulirschStoerDenseDriver::AccurateAdvance(G4FieldTrack&  track,
                                               G4double hstep,
                                               G4double eps){


    bulirsch_stoer_dense_out<state_type> stepper(0,eps,1,0);
    G4double arrayY[ncomp];
    state_type y;
    track.DumpToArray(arrayY);
    memcpy(y.data(),arrayY,sizeof(G4double)*nvar);
    G4double curveLength = track.GetCurveLength();
    G4double clEnd = curveLength + hstep;
    G4double hrest = hstep;

    G4cout<<"in: ";
    print(y);

    while(hrest > eps){
        stepper.initialize(y,curveLength,hrest);
        stepper.do_step(system);
        curveLength = stepper.current_time();
        hrest = clEnd - curveLength;
        y = stepper.current_state();
    }
    G4cout<<" out: ";
    print(y);

    track.SetCurveLength(curveLength);
    memcpy(arrayY,y.data(),sizeof(G4double)*nvar);
    track.LoadFromArray(arrayY,ncomp);
}

