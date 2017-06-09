#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"

using namespace CLHEP;


SteppingAction::SteppingAction(): G4UserSteppingAction()
{}



SteppingAction::~SteppingAction()
{}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
    static G4double time = 0;
    G4double BField = 1*tesla;

    const G4ParticleDefinition* particle = step->GetTrack()->GetParticleDefinition();
    G4double mass = particle->GetPDGMass();

    G4double charge = particle->GetPDGCharge();

    G4ThreeVector momentum = step->GetPreStepPoint()->GetMomentum();
    G4double mom2 = momentum.mag2();
    G4ThreeVector momDirection = step->GetPreStepPoint()->GetMomentumDirection();


    G4double betta2 = mom2/(mass*mass + mom2);
    G4double gamma = 1./sqrt(1-betta2);

    G4double velocity = sqrt(betta2)*c_light;

    //G4cout<<"betta: "<<sqrt(betta2)<<G4endl;

    G4double stepLength = step->GetStepLength();

    G4double timeStep = stepLength/velocity;

    time += timeStep;

    G4double omega = charge*BField/(gamma*mass/c_squared);

    G4double momXZ = sqrt(sqr(momentum.x()) + sqr(momentum.z()))*c_light;

    G4double radius = momXZ/(charge*BField*c_squared);

    //G4cout<<"radius: "<<radius<<G4endl;

    G4ThreeVector pre = step->GetPreStepPoint()->GetPosition();
    G4ThreeVector post = step->GetPostStepPoint()->GetPosition();
    G4ThreeVector calc;


    //calc.setY(pre.y() + velocityVector.y()*timeStep);
    calc.setY((pre.y() + momDirection.y()*stepLength));
    G4double y_error2 = sqr(calc.y() - post.y());


    //calc.setX(pre.x() + radius*sin(omega*timeStep));
    calc.setX(radius*sin(omega*time));
    G4double x_error2 = sqr(calc.x() - post.x());

    //calc.setZ(pre.z() - radius*cos(omega*timeStep));
    calc.setZ(-radius*cos(omega*time));
    G4double z_error2 = sqr(calc.z() - post.z());

    //G4cout<<2*pi/omega*velocityVector.y()<<G4endl;

    G4double error = sqrt((x_error2 + y_error2 + z_error2)/sqr(stepLength));
    //G4cout<<"hstep: "<<stepLength<<" error: "<<error<<G4endl;
    if (error > 1e-5){
        G4cout<<"pre: "<<pre<<"\npost: "<<post<<"\ncalc: "<<calc<<G4endl;
    }


    if (/*error > perMillion*/ false){
        G4cout<<"ycalc: "<<calc.y()<<" yget: "<<post.y()<<" y_error2: "<<y_error2<<G4endl;
        G4cout<<"xcalc: "<<calc.x()<<" xget: "<<post.x()<<" x_error2: "<<x_error2<<G4endl;
        G4cout<<"zcalc: "<<calc.z()<<" zget: "<<post.z()<<" z_error2: "<<z_error2<<G4endl;
    }


    //G4double Ekin = step->GetPreStepPoint()->GetKineticEnergy() * e_SI; //in Dj
    //G4double pMass = step->GetTrack()->GetParticleDefinition()->GetPDGMass()/c_squared*e_SI;
    //G4cout<<"pMass: "<<pMass<<G4endl;
    //electron_charge
}



