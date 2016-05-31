#include "RKTestTrack.hh"
#include "G4SystemOfUnits.hh"

G4FieldTrack CreateTrack(){
    
    
    //Setup a field
    G4double
    x_mom = 0.,                 //mom - momentum
    y_mom = 10.,//0.5*proton_mass_c2,
    z_mom = 0.;
    
    G4double mass = 1;
    
    
    G4double momentum = G4ThreeVector(x_mom,y_mom, z_mom).mag();
    //for the G4FieldTrack - momentum, position etc.
    
    G4ThreeVector position(0.,0.,0.);
    
    G4ThreeVector UnitMomentum = G4ThreeVector(x_mom,y_mom, z_mom).unit();
    
    //    G4double mass = 1.67262178e-27*kilogram;
    G4double kineticEnergy = sqrt(sqr(momentum) + sqr(mass)) - mass;
    G4double velocity = momentum/(mass + kineticEnergy);
    G4double labTof = 10.0*ns ;
    G4double properTof = 0.1*ns;
    
    
    //Declare and construct the fieldtrack :
    G4FieldTrack f( position,
                   UnitMomentum,
                   0.0, //for S curve len - starting point
                   kineticEnergy,
                   mass,
                   velocity,
                   labTof,
                   properTof,
                   0    //for Spin
                   );
    
    return f;
    
    
}
