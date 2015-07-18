#include "RKTestField.hh"
using namespace std;
using namespace CLHEP;


G4QuadrupoleMagField SetupQMF(G4double gradient ){
    
    G4QuadrupoleMagField myField(gradient);
    return myField;
}

G4UniformMagField SetupUMF(G4double x_field,
                           G4double y_field,
                           G4double z_field){
    
    G4UniformMagField myField(G4ThreeVector(x_field, y_field, z_field));
    return myField;
}

G4double toStepLen(G4double theta, G4double R){
    if(R == 0){
        std::cerr<<"Error";
        return 0;
    }
    
    return theta*R;
    
}

G4double findRFrom2pt(G4double yOut2[], G4double yOut1[], G4double step_len){
    G4double dy, dx, phi, R, theta;
    
    dx = yOut2[0] - yOut1[0];
    dy = yOut2[1] - yOut1[1];
    phi = atan(dy/dx);
    theta = pi - 2.0*phi;
    R = step_len/theta;
    return R;
}