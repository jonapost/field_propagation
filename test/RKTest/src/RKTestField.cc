#include "RKTestField.hh"

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
