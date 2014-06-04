#ifndef TESTUNIFORMELECTRICFIELD_HH
#define TESTUNIFORMELECTRICFIELD_HH

#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include "G4Field.hh"
template
<class Vector>
class TestUniformElectricField 
{
  public: 
    TestUniformElectricField(const Vector FieldVector);
      // A field with value equal to FieldVector.

    TestUniformElectricField(G4double vField,
                           G4double vTheta,
                           G4double vPhi     ) ;
       
    ~TestUniformElectricField() ;
    G4bool DoesFieldChangeEnergy() const {return true;}
    TestUniformElectricField(const TestUniformElectricField &p);
    TestUniformElectricField& operator = (const TestUniformElectricField &p);
      // Copy constructor and assignment operator

    void GetFieldValue(const G4double pos[4], G4double *field) const;

    //TestUniformElectricField* Clone() const;
  private:
  
    G4double fFieldComponents[6] ;
};

#endif
