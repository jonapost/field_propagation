#ifndef TUniformMagField_HH
#define TUniformMagField_HH

#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include "G4MagneticField.hh"
#include "G4PhysicalConstants.hh"

class TUniformMagField : public G4MagneticField
{
    public:  // with description

        TUniformMagField(const G4ThreeVector& FieldVector ) 
            // A field with value equal to FieldVector.
        {
            fFieldComponents[0] = FieldVector.x();
            fFieldComponents[1] = FieldVector.y();
            fFieldComponents[2] = FieldVector.z();
        }


        TUniformMagField(G4double vField,
                G4double vTheta,
                G4double vPhi     ) 
        {
            if ( (vField<0) || (vTheta<0) || (vTheta>pi) || (vPhi<0) || (vPhi>twopi) )
            {
                G4Exception("TUniformMagField::TUniformMagField()",
                        "GeomField0002", FatalException, "Invalid parameters.") ;
            }
            fFieldComponents[0] = vField*std::sin(vTheta)*std::cos(vPhi) ;
            fFieldComponents[1] = vField*std::sin(vTheta)*std::sin(vPhi) ;
            fFieldComponents[2] = vField*std::cos(vTheta) ;
        }



        virtual ~TUniformMagField() {;}


        TUniformMagField(const TUniformMagField &p)
            : G4MagneticField(p)
        {
            for (G4int i=0; i<3; i++)
                fFieldComponents[i] = p.fFieldComponents[i];
        }

        TUniformMagField& operator = (const TUniformMagField &p)
            // Copy constructor and assignment operator.
        {
            if (&p == this) return *this;
            for (G4int i=0; i<3; i++)
                fFieldComponents[i] = p.fFieldComponents[i];
            return *this;
        }

        inline void GetFieldValue(const G4double yTrack[4],
                G4double *B) const 
        {
            B[0]= fFieldComponents[0] ;
            B[1]= fFieldComponents[1] ;
            B[2]= fFieldComponents[2] ;
        }

        void SetFieldValue(const G4ThreeVector& newFieldVector)
        {
            fFieldComponents[0] = newFieldVector.x();
            fFieldComponents[1] = newFieldVector.y();
            fFieldComponents[2] = newFieldVector.z();
        }

        G4ThreeVector GetConstantFieldValue() const
        {
            G4ThreeVector B(fFieldComponents[0],
                    fFieldComponents[1],
                    fFieldComponents[2]);
            return B;
        }
        // Return the field value

        virtual TUniformMagField* Clone() const
        { 
            return new TUniformMagField( G4ThreeVector(this->fFieldComponents[0],
                        this->fFieldComponents[1],
                        this->fFieldComponents[2]) );
        }

    private:

        G4double fFieldComponents[3] ;
};

#endif
