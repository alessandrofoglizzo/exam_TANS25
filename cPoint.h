#ifndef CPOINT_H
#define CPOINT_H

#include "TObject.h"

class cPoint: public TObject{
    // class for cPoints in space, in cylindrical coord.

    public:
    //constructors
    cPoint();
    cPoint(double R, double Z, double Phi);
    cPoint(const cPoint& Source);
    virtual ~cPoint();

    //assignment operator
    cPoint& operator=(const cPoint& source);

    //member functions
    double GetX() const;
    double GetY() const;
    double GetZ() const;
    double GetR() const;
    double GetPhi() const;
    double GetTheta() const;
    double GetPseudorapidity() const;
    cPoint& Set(double X, double Y, double Z);
        //smearing functions are in-place: we don't want to save MC truth (HITS)
    cPoint& SmearingZ(double sigmaZ); 
    cPoint& SmearingPhi(double sigmaAR);

    //data members
    private:
    double fR;
    double fZ;
    double fPhi;

    ClassDef(cPoint,1)

};

#endif