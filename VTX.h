#ifndef VTX_H
#define VTX_H

#include "TObject.h"
#include "TH1.h"

class VTX: public TObject{
    // class for vertex

    public:
    VTX();
    VTX(double X, double Y, double Z, int mean_mult);
    //VTX(double X, double Y, double Z, int mult, const char* distr);
    VTX(const VTX& source);
    virtual ~VTX();

    VTX& operator=(const VTX& source);

    //member functions
    double GetX() const;
    double GetY() const;
    double GetZ() const;
    double GetR() const;
    double GetPhi() const;
    double GetTheta() const;
    int GetMult() const;
    VTX& Initialize(int mean_mult);
    VTX& Initialize(int mean_mult, const TH1* mult_dist);

    //
    private:
    double fX;
    double fY;
    double fZ;
    int fMult;   //multiplicity of event

    //
    ClassDef(VTX, 1)

};

#endif