#ifndef CYL_H
#define CYL_H

#include "TObject.h"
#include "TString.h"
//#include "TGeoVolume.h"
#include "TGeoMaterial.h"

class CYL : public TObject{
    //   Class for cylindrical objects   //

    public:
    CYL();        //default constructor
    CYL(const char* cylname, double R, double H, double t, const char* element);
    CYL(const CYL& source);
    virtual ~CYL();        //destructor

    CYL& operator=(const CYL& source);

    //  member functions
    double GetR() const;
    double GetH() const;
    int GetZ() const;
    double Gett() const;
    double GetX0() const; //get radiation lenght in cm
    double GetThetaRMS(double pc) const; // get ThetaRMS of MS

    //
    private:
    //
    TString fName;
    double fR;
    double fH;
    double ft;
    TGeoMaterial* fMat;


    ClassDef(CYL,1)

};


#endif