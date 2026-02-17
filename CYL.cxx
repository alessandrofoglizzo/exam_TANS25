#include <Riostream.h>

#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include <cmath>
//#include "TGeoTube.h"
#include "CYL.h" 

ClassImp(CYL);

//__________________________
CYL::CYL(): TObject(),
            fName(""),
            fR(0.),
            fH(0.),
            ft(0.),
            fMat(nullptr)
{
    //cout << "DEFAULT CONSTRUCTOR - this = " << this << endl;
}

//________________________________
CYL::CYL(const char* cylname, double R, double H, double t, const char* element):
    TObject(),
    fName(cylname),
    fR(R),
    fH(H),
    ft(t)
{
    //cout << "STANDARD CONSTRUCTOR - this = " << this << endl;
    //create a gometry manager if it dosent exist
    if(!gGeoManager) new TGeoManager("manager", "cylinders geometry");

    //select and searhces for material in the manager
    fMat = gGeoManager->GetMaterial(element);
    //if the material has not been created yet, it creates it
    if(!fMat){
        if (strcmp(element, "Si")==0) fMat = new TGeoMaterial("Si", 28.08, 14, 2.33);
        else if (strcmp(element, "Be")==0) fMat = new TGeoMaterial("Be", 9.01, 4, 1.85);
        else {
            cerr << "ERROR! Invalid Material: select Si or Be!" << endl;
            fMat = nullptr;
        }
    }
}

//______________________________________________________
CYL::CYL(const CYL& source):
    TObject(source),
    fName(source.fName),
    fR(source.fR),
    fH(source.fH),
    ft(source.ft),
    fMat(source.fMat)
{
   // cout << "COPY CONSTRUCTOR - this = " << this << endl;
}


CYL::~CYL() {
    //cout << "DESTRUCTOR - this = " << this << endl;
}

//__________________________________________
CYL& CYL::operator=(const CYL& source){
    //Assignment operator
    if(this == &source) return *this;

    TObject::operator=(source);
    fName = source.fName;
    fR = source.fR;
    fH = source.fH;
    ft = source.ft;
    fMat = source.fMat;
    return *this;
}
//_________________________________
double CYL::GetR() const {
    return fR;
}
//_________________________________
double CYL::GetH() const {
    return fH;
}
//_________________________________
int CYL::GetZ() const {
    if (fMat) return fMat->GetZ();
    return 0;
}
//_________________________________
double CYL::Gett() const {
    return ft;
}
double CYL::GetX0() const{
    if (fMat) return fMat->GetRadLen(); //radiation lenght in cm
    return 0;
}
double CYL::GetThetaRMS(double pc) const{
    double z=1; //unitary charged particle
    double RadLen = fMat->GetRadLen();
    double x = ft;
    return ((13.6)/pc)*z*std::sqrt(x/RadLen)*(1+0.038*std::log(x/RadLen));
}


