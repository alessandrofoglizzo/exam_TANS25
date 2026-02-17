#include "VTX.h"
#include <cmath>
#include "TRandom3.h"
#include <Riostream.h>

namespace {
    const double kPiVTX = 3.14159265358979323846;
}

ClassImp(VTX);

//___________________________
VTX::VTX():
    TObject(),
    fX(0.),
    fY(0.),
    fZ(0.),
    fMult(1)
    {
        // default constructor
    }
//___________________________
VTX::VTX(double X, double Y, double Z, int mean_mult):
    TObject(),
    fX(X),
    fY(Y),
    fZ(Z),
    fMult(gRandom->Exp(mean_mult))
    {
        // std constructor
    }
//____________________________________________________
VTX::VTX(const VTX& source):
    TObject(source),
    fX(source.fX),
    fY(source.fY),
    fZ(source.fZ),
    fMult(source.fMult)
    {
        //copy constructor
    }
//_____________________________________________________
VTX::~VTX(){
    //destructor
}
//_____________________________________________________
VTX& VTX::operator=(const VTX& source){
    if(this==&source) return *this;

    TObject::operator=(source);
    fX = source.fX;
    fY = source.fY;
    fZ = source.fZ;
    fMult = source.fMult;
    return *this;
}
//______________________________________________________
//______________________MEMBER FUNCTIONS________________
//______________________________________________________
double VTX::GetX() const{
    return fX;
}
//______________________________________________________
double VTX::GetY() const{
    return fY;
}
//______________________________________________________
double VTX::GetZ() const{
    return fZ;
}
//______________________________________________________
double VTX::GetR() const{
    return std::sqrt(fX*fX + fY*fY);
}
//______________________________________________________
double VTX::GetPhi() const{
    double Phi = std::atan2(fY, fX);
    if (Phi<0) Phi += 2*kPiVTX;
    return Phi;
}
//______________________________________________________
double VTX::GetTheta() const{
    double r = std::sqrt(fX*fX + fY*fY + fZ*fZ);
    return std::acos(fZ/r);
}
//______________________________________________________
int VTX::GetMult() const{
    return fMult;
}
//______________________________________________________
VTX& VTX::Initialize(int mean_mult){
    fX = gRandom->Gaus(0, 0.01);
    fY = gRandom->Gaus(0, 0.01);
    fZ = gRandom->Gaus(0, 5.3);
    fMult = gRandom->Poisson(mean_mult);
    return *this;
}
//______________________________________________________
VTX& VTX::Initialize(int mean_mult, const TH1* mult_dist){
    fX = gRandom->Gaus(0, 0.01);
    fY = gRandom->Gaus(0, 0.01);
    fZ = gRandom->Gaus(0, 5.3);
    if (mult_dist){
        fMult = std::round(2*mult_dist->GetRandom());    //histo is given in 2 units of rapidity, we use 4: eta(-2, 2)
    } else fMult = gRandom->Poisson(mean_mult); //if the histo is null
    return *this;
}

