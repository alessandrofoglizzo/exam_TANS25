#include "cPoint.h"
#include <cmath>
#include "TRandom3.h"
#include <Riostream.h>

ClassImp(cPoint);

//________________________________
cPoint::cPoint():
    TObject(),
    fR(0.),
    fZ(0.),
    fPhi(0.)
    {
        //Default constructor
    }
//________________________________
cPoint::cPoint(double R, double Z, double Phi):
    TObject(),
    fR(R),
    fZ(Z),
    fPhi(Phi)
    {
        //STD constructor
    }
//________________________________
cPoint::cPoint(const cPoint& source):
    TObject(source),
    fR(source.fR),
    fZ(source.fZ),
    fPhi(source.fPhi)
    {
        //copy constructor
    }
//________________________________
cPoint::~cPoint(){
    //destructor
}
//________________________________
cPoint& cPoint::operator=(const cPoint& source){
    if(this==&source) return *this;

    TObject::operator=(source);
    fR = source.fR;
    fZ = source.fZ;
    fPhi = source.fPhi;
    return *this;
}
//________________________________
//__MEMBER FUNCTIONS______________
//________________________________
double cPoint::GetX() const{
    return fR*std::cos(fPhi);
}
//________________________________
double cPoint::GetY() const{
    return fR*std::sin(fPhi);
}
//________________________________
double cPoint::GetZ() const{
    return fZ;
}
//________________________________
double cPoint::GetR() const{
    return fR;
}
//________________________________
double cPoint::GetPhi() const{
    return fPhi;
}
//________________________________
double cPoint::GetTheta() const{
    return fR*std::atan(fR/fZ);
}
//________________________________
double cPoint::GetPseudorapidity() const{
    double r = std::sqrt(fR*fR + fZ*fZ); //distance from the origin
    double sinTh = fR/r;
    double cosTh = fZ/r;
    return -std::log(sinTh/(1+cosTh)); //from bisection formulas
}
cPoint& cPoint::Set(double R, double Z, double Phi){
    fR = R;
    fZ = Z;
    fPhi = Phi;
    return *this;
}
//________________________________
//____SMEARING FUNCTIONS__________
//________________________________
cPoint& cPoint::SmearingZ(double sigmaZ){
    //in-place function for smearing in Z from a gaussian distribution
    fZ += gRandom->Gaus(0., sigmaZ);
    return *this;
}
//________________________________
cPoint& cPoint::SmearingPhi(double sigmaAR){
    //in-place function for smearing in RPhi direction, from a gaussian distribution of arches
    double AR = gRandom->Gaus(0., sigmaAR);  //arch: noise on RPhi direction
    fPhi += AR/fR;
    return *this;
}



