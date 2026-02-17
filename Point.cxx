#include "Point.h"
#include <cmath>
#include "TRandom3.h"
#include <Riostream.h>

ClassImp(Point);

//________________________________
Point::Point():
    TObject(),
    fR(0.),
    fZ(0.),
    fPhi(0.)
    {
        //Default constructor
    }
//________________________________
Point::Point(double R, double Z, double Phi):
    TObject(),
    fR(R),
    fZ(Z),
    fPhi(Phi)
    {
        //STD constructor
    }
//________________________________
Point::Point(const Point& source):
    TObject(source),
    fR(source.fR),
    fZ(source.fZ),
    fPhi(source.fPhi)
    {
        //copy constructor
    }
//________________________________
Point::~Point(){
    //destructor
}
//________________________________
Point& Point::operator=(const Point& source){
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
double Point::GetX() const{
    return fR*std::cos(fPhi);
}
//________________________________
double Point::GetY() const{
    return fR*std::sin(fPhi);
}
//________________________________
double Point::GetZ() const{
    return fZ;
}
//________________________________
double Point::GetR() const{
    return fR;
}
//________________________________
double Point::GetPhi() const{
    return fPhi;
}
//________________________________
double Point::GetTheta() const{
    return fR*std::atan(fR/fZ);
}
//________________________________
double Point::GetPseudorapidity() const{
    double r = std::sqrt(fR*fR + fZ*fZ); //distance from the origin
    double sinTh = fR/r;
    double cosTh = fZ/r;
    return -std::log(sinTh/(1+cosTh)); //from bisection formulas
}
Point& Point::Set(double R, double Z, double Phi){
    fR = R;
    fZ = Z;
    fPhi = Phi;
    return *this;
}
//________________________________
//____SMEARING FUNCTIONS__________
//________________________________
Point& Point::SmearingZ(double sigmaZ){
    //in-place function for smearing in Z from a gaussian distribution
    fZ += gRandom->Gaus(0., sigmaZ);
    return *this;
}
//________________________________
Point& Point::SmearingPhi(double sigmaAR){
    //in-place function for smearing in RPhi direction, from a gaussian distribution of arches
    double AR = gRandom->Gaus(0., sigmaAR);  //arch: noise on RPhi direction
    fPhi += AR/fR;
    return *this;
}



