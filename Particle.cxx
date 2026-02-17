#include "Particle.h"
#include "TRandom3.h"
#include <cmath>

namespace {
    const double kPiPART = 3.14159265358979323846;
}

ClassImp(Particle);

//_________________________________
Particle::Particle():
    TObject(),
    fPhi(0.),
    fEta(0.),
    fTheta(kPiPART/2),
    fX0(0.),
    fY0(0.),
    fZ0(0.),
    fVisible(true)
    {
        //default constructor
    }
//_________________________________
Particle::Particle(double Phi, double Eta):
    TObject(),
    fPhi(Phi),
    fEta(Eta),
    fTheta(2*std::atan(std::exp(-fEta))),
    fX0(0.),
    fY0(0.),
    fZ0(0.),
    fVisible(true)
    {
        //standard constructor
    }
//_________________________________
Particle::Particle(double Phi, double Eta, double Z0):
    TObject(),
    fPhi(Phi),
    fEta(Eta),
    fTheta(2*std::atan(std::exp(-fEta))),
    fX0(0.),
    fY0(0.),
    fZ0(Z0),
    fVisible(true)
    {
        //tracklets constructor
    }
//_________________________________
Particle::Particle(double Phi, double Eta, const VTX& vertex):
TObject(),
    fPhi(Phi),
    fEta(Eta),
    fTheta(2*std::atan(std::exp(-fEta))),
    fX0(vertex.GetX()),
    fY0(vertex.GetY()),
    fZ0(vertex.GetZ()),
    fVisible(true)
    {
        //constructor by vertex
    }
//_________________________________
Particle::Particle(double Phi, double Eta, const Point& point):
    TObject(),
    fPhi(Phi),
    fEta(Eta),
    fTheta(2*std::atan(std::exp(-fEta))),
    fX0(point.GetX()),
    fY0(point.GetY()),
    fZ0(point.GetZ()),
    fVisible(true)
    {
        //constructor by point
    }
//_________________________________
Particle::Particle(const Particle& source):
    TObject(source),
    fPhi(source.fPhi),
    fEta(source.fEta),
    fX0(source.fX0),
    fY0(source.fY0),
    fZ0(source.fZ0),
    fVisible(source.fVisible)
    {
        //copy constructor
    }
//_________________________________
Particle::~Particle(){
    //destructor
}
//_________________________________
Particle& Particle::operator=(const Particle& source){
    if(this==&source) return *this;

    TObject::operator=(source);
    fPhi = source.fPhi;
    fEta = source.fEta;
    fTheta = source.fTheta;
    fX0 = source.fX0;
    fY0 = source.fY0;
    fZ0 = source.fZ0;
    fVisible = source.fVisible;
    return *this;
}
//_________________________________
//___________MEMBER FUNCTIONS______
//_________________________________
double Particle::GetPhi() const{
    return fPhi;
}
//_________________________________
double Particle::GetPseudorapidity() const{
    return fEta;
}
//_________________________________
double Particle::GetTheta() const{
    return fTheta;
}
//_________________________________
double Particle::GetX0() const{
    return fX0;
}
//_________________________________
double Particle::GetY0() const{
    return fY0;
}
//_________________________________
double Particle::GetZ0() const{
    return fZ0;
}
//_________________________________
bool Particle::GetVisible() const{
    return fVisible;
}
//___________________________________
Particle& Particle::Initialize(const VTX& vertex){
    //function to initialize particles efficiently
    fX0 = vertex.GetX();
    fY0 = vertex.GetY();
    fZ0 = vertex.GetZ();
    fPhi = gRandom->Uniform(2*kPiPART);
    fEta = gRandom->Uniform(-2, 2);
    fTheta = 2*std::atan(std::exp(-fEta));
    fVisible = true;
    return *this;
}
//_______________________________________________
Particle& Particle::Initialize(const VTX& vertex, const TH1* eta_dist){
    //function to initialize particles efficiently with assigned distribution
    fX0 = vertex.GetX();
    fY0 = vertex.GetY();
    fZ0 = vertex.GetZ();
    fPhi = gRandom->Uniform(2*kPiPART);
    if (eta_dist){
        fEta = eta_dist->GetRandom();
    } else fEta = gRandom->Uniform(-2, 2); //if the histo is null, use uniform
    fTheta = 2*std::atan(std::exp(-fEta));
    fVisible = true;
    return *this;
}
//_________________________________
Particle& Particle::MultipleScattering(double pc, const CYL& cylinder, bool MS){ //pc momentum in MeV
    if(fVisible && MS){  //make MS if the particle hits the layer (fVisible) and if we turn on MS (MS=1)
    double z=1; //unitary charged particle
    double RadLen = cylinder.GetX0();
    double t = cylinder.Gett();
    double ThetaRMS = ((13.6)/pc)*z*std::sqrt(t/RadLen)*(1+0.038*std::log(t/RadLen));
    //extract multiple scattering angles in 
    double Thp = gRandom->Gaus(0., ThetaRMS);
    double Php = gRandom->Uniform(2*kPiPART);

    double RM[3][3];   //rotation matrix
    double up[3]; // versor in particle RS
    double u[3];  // versor in LAB after MS

    RM[0][0] = -std::sin(fPhi);
    RM[0][1] = -std::cos(fPhi)*std::cos(fTheta);
    RM[0][2] = std::cos(fPhi)*std::sin(fTheta);
    
    RM[1][0] = std::cos(fPhi);
    RM[1][1] = -std::sin(fPhi)*std::cos(fTheta);
    RM[1][2] = std::sin(fPhi)*std::sin(fTheta);

    RM[2][0] = 0.;
    RM[2][1] = std::sin(fTheta);
    RM[2][2] = std::cos(fTheta);


    up[0] = std::sin(Thp)*std::cos(Php);
    up[1] = std::sin(Thp)*std::sin(Php);
    up[2] = std::cos(Thp);

    for(int i=0; i<3; i++){
        u[i]=0.;
        for(int j=0; j<3; j++){
            u[i] += RM[i][j]*up[j];
        }
    }

    //get final angles
    double Phi_final = std::atan2(u[1],u[0]);
    if(Phi_final < 0) Phi_final += 2*kPiPART;

    fPhi = Phi_final;
    fTheta = std::acos(u[2]);
    fEta = -std::log(std::tan(fTheta/2));
    }
    return *this;
}
//_________________________________
Particle& Particle::MultipleScattering(double pc, double RadLen, double t, bool MS){ //pc momentum in MeV
    if(fVisible && MS){  //make MS if the particle hits the layer (fVisible) and if we turn on MS (MS=1)
    double z=1; //unitary charged particle
    double ThetaRMS = ((13.6)/pc)*z*std::sqrt(t/RadLen)*(1+0.038*std::log(t/RadLen));
    //extract multiple scattering angles in 
    double Thp = gRandom->Gaus(0., ThetaRMS);
    double Php = gRandom->Uniform(2*kPiPART);

    double RM[3][3];   //rotation matrix
    double up[3]; // versor in particle RS
    double u[3];  // versor in LAB after MS

    RM[0][0] = -std::sin(fPhi);
    RM[0][1] = -std::cos(fPhi)*std::cos(fTheta);
    RM[0][2] = std::cos(fPhi)*std::sin(fTheta);
    
    RM[1][0] = std::cos(fPhi);
    RM[1][1] = -std::sin(fPhi)*std::cos(fTheta);
    RM[1][2] = std::sin(fPhi)*std::sin(fTheta);

    RM[2][0] = 0.;
    RM[2][1] = std::sin(fTheta);
    RM[2][2] = std::cos(fTheta);


    up[0] = std::sin(Thp)*std::cos(Php);
    up[1] = std::sin(Thp)*std::sin(Php);
    up[2] = std::cos(Thp);

    for(int i=0; i<3; i++){
        u[i]=0.;
        for(int j=0; j<3; j++){
            u[i] += RM[i][j]*up[j];
        }
    }

    //get final angles
    double Phi_final = std::atan2(u[1],u[0]);
    if(Phi_final < 0) Phi_final += 2*kPiPART;

    fPhi = Phi_final;
    fTheta = std::acos(u[2]);
    fEta = -std::log(std::tan(fTheta/2));
    }
    return *this;
}
//_____________________________________________________
Particle& Particle::Transport(const CYL& cylinder, Point& hit_point){
    double c1 = std::sin(fTheta)*std::cos(fPhi);
    double c2 = std::sin(fTheta)*std::sin(fPhi);
    double c3 = std::cos(fTheta);
    double R = cylinder.GetR();
    double SqrtDelta = std::sqrt((fX0*c1 + fY0*c2)*(fX0*c1 + fY0*c2) - (c1*c1+c2*c2)*(fX0*fX0+fY0*fY0 - R*R));
    double t = (-(fX0*c1 + fY0*c2) + SqrtDelta)/(c1*c1 + c2*c2);

    double H = cylinder.GetH();
    //check that particle has an intersection with cylinder and
    //assign it to the point
    if((fZ0+c3*t >= -H/2) && (fZ0+c3*t <= H/2)){
        double XHit = fX0+c1*t;
        double YHit = fY0+c2*t;
        double ZHit = fZ0+c3*t;
        double RHit = std::sqrt(XHit*XHit + YHit*YHit);
        double PhiHit = std::atan2(YHit, XHit);
        if (PhiHit<0) PhiHit += 2*kPiPART;

        hit_point.Set(RHit, ZHit, PhiHit); //save hit point
        // update new initial point for next transport
        fX0 = fX0+c1*t;
        fY0 = fY0+c2*t;
        fZ0 = fZ0+c3*t;
    } else fVisible = false;
    return *this;
}
//_____________________________________________________
Particle& Particle::Transport(double R, double H, Point& hit_point){
    double c1 = std::sin(fTheta)*std::cos(fPhi);
    double c2 = std::sin(fTheta)*std::sin(fPhi);
    double c3 = std::cos(fTheta);
    double SqrtDelta = std::sqrt((fX0*c1 + fY0*c2)*(fX0*c1 + fY0*c2) - (c1*c1+c2*c2)*(fX0*fX0+fY0*fY0 - R*R));
    double t = (-(fX0*c1 + fY0*c2) + SqrtDelta)/(c1*c1 + c2*c2);

    //check that particle has an intersection with cylinder and
    //assign it to the point
    if((fZ0+c3*t >= -H/2) && (fZ0+c3*t <= H/2)){
        double XHit = fX0+c1*t;
        double YHit = fY0+c2*t;
        double ZHit = fZ0+c3*t;
        double RHit = std::sqrt(XHit*XHit + YHit*YHit);
        double PhiHit = std::atan2(YHit, XHit);
        if (PhiHit<0) PhiHit += 2*kPiPART;

        hit_point.Set(RHit, ZHit, PhiHit); //save hit point
        // update new initial point for next transport
        fX0 = fX0+c1*t;
        fY0 = fY0+c2*t;
        fZ0 = fZ0+c3*t;
    } else fVisible = false;
    return *this;
}