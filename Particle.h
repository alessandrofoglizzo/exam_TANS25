#ifndef PARTICLE_H
#define PARTICLE_H

#include "TObject.h"
#include "CYL.h"
#include "VTX.h"
#include "Point.h"
#include "TH1.h"

class Particle: public TObject{
    // class for particles generated in events

    public:
    //constructors
    Particle();
    Particle(double Phi, double Eta);
    Particle(double Phi, double Eta, double Z0);
    Particle(double Phi, double Eta, const VTX& vertex);
    Particle(double Phi, double Eta, const Point& point);
    Particle(const Particle& source);
    virtual ~Particle();

    //assignment operator
    Particle& operator=(const Particle& source);

    //member functions
    double GetPhi() const;
    double GetTheta() const;
    double GetPseudorapidity() const;
    double GetX0() const;
    double GetY0() const;
    double GetZ0() const;
    bool GetVisible() const;
    Particle& Initialize(const VTX& vertex);
    Particle& Initialize(const VTX& vertex, const TH1* eta_dist);
    Particle& MultipleScattering(double pc, const CYL& cylinder, bool MS); //pc is momentum in MeV, beta=1
    Particle& MultipleScattering(double pc, double RadLen, double t, bool MS);
    Particle& Transport(const CYL& cylinder, Point& hit_point);   //function for transport of particle between a point and cylinder
    Particle& Transport(double R, double H, Point& hit_point);

    private:
    //data members
    double fPhi;
    double fEta;
    double fTheta;
    double fX0;
    double fY0;
    double fZ0;
    bool fVisible;  //1 if the particle has intersection, 0 if not

    ClassDef(Particle, 1)

};

#endif