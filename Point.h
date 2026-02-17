#ifndef POINT_H
#define POINT_H

#include "TObject.h"

class Point: public TObject{
    // class for points in space, in cylindrical coord.

    public:
    //constructors
    Point();
    Point(double R, double Z, double Phi);
    Point(const Point& Source);
    virtual ~Point();

    //assignment operator
    Point& operator=(const Point& source);

    //member functions
    double GetX() const;
    double GetY() const;
    double GetZ() const;
    double GetR() const;
    double GetPhi() const;
    double GetTheta() const;
    double GetPseudorapidity() const;
    Point& Set(double X, double Y, double Z);
        //smearing functions are in-place: we don't want to save MC truth (HITS)
    Point& SmearingZ(double sigmaZ); 
    Point& SmearingPhi(double sigmaAR);

    //data members
    private:
    double fR;
    double fZ;
    double fPhi;

    ClassDef(Point,1)

};

#endif