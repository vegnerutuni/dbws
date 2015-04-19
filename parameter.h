#ifndef PARAMETER_H
#define PARAMETER_H

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;

class ParameterReal
{
public:
    ParameterReal();
    virtual void set(double);
    virtual double get(void);
    ParameterReal& operator=(ParameterReal);
    ParameterReal& operator=(double);
    ParameterReal& operator=(string);
    operator double();
private:
    double value;
};

class ParameterInteger
{
public:
    ParameterInteger();
    virtual void set(int);
    virtual int get(void);
    ParameterInteger& operator=(ParameterInteger);
    ParameterInteger& operator=(int);
    ParameterInteger& operator=(string);
    friend ostream& operator<<(ostream&, ParameterInteger&);
    operator int();
private:
    int value;
};

class Parameter:public ParameterReal
{
public:
    Parameter();
    Parameter& operator=(Parameter);
    Parameter& operator=(double);
    Parameter& operator=(string);
public:
    double codeword;
    int L;
};


// Particle absorption factor for phase

class ParticleAbsorption:public ParameterReal
{
public:
    ParticleAbsorption();
    void set(double);
    ParticleAbsorption& operator=(ParticleAbsorption);
    ParticleAbsorption& operator=(double);
    ParticleAbsorption& operator=(string);
public:
    double A;
    int L;
};

#endif // PARAMETER_H
