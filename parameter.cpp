#include "parameter.h"


ParameterReal::ParameterReal()
{
    value = 0.0;
}

void ParameterReal::set(double p)
{
    value = p;
}

double ParameterReal::get(void)
{
    return value;
}

ParameterReal& ParameterReal::operator=(ParameterReal p)
{
    set(p.get());

    return *this;
}

ParameterReal& ParameterReal::operator=(double p)
{
    set(p);

    return *this;
}

ParameterReal& ParameterReal::operator=(string p)
{
    double tmp;

    stringstream(p) >> tmp;
    set(tmp);

    return *this;
}

ParameterReal::operator double()
{
    return value;
}









ParameterInteger::ParameterInteger()
{
    value = 0;
}

void ParameterInteger::set(int p)
{
    value = p;
}

int ParameterInteger::get(void)
{
    return value;
}

ParameterInteger& ParameterInteger::operator=(ParameterInteger p)
{
    set(p.get());

    return *this;
}

ParameterInteger& ParameterInteger::operator=(int p)
{
    set(p);

    return *this;
}

ParameterInteger& ParameterInteger::operator=(string p)
{
    int tmp;

    stringstream(p) >> tmp;
    set(tmp);

    return *this;
}

ParameterInteger::operator int()
{
    return value;
}

ostream& operator<<(ostream& os, ParameterInteger& p)
{
    os << setw(4) << p.get();
    return os;
}

















Parameter::Parameter()
{
    codeword = 0.0;
    L = 0;
}

Parameter& Parameter::operator=(Parameter p)
{
    set(p.get());

    return *this;
}

Parameter& Parameter::operator=(double p)
{
    set(p);

    return *this;
}

Parameter& Parameter::operator=(string p)
{
    double tmp;

    stringstream(p) >> tmp;
    set(tmp);

    return *this;
}














ParticleAbsorption::ParticleAbsorption()
{
}

void ParticleAbsorption::set(double p)
{
    double tmp;

    tmp = p;

    if (tmp <= 0.0)
    {
        tmp = 1.0;
        cout << "Particle absorption factor is now 1." << endl;
    }
    ParameterReal::set(tmp);
}

ParticleAbsorption& ParticleAbsorption::operator=(ParticleAbsorption p)
{
    set(p.get());

    return *this;
}

ParticleAbsorption& ParticleAbsorption::operator=(double p)
{
    set(p);

    return *this;
}

ParticleAbsorption& ParticleAbsorption::operator=(string p)
{
    double tmp;

    stringstream(p) >> tmp;
    set(tmp);

    return *this;
}

