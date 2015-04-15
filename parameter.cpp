#include "parameter.h"

Parameter::Parameter()
{
    value = 0.0;
    A = 0.0;
    L = 0;
}

void Parameter::set(double p)
{
    value = p;
}

double Parameter::get(void)
{
    return value;
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

Parameter::operator double()
{
    return value;
}
