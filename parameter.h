#ifndef PARAMETER_H
#define PARAMETER_H

#include <string>
#include <sstream>


using namespace std;

class Parameter
{
public:
    Parameter();
    virtual void set(double);
    virtual double get(void);
    Parameter& operator=(Parameter);
    Parameter& operator=(double);
    Parameter& operator=(string);
    operator double();
public:
    double A;
    int L;
private:
    double value;
};

#endif // PARAMETER_H
