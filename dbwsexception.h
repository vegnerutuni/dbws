#ifndef DBWSEXCEPTION_H
#define DBWSEXCEPTION_H

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <iomanip>
#include <fstream>

using namespace std;

class DBWSException {
public:
    DBWSException(string m);
    ~DBWSException();
    const string Show() const
    {
        return msg;
    }
private:
    string msg;
};

#endif // DBWSEXCEPTION_H
