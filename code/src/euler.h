#ifndef EULER_H
#define EULER_H

#include "problem.h"

class Euler : public Problem {
private:
public:
    Euler(std::string inputFile);
    virtual void solve();
    virtual void plot() const;
    virtual void print() const;
    virtual void writeToFile(std::string filename, int filetype) const;
};

#endif // EULER_H
