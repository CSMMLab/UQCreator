#ifndef EULER_H
#define EULER_H

#include "problem.h"

class Euler : public Problem {
private:
public:
    Euler(std::string inputFile);
    virtual void Solve();
    virtual void Plot() const;
    virtual void Print() const;
    virtual void WriteToFile(std::string filename, int filetype) const;
};

#endif // EULER_H
