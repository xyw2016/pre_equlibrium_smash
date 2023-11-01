#ifndef READER_H
#define READER_H
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <chrono>
#include <random>
#include <map>
#include <iomanip>

#include "parameter.h"




class Events{

    public:

    Events();

    ~Events();

    void Readprestage();
    void Readfinal();

    
    TEPaticlelist NIsotauptc;
    EventsofTEPaticlelist allevents;


    private:
    
};


#endif
