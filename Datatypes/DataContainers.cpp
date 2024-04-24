#include "DataContainers.h"

std::string toStr(Metric metric)
{   
    switch (metric){

    case Metric::MSD:
        return("MSD");
    case Metric::BUB:
        return("BUB");
    case Metric::FAI:
        return("FAI");
    case Metric::GLE:
        return("GLE");
    case Metric::JA:
        return("JA");
    case Metric::JT:
        return("JT");
    case Metric::RT:
        return("RT");
    case Metric::RR:
        return("RR");
    case Metric::SM:
        return("SM");
    case Metric::SS1:
        return("SS1");
    case Metric::SS2:
        return("SS2");
    default:
        return("Unknown Metric");
    }

    return("Unknown Metric");
}


