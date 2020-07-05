//
// Created by cc on 5/18/20.
//

#ifndef GNSS_AMBMODEL_H
#define GNSS_AMBMODEL_H

#include "PPPLibGlo.h"

typedef struct {
    double amb[NFREQ];
    double LC[NFREQ];  //L1_L2 L1_L3 L2_L3
    double MW;         //L1_L2
}Amb_t;


class AmbModel {
public:
    AmbModel();
    ~AmbModel();

public:
    double LcAmb(ObsData& obs);


};


#endif //GNSS_AMBMODEL_H
