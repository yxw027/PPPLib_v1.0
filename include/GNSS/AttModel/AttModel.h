//
// Created by cc on 5/21/20.
//

#ifndef GNSS_ATTMODEL_H
#define GNSS_ATTMODEL_H

#include "CmnFunc.h"
#include "GNSS.h"

class AttModel {
public:
    AttModel();
    ~AttModel();

protected:
    int SatYaw(ObsData obs);

public:
    double SatPhw(ObsData& obs, double *rec_pos);

public:
    double exs_[3],eys_[3];
    double phw_;
};


#endif //GNSS_ATTMODEL_H
