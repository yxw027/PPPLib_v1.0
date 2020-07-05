//
// Created by cc on 4/22/20.
//

#ifndef MSNAVS_BIAMODEL_H
#define MSNAVS_BIAMODEL_H

#include "CmnFunc.h"
#include "GNSS.h"

typedef struct {
    int sat;
    Time ts,te;
    int type;
    int code;
    double val;
}Osb_t;

typedef struct {
    int nb,nb_max;
    Osb_t *osbs;
}SatOsb_t;

typedef struct {
    double code[MAXSAT][MAXCODE];
    double phase[MAXSAT][MAXCODE];
}Bias_t;

typedef struct {
    Time tmin,tmax;
    double dt;
    Bias_t *bias;
}SatBias_t;

class BiaModel {
public:
    BiaModel();
    ~BiaModel();

public:
    double BDMultiPathCorr();
    int MatchOsb2Sat();
    int AlignOsb2SatBias(SatOsb_t sat_osb);
    double GetTGD();
    double DCBCorr(PrcOpt popt,int f);
    int OSBCorr();
public:
    void InitBiaCorr(ObsData& obs);
    double CBiaCorr(PrcOpt popt,int f,ObsData& obs);
    double PBiaCorr(int f,ObsData& obs);

public:
    ObsData *sat_data_;

public:
    double dcb_[MAXSAT][MAXCBIASPAIR];
    SatBias_t sat_bias_;
    int ac_,type_;
    double cbias_;
};



#endif //MSNAVS_BIAMODEL_H
