//
// Created by cc on 4/9/20.
//

#ifndef MSNAVS_TRPMODEL_H
#define MSNAVS_TRPMODEL_H

#include "GNSS.h"
#include "ParSetting.h"

class TrpModel {
public:
    TrpModel();
    ~TrpModel();

protected:
    void StandAtmoPara(Vec3 blh,double hgt,double humi);
    void TrpMapEl(double z);
    void TrpMapGMF(Time t,Vec3 pos,double el);
    void TrpMapNeil(Time t,Vec3 pos,double el);

public:
    double SaasModel(Vec3 blh,double humi);
    virtual double TrpCorr(const Vec3 pos,const double humi);

protected:
    // standard atmosphere parameters
    double P_,T_,e_;

public:
    ObsData *sat_data_;
    Isat_t *isat_;
    Parameter para_;

public:
    double z_wet_,z_hyd_;
    double map_w_,map_h_;
    vector<double> par_X_;
};

class SaasTrp:public TrpModel {
public:
    SaasTrp();
    ~SaasTrp();

public:
    double TrpCorr(Vec3 pos,double humi) override;
};

class EstTrp:public TrpModel {
public:
    EstTrp();
    ~EstTrp();

public:
    double TrpCorr(const Vec3 pos,const double humi);
};


#endif //MSNAVS_TRPMODEL_H
