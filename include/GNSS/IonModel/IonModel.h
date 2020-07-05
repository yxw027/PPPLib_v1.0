//
// Created by cc on 4/9/20.
//

#ifndef MSNAVS_IONMODEL_H
#define MSNAVS_IONMODEL_H

#include "GNSS.h"

class Tec {
public:
    Tec();
    ~Tec();

public:
    Time t_;            // epoch time
    int ndata_[3]{};      // TEC grid data size {nlat,nlon,nhgt}
    double re_;         // earth radius
    double lats_[3]{};    // latitude start/interval (deg)
    double lons_[3]{};
    double hgts_[3]{};
    vector<double> data_;// TEC grid data(tecu)
    vector<float> rms_;  // rms value(tecu)
};

class IonModel {
public:
    IonModel();
    ~IonModel();

public:
    double IonMapFun(const Vec3& blh,double el);
    void GetBrdIonPara(double *paras);
    int KlobModel(const Vec3& blh);
    virtual int IonCorr(Vec3 blh);

public:
    double ion_paras_[8];
    ObsData *sat_data_;
    Isat_t *isat_;

public:
    double map_;
    double ion_;
    double ion_var_;
    vector<Tec>tecs_;
};

class KlobIon:public IonModel {
public:
    KlobIon();
    ~KlobIon();

public:
     int IonCorr(Vec3 blh) override;
};

class GIMIon:public IonModel {
public:
    GIMIon();
    ~GIMIon();

protected:
    int DataIndex(int i, int j, int k, const int *ndata);
    int InterpTec(Tec& tec,int k,const double *posp,double& val,double& rms);
    double IonPpp(const Vec3& blh,double re,double hion,double *posp);
    int GIMModel(const Vec3& blh,Tec& tec,double *delay,double *var);

public:
    int IonCorr(Vec3 blh) override;
};

class IFIon:public IonModel {
public:
    IFIon();
    ~IFIon();

public:
    int IonCorr(Vec3 blh) override;
};

#endif //MSNAVS_IONMODEL_H
