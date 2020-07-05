//
// Created by cc on 4/9/20.
//

#ifndef MSNAVS_TIDMODEL_H
#define MSNAVS_TIDMODEL_H

typedef struct {
    double mjd;
    double xp,yp;
    double xpr,ypr;
    double ut1_utc;
    double lod;
}Erpd_t;

class TidModel {
public:
    TidModel();
    ~TidModel();

protected:
    Time tut_;
    double re_;
    Vec3 blh_;
    double E_[9];
    double tid_dr_[3][3];
    double denu_[2][3];
    double sun_pos_[3];
    double moon_pos_[3];
    double gmst_;
    double erp_val_[5];

protected:
    int GetErpVal(Time t);
    void TideSl(const double *eu,const double *rp,double GMp,Vec3 blh,double *dr);
    void TidSolid();
    void TidOcean(int rov_bas);
    void IersMeanPole(double *xp,double *yp);
    void TidPole();

public:
    void TidCorr(Time time,int rov_bas,double *xyz,double *dr);

public:
    int tid_opt_;

    vector<Erpd_t> erp_;
    double ocean_par_[2][6*11];
};


#endif //MSNAVS_TIDMODEL_H
