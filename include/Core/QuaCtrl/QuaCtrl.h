//
// Created by cc on 4/9/20.
//

#ifndef MSNAVS_QUACTRL_H
#define MSNAVS_QUACTRL_H

#include "PPPLibGlo.h"
#include "GNSS.h"
class QuaCtrl {
public:
    QuaCtrl();
    ~QuaCtrl();

public:
    double MWMeas(Isat_t &isat,ObsData data,double dcb[MAXSAT][MAXCBIASPAIR], int flag);
    double GFMeas(Isat_t &isat,ObsData data, int flag);
    void GFCycSlip(PrcOpt popt,Obss& obss,double dt);
    void MWCycSlip(PrcOpt popt,Obss& obss,double dt,double dcb[MAXSAT][MAXCBIASPAIR]);
    void PriResCheck(double omc_CP, double omc_PR,double& var_PR);
    int PostResCheck(vector<double>res,vector<int>sats,vector<int>frqs,vector<int>types,Obss& obss);
    int PostResQC(PrcOpt popt, Obss& obss,vector<double>v,vector<double>norm_v,int nv,int type);
    void ExcludeSat(ObsData& obs);

public:
    int epoch_;
    vector<int> v_info_;
    vector<int> PRv_info_;
    vector<int> CPv_info_;
    Isat_t *isat_;
};

#endif //MSNAVS_QUACTRL_H
