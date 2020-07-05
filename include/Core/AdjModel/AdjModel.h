//
// Created by cc on 4/17/20.
//

#ifndef MSNAVS_ADJMODEL_H
#define MSNAVS_ADJMODEL_H

#include "CmnFunc.h"

class AdjModel {
public:
    AdjModel();
    ~AdjModel();

public:
    int Lsq(const vector<double>& L,const vector<double>& A,int nl,int nx,
            const vector<double>& R,vector<double>& X,vector<double>& Px);
    virtual int Adjustment(const vector<double>& L,const vector<double>& A,int nl,int nx,
                           const vector<double>& R,vector<double>& X,vector<double>& Px);

public:
    int epoch_;
    int num_PR_;
    int num_CP_;
    vector<double> v_;
    vector<double> v_PR_;
    vector<double> v_CP_;
    vector<double> norm_PRv_;
    vector<double> norm_CPv_;
    vector<double> dx_;
    vector<double> Qx_;
    vector<double> Qvv_;
    vector<int>v_info_;
    vector<int>PRv_info_;
    vector<int>CPv_info_;
    double sigma0_;
    // Helmert component covariance estimate for multi-GNSS observation
    int use_sys_;
    int num_obs_[NSYS][NFREQ*2];
    vector<double> hel_R_;
    int hel_iter_;
    vector<double>sgm2;

};

class LSQAdj: public AdjModel {
public:
    LSQAdj();
    ~LSQAdj();

public:
    int Adjustment(const vector<double>& L,const vector<double>& A,int nl,int nx,
                   const vector<double>& R,vector<double>& X,vector<double>& Px);
};

class KFAdj: public AdjModel {
public:
    KFAdj();
    ~KFAdj();

public:
    int Adjustment(const vector<double>& L,const vector<double>& A,int nl,int nx,
                   const vector<double>& R,vector<double>& X,vector<double>& Px);

};

class HelmertAdj:public KFAdj {
public:
    HelmertAdj();
    ~HelmertAdj();

protected:
    int HelmertEst(const vector<double>&A,const vector<double>& Px,int ml,int nx);

public:
    int Adjustment(const vector<double>& L,const vector<double>& A,int nl,int nx,
                   const vector<double>& R,vector<double>& X,vector<double>& Px);
};


#endif //MSNAVS_ADJMODEL_H
