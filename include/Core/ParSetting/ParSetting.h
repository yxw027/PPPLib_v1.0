//
// Created by cc on 4/9/20.
//

#ifndef MSNAVS_PARSETTING_H
#define MSNAVS_PARSETTING_H

#include "OutSol.h"
#include "GNSS.h"
#include "IonModel.h"

class Parameter {
public:
    Parameter();
    ~Parameter();

public:
    int NumUseFrq();
    int NumPosPar();
    int NumClkPar();
    int NumSPPPar();
    int NumDcbPar();
    int NumIfbPar();
    int NumGloPar();
    int NumTrpPar();
    int NumIonPar();
    int NumRPar();
    int NumAmbPar();

    int IdxClkPar(int i);
    int IdxDcbPar();
    int IdxIfbPar(int i);
    int IdxGloPar();
    int IdxTrpPar();
    int IdxIonPar(int sat);
    int IdxAmbPar(int sat,int f);

    void CoePosPar(vector<double>& A_coe);
    void CoeClkPar(vector<double>& A_coe);
    void CoeDcbPar(vector<double>& A_coe);
    void CoeIfbPar(vector<double>& A_coe);
    void CoeGloPar(vector<double>& A_coe);
    void CoeTrpPar(vector<double>& A_coe);
    void CoeIonPar(vector<double>& A_coe,int f,int type);
    void CoeAmbPar(int f,vector<double>& A_coe);

    int InitX(double xi,double var,int i);
    void PosParUpdate();
    void ClkParUpdate();
    void DcbParUpdate();
    void IfbParUpdate();
    void GloParUpdate();
    void TrpParUpdate();
    void IonParUpdate();
    void AmbParUpdate();

public:
    PrcOpt popt_;
    int nx_;
    int num_L_;
    int epoch_;
    int exist_sys_mark_[NSYS];
    ObsData sat_info_;
    Obss obss_;
    Isat_t *isat_;
    GIMIon gim_model_;
    int *glo_fcn_;
    double tt_;
    double *x_,*Px_;
    SolStat last_sol_;
};

#endif //MSNAVS_PARSETTING_H
