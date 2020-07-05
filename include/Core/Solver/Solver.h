//
// Created by cc on 4/9/20.
//

#ifndef MSNAVS_SOLVER_H
#define MSNAVS_SOLVER_H


#include "EphModel.h"
#include "PPPLibGlo.h"
#include "ReadFiles.h"
#include "GNSS.h"
#include "TrpModel.h"
#include "IonModel.h"
#include "AttModel.h"
#include "BiaModel.h"
#include "ParSetting.h"
#include "AdjModel.h"
#include "QuaCtrl.h"
#include "OutSol.h"

class Solver {
public:
    Solver();
    ~Solver();

public:
    void GetUseFrq(PrcOpt popt,ObsData& obs,double *lam);
    void ReSetIsat(int sat);
    void UpdateSatInfo(PrcOpt popt);
    void SetSysMask();
    double Dops();
    void UpdateSol(vector<double>x,vector<double>Px,double zhtrp,double zwtrp);
    int Fullx2Zipx(vector<double>full_x,vector<double>full_Px,vector<double>&zip_x,vector<double>&zip_Px,int full_nx);
    void Zipx2Fullx(vector<double>zip_x,vector<double>zip_Px,vector<double>&full_x,vector<double>&full_Px,int full_nx);
    double MeasVar(PrcOpt popt,int type);
    void NormalEqLine(PrcOpt popt,int use_f,int idx_f,int nl,int nx,int type);
    double GetCorrObs(PrcOpt popt,int f, int type);
    void CorrPRMeas(PrcOpt popt,ObsData &obs,double dantr[NFREQ],double dants[NFREQ]);
    void CorrCPMeas(PrcOpt popt,ObsData &obs,double dantr[NFREQ],double dants[NFREQ],double phw);
    double RecClkCorr(vector<double> est_x);
    double RecIfbCorr(vector<double> est_x);
    double GloIcbCorr(PrcOpt popt,vector<double> est_x,int type);
    int SortObs(Obss &obs);
    virtual int Algorithm();

public:
    Nav nav_;
    Obss epoch_rover_;
    Obss epoch_base_;
    Parameter  para_;
    BiaModel   bia_corr_;
    AntModel   ant_corr_;
    TidModel   tid_corr_;
    ReadGIM    reader_gim_;
    QuaCtrl    qc_;
    SolStat    epoch_sol_;

public:
    int site_count_;
    int epoch_count_;
    int exist_sys_mask_[NSYS];
    double dt_;
    ObsData *sig_sat_data_;
    Isat_t  isat_[MAXSAT];
    int num_L_;
    int num_PR_,num_CP_;
    int vaild_sat_num_;
    vector<double>A_coe_,OMCs_,R_meas_,R_vec_;
    vector<int>L_info_;
};

class SPPSolver: public Solver {
public:
    SPPSolver();
    ~SPPSolver();

private:
    double ModPRMeas(int post);
    int PRErrEq(PrcOpt popt,int use_f,int idx_f,int iter,int post);

protected:
    int SPPPostRes(PrcOpt popt,AdjModel adj_model);
    void MakeSPPEqs(int post, int iter,PrcOpt popt);
    int SolverSPP(PrcOpt popt);
    int InitSPPSolver(PrcOpt popt);
    int SPPAlgorithm(PrcOpt popt);

public:
    virtual int Algorithm();

private:
    SatEph   *brd_eph_;
    TrpModel *spp_trp_;
    IonModel *spp_ion_;
    AdjModel *spp_adj_;

private:
    PrcOpt popt_;
    int full_nx_,zip_nx_;
    vector<double> full_x_,full_Px_;
    vector<double> zip_x_,zip_Px_;
};

class PPPSolver: public SPPSolver {
public:
    PPPSolver();
    ~PPPSolver();

private:
    SatEph   *pre_eph_;
    TrpModel *trp_corr_;
    IonModel *ion_corr_;
    AttModel  att_corr_;
    AdjModel *ppp_adj_;

private:
    int PPPPostRes(AdjModel adj_model,vector<double>x);
    void UpdateIsat();
    double ModPRMeas(int use_f,int idx_f,vector<double> x);
    double ModCPMeas(int use_f,int idx_f,vector<double> x);
    int PRErrEq(PrcOpt popt,int use_f,int idx_f,int iter,int post,vector<double> x);
    int CPErrEq(PrcOpt popt,int use_f,int idx_f,int iter,int post,vector<double> x);
    int MakePPPEqs(int post,int iter,vector<double> x);
    void ParTimeUpdate();
    void PPPCycSlip();
    int SolverPPP();
    int P3SPP();
    int InitPPPSolver();
    int PPPAlgorithm();

public:
    int Algorithm();

private:
    double rr_[3];
    PrcOpt popt_;
    int full_nx_,zip_nx_;
    vector<double> full_x_,full_Px_;
    vector<double> zip_x_,zip_Px_;
    vector<int> lar_res_sat_;
    vector<int> lar_res_frq_;
    vector<int> lar_res_type_;
    vector<double> lar_res_;
};

class MainSolver {
public:
    MainSolver();
    ~MainSolver();

public:
    int InitMainSolver();
    int InitReader();
    int InitSolver();
    int EpochSol();

public:
    ReadRnxNav reader_nav_;
    ReadPreEph reader_pre_;
    ReadBias   reader_bia_;
    ReadAnt    reader_ant_;
    ReadBlq    reader_blq_;
    ReadErp    reader_erp_;
    ReadRnxObs *reader_rov_;
    ReadRnxObs *reader_bas_;

public:
    AntModel ant_corr_;
    BiaModel bia_corr_;
    TidModel tid_corr_;
    Nav nav_;

public:
    int site_count_;
    Time ts_,te_;
    double tint_;
    Solver *solver_;
    OutSol *out_sol_;
};

#endif //MSNAVS_SOLVER_H
