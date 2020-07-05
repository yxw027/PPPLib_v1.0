//
// Created by cc on 4/9/20.
//

#ifndef MSNAVS_MSNAVSGLO_H
#define MSNAVS_MSNAVSGLO_H

#include "CmnFunc.h"
#include "GNSS.h"
#include "Config.h"

// Core options class
enum kAdjModOpt{ADJ_LSQ,ADJ_KF,ADJ_HEL};
enum kFilterType{FORWARD,BACKWARD,COMBINATION};
enum kGNSSSys{NONE=0x00,GPS=0x01,BD2=0x02,BD3=0x04,GAL=0x08,GLO=0x10,QZS=0x20,IRN=0x40,SBS=0x80,ALL=0xFF};
enum kPosMode{SPP,PPP,DGPS,PPK,FIX,IGLC,IGTC};
enum kPosModeOpt{MDOPT_STATIC,MDOPT_KINE_SIM,MDOPT_KINE,MDOPT_SPP,MDOPT_PPP,MDOPT_DGPS,MDOPT_PPK};
enum kSlipOpt{SLIP_MW=0x01,SLIP_GF=0x01};
enum kProdsAC{BRD,COD,WUM,GFZ,GBM,GRG};
enum kUseFrqNum{SINGLE_FRQ=1,DUAL_FRQ=2,TRIPLE_FRQ=3,QUAD_FRQ=4};
enum kWeightOpt{ELEVATION,SNR};
enum kSatEph{BRDC,PRECISE};
enum kCodeBias{CBC_OFF,CBC_TGD,CBC_DCB,CBC_OSB};
enum kBiaAC{BIA_CODE,BIA_CAS,BIA_CNES};
enum kBiasType{CODE_BIAS,PHASE_BIAS};
enum kIonOpt{ION_OFF,ION_BRDC,ION_TEC,ION_IF,ION_UC,ION_CONST};
enum kIF3Opt{IF3_SINGLE,IF3_DUAL};
enum kTrpOpt{TRP_OFF,TRP_SAAS,TRP_EST_WET,TRP_EST_GRAD};
enum kMapOpt{TRPMAP_EL,TRPMAP_GMF,TRPMAP_VMF1,TRPMAP_NEIL};
enum kTidOpt{TIDE_OFF=0x00,TIDE_SOLID=0x01,TIDE_OCEAN=0x02,TIDE_POLE=0x04};
enum kAROpt{AR_OFF,AR_CONTINUOUS,AR_INSTANTANEOUS,AR_FIX_HOLD};
enum kGloAROpt{GLOAR_OFF,GLOAR_ON,GLOAR_AUTOCAL,GLOAR_FIX_HOLD};
enum kStoOpt{STO_OFF,STO_CT,STO_PWC,STO_RW,STO_WN};
enum kDymOpt{DYM_OFF,DYM_VEL,DYM_ACC};
enum kSolStat{SOL_NONE,SOL_SPP,SOL_P3SPP,SOL_PPP,SOL_PPP_FIX,SOL_DGPS,SOL_PPK,SOL_PPK_FIX};
enum kGloICB{GLOICB_OFF,GLOICB_LNF,GLOICB_QUAD};
enum kPosFmt{POSFMT_XYZ,POSFMT_ENU};
enum kSatExOpt{NO_USE=0,USED=1,EX_MAN=2,EX_NOPOBS=3,EX_NOLOBS=4,EX_SVH=5,EX_URA=6,EX_NOPRED=7,EX_EL=8,EX_PR_PRI_RES=9,EX_PR_NOR_RES=10,EX_CP_PRI_RES=11,
             EX_CP_NOR_RES=12,EX_PR_POST_RES=13,EX_CP_POST_RES=14,EX_SLIP=15};
const string kListACsStr[]={
        "","cod","wum","gfz","gbm","grg"
};

const string kListFrqsStr[]={
        "","SF","DF","TF","QF"
};

const string kListModesStr[]={
        "SPP","PPP","DGPS","PPK","FIXED"
};
const string kListModOptsStr[]={
        "STATIC","KINE","KINE-SIM","SPP","PPP","DGPS","PPK"
};
const string kListSolStatStr[]={
        "NONE","SPP","P3SPP","PPP","PPP-FIX","DGPS","PPK","PPK-FIX"
};
const string kListTropStr[]={
        "OFF","SAAS","EST_WET","EST_GRAD"
};
const string kListIonsStr[]={
        "OFF","BRD","TEC","IF","UC","CONST"
};
const string kListEphStr[]={
        "BRDC","PRECISE"
};
const string kListPosFmtStr[]={
        "XYZ","ENU"
};

const string kListSatEXOptStr[]={
       "","USED","MAN-MADE","NOPOBS","NOLOBS","SVH","URA","NO_PROD","LOW_EL","PRI_PRES","P_NRES","PRI_LRES","L_NRES",
       "POS_PRES","POS_LRES","SLIP"
};

class PosOpt {
public:
    PosOpt();
    ~PosOpt();

public:
    void SetPosOpt();

public:
    Vec3 rover_;
    Vec3 base_;
};

class AROpt {
public:
    AROpt();
    ~AROpt();

public:
    void SetAROpt();

public:
    int mode_;
    bool ar_gps_;
    bool ar_bds_;
    int ar_glo_;
    bool filter_;
    int min_lock2fix_;
    int min_sat2fix_;
    int min_sat2drop_;
    double min_ele2fix_;
    int max_iter_;
    int min_fix2hold_;
    int min_sat2hold_;
    double min_ele2hold_;
    double pseudo_var2hold_;
    double thres_ar_[6];
};

class PPKOpt {
public:
    PPKOpt();
    ~PPKOpt();

public:
    void SetPPKOpt();

public:
    bool base_obs_interp_;
};

class PPPOpt {
public:
    PPPOpt();
    ~PPPOpt();

public:
    void SetPPPOpt();

public:
    bool phw_;
    bool clk_jump_repair_;
    bool eclips_sat_test;
    bool bd2_multipath_corr;
};

class GNSSOpt{
public:
    GNSSOpt();
    ~GNSSOpt();

public:
    void SetGNSSOpt();

public:
    double obs_sample_rate_;
    int nav_sys_;
    int prods_ac_;
    int use_frq_;
    int use_GPS_frq_;
    int use_BD2_frq_;
    int use_BD3_frq_;
    int use_GAL_frq_;
    int use_GLO_frq_;
    int use_QZS_frq_;
    double p_l_err_ratio_[3];
    double meas_err_factor_[5];
    double ele_min_;
    int weight_mode_;
    bool code_smooth_;
    int sat_eph_;
    int cb_prd_type_;   // code bias correction (DCB or OSB)
    int cb_prd_ac_;
    int ion_opt_;
    int IF3_format_;
    int trp_opt_;
    int trp_map_opt_;
    int tid_opt_;
    bool rec_pcv_;
    bool sat_pcv_;
    bool dynamic_modeling_;
    PosOpt pos_opt_;
    PPPOpt ppp_opt_;
    AROpt ar_opt_;
};

class AdjOpt {
public:
    AdjOpt();
    ~AdjOpt();

public:
    void SetAdjOpt();

public:
    int mode_;
    double gnss_std_[3];     // bias, ion, trp
    double
    gnss_psd_[5];     // amb, ion, trp, ifb, dcb
    int gnss_isb_modeling_;
    int glo_icb_modeling_;
    int bds_isb_;
};

class SlipOpt {
public:
    SlipOpt();
    ~SlipOpt();

public:
    void SetSlipOpt();
    int method_;
    double thres_mw_;
    double thres_gf_;
};

class QCOpt {
public:
    QCOpt();
    ~QCOpt();

public:
    void SetQCOpt();
    bool raim_;
    double max_rej_inno_;
    double max_rej_dop_;
    int max_outage2reset_amb_;
    unsigned char exclude_sat_[MAXSAT];
    bool adj_meas_var;
    SlipOpt slip_opt_;
};


class PrcOpt {
public:
    PrcOpt();
    ~PrcOpt();

public:
    void SetPrcOpt();
    void SetDefaultPPPOpt();
    void SetDefaultPPKOpt();

public:
    string dir_;
    Time prc_date_;
    Time ts_,te_;
    double tu_;
    double ti_;

    bool use_def_opt_;
    int debug_level_;
    int mode_;
    int mode_opt_;
    int filter_type_;
    GNSSOpt GNSS_opt_;
    AdjOpt  adj_opt_;
    QCOpt  QC_opt_;
};

class FileOpt {
public:
    FileOpt();
    ~FileOpt();

public:
    void SetFileOpt();

public:
    int customize_file_;
    string customize_dir_;
    string imu_;
    string rover_;
    string base_;
    string brdc_;
    string clk_;
    string sp3_[3];
    string bia_;
    string erp_;
    string atx_;
    string blq_;
    string ion_;
    string sol_;
    string ref_;
};

typedef struct {
    Vec3 pos;
    Vec3 vel;
    Vec3 rpy;
}RefSol_t;

class SolOpt {
public:
    SolOpt();
    ~SolOpt();

public:
    void SetSolOpt();

public:
    bool err_fmt_;
    bool out_head_;
    int  pos_fmt_;
    RefSol_t ref_sol_;
    bool out_sat_;
    bool out_bug_;
};

//// global class declaration
extern el::Logger* defaultLogger;
extern PrcOpt   kPrcOpt;
extern FileOpt  kFileOpt;
extern SolOpt   kSolOpt;
extern string   kSite;

#endif //MSNAVS_MSNAVSGLO_H
