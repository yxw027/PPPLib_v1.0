/**************************************************************************

Copyright(c): 2020-, by Chao Chen, All rights reserved.
              China University of Mining and Technology, XuZhou, P.R. China

Author: Chao Chen(cchen@cumt.edu.cn)

Date: 2020-05-30

Description:

**************************************************************************/
#include "PPPLibGlo.h"

// Core options class
PosOpt::PosOpt() {
    rover_=Vec3(0.0);
    base_=Vec3(0.0);
}

PosOpt::~PosOpt() {

}

void PosOpt::SetPosOpt() {
    Config::Ptr_ config=Config::GetInstance();
    vector<double>pos;
    pos=config->GetArray<double>("rover_pos");
    rover_.i_=pos[0];rover_.j_=pos[1];rover_.k_=pos[2];
    pos.clear();
    pos=config->GetArray<double>("base_pos");
    base_.i_=pos[0];rover_.j_=pos[1];rover_.k_=pos[2];
}

AROpt::AROpt() {
    mode_=AR_OFF;
    ar_gps_= false;
    ar_bds_= false;
    ar_glo_=AR_OFF;
    min_lock2fix_=5;
    min_sat2fix_=5;
    min_sat2drop_=5;
    min_ele2fix_=10.0;
    max_iter_=1;
    min_fix2hold_=5;
    min_sat2hold_=5;
    min_ele2hold_=10.0;
    pseudo_var2hold_=1.0;
    for(double & i : thres_ar_) i=0.0;
}

AROpt::~AROpt() {

}

void AROpt::SetAROpt() {
    Config::Ptr_ config=Config::GetInstance();
}

PPKOpt::PPKOpt() {
    base_obs_interp_= false;
}

PPKOpt::~PPKOpt() {

}

void PPKOpt::SetPPKOpt() {
    Config::Ptr_ config=Config::GetInstance();

}

PPPOpt::PPPOpt() {
    phw_= true;
    clk_jump_repair_= false;
    eclips_sat_test= false;
    bd2_multipath_corr= false;
}

PPPOpt::~PPPOpt() {

}

void PPPOpt::SetPPPOpt() {
    Config::Ptr_ config=Config::GetInstance();
    phw_=config->Get<int>("phw_corr");
    clk_jump_repair_=config->Get<int>("clk_jump_repair");
    eclips_sat_test=config->Get<int>("eclips_sat_test");
    bd2_multipath_corr=config->Get<int>("bd2_multipath_corr");
}

AdjOpt::AdjOpt() {
    mode_=ADJ_LSQ;
    gnss_std_[0]=5.0;gnss_std_[1]=0.03;gnss_std_[2]=0.3;
    gnss_psd_[0]=10e-11; gnss_psd_[1]=10e-3;gnss_psd_[2]=10e-8;gnss_psd_[3]=10e-6;gnss_psd_[4]=0.0;
    glo_icb_modeling_= GLOICB_LNF;
    gnss_isb_modeling_=STO_RW;
    bds_isb_= false;
}

AdjOpt::~AdjOpt() {

}

void AdjOpt::SetAdjOpt() {
    Config::Ptr_ config=Config::GetInstance();

    mode_=config->Get<int>("adj_mode");
    vector<double>std;
    std=config->GetArray<double>("gnss_prn");
    for(int i=0;i<std.size();i++) gnss_psd_[i]=std[i];
    std.clear();
    std=config->GetArray<double>("gnss_std");
    for(int i=0;i<std.size();i++) gnss_std_[i]=std[i];
    gnss_isb_modeling_=config->Get<int>("gnss_isb_modeling");
    glo_icb_modeling_=config->Get<int>("glo_icb_modeling");
    bds_isb_=config->Get<int>("bd2_bd3_isb");
}

SlipOpt::SlipOpt() {
    method_=SLIP_MW;
    thres_mw_=5.0;
    thres_gf_=0.15;
}

SlipOpt::~SlipOpt() {

}

void SlipOpt::SetSlipOpt() {
    Config::Ptr_ config=Config::GetInstance();

//    method_=config->Get<int>("cs_method");
//    thres_mw_=config->Get<double>("cs_mw_thres");
//    thres_gf_=config->Get<double>("cs_gf_thres");
}

QCOpt::QCOpt() {
    raim_= true;
    max_rej_inno_=500.0;
    max_rej_dop_=30.0;
    max_outage2reset_amb_=5;
    adj_meas_var= false;
    for(auto &i:exclude_sat_) i=USED;
    exclude_sat_[3]=EX_MAN;
//    exclude_sat_[49]=EX_MAN;
//    exclude_sat_[87]=EX_MAN;
}

QCOpt::~QCOpt() {

}

void QCOpt::SetQCOpt() {
    Config::Ptr_ config=Config::GetInstance();
    max_rej_inno_=config->Get<double>("max_rej_inno");
    max_rej_dop_=config->Get<double>("max_rej_dop");
    raim_=config->Get<int>("spp_raim");
    vector<string>exc_sat;
    exc_sat=config->GetArray<string>("exclude_sat");

    slip_opt_.SetSlipOpt();
}

GNSSOpt::GNSSOpt() {
    obs_sample_rate_=0.0;nav_sys_=GPS;
    prods_ac_=BRD;
    use_frq_=SINGLE_FRQ;
    use_GPS_frq_=(GPS_L5<<8)|(GPS_L2<<4)|(GPS_L1);
    use_BD2_frq_=(BD2_B2I<<8)|(BD2_B3I<<4)|(BD2_B1I);
    use_BD3_frq_=(BD3_B2a<<8)|(BD3_B3I<<4)|(BD3_B1I);
    use_GAL_frq_=(GAL_E5b<<8)|(GAL_E5a<<4)|(GAL_E1);
    use_GLO_frq_=(GLO_G2<<4)|(GLO_G1);
    use_QZS_frq_=(QZS_L5<<8)|(QZS_L2<<4)|(QZS_L1);
    for(double & i : p_l_err_ratio_) i=100.0;
    meas_err_factor_[0]=100.0;meas_err_factor_[1]=0.003;meas_err_factor_[2]=0.003;
    meas_err_factor_[3]=0.0;meas_err_factor_[4]=0.0;
    ele_min_=10.0;
    weight_mode_=ELEVATION;
    code_smooth_= false;
    sat_eph_=BRDC;
    cb_prd_type_=CBC_DCB;
    cb_prd_ac_=BIA_CAS;
    ion_opt_=ION_BRDC;
    IF3_format_=IF3_DUAL;
    trp_opt_=TRP_SAAS;
    trp_map_opt_=TRPMAP_EL;
    tid_opt_=TIDE_OFF;
    rec_pcv_= false;
    sat_pcv_= false;
    dynamic_modeling_= false;
}

GNSSOpt::~GNSSOpt() {

}

void GNSSOpt::SetGNSSOpt() {
    Config::Ptr_ config=Config::GetInstance();
    obs_sample_rate_=config->Get<double>("obs_sample");
//    nav_sys_=config->Get<int>("nav_sys");
    prods_ac_=config->Get<int>("prod_ac");
    use_frq_=config->Get<int>("num_use_frqs");
    use_GPS_frq_=config->Get<int>("use_GPS_frq");
    use_BD2_frq_=config->Get<int>("use_BD2_frq");
    use_BD3_frq_=config->Get<int>("use_BD3_frq");
    use_GAL_frq_=config->Get<int>("use_GAL_frq");
    use_GLO_frq_=config->Get<int>("use_GLO_frq");
    use_QZS_frq_=config->Get<int>("use_QZS_frq");
    vector<double>ratio;
    ratio=config->GetArray<double>("p_h_err_ratio");
    for(int i=0;i<ratio.size();i++) p_l_err_ratio_[i]=ratio[i];
    ratio.clear();
    ratio=config->GetArray<double>("meas_err_ratio");
    for(int i=0;i<ratio.size();i++) meas_err_factor_[i]=ratio[i];
    ele_min_=config->Get<double>("ele_min");
    weight_mode_=config->Get<int>("weight_mode");
    code_smooth_=config->Get<int>("code_smooth");
    sat_eph_=config->Get<int>("sat_eph");
    cb_prd_type_=config->Get<int>("cbias_prd_type");
    cb_prd_ac_=config->Get<int>("cbias_prd_ac");
    ion_opt_=config->Get<int>("ion_opt");
    IF3_format_=config->Get<int>("if3_opt");
    trp_opt_=config->Get<int>("trp_opt");
    trp_map_opt_=config->Get<int>("trp_map_opt");
    tid_opt_=config->Get<int>("tid_opt");
    rec_pcv_=config->Get<int>("rec_pcv");
    sat_pcv_=config->Get<int>("sat_pcv");
    dynamic_modeling_=config->Get<int>("dyc_modeling");
}

PrcOpt::PrcOpt() {
    dir_=DEFPRC_DIR;
    tu_=0.0;
    ti_=30.0;
    use_def_opt_= true;
    debug_level_=128;
    mode_=PPP;mode_opt_=MDOPT_STATIC;  // default SPP-STATIC processing settings
    filter_type_=FORWARD;
};

PrcOpt::~PrcOpt() {

}

void PrcOpt:: SetPrcOpt() {
    Config::Ptr_ config=Config::GetInstance();
    dir_=config->Get<string>("prc_data_dir");
    vector<int>date;
    date=config->GetArray<int>("prc_date");
    if(!date.empty()){
        prc_date_.epoch_[0]=date[0];prc_date_.epoch_[1]=date[1];prc_date_.epoch_[2]=date[2];
        prc_date_.Epoch2Time(prc_date_.epoch_);
        ts_.epoch_[0]=te_.epoch_[0]=date[0];ts_.epoch_[1]=te_.epoch_[1]=date[1];ts_.epoch_[2]=te_.epoch_[2]=date[2];
    }
    date.clear();
    date=config->GetArray<int>("prc_ts");
    if(!date.empty()){
        ts_.epoch_[3]=date[0];ts_.epoch_[4]=date[1];ts_.epoch_[5]=date[2];
        ts_.Epoch2Time(ts_.epoch_);
    }
    date.clear();
    date=config->GetArray<int>("prc_te");
    if(!date.empty()){
        te_.epoch_[3]=date[0];te_.epoch_[4]=date[1];te_.epoch_[5]=date[2];
        te_.Epoch2Time(te_.epoch_);
    }
    tu_=config->Get<double>("prc_unit");
//    mode_=config->Get<int>("prc_mode");
//    mode_opt_=config->Get<int>("prc_mode_opt");
    filter_type_=config->Get<int>("filter_type");

    GNSS_opt_.SetGNSSOpt();
    GNSS_opt_.ppp_opt_.SetPPPOpt();
    adj_opt_.SetAdjOpt();
    QC_opt_.SetQCOpt();
}

void PrcOpt::SetDefaultPPPOpt() {
    GNSS_opt_.prods_ac_=WUM;
//    GNSS_opt_.use_frq_=SINGLE_FRQ;
    GNSS_opt_.use_GPS_frq_=(GPS_L5<<8)|(GPS_L1<<4)|(GPS_L1);
    GNSS_opt_.use_BD2_frq_=(BD2_B3I<<8)|(BD2_B2I<<4)|(BD2_B1I);
    GNSS_opt_.use_BD3_frq_=(BD3_B2a<<8)|(BD3_B2a<<4)|(BD3_B1I);
    GNSS_opt_.use_GAL_frq_=(GAL_E5b<<8)|(GAL_E5a<<4)|(GAL_E1);
    GNSS_opt_.use_GLO_frq_=(GLO_G2<<4)|(GLO_G1);
    GNSS_opt_.use_QZS_frq_=(QZS_L5<<8)|(QZS_L2<<4)|(QZS_L1);
    GNSS_opt_.sat_eph_=PRECISE;
//    GNSS_opt_.ion_opt_=ION_UC;
    GNSS_opt_.IF3_format_=IF3_DUAL;
    GNSS_opt_.trp_opt_=TRP_EST_WET;
    GNSS_opt_.trp_map_opt_=TRPMAP_GMF;
    GNSS_opt_.rec_pcv_= true;
    GNSS_opt_.sat_pcv_= true;
    GNSS_opt_.tid_opt_=TIDE_SOLID|TIDE_OCEAN|TIDE_POLE;
    adj_opt_.mode_=ADJ_KF;
    GNSS_opt_.ppp_opt_.phw_=true;
}

void PrcOpt::SetDefaultPPKOpt() {

}

FileOpt::FileOpt() {
    customize_file_= false;
}

FileOpt::~FileOpt() {

}

void FileOpt::SetFileOpt() {
    Config::Ptr_ config=Config::GetInstance();
    customize_file_=config->Get<int>("customize_file");
    if(customize_file_!=0)
        customize_dir_=config->Get<string>("customize_dir");
}

SolOpt::SolOpt() {
    err_fmt_=false;
    out_head_=true;
    pos_fmt_=POSFMT_ENU;
    out_sat_= true;
    out_bug_= false;
    ref_sol_.pos=Vec3(0.0);
    ref_sol_.vel=Vec3(0.0);
    ref_sol_.rpy=Vec3(0.0);
}

SolOpt::~SolOpt() {

}

// global
el::Logger* defaultLogger;
PrcOpt    kPrcOpt;
FileOpt   kFileOpt;
SolOpt    kSolOpt;
string    kSite="";
