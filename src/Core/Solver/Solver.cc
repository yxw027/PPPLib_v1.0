//
// Created by cc on 4/9/20.
//

#include "Solver.h"
#include "EphModel.h"
#include <algorithm>

// 初始化一些参数
Solver::Solver() {
    Isat_t *isat= nullptr;
    epoch_count_=site_count_=dt_=vaild_sat_num_=0;
    for(auto &i:exist_sys_mask_) i=0;
    for(int i=0;i<MAXSAT;i++){
        isat=&isat_[i];
        for(int j=0;j<NFREQ;j++){
            isat->corr_P[j]=isat->corr_L[j]=isat->IF_P[j]=isat->IF_L[j]=0.0;
            isat->slip[j]=isat->rejc[j]=isat->lock[j]=isat->rejc[j]=0;

            isat->amb[j]=isat->cbias[j]=isat->pbias[j]=isat->dants[j]=isat->dantr[j]=0.0;
            for(int k=0;k<2;k++) isat->post_res[k][j]=isat->pri_res[k][j]=0.0;
        }
        isat->gf[0]=isat->gf[1]=isat->mw[0]=isat->mw[1]=isat->smw[0]=isat->smw[1]=isat->mw_idx[0]=isat->mw_idx[1]=isat->mw_var[0]=isat->mw_var[1]=0.0;
        isat->strp=isat->strp_w=isat->strp_h=0.0;
        isat->map_trp_w=isat->map_trp_h=isat->map_grad_e=isat->map_grad_n=0.0;
        isat->var_trp=isat->var_ion=0.0;
        isat->dist=isat->tide=isat->shapiro=isat->sagnac=0.0;
        isat->sion=isat->map_ion=0.0;
        isat->stat=USED;isat->clk_err[0]=isat->clk_err[1]=0.0;
        isat->last_ep=1;
    }
}

Solver::~Solver() {
    sig_sat_data_= nullptr;
    A_coe_.clear();OMCs_.clear();R_meas_.clear();R_vec_.clear();
}

void Solver::GetUseFrq(PrcOpt popt, ObsData &obs, double *lam) {
    int i,sys=obs.sat_.sat_sys_,idx_sys=obs.sat_.sat_sys_idx_;
    if(obs.sat_.sat_sys_==BD2){
        if(obs.sat_.sat_prn_>=19) {
            idx_sys=INDEXBD3;
            sys=BD3;
        }
    }

    for(i=0;i<NFREQ;i++) {
        obs.lam_[i]=lam[i];
        obs.frq_[i]=kGNSSFreqs[idx_sys][i];
    }

    int np=11;
    if(sys==GPS){
        LOG_N_TIMES(1,DEBUG)<<"GPS frequency: "<<"L1  "<<setprecision(np)<<obs.frq_[0]<<" L2  "<<obs.frq_[1]\
                               <<" L5  "<<obs.frq_[2];
    }
    else if(sys==BD2){
        LOG_N_TIMES(1,DEBUG)<<"BD2 frequency: "<<"B1I "<<setprecision(np)<<obs.frq_[0]<<" B2I "<<obs.frq_[1]\
                               <<" B3I "<<obs.frq_[2];
    }
    else if(sys==BD3){
        LOG_N_TIMES(1,DEBUG)<<"BD3 frequency: "<<"B1I "<<setprecision(np)<<obs.frq_[0]<<" B2a "<<obs.frq_[1]\
                               <<" B3I "<<obs.frq_[2];
    }
    else if(sys==GAL){
        LOG_N_TIMES(1,DEBUG)<<"GAL frequency: "<<"E1  "<<setprecision(np)<<obs.frq_[0]<<" E5a "<<obs.frq_[1]\
                               <<" E5b "<<obs.frq_[2];
    }
    else if(sys==GLO){
        LOG_N_TIMES(1,DEBUG)<<"GLO frequency: "<<"G1  "<<setprecision(np)<<obs.frq_[0]<<" G2  "<<obs.frq_[1]\
                               <<" G3  "<<obs.frq_[2];
    }
    else if(sys==QZS){
        LOG_N_TIMES(1,DEBUG)<<"QZS frequency: "<<"L1  "<<setprecision(np)<<obs.frq_[0]<<" L2  "<<obs.frq_[1]\
                               <<" L5  "<<obs.frq_[2];
    }

    GNSSOpt gopt=popt.GNSS_opt_;
    int f=0;
    switch(sys){
        case GPS: f=gopt.use_GPS_frq_;break;
        case BD2: f=gopt.use_BD2_frq_;break;
        case BD3: f=gopt.use_BD3_frq_;break;
        case GAL: f=gopt.use_GAL_frq_;break;
        case GLO: f=gopt.use_GLO_frq_;break;
        case QZS: f=gopt.use_QZS_frq_;break;
    }
    if(gopt.use_frq_==SINGLE_FRQ){
        // 单频用户不考虑使用双频进行数据处理
        obs.use_f1_=f&0x0F;
        if(obs.P_[obs.use_f1_]==0.0||obs.L_[obs.use_f1_]==0.0) obs.stat_=EX_NOPOBS;
    }
    else if(gopt.use_frq_==DUAL_FRQ){
        obs.use_f1_=f&0x0F;
        obs.use_f2_=(f>>4)&0x0F;
         if(obs.P_[obs.use_f1_]==0.0||obs.P_[obs.use_f2_]==0.0||obs.L_[obs.use_f1_]==0.0||obs.L_[obs.use_f2_]==0.0)
            obs.stat_=EX_NOPOBS;
    }
    else if(gopt.use_frq_==TRIPLE_FRQ){
        obs.use_f3_=f>>8&0xFF;obs.use_f2_=f>>4&0x0F;obs.use_f1_=f&0x0F;
        if((obs.P_[obs.use_f1_]==0.0||obs.P_[obs.use_f2_]==0.0)||(obs.L_[obs.use_f1_]==0.0||obs.L_[obs.use_f2_]==0.0))
            obs.stat_=EX_NOPOBS;
    }
    else if(gopt.use_frq_==QUAD_FRQ){
        //only BD3 supported quad frequency
        obs.use_f4_<<16&0xFF;obs.use_f3_=f>>8&0xFF;obs.use_f2_=f>>4&0x0F;obs.use_f1_=f&0x0F;
    }
}

void Solver::ReSetIsat(int sat) {
    Isat_t *isat;
    isat=&isat_[sat-1];
    for(int j=0;j<NFREQ;j++){
        isat->corr_P[j]=isat->corr_L[j]=isat->IF_P[j]=isat->IF_L[j]=0.0;
        isat->slip[j]=isat->rejc[j]=isat->lock[j]=isat->rejc[j]=0;

        isat->amb[j]=isat->cbias[j]=isat->pbias[j]=isat->dants[j]=isat->dantr[j]=0.0;
        for(int k=0;k<2;k++) isat->post_res[k][j]=isat->pri_res[k][j]=0.0;
    }
    isat->strp=isat->strp_w=isat->strp_h=0.0;
    isat->map_trp_w=isat->map_trp_h=isat->map_grad_e=isat->map_grad_n=0.0;
    isat->var_trp=isat->var_ion=0.0;
    isat->dist=isat->tide=isat->shapiro=isat->sagnac=0.0;
    isat->sion=isat->map_ion=0.0;
    isat->clk_err[0]=isat->clk_err[1]=0.0;
}

void Solver::UpdateSatInfo(PrcOpt popt) {
    int n=epoch_rover_.num_;
    ObsData *data;
    double *lam;
    for(int i=0;i<n;i++){
        data=&epoch_rover_.sat_infos_[i];
        data->epoch_=epoch_count_;
        qc_.ExcludeSat(epoch_rover_.sat_infos_[i]);
        if(!popt.adj_opt_.bds_isb_&&data->sat_.sat_sys_==BD3){
            data->sat_.sat_sys_=BD2;data->sat_.sat_sys_idx_=INDEXBD2;
        }
        if(data->stat_>USED&&data->stat_!=EX_SLIP) data->ReSetSat();
        lam=nav_.lam_[data->sat_.sat_no_-1];
        GetUseFrq(popt,*data,lam);
        if(popt.GNSS_opt_.cb_prd_type_>CBC_OFF){
            bia_corr_.sat_data_=data;
            bia_corr_.MatchOsb2Sat();
        }
    }
}

void Solver::SetSysMask() {
    for(int i=0;i<NSYS;i++) exist_sys_mask_[i]=0;
    for(int i=0;i<epoch_rover_.num_;i++){
        int sys=epoch_rover_.sat_infos_[i].sat_.sat_sys_;
        if(epoch_rover_.sat_infos_[i].stat_>USED) continue;
        switch(sys){
            case GPS: exist_sys_mask_[INDEXGPS]=1;break;
            case BD2: exist_sys_mask_[INDEXBD2]=1;break;
            case BD3: exist_sys_mask_[INDEXBD3]=1;break;
            case GAL: exist_sys_mask_[INDEXGAL]=1;break;
            case GLO: exist_sys_mask_[INDEXGLO]=1;break;
            case QZS: exist_sys_mask_[INDEXQZS]=1;break;
        }
    }
}

double Solver::Dops() {
    double dops[4]={0};
    ObsData *data= nullptr;
    vector<double> H;
    H.assign(4*epoch_rover_.num_,0.0);
    int n=0;

    for(int i=0;i<epoch_rover_.num_;i++){
        data=&epoch_rover_.sat_infos_[i];
        if(data->stat_>USED&&data->stat_!=EX_SLIP) continue;
        if(data->azel_[1]*R2D<kPrcOpt.GNSS_opt_.ele_min_||data->azel_[1]<=0.0) continue;
        double cosel=cos(data->azel_[1]);
        double sinel=sin(data->azel_[1]);
        H[4*n]=cosel*sin(data->azel_[0]);
        H[1+4*n]=cosel*cos(data->azel_[0]);
        H[2+4*n]=sinel;
        H[3+4*n++]=1.0;
    }

    if(n<4) return 0.0;
    vector<double>Q;
    Q.assign(4*4,0.0);
    MatMulVec("NT",4,4,n,1.0,H,H,0.0,Q);
    if(!MatInv(Q,4)){
        dops[0]=SQRT(Q[0]+Q[5]+Q[10]+Q[15]);  // GDOP
        dops[1]=SQRT(Q[0]+Q[5]+Q[10]);        // PDOP
        dops[2]=SQRT(Q[0]+Q[5]);              // HDOP
        dops[3]=SQRT(Q[10]);                  // VDOP
    }
    Q.clear();H.clear();
    return dops[1];
}

void Solver::UpdateSol(vector<double> x, vector<double> Px,double zhtrp,double zwtrp) {
    int npos=para_.NumPosPar();
    epoch_sol_.i_epoch_=epoch_count_;
    epoch_sol_.num_used_sat_=vaild_sat_num_;

    if(epoch_sol_.stat_){
        // sol time
        if(epoch_sol_.stat_!=SOL_P3SPP) epoch_sol_.sol_time_=epoch_rover_.sat_infos_[0].sig_recep_;
        // position
        epoch_sol_.pos_.i_=x[0];
        epoch_sol_.pos_.j_=x[1];
        epoch_sol_.pos_.k_=x[2];
        // clock
        epoch_sol_.clk_G_=x[npos+INDEXGPS];
        epoch_sol_.B2_isb_=x[npos+INDEXBD2];
        epoch_sol_.B3_isb_=x[npos+INDEXBD3];
        epoch_sol_.E_isb_=x[npos+INDEXGAL];
        epoch_sol_.R_isb_=x[npos+INDEXGLO];
        epoch_sol_.J_isb_=x[npos+INDEXQZS];
        // wet zenith trp
        epoch_sol_.ztrp_dry_=zhtrp;
        epoch_sol_.ztrp_wet_=zwtrp;
    }
}

int Solver::Fullx2Zipx(vector<double> full_x,vector<double> full_Px, vector<double>&zip_x,vector<double>&zip_Px,
        int full_nx) {
    int i=0,k=0;
    vector<int>ix;vector<double>H;
    H.assign(A_coe_.begin(),A_coe_.end());

#if 0
    // pos
    for(i=0;i<para_.NumPosPar();i++) ix.push_back(i);
    // clk
    for(i=para_.NumPosPar();i<para_.NumSPPPar();i++){
        if(exist_sys_mask_[i-para_.NumPosPar()]) ix.push_back(i);
    }

    // dcb/ifb/trp/ion/amb
    if(full_nx>para_.NumSPPPar()){
        for(i=para_.NumSPPPar();i<full_nx;i++){
            if(full_x[i]!=0.0&&full_Px[i+i*full_nx]>0.0){
                ix.push_back(i);
            }
        }
    }
#else
    for(i=0;i<para_.NumPosPar();i++) ix.push_back(i);
    for(i=para_.NumPosPar();i<para_.NumSPPPar();i++){
        if(exist_sys_mask_[i-para_.NumPosPar()]) ix.push_back(i);
    }
    if(para_.NumDcbPar()){
        for(i=para_.IdxDcbPar();i<para_.IdxDcbPar()+para_.NumDcbPar();i++){
            if(exist_sys_mask_[i-para_.IdxDcbPar()]) ix.push_back(i);
        }
    }
    if(para_.NumIfbPar()){
        for(i=para_.IdxIfbPar(INDEXGPS);i<para_.IdxIfbPar(INDEXGPS)+para_.NumIfbPar();i++){
            if(exist_sys_mask_[i-para_.IdxIfbPar(INDEXGPS)]) ix.push_back(i);
        }
    }
    if(para_.NumGloPar()){
        for(i=para_.IdxGloPar();i<para_.IdxGloPar()+para_.NumGloPar();i++){
            if(full_x[i]!=0.0&&full_Px[i+i*full_nx]>0.0) ix.push_back(i);
        }
    }
    if(para_.NumTrpPar()){
        for(i=para_.IdxTrpPar();i<para_.IdxTrpPar()+para_.NumTrpPar();i++){
            if(full_x[i]!=0.0&&full_Px[i+i*full_nx]>0.0) ix.push_back(i);
        }
    }

    if(para_.NumIonPar()){
        for(k=0;k<epoch_rover_.num_;k++){
            if(epoch_rover_.sat_infos_[k].stat_>USED&&epoch_rover_.sat_infos_[k].stat_!=EX_SLIP) continue;
            int sat=epoch_rover_.sat_infos_[k].sat_.sat_no_;
            int iI=para_.IdxIonPar(sat);
            if(full_x[iI]!=0.0&&full_Px[iI+iI*full_nx]>0.0) ix.push_back(iI);
        }
    }
    if(para_.NumAmbPar()){
        for(k=0;k<para_.NumUseFrq();k++){
            for(i=0;i<epoch_rover_.num_;i++){
                if(epoch_rover_.sat_infos_[i].stat_>USED&&epoch_rover_.sat_infos_[i].stat_!=EX_SLIP) continue;
                int sat=epoch_rover_.sat_infos_[i].sat_.sat_no_;
                int iA=para_.IdxAmbPar(sat,k);
                if(full_x[iA]!=0.0&&full_Px[iA+iA*full_nx]>0.0) ix.push_back(iA);
            }
        }
    }


#endif

    int zip_nx=ix.size();
    zip_x.assign(zip_nx,0.0);
    zip_Px.assign(zip_nx*zip_nx,0.0);
    A_coe_.assign(num_L_*zip_nx,0.0);
    for(i=0;i<zip_nx;i++){
        zip_x[i]=full_x[ix[i]];
        for(int j=0;j<zip_nx;j++) zip_Px[i+j*zip_nx]=full_Px[ix[i]+ix[j]*full_nx];
        for(int j=0;j<num_L_;j++)  A_coe_[i+j*zip_nx]=H[ix[i]+j*full_nx];
    }
    ix.clear();H.clear();
    return zip_nx;
}

void Solver::Zipx2Fullx(vector<double>zip_x, vector<double>zip_Px,vector<double>&full_x,vector<double>&full_Px,int full_nx) {
    int i=0,k=0;vector<int>ix;
#if 0
    // pos
    for(i=0;i<para_.NumPosPar();i++) ix.push_back(i);
    // clk
    for(i=para_.NumPosPar();i<para_.NumSPPPar();i++){
        if(exist_sys_mask_[i-para_.NumPosPar()]) ix.push_back(i);
    }

    // dcb/ifb/trp/ion/amb
    if(full_nx>para_.NumSPPPar()){
        for(i=para_.NumSPPPar();i<full_nx;i++){
            if(full_x[i]!=0.0&&full_Px[i+i*full_nx]>0.0){
                ix.push_back(i);
            }
        }
    }
#else
    for(i=0;i<para_.NumPosPar();i++) ix.push_back(i);
    for(i=para_.NumPosPar();i<para_.NumSPPPar();i++){
        if(exist_sys_mask_[i-para_.NumPosPar()]) ix.push_back(i);
    }
    if(para_.NumDcbPar()){
        for(i=para_.IdxDcbPar();i<para_.IdxDcbPar()+para_.NumDcbPar();i++){
            if(exist_sys_mask_[i-para_.IdxDcbPar()]) ix.push_back(i);
        }
    }
    if(para_.NumIfbPar()){
        for(i=para_.IdxIfbPar(INDEXGPS);i<para_.IdxIfbPar(INDEXGPS)+para_.NumIfbPar();i++){
            if(exist_sys_mask_[i-para_.IdxIfbPar(INDEXGPS)]) ix.push_back(i);
        }
    }
    if(para_.NumGloPar()){
        for(i=para_.IdxGloPar();i<para_.IdxGloPar()+para_.NumGloPar();i++){
            if(full_x[i]!=0.0&&full_Px[i+i*full_nx]>0.0) ix.push_back(i);
        }
    }
    if(para_.NumTrpPar()){
        for(i=para_.IdxTrpPar();i<para_.IdxTrpPar()+para_.NumTrpPar();i++){
            if(full_x[i]!=0.0&&full_Px[i+i*full_nx]>0.0) ix.push_back(i);
        }
    }

    if(para_.NumIonPar()){
        for(k=0;k<epoch_rover_.num_;k++){
            if(epoch_rover_.sat_infos_[k].stat_>USED&&epoch_rover_.sat_infos_[k].stat_!=EX_SLIP) continue;
            int sat=epoch_rover_.sat_infos_[k].sat_.sat_no_;
            int iI=para_.IdxIonPar(sat);
            if(full_x[iI]!=0.0&&full_Px[iI+iI*full_nx]>0.0) ix.push_back(iI);
        }
    }
    if(para_.NumAmbPar()){
        for(k=0;k<para_.NumUseFrq();k++){
            for(i=0;i<epoch_rover_.num_;i++){
                if(epoch_rover_.sat_infos_[i].stat_>USED&&epoch_rover_.sat_infos_[i].stat_!=EX_SLIP) continue;
                int sat=epoch_rover_.sat_infos_[i].sat_.sat_no_;
                int iA=para_.IdxAmbPar(sat,k);
                if(full_x[iA]!=0.0&&full_Px[iA+iA*full_nx]>0.0) ix.push_back(iA);
            }
        }
    }
#endif

    for(i=0;i<ix.size();i++){
        full_x[ix[i]]=zip_x[i];
        for(int j=0;j<ix.size();j++) full_Px[ix[i]+ix[j]*full_nx]=zip_Px[i+j*ix.size()];
    }
    ix.clear();
}

#define EFACT_GPS     1.0
#define EFACT_BD2     1.0
#define EFACT_BD2_GEO 5.0
#define EFACT_BD3     1.0
#define EFACT_GAL     1.0
#define EFACT_GLO     1.0
#define EFACT_QZS     1.0
double Solver::MeasVar(PrcOpt popt,int type) {
    GNSSOpt gopt=popt.GNSS_opt_;
    double a,b;
    double fact=type==OBS_CP?1.0:gopt.p_l_err_ratio_[0];
    double sinel=sin(sig_sat_data_->azel_[1]);
    int sys=sig_sat_data_->sat_.sat_sys_,prn=sig_sat_data_->sat_.sat_prn_;

    switch(sys){
        case GPS: fact*=EFACT_GPS;break;
        case BD3: fact*=EFACT_BD3;break;
        case GAL: fact*=EFACT_GAL;break;
        case GLO: fact*=EFACT_GLO;break;
        case QZS: fact*=EFACT_QZS;break;
    }

    if(sys==BD2){
        if(sig_sat_data_->sat_.sat_prn_>=19) fact*=EFACT_BD3;
        else if(std::binary_search(kBDS2_GEO,kBDS2_GEO+NUM_BDS2_GEO,prn)){
            fact*=EFACT_BD2_GEO;
        }else fact*=EFACT_BD2;
    }

    a=gopt.meas_err_factor_[1];
    b=gopt.meas_err_factor_[2];

    if(popt.GNSS_opt_.use_frq_==SINGLE_FRQ&&popt.GNSS_opt_.ion_opt_==ION_IF){
        fact*=0.5*100;
    } else fact*=(gopt.ion_opt_==ION_IF)?3.0:1.0;

    switch(gopt.weight_mode_){
        case ELEVATION: return SQR(fact)*(SQR(a)+SQR(b/sinel));
    }
}

void Solver::NormalEqLine(PrcOpt popt,int use_f,int idx_f,int nl,int nx,int type) {
    int idx=0;

    if(popt.GNSS_opt_.use_frq_==TRIPLE_FRQ){
        if(popt.GNSS_opt_.ion_opt_==ION_IF) idx=1;
        else idx=2;
    }

    para_.popt_=popt;
    para_.num_L_=nl;
    para_.nx_=nx;
    for(int i=0;i<NSYS;i++) para_.exist_sys_mark_[i]=exist_sys_mask_[i];
    para_.sat_info_=*sig_sat_data_;
    para_.CoePosPar(A_coe_);
    para_.CoeClkPar(A_coe_);
    para_.CoeDcbPar(A_coe_);
    if(!type&&idx_f==idx) para_.CoeIfbPar(A_coe_);
    if((!type||(popt.GNSS_opt_.use_frq_==SINGLE_FRQ&&popt.GNSS_opt_.ion_opt_==ION_IF))&&sig_sat_data_->sat_.sat_sys_==GLO)
        para_.CoeGloPar(A_coe_);
    para_.CoeTrpPar(A_coe_);
    para_.CoeIonPar(A_coe_,use_f,type);
    if(type) para_.CoeAmbPar(idx_f,A_coe_);
}

double Solver::GetCorrObs(PrcOpt popt,int f,int type) {
    //type　0:伪距　1:相位
    double Pc=0.0;
    int sat=sig_sat_data_->sat_.sat_no_;
    int flag=0;
    GNSSOpt gopt=popt.GNSS_opt_;
    if(gopt.ion_opt_==ION_IF){
        if(gopt.use_frq_==SINGLE_FRQ){
            if(isat_[sat-1].corr_L[f]!=0.0&&isat_[sat-1].corr_P[f]!=0.0)
                Pc=0.5*isat_[sat-1].corr_L[f]+0.5*isat_[sat-1].corr_P[f];
        }
        if(gopt.use_frq_==DUAL_FRQ){
            // L1+L2 0+1-1=0
            // L1+L3 0+2-1=1
            // L2+L3 1+2-1=2
            // L1+L2+L3 0+1+2=3
            flag=sig_sat_data_->use_f1_+sig_sat_data_->use_f2_-1;
            if(type) Pc=isat_[sat-1].IF_L[flag];
            else     Pc=isat_[sat-1].IF_P[flag];
        }
        if(gopt.use_frq_==TRIPLE_FRQ){

            if(gopt.IF3_format_==IF3_DUAL){
                // f=0 L1_L2
                // f=1 L1_L3
                // BDS 存的是B1_B3 B1_B2
                if(type) Pc=isat_[sat-1].IF_L[f];
                else     Pc=isat_[sat-1].IF_P[f];
            }
            else if(gopt.IF3_format_==IF3_SINGLE){
                // f=0 L1_L2
                // GPS 并不是所有卫星都支持三频，对于非三频卫星或者三频数据缺失应使用L1_L2
                // f=1 L1_L2_L3
                if(f==0&&isat_[sat-1].IF_P[3]!=0.0&&isat_[sat-1].IF_L[3]!=0.0){
                    if(type) Pc=isat_[sat-1].IF_L[3];
                    else     Pc=isat_[sat-1].IF_P[3];
                }
                else{
                    if(type) Pc=isat_[sat-1].IF_L[0];
                    else     Pc=isat_[sat-1].IF_P[0];
                }
            }
        }
    }
    else{
        if(type) Pc=isat_[sat-1].corr_L[f];
        else Pc=isat_[sat-1].corr_P[f];
    }

    return Pc;
}

void Solver::CorrPRMeas(PrcOpt popt,ObsData &obs,double dantr[NFREQ],double dants[NFREQ]) {
    int sat=obs.sat_.sat_no_;
    double cbias=0.0;
    int bd3_flag=0;
    if(obs.sat_.sat_sys_==BD2&&obs.sat_.sat_prn_>=19) bd3_flag=1;

    for(int i=0;i<NFREQ;i++){
        isat_[sat-1].corr_P[i]=isat_[sat-1].corr_L[i]=0.0;
        if(obs.P_[i]==0.0) continue;
        cbias=isat_[sat-1].cbias[i]=bia_corr_.CBiaCorr(popt,i,obs);
        isat_[sat-1].corr_P[i]=obs.P_[i]-dants[i]-dantr[i]-cbias;
    }

    // L1_L2 L1_L3 L2_L3 L1_L2_L3
    double alpha=0.0,beta=0.0;
    double f1=0.0,f2=0.0,f3=0.0;
    double P1=0.0,P2=0.0,P3=0.0;
    double Pc=0.0;
    switch(obs.sat_.sat_sys_){
        case GPS:
            f1=obs.frq_[GPS_L1];f2=obs.frq_[GPS_L2];f3=obs.frq_[GPS_L5];
            P1=isat_[sat-1].corr_P[GPS_L1];P2=isat_[sat-1].corr_P[GPS_L2];P3=isat_[sat-1].corr_P[GPS_L5];
            break;
        case BD2:
            if(bd3_flag){
                f1=obs.frq_[BD3_B1I];f2=obs.frq_[BD3_B2a];f3=obs.frq_[BD3_B3I];
                P1=isat_[sat-1].corr_P[BD3_B1I];P2=isat_[sat-1].corr_P[BD3_B2a];P3=isat_[sat-1].corr_P[BD3_B3I];
                break;
            }
            else{
                f1=obs.frq_[BD2_B1I];f2=obs.frq_[BD2_B2I];f3=obs.frq_[BD2_B3I];
                P1=isat_[sat-1].corr_P[BD2_B1I];P2=isat_[sat-1].corr_P[BD2_B2I];P3=isat_[sat-1].corr_P[BD2_B3I];
                break;
            }
        case BD3:
            f1=obs.frq_[BD3_B1I];f2=obs.frq_[BD3_B2a];f3=obs.frq_[BD3_B3I];
            P1=isat_[sat-1].corr_P[BD3_B1I];P2=isat_[sat-1].corr_P[BD3_B2a];P3=isat_[sat-1].corr_P[BD3_B3I];
            break;
        case GAL:
            f1=obs.frq_[GAL_E1];f2=obs.frq_[GAL_E5a];f3=obs.frq_[GAL_E5b];
            P1=isat_[sat-1].corr_P[GAL_E1];P2=isat_[sat-1].corr_P[GAL_E5a];P3=isat_[sat-1].corr_P[GAL_E5b];
            break;
        case GLO:
            f1=obs.frq_[GLO_G1];f2=obs.frq_[GLO_G2];
            P1=isat_[sat-1].corr_P[GLO_G1];P2=isat_[sat-1].corr_P[GLO_G2];
            break;
        case QZS:
            f1=obs.frq_[QZS_L1];f2=obs.frq_[QZS_L2];f3=obs.frq_[QZS_L5];
            P1=isat_[sat-1].corr_P[QZS_L1];P2=isat_[sat-1].corr_P[QZS_L2];P3=isat_[sat-1].corr_P[QZS_L5];
            break;
    }
    if(P1!=0.0&&P2!=0.0&&P3!=0.0) obs.frq3_=1;
    else obs.frq3_=0;

    alpha=SQR(f1)/(SQR(f1)-SQR(f2));beta=-SQR(f2)/(SQR(f1)-SQR(f2));
    if(P1!=0.0&&P2!=0.0) Pc=alpha*P1+beta*P2;
    isat_[obs.sat_.sat_no_-1].IF_P[0]=Pc;
    isat_[obs.sat_.sat_no_-1].IF_L[0]=0.0;

    //1-3 0+2 -1= 1
    Pc=0.0;
    alpha=SQR(f1)/(SQR(f1)-SQR(f3));beta=-SQR(f3)/(SQR(f1)-SQR(f3));
    if(P1!=0.0&&P3!=0.0) Pc=alpha*P1+beta*P3;
    isat_[obs.sat_.sat_no_-1].IF_P[1]=Pc;
    isat_[obs.sat_.sat_no_-1].IF_L[1]=0.0;

    //2-3 1+2-1 = 2
    Pc=0.0;
    alpha=SQR(f2)/(SQR(f2)-SQR(f3));beta=-SQR(f3)/(SQR(f2)-SQR(f3));
    if(P2!=0.0&&P3!=0.0) Pc=alpha*P2+beta*P3;
    isat_[obs.sat_.sat_no_-1].IF_P[2]=Pc;
    isat_[obs.sat_.sat_no_-1].IF_L[2]=0.0;

    //1-2-3 0+1+2 = 3
    Pc=0.0;
    if(obs.sat_.sat_sys_==GLO) Pc=0.0;
    else{
        double gam1=SQR(f1)/SQR(f1),gam2=SQR(f1)/SQR(f2),gam3=SQR(f1)/SQR(f3);
        double e=2*(SQR(gam2)+SQR(gam3)-gam2*gam3-gam2-gam3+1.0);
        double e1=(SQR(gam2)+SQR(gam3)-gam2-gam3)/e;
        double e2=(SQR(gam3)-gam2*gam3-gam2+1.0)/e;
        double e3=(SQR(gam2)-gam2*gam3-gam3+1.0)/e;
        if(P1!=0.0&&P2!=0.0&&P3!=0.0) Pc=e1*P1+e2*P2+e3*P3;
        isat_[obs.sat_.sat_no_-1].IF_P[3]=Pc;
        isat_[obs.sat_.sat_no_-1].IF_L[3]=0.0;
    }

    // obs scan
    if(popt.GNSS_opt_.use_frq_==SINGLE_FRQ){
        if(obs.P_[obs.use_f1_]==0.0) obs.stat_=EX_NOPOBS;
    }
    else if(popt.GNSS_opt_.use_frq_==DUAL_FRQ){
        if(popt.GNSS_opt_.ion_opt_==ION_IF){

        }
        else{
            if(obs.P_[obs.use_f1_]==0.0||obs.P_[obs.use_f2_]==0.0) obs.stat_=EX_NOPOBS;
        }
    }
    else if(popt.GNSS_opt_.ion_opt_==TRIPLE_FRQ){
        if(popt.GNSS_opt_.ion_opt_==ION_IF){
            if(isat_[sat-1].IF_P[0]==0.0&&isat_[sat-1].IF_P[1]==0.0) obs.stat_=EX_NOPOBS;
        }
        else{
            if(obs.P_[obs.use_f1_]==0.0&&obs.P_[obs.use_f2_]==0.0&&obs.P_[obs.use_f3_]==0.0) obs.stat_=EX_NOPOBS;
        }
    }
}

void Solver::CorrCPMeas(PrcOpt popt,ObsData &obs,double dantr[NFREQ],double dants[NFREQ],double phw) {
    int sat=obs.sat_.sat_no_;
    int bd3_flag=0;

    if(obs.sat_.sat_sys_==BD2&&obs.sat_.sat_prn_>=19) bd3_flag=1;

    for(int i=0;i<NFREQ;i++){
        isat_[sat-1].corr_L[i]=0.0;
        if(obs.L_[i]==0.0) continue;
        isat_[sat-1].corr_L[i]=obs.L_[i]*obs.lam_[i]-dants[i]-dantr[i]-phw*obs.lam_[i];
    }

    // L1_L2 L1_L3 L2_L3 L1_L2_L3
    double alpha=0.0,beta=0.0;
    double f1=0.0,f2=0.0,f3=0.0;
    double l1=0.0,l2=0.0,l3=0.0;
    double L1=0.0,L2=0.0,L3=0.0;
    double Lc=0.0;
    switch(obs.sat_.sat_sys_){
        case GPS:
            f1=obs.frq_[GPS_L1];f2=obs.frq_[GPS_L2];f3=obs.frq_[GPS_L5];
            L1=isat_[sat-1].corr_L[GPS_L1];L2=isat_[sat-1].corr_L[GPS_L2];L3=isat_[sat-1].corr_L[GPS_L5];
            break;
        case BD2:
            if(bd3_flag){
                f1=obs.frq_[BD2_B1I];f2=obs.frq_[BD2_B2I];f3=obs.frq_[BD2_B3I];
                L1=isat_[sat-1].corr_L[BD2_B1I];L2=isat_[sat-1].corr_L[BD2_B2I];L3=isat_[sat-1].corr_L[BD2_B3I];
                break;
            }
            else{
                f1=obs.frq_[BD3_B1I];f2=obs.frq_[BD3_B2a];f3=obs.frq_[BD3_B3I];
                L1=isat_[sat-1].corr_L[BD3_B1I];L2=isat_[sat-1].corr_L[BD3_B2a];L3=isat_[sat-1].corr_L[BD3_B3I];
                break;
            }
        case BD3:
            f1=obs.frq_[BD3_B1I];f2=obs.frq_[BD3_B2a];f3=obs.frq_[BD3_B3I];
            L1=isat_[sat-1].corr_L[BD3_B1I];L2=isat_[sat-1].corr_L[BD3_B2a];L3=isat_[sat-1].corr_L[BD3_B3I];
            break;
        case GAL:
            f1=obs.frq_[GAL_E1];f2=obs.frq_[GAL_E5a];f3=obs.frq_[GAL_E5b];
            L1=isat_[sat-1].corr_L[GAL_E1];L2=isat_[sat-1].corr_L[GAL_E5a];L3=isat_[sat-1].corr_L[GAL_E5b];
            break;
        case GLO:
            f1=obs.frq_[GLO_G1];f2=obs.frq_[GLO_G2];
            L1=isat_[sat-1].corr_L[GLO_G1];L2=isat_[sat-1].corr_L[GLO_G2];
            break;
        case QZS:
            f1=obs.frq_[QZS_L1];f2=obs.frq_[QZS_L2];f3=obs.frq_[QZS_L5];
            L1=isat_[sat-1].corr_L[QZS_L1];L2=isat_[sat-1].corr_L[QZS_L2];L3=isat_[sat-1].corr_L[QZS_L5];
            break;
    }

    if(L1!=0.0&&L2!=0.0&&L3!=0.0) obs.frq3_=1;
    else obs.frq3_=0;

    alpha=SQR(f1)/(SQR(f1)-SQR(f2));beta=-SQR(f2)/(SQR(f1)-SQR(f2));
    if(L1!=0.0&&L2!=0.0) Lc=alpha*L1+beta*L2;
    isat_[obs.sat_.sat_no_-1].IF_L[0]=Lc;

    //1-3 0+2 -1= 1
    Lc=0.0;
    alpha=SQR(f1)/(SQR(f1)-SQR(f3));beta=-SQR(f3)/(SQR(f1)-SQR(f3));
    if(L1!=0.0&&L3!=0.0) Lc=alpha*L1+beta*L3;
    isat_[obs.sat_.sat_no_-1].IF_L[1]=Lc;

    //2-3 1+2-1 = 2
    Lc=0.0;
    alpha=SQR(f2)/(SQR(f2)-SQR(f3));beta=-SQR(f3)/(SQR(f2)-SQR(f3));
    if(L2!=0.0&&L3!=0.0) Lc=alpha*L2+beta*L3;
    isat_[obs.sat_.sat_no_-1].IF_L[2]=Lc;

    //1-2-3 0+1+2 = 3
    Lc=0.0;
    if(obs.sat_.sat_sys_==GLO) Lc=0.0;
    else{
        double gam1=SQR(f1)/SQR(f1),gam2=SQR(f1)/SQR(f2),gam3=SQR(f1)/SQR(f3);
        double e=2*(SQR(gam2)+SQR(gam3)-gam2*gam3-gam2-gam3+1.0);
        double e1=(SQR(gam2)+SQR(gam3)-gam2-gam3)/e;
        double e2=(SQR(gam3)-gam2*gam3-gam2+1.0)/e;
        double e3=(SQR(gam2)-gam2*gam3-gam3+1.0)/e;
        if(L1!=0.0&&L2!=0.0&&L3!=0.0) Lc=e1*L1+e2*L2+e3*L3;
        isat_[obs.sat_.sat_no_-1].IF_L[3]=Lc;
    }

    if(popt.GNSS_opt_.use_frq_==SINGLE_FRQ){
        if(obs.L_[obs.use_f1_]==0.0) obs.stat_=EX_NOLOBS;
    }
    else if(popt.GNSS_opt_.use_frq_==DUAL_FRQ){
        if(popt.GNSS_opt_.ion_opt_==ION_IF){

        }
        else{
            if(obs.L_[obs.use_f1_]==0.0||obs.L_[obs.use_f2_]==0.0) obs.stat_=EX_NOLOBS;
        }
    }
    else if(popt.GNSS_opt_.ion_opt_==TRIPLE_FRQ){
        if(popt.GNSS_opt_.ion_opt_==ION_IF){
            if(isat_[sat-1].IF_L[0]==0.0&&isat_[sat-1].IF_L[1]==0.0) obs.stat_=EX_NOLOBS;
        }
        else{
            if(obs.L_[obs.use_f1_]==0.0&&obs.L_[obs.use_f2_]==0.0) obs.stat_=EX_NOLOBS;
        }
    }
}

double Solver::RecClkCorr(vector<double> est_x) {
    int sys=sig_sat_data_->sat_.sat_sys_idx_;
    int nspp=para_.NumSPPPar();
    vector<int>idx_sys;
    vector<double>x;
    x.assign(nspp,0.0);
    int npos=para_.NumPosPar();

    for(int i=0,k=0;i<nspp;i++){
        if(est_x[i]!=0.0) x[k++]=est_x[i];
    }

    for(int i=0;i<NSYS;i++){
        if(exist_sys_mask_[i]!=0){
            idx_sys.push_back(i);
        }
    }
    double clock=x[npos];
    double isb=0;
    for(int i=0;i<idx_sys.size();i++){
        if(sys==idx_sys[i]){
            if(i==0) isb=0.0;
            else isb=x[i+npos];
        }
    }

    return clock+isb;
}

double Solver::RecIfbCorr(vector<double> est_x) {
    int sys=sig_sat_data_->sat_.sat_sys_idx_;
    vector<int>idx_sys;

    int ifb=para_.IdxIfbPar(INDEXGPS);

    for(int i=0;i<NSYS;i++){
        if(exist_sys_mask_[i]!=0){
            idx_sys.push_back(i);
        }
    }

    for(int i=0;i<idx_sys.size();i++){
        if(sys==idx_sys[i]){
            return est_x[ifb+idx_sys[i]];
        }
    }
}

double Solver::GloIcbCorr(PrcOpt popt,vector<double> est_x,int type) {
    int sys=sig_sat_data_->sat_.sat_sys_;

    if(sys!=GLO||type!=OBS_PR) return 0.0;

    int idx=para_.IdxGloPar();
    int frq=nav_.glo_frq_num_[sig_sat_data_->sat_.sat_prn_-1];
    double icb=0.0;
    icb=est_x[idx]*frq;
    if(popt.adj_opt_.glo_icb_modeling_==GLOICB_QUAD) icb+=frq*frq*est_x[idx+1];
    return icb;
}

bool CmpObs(const ObsData& p1,const ObsData& p2){
    return p1.sat_.sat_no_<p2.sat_.sat_no_;
}

int Solver::SortObs(Obss &obs) {
    if (obs.num_<=0) return 0;
    sort(obs.sat_infos_.begin(),obs.sat_infos_.end(),CmpObs);
    return 1;
}

int Solver::Algorithm() {
    return 1;
}

SPPSolver::SPPSolver() {
}

SPPSolver::~SPPSolver() {
    spp_ion_= nullptr;spp_trp_= nullptr;spp_adj_= nullptr;
    full_x_.clear();full_Px_.clear();zip_x_.clear();zip_Px_.clear();

}

#define MAXSPP_ITER   15
int SPPSolver::SolverSPP(PrcOpt popt) {
    int stat=0,iter=0,qc_flag=0;
    char buff[MAXBUFF]={'\0'};


    do{
        iter++;
        SetSysMask();

        para_.popt_=popt;
        MakeSPPEqs(0,iter,popt);
        zip_nx_=Fullx2Zipx(full_x_,full_Px_,zip_x_,zip_Px_,full_nx_);
        if(num_L_<zip_nx_||vaild_sat_num_<=4){
            LOG(WARNING)<<epoch_count_<<"th epoch insufficient number of necessary observations";
            return 0;
        }

        R_meas_.assign(num_L_*num_L_,0.0);
        for(int j=0;j<num_L_;j++) R_meas_[j*num_L_+j]=R_vec_[j];

        if((stat=spp_adj_->Adjustment(OMCs_,A_coe_,num_L_,zip_nx_,R_meas_,zip_x_,zip_Px_))==-1){
            LOG(ERROR)<<epoch_count_<<"th epoch least square estimator error";
            return 0;
        }
        Zipx2Fullx(zip_x_,zip_Px_,full_x_,full_Px_,full_nx_);

        if(iter>MAXSPP_ITER){
            sprintf(buff,"%4d %20s:  LSQ failed to converge",epoch_count_,epoch_rover_.sat_infos_[0].sig_recep_.time_str_.c_str());
            LOG(WARNING)<<buff;
            return 0;
        }

        if(stat==1) qc_flag=SPPPostRes(popt,*spp_adj_);

    }while(stat<1||qc_flag);

    epoch_sol_.sigma0_=spp_adj_->sigma0_;
    // 验证spp解是否可靠

    return 1;
}

int SPPSolver::SPPAlgorithm(PrcOpt popt) {
    epoch_sol_.stat_=SOL_NONE;
    if(epoch_rover_.num_<0){
        LOG(WARNING)<<epoch_count_<<"th no observation";
        return 0;
    }

    UpdateSatInfo(popt_);

    if(brd_eph_->BrdcPosClk(epoch_rover_)<=4){
        LOG(WARNING)<<epoch_count_<<"th no enough valid satellite";
        return 0;
    }

    if(SolverSPP(popt)){
        epoch_sol_.stat_=popt_.mode_==SPP?SOL_SPP:SOL_P3SPP;
        epoch_sol_.pdop_=Dops();
        UpdateSol(full_x_,full_Px_,spp_trp_->z_hyd_,spp_trp_->z_wet_);
    }
    else{
        epoch_sol_.stat_=SOL_NONE;
        LOG(WARNING)<<epoch_count_<<"th spp solver fail";
    }


    return epoch_sol_.stat_;
}

int SPPSolver::SPPPostRes(PrcOpt popt, AdjModel adj_model) {
    int qc_PR=0;
    // 验后残差
    MakeSPPEqs(1,0,popt);
    // 伪距验后残差检验
    qc_.epoch_=epoch_count_;
    qc_.PRv_info_=L_info_;
    qc_PR=qc_.PostResQC(popt,epoch_rover_,adj_model.v_PR_,adj_model.norm_PRv_,num_PR_,OBS_PR);
    return qc_PR;
}

int SPPSolver::PRErrEq(PrcOpt popt,int use_f,int idx_f,int iter,int post) {
    int sat=sig_sat_data_->sat_.sat_no_;
    double mod_pr=0.0,corr_pr=0.0,omc=0.0,q_meas=0.0;
    if((corr_pr=GetCorrObs(popt,use_f,OBS_PR))==0.0){
        LOG(WARNING)<<epoch_count_<<" "<<sig_sat_data_->sig_recep_.time_str_<<": "<<sig_sat_data_->sat_.sat_id_<<" corr_pr=0.0";
        return 0;
    }
    if((mod_pr=ModPRMeas(post))==0.0){
        LOG(DEBUG)<<epoch_count_<<" "<<sig_sat_data_->sig_recep_.time_str_<<": "<<sig_sat_data_->sat_.sat_id_<<" mod_pr=0.0";
        return 0;
    }
    omc=corr_pr-mod_pr;
    if(epoch_count_>1&&fabs(omc)>500000) {
        LOG(WARNING)<<epoch_count_<<" "<<sig_sat_data_->sig_recep_.time_str_<<": "<<sig_sat_data_->sat_.sat_id_<<" omc abnormal";
        return 0;
    }
    OMCs_.push_back(omc);
    if(!post&&iter==1) isat_[sat-1].pri_res[0][use_f]=omc;
    else if(post) isat_[sat-1].post_res[0][use_f]=omc;
    q_meas=MeasVar(popt,OBS_PR)+sig_sat_data_->sat_var_+isat_[sat-1].var_trp+isat_[sat-1].var_ion;
    R_vec_.push_back(q_meas);
    NormalEqLine(popt,use_f,idx_f,num_PR_,full_nx_,OBS_PR);
    char sat_info[MAXBUFF]={'\0'};
    sprintf(sat_info,"%5d %s P%d omc=%10.4f var=%7.3f el=%3.1f obs=%10.3f r=%10.3f dtr=%10.3f dts=%10.3f trp=%6.3f ion=%6.3f cbias=%6.3f shapiro=%4.3f",
            epoch_count_,sig_sat_data_->sat_.sat_id_.c_str(),use_f+1,omc,q_meas,sig_sat_data_->azel_[1]*R2D,corr_pr,
            isat_[sat-1].dist+isat_[sat-1].sagnac,RecClkCorr(full_x_),sig_sat_data_->clk_err_[0]*CLIGHT,isat_[sat-1].strp,
            isat_[sat-1].sion,popt.GNSS_opt_.ion_opt_==ION_IF?0.0:isat_[sat-1].cbias[use_f],isat_[sat-1].shapiro);
    LOG(DEBUG)<<sat_info;
    num_PR_++;
    return num_PR_;
}

void SPPSolver::MakeSPPEqs(int post,int iter,PrcOpt popt) {
    int n=epoch_rover_.num_,nf=popt.GNSS_opt_.use_frq_,frq=0;
    int info_flag=0,use_flag=0;
    double dants[NFREQ]={0},dantr[NFREQ]={0};
    num_PR_=vaild_sat_num_=0;
    A_coe_.clear();OMCs_.clear();R_meas_.clear();R_vec_.clear();L_info_.clear();
    double rr[3];
    for(int i=0;i<3;i++) rr[i]=full_x_[i];

    for(int i=0;i<n;i++){
        sig_sat_data_=&epoch_rover_.sat_infos_[i];
        sig_sat_data_=&epoch_rover_.sat_infos_[i];
        CorrPRMeas(popt,*sig_sat_data_,dantr,dants);
        if(sig_sat_data_->stat_>USED) continue;
        for(int f=0;f<nf;f++){
            switch(f){
                case 0: frq=sig_sat_data_->use_f1_;break;
                case 1: frq=sig_sat_data_->use_f2_;break;
                case 2: frq=sig_sat_data_->use_f3_;break;
                case 3: frq=sig_sat_data_->use_f4_;break;
            }
            char buff[MAXBUFF]={'\0'};
            sprintf(buff,"prior position %10.3f %10.3f %10.3f",rr[0],rr[1],rr[2]);
            LOG_IF(i==0&&f==0,DEBUG)<<"SPP iter("<<iter<<") "<<buff;
            if(!(use_flag=PRErrEq(popt,frq,f,iter,post))) continue;
            info_flag=(sig_sat_data_->sat_.sat_no_<<8)|(OBS_PR<<4)|(f);
            L_info_.push_back(info_flag);
        }
        if(use_flag){
            sig_sat_data_->stat_=USED;
            vaild_sat_num_++;
        }
    }
    num_L_=num_PR_;
}

double SPPSolver::ModPRMeas(int post) {
    int sat=sig_sat_data_->sat_.sat_no_;
    double mod_pr=0.0;
    double rr[3]={0},tide_dr[3]={0};
    for(int i=0;i<3;i++) rr[i]=zip_x_[i];

    Vec3 rec_pos(rr[0],rr[1],rr[2]);
    Vec3 blh=Xyz2Blh(rec_pos,WGS84);
    isat_[sat-1].dist=GeoDist(sig_sat_data_->pos_,rec_pos,sig_sat_data_->sig_vec_);
    isat_[sat-1].sagnac=sig_sat_data_->SagnacCorr(rec_pos);
    double el=SatAzEl(blh,sig_sat_data_->sig_vec_,sig_sat_data_->azel_);
    if(el<popt_.GNSS_opt_.ele_min_*D2R||isat_[sat-1].dist<=0.0){
        if(post) sig_sat_data_->stat_=EX_EL;
        return 0.0;
    }

    spp_trp_->isat_=isat_;
    spp_trp_->sat_data_=sig_sat_data_;
    spp_trp_->TrpCorr(blh,0.7);
    spp_ion_->isat_=isat_;
    spp_ion_->sat_data_=sig_sat_data_;
    spp_ion_->IonCorr(blh);

    mod_pr=isat_[sat-1].dist+isat_[sat-1].sagnac+RecClkCorr(zip_x_)-sig_sat_data_->clk_err_[0]*CLIGHT
            +isat_[sat-1].strp+isat_[sat-1].sion;
    return mod_pr;
}

int SPPSolver::InitSPPSolver(PrcOpt popt) {
    popt_=popt;
    brd_eph_=new BrdcEph;
    brd_eph_->InitSatEph(nav_);

    spp_trp_=new SaasTrp;
    spp_ion_=new KlobIon;
    spp_adj_=new LSQAdj;
    spp_ion_->GetBrdIonPara(nav_.ion_para_[INDEXGPS].ion);

    full_nx_=para_.NumSPPPar();
    full_x_.assign(full_nx_,0.0);
    full_Px_.assign(full_nx_*full_nx_,0.0);
    zip_x_.assign(full_x_.begin(),full_x_.end());
    zip_Px_.assign(full_Px_.begin(),full_Px_.end());
}

int SPPSolver::Algorithm() {
    if(epoch_count_==1) InitSPPSolver(kPrcOpt);
    return SPPAlgorithm(popt_);
}

PPPSolver::PPPSolver() {
    for(auto &i:rr_) i=0.0;

}

PPPSolver::~PPPSolver() {
    full_x_.clear();full_Px_.clear();
    zip_x_.clear();zip_Px_.clear();
    lar_res_.clear();lar_res_sat_.clear();
    lar_res_frq_.clear();lar_res_type_.clear();

}

int PPPSolver::P3SPP() {
    PrcOpt popt=kPrcOpt;
    popt.mode_=SPP;
    popt.GNSS_opt_.nav_sys_=popt_.GNSS_opt_.nav_sys_;
    popt.GNSS_opt_.ion_opt_=ION_BRDC;
    popt.GNSS_opt_.trp_opt_=TRP_SAAS;
    popt.GNSS_opt_.sat_eph_=BRDC;
    popt.GNSS_opt_.use_frq_=SINGLE_FRQ;
    popt.GNSS_opt_.use_GPS_frq_=GPS_L1;
    popt.GNSS_opt_.use_BD2_frq_=BD2_B1I;
    popt.GNSS_opt_.use_BD3_frq_=BD3_B1I;
    popt.GNSS_opt_.use_GAL_frq_=GAL_E1;
    popt.GNSS_opt_.use_GLO_frq_=GLO_G1;
    popt.GNSS_opt_.use_QZS_frq_=QZS_L1;
    popt.GNSS_opt_.sat_pcv_=false;
    popt.GNSS_opt_.rec_pcv_=false;
    popt.GNSS_opt_.cb_prd_type_=BIA_CAS;
    if(epoch_count_==1) InitSPPSolver(popt);
    else{
        full_x_[0]=epoch_sol_.pos_.i_;
        full_x_[1]=epoch_sol_.pos_.j_;
        full_x_[2]=epoch_sol_.pos_.k_;
    }
    return SPPAlgorithm(popt);
}

#define MAXPPP_ITER  8
int PPPSolver::SolverPPP() {
    int stat=0,iter=0;
    vector<double>x,Px;
    x.assign(full_nx_,0.0);Px.assign(full_nx_*full_nx_,0.0);
    epoch_sol_.stat_=SOL_SPP;

    do{
        x.assign(full_x_.begin(),full_x_.end());
        Px.assign(full_Px_.begin(),full_Px_.end());

        iter++;
        SetSysMask();

        if(!MakePPPEqs(0,iter,x)) return 0;
        if(vaild_sat_num_<=4) return 0;

        R_meas_.assign(num_L_*num_L_,0.0);
        for(int j=0;j<num_L_;j++) R_meas_[j*num_L_+j]=R_vec_[j];

        para_.popt_=popt_;
        zip_nx_=Fullx2Zipx(x,Px,zip_x_,zip_Px_,full_nx_);

#if 0
        Mat X(zip_nx_,1,zip_x_);
        X.TraceMat(10,3);

        Mat A(num_L_,zip_nx_,A_coe_);
        A.TraceMat(6,3);
        Mat Pp(zip_nx_,zip_nx_,zip_Px_);
        Pp.TraceMat(10,6);

        Mat R(num_L_,num_L_,R_meas_);
        R.TraceMat(6,3);
#endif

        ppp_adj_->num_PR_=num_PR_;ppp_adj_->num_CP_=num_CP_;ppp_adj_->v_info_=L_info_;
        ppp_adj_->epoch_=epoch_count_;
        stat=ppp_adj_->Adjustment(OMCs_,A_coe_,num_L_,zip_nx_,R_meas_,zip_x_,zip_Px_);
        Zipx2Fullx(zip_x_,zip_Px_,x,Px,full_nx_);

        if(iter>=MAXPPP_ITER){
            LOG(WARNING)<<epoch_count_<<"th epoch ppp filter no coverage";
            break;
        }

        if(!PPPPostRes(*ppp_adj_,x)){
            full_x_.assign(x.begin(),x.end());
            full_Px_.assign(Px.begin(),Px.end());
            break;
        }

    }while(stat<1);

    x.clear();Px.clear();
    epoch_sol_.sigma0_=ppp_adj_->sigma0_;
    return 1;
}

void PPPSolver::PPPCycSlip() {
    qc_.epoch_=epoch_count_;
    qc_.isat_=isat_;
    qc_.MWCycSlip(popt_,epoch_rover_,dt_,bia_corr_.dcb_);
    qc_.GFCycSlip(popt_,epoch_rover_,dt_);

    Isat_t *isat= nullptr;
    ObsData *data= nullptr;
    for(int i=0;i<MAXSAT;i++){
        isat=&isat_[i];
        isat->mw[0]=isat->mw[1]=isat->gf[0]=isat->gf[1]=0.0;
    }
    int sat=0;

    for(int i=0;i<epoch_rover_.num_;i++){
        data=&epoch_rover_.sat_infos_[i];
        sat=data->sat_.sat_no_;
        isat=&isat_[sat-1];
        isat->last_time=epoch_rover_.sat_infos_[0].sig_recep_;
        double gf=0.0,w1=0.0;

#if 0
        if((gf=qc_.GFMeas(*isat,*data,0))!=0.0) isat->gf=gf;
        if((w1=qc_.MWMeas(*isat,*data,bia_corr_.dcb_,0))==0.0) continue;

        double w0=isat->smw;
        isat->mw=w1;
        if(isat->mw_idx>0){
            int j=isat->mw_idx;
            double var0=isat->mw_var;
            double var1=(w1-w0)*(w1-w0)-var0;
            var1=var0+var1/j;

            isat->smw=(w0*j+w1)/(j+1);
            isat->mw_idx++;
            isat->mw_var=var1;
        }
        else{
            isat->smw=w1;
            isat->mw_idx++;
            isat->mw_var=0.25;
        }
#else
        int f=data->frq3_?2:1;
        for(int k=0;k<f;k++){
            if((gf=qc_.GFMeas(*isat,*data,k))!=0.0) isat->gf[k]=gf;
            if((w1=qc_.MWMeas(*isat,*data,bia_corr_.dcb_,k))==0.0) continue;

            double w0=isat->smw[k];
            isat->mw[k]=w1;
            if(isat->mw_idx[k]>0){
                int j=isat->mw_idx[k];
                double var0=isat->mw_var[k];
                double var1=(w1-w0)*(w1-w0)-var0;
                var1=var0+var1/j;

                isat->smw[k]=(w0*j+w1)/(j+1);
                isat->mw_idx[k]++;
                isat->mw_var[k]=var1;
            }
            else{
                isat->smw[k]=w1;
                isat->mw_idx[k]++;
                isat->mw_var[k]=0.25;
            }
        }
#endif

    }
}

void PPPSolver::ParTimeUpdate() {
    para_.popt_=popt_;
    para_.nx_=full_nx_;
    para_.x_=&full_x_[0];
    para_.Px_=&full_Px_[0];
    para_.last_sol_=epoch_sol_;
    para_.epoch_=epoch_count_;
    para_.isat_=isat_;
    para_.obss_=epoch_rover_;
    para_.glo_fcn_=nav_.glo_frq_num_;
    if(epoch_count_==1) para_.tt_=0.0;
    else para_.tt_=fabs(dt_);
    if(popt_.GNSS_opt_.ion_opt_==ION_CONST){
        para_.gim_model_.tecs_=ion_corr_->tecs_;
        para_.gim_model_.isat_=isat_;
    }
    para_.PosParUpdate();
    para_.ClkParUpdate();
    para_.DcbParUpdate();
    para_.IfbParUpdate();
    para_.GloParUpdate();
    para_.TrpParUpdate();
    para_.IonParUpdate();
    para_.AmbParUpdate();
}

int PPPSolver::PPPAlgorithm() {
    Time t=epoch_sol_.sol_time_;

    if(!P3SPP()){
        epoch_sol_.sol_time_=epoch_rover_.sat_infos_[0].sig_recep_;
        LOG(WARNING)<<"Precise point positioning pre-process error";
        return 0;
    }

    if(t.time_!=0.0) dt_=epoch_rover_.sat_infos_[0].sig_recep_.TimeDiff(t);

    if(pre_eph_->PrecPosClk(epoch_rover_)<=4){
        LOG(WARNING)<<epoch_count_<<"th no enough vaild satellite";
        epoch_sol_.sol_time_=epoch_rover_.sat_infos_[0].sig_recep_;
        return 0;
    }

    UpdateSatInfo(popt_);

    double rr[3]={0};
    for(int i=0;i<3;i++) rr[i]=full_x_[i];
    double dantr[NFREQ]={0},dants[NFREQ]={0};


    for(int i=0;i<epoch_rover_.num_;i++){
        sig_sat_data_=&epoch_rover_.sat_infos_[i];
        ReSetIsat(sig_sat_data_->sat_.sat_no_);
        if(popt_.GNSS_opt_.sat_pcv_)
            ant_corr_.SatPcvCorr(*sig_sat_data_,rr,dants);
        if(popt_.GNSS_opt_.rec_pcv_)
            ant_corr_.RecAntCorr(0,*sig_sat_data_,dantr);
        if(popt_.GNSS_opt_.ppp_opt_.phw_){
            att_corr_.phw_=isat_[sig_sat_data_->sat_.sat_no_-1].phw;
            isat_[sig_sat_data_->sat_.sat_no_-1].phw=att_corr_.SatPhw(*sig_sat_data_,rr);
        }
        CorrPRMeas(popt_,*sig_sat_data_,dantr,dants);
        CorrCPMeas(popt_,*sig_sat_data_,dantr,dants,att_corr_.phw_);
    }
    sig_sat_data_=nullptr;

    PPPCycSlip();
    ParTimeUpdate();

    if(SolverPPP()) {
        epoch_sol_.ppp_flag_=1;
        epoch_sol_.stat_=SOL_PPP;
    }
    else{
        epoch_sol_.ppp_flag_=0;
        epoch_sol_.stat_=SOL_P3SPP;
        epoch_sol_.sol_time_=epoch_rover_.sat_infos_[0].sig_recep_;
        LOG(WARNING)<<epoch_count_<<"th ppp solver fail";
    }
    UpdateIsat();
    int it=para_.IdxTrpPar();
    UpdateSol(full_x_,full_Px_,trp_corr_->z_hyd_,full_x_[it]);

    return epoch_sol_.stat_;
}

void PPPSolver::UpdateIsat() {

    for(int i=0;i<MAXSAT;i++) isat_[i].stat=NO_USE;

    for(int i=0;i<epoch_rover_.num_;i++){
        int sat=epoch_rover_.sat_infos_[i].sat_.sat_no_;
        for(int j=0;j<para_.NumUseFrq();j++){
            int ia=para_.IdxAmbPar(sat,j);
            if(epoch_rover_.sat_infos_[i].stat_>USED) continue;
            isat_[sat-1].stat=USED;
            isat_[sat-1].last_ep=epoch_count_;
            isat_[sat-1].last_time=epoch_rover_.sat_infos_[0].sig_recep_;
            isat_[sat-1].amb[j]=full_x_[ia];
            isat_[sat-1].lock[j]++;
            isat_[sat-1].outc[j]=0;
        }
    }
}

double PPPSolver::ModPRMeas(int use_f,int idx_f,vector<double> x) {
    int sat=sig_sat_data_->sat_.sat_no_;
    double mod_pr=0.0,C=0.0;
    double rr[3]={0},tide_dr[3]={0};
    Time t=sig_sat_data_->sig_recep_;
    for(int i=0;i<3;i++) rr[i]=x[i];

    if(popt_.GNSS_opt_.tid_opt_>TIDE_OFF){
        tid_corr_.TidCorr(*t.GPST2UTC(),0,rr,tide_dr);
    }
    Vec3 rec_pos(rr[0],rr[1],rr[2]);
    Vec3 blh=Xyz2Blh(rec_pos,WGS84);
    isat_[sat-1].dist=GeoDist(sig_sat_data_->pos_,rec_pos,sig_sat_data_->sig_vec_);
    isat_[sat-1].sagnac=sig_sat_data_->SagnacCorr(rec_pos);
    double el=SatAzEl(blh,sig_sat_data_->sig_vec_,sig_sat_data_->azel_);
    if(el<popt_.GNSS_opt_.ele_min_*D2R||isat_[sat-1].dist<=0.0) return 0.0;

    trp_corr_->para_.popt_=popt_;
    trp_corr_->par_X_=x;
    trp_corr_->isat_=isat_;
    trp_corr_->sat_data_=sig_sat_data_;
    trp_corr_->TrpCorr(blh,0.0);
    C=SQR(sig_sat_data_->lam_[use_f]/sig_sat_data_->lam_[sig_sat_data_->use_f1_]);

    if(popt_.GNSS_opt_.ion_opt_==ION_IF){
        isat_[sat-1].sion=0.0;
    }
    else if(popt_.GNSS_opt_.ion_opt_==ION_TEC){
        ion_corr_->sat_data_=sig_sat_data_;
        ion_corr_->isat_=isat_;
        ion_corr_->IonCorr(blh);
    }
    else if(popt_.GNSS_opt_.ion_opt_==ION_CONST){
        int ii=para_.IdxIonPar(sat);
        isat_[sat-1].sion=x[ii];
    }
    else{
        int ii=para_.IdxIonPar(sat);
        isat_[sat-1].sion=x[ii];
    }

    double bd2_multipath=0.0;
    if(sig_sat_data_->sat_.sat_sys_==BD2){
        bia_corr_.sat_data_=sig_sat_data_;
        bd2_multipath=bia_corr_.BDMultiPathCorr();
    }

    isat_[sat-1].shapiro=sig_sat_data_->ShapiroCorr(rec_pos);

    mod_pr=isat_[sat-1].dist+isat_[sat-1].sagnac+RecClkCorr(x)-sig_sat_data_->clk_err_[0]*CLIGHT
           +isat_[sat-1].strp+C*isat_[sat-1].sion+isat_[sat-1].shapiro+bd2_multipath;


    if(para_.NumIfbPar()){
        int idx=0;
        if(popt_.GNSS_opt_.use_frq_==TRIPLE_FRQ){
            if(popt_.GNSS_opt_.ion_opt_==ION_IF) idx=1;
            else idx=2;
        }
        if(idx_f==idx){
            double a = RecIfbCorr(x);
            mod_pr+=RecIfbCorr(x);
        }
    }

    if(sig_sat_data_->sat_.sat_sys_==GLO)
        mod_pr+=GloIcbCorr(popt_,x,OBS_PR);

    return mod_pr;
}

double PPPSolver::ModCPMeas(int use_f,int idx_f,vector<double> x) {
    int sat=sig_sat_data_->sat_.sat_no_;
    double mod_cp=0.0;
    double C=0.0,fact=1.0;
    double rr[3]={0},tide_dr[3]={0};
    for(int i=0;i<3;i++) rr[i]=x[i];
    Time t=sig_sat_data_->sig_recep_;
    if(popt_.GNSS_opt_.tid_opt_>TIDE_OFF){
        tid_corr_.TidCorr(*t.GPST2UTC(),0,rr,tide_dr);
        for(int i=0;i<3;i++) rr[i]+=tide_dr[i];
    }
    Vec3 rec_pos(rr[0],rr[1],rr[2]);
    Vec3 blh=Xyz2Blh(rec_pos,WGS84);
    isat_[sat-1].dist=GeoDist(sig_sat_data_->pos_,rec_pos,sig_sat_data_->sig_vec_);
    isat_[sat-1].sagnac=sig_sat_data_->SagnacCorr(rec_pos);
    double el=SatAzEl(blh,sig_sat_data_->sig_vec_,sig_sat_data_->azel_);
    if(el<popt_.GNSS_opt_.ele_min_*D2R||isat_[sat-1].dist<=0.0) return 0.0;
    isat_[sat-1].shapiro=sig_sat_data_->ShapiroCorr(rec_pos);

    trp_corr_->para_=para_;
    trp_corr_->para_.popt_=popt_;
    trp_corr_->par_X_=x;
    trp_corr_->isat_=isat_;
    trp_corr_->sat_data_=sig_sat_data_;
    trp_corr_->TrpCorr(blh,0.0);
    C=SQR(sig_sat_data_->lam_[use_f]/sig_sat_data_->lam_[sig_sat_data_->use_f1_]);

    if(popt_.GNSS_opt_.ion_opt_==ION_IF){
        isat_[sat-1].sion=0.0;
    }
    else if(popt_.GNSS_opt_.ion_opt_==ION_TEC){
        ion_corr_->sat_data_=sig_sat_data_;
        ion_corr_->isat_=isat_;
        ion_corr_->IonCorr(blh);
    }
    else if(popt_.GNSS_opt_.ion_opt_==ION_CONST){
        int ii=para_.IdxIonPar(sat);
        isat_[sat-1].sion=x[ii];
    }else{
        int ii=para_.IdxIonPar(sat);
        isat_[sat-1].sion=x[ii];
    }

    if(popt_.GNSS_opt_.use_frq_==SINGLE_FRQ&&popt_.GNSS_opt_.ion_opt_==ION_IF) fact=0.5;
    int ia=para_.IdxAmbPar(sat,idx_f);
    mod_cp=isat_[sat-1].dist+isat_[sat-1].sagnac+RecClkCorr(x)-sig_sat_data_->clk_err_[0]*CLIGHT
           +isat_[sat-1].strp-C*isat_[sat-1].sion+isat_[sat-1].shapiro+x[ia]*fact;
    if(popt_.GNSS_opt_.use_frq_==SINGLE_FRQ&&popt_.GNSS_opt_.ion_opt_==ION_IF)
        mod_cp+=GloIcbCorr(popt_,x,OBS_PR);

    return mod_cp;
}

#define THRES_REJECT 5.0
int PPPSolver::PRErrEq(PrcOpt popt, int use_f,int idx_f, int iter, int post,vector<double>x) {
    int sat=sig_sat_data_->sat_.sat_no_;
    double mod_pr=0.0,corr_pr=0.0,omc=0.0,q_meas=0.0;
    if((corr_pr=GetCorrObs(popt,use_f,OBS_PR))==0.0) return 0;
    if((mod_pr=ModPRMeas(use_f,idx_f,x))==0.0) return 0;
    omc=corr_pr-mod_pr;
    q_meas=MeasVar(popt,OBS_PR)+isat_[sat-1].var_trp+isat_[sat-1].var_ion;

    // record large post residual
    if(post&&fabs(omc)>sqrt(q_meas)*THRES_REJECT){
        lar_res_.push_back(omc);
        lar_res_sat_.push_back(sig_sat_data_->sat_.sat_no_);
        lar_res_frq_.push_back(use_f);
        lar_res_type_.push_back(OBS_PR);
    }

    OMCs_.push_back(omc);
    R_vec_.push_back(q_meas);
    NormalEqLine(popt,use_f,idx_f,num_L_,full_nx_,OBS_PR);

    if(!post&&iter==1) isat_[sat-1].pri_res[0][use_f]=omc;
    else if(post) isat_[sat-1].post_res[0][use_f]=omc;

    char sat_info[MAXBUFF]={'\0'};
    double C=SQR(sig_sat_data_->lam_[use_f]/sig_sat_data_->lam_[sig_sat_data_->use_f1_]);
    sprintf(sat_info,"%5d %s P%d omc=%10.4f var=%7.4f el=%3.1f obs=%10.4f r=%10.4f dtr=%10.4f dts=%10.4f trp=%7.4f ion=%6.3f cbias=%7.4f shapiro=%4.3f",
            epoch_count_,sig_sat_data_->sat_.sat_id_.c_str(),use_f+1,omc,q_meas,sig_sat_data_->azel_[1]*R2D,corr_pr,
            isat_[sat-1].dist+isat_[sat-1].sagnac,RecClkCorr(x),sig_sat_data_->clk_err_[0]*CLIGHT,isat_[sat-1].strp,
            C*isat_[sat-1].sion,popt.GNSS_opt_.ion_opt_==ION_IF?0.0:isat_[sat-1].cbias[use_f],isat_[sat-1].shapiro);
    LOG(DEBUG)<<sat_info;

    num_PR_++;
    num_L_++;
    return num_PR_;
}

int PPPSolver::CPErrEq(PrcOpt popt,int use_f,int idx_f, int iter, int post,vector<double> x) {
    int sat=sig_sat_data_->sat_.sat_no_;
    double mod_cp=0.0,corr_cp=0.0,omc=0.0,q_meas=0.0;

    if((corr_cp=GetCorrObs(popt_,use_f,OBS_CP))==0.0){
//        sig_sat_data_->stat_=NO_USE;
        return 0;
    }
    if((mod_cp=ModCPMeas(use_f,idx_f,x))==0.0){
//        sig_sat_data_->stat_=NO_USE;
        return 0;
    }
    omc=corr_cp-mod_cp;
    q_meas=MeasVar(popt_,OBS_CP)+isat_[sat-1].var_trp;

    // record large post residual
    if(post&&fabs(omc)>sqrt(q_meas)*THRES_REJECT){
        lar_res_.push_back(omc);
        lar_res_sat_.push_back(sig_sat_data_->sat_.sat_no_);
        lar_res_frq_.push_back(use_f);
        lar_res_type_.push_back(OBS_CP);
    }

    OMCs_.push_back(omc);
    R_vec_.push_back(q_meas);
    NormalEqLine(popt,use_f,idx_f,num_L_,full_nx_,OBS_CP);

    if(!post&&iter==1) isat_[sat-1].pri_res[OBS_CP][use_f]=omc;
    else if(post) isat_[sat-1].post_res[OBS_CP][use_f]=omc;

    char sat_info[MAXBUFF]={'\0'};
    int ia=para_.IdxAmbPar(sat,idx_f);
    double C=SQR(sig_sat_data_->lam_[use_f]/sig_sat_data_->lam_[sig_sat_data_->use_f1_]);
    sprintf(sat_info,"%5d %s L%d omc=%10.4f var=%9.4f el=%3.1f obs=%10.4f r=%10.4f dtr=%10.4f dts=%10.4f trp=%7.4f ion=%6.3f amb=%7.4f shapiro=%4.3f",
            epoch_count_,sig_sat_data_->sat_.sat_id_.c_str(),use_f+1,omc,q_meas,sig_sat_data_->azel_[1]*R2D,corr_cp,
            isat_[sat-1].dist+isat_[sat-1].sagnac,RecClkCorr(x),sig_sat_data_->clk_err_[0]*CLIGHT,isat_[sat-1].strp,
            -C*isat_[sat-1].sion,x[ia],isat_[sat-1].shapiro);
    LOG(DEBUG)<<sat_info;

    num_CP_++;
    num_L_++;
    return num_CP_;
}

int PPPSolver::PPPPostRes(AdjModel adj_model,vector<double>x) {
    int qc_PR=0,qc_CP=0,qc=0;
    if((qc=MakePPPEqs(1,0,x))) return qc;
#if 1
    qc_.v_info_=adj_model.v_info_;
    qc_.PRv_info_=adj_model.PRv_info_;
    qc_.CPv_info_=adj_model.CPv_info_;
    if(!(popt_.GNSS_opt_.use_frq_==SINGLE_FRQ&&popt_.GNSS_opt_.ion_opt_==ION_IF)){
        qc_PR=qc_.PostResQC(popt_,epoch_rover_,adj_model.v_PR_,adj_model.norm_PRv_,num_PR_,OBS_PR);
        if(qc_PR) return qc_PR;
    }

    qc_CP=qc_.PostResQC(popt_,epoch_rover_,adj_model.v_CP_,adj_model.norm_CPv_,num_CP_,OBS_CP);
#endif
    return qc_CP;
}

int PPPSolver::MakePPPEqs(int post, int iter, vector<double> x) {
    int qc_flag=0;
    int n=epoch_rover_.num_,nf=para_.NumUseFrq(),frq=0,flag=0;
    num_PR_=num_CP_=num_L_=vaild_sat_num_=0;
    A_coe_.clear();OMCs_.clear();R_meas_.clear();R_vec_.clear();L_info_.clear();
    double rr[3];
    for(int i=0;i<3;i++) rr[i]=x[i];

    char buff[MAXBUFF]={'\0'};
    sprintf(buff,"%s position %10.4f %10.4f %10.4f",(post?"post":"prior"),rr[0],rr[1],rr[2]);
    LOG(DEBUG)<<"PPP iter("<<iter<<") "<<buff;
    for(int i=0;i<n;i++){
        sig_sat_data_=&epoch_rover_.sat_infos_[i];
        if(sig_sat_data_->stat_>USED&&sig_sat_data_->stat_!=EX_SLIP){
            ReSetIsat(sig_sat_data_->sat_.sat_no_);
            continue;
        }
        for(int f=0;f<nf;f++){
            switch(f){
                case 0: frq=sig_sat_data_->use_f1_;break;
                case 1: frq=sig_sat_data_->use_f2_;break;
                case 2: frq=sig_sat_data_->use_f3_;break;
                case 3: frq=sig_sat_data_->use_f4_;break;
            }

            if(!CPErrEq(popt_,frq,f,iter,post,x)) continue;
            flag=(sig_sat_data_->sat_.sat_no_<<8)|(OBS_CP<<4)|(f);
            L_info_.push_back(flag);

            if(popt_.GNSS_opt_.use_frq_==SINGLE_FRQ&&popt_.GNSS_opt_.ion_opt_==ION_IF){
                if(f==0) vaild_sat_num_++;
                if(sig_sat_data_->stat_!=EX_SLIP) sig_sat_data_->stat_=USED;
                continue;
            }

            if(!PRErrEq(popt_,frq,f,iter,post,x)) continue;
            flag=(sig_sat_data_->sat_.sat_no_<<8)|(OBS_PR<<4)|(f);
            L_info_.push_back(flag);
//            if(!post&&popt_.GNSS_opt_.use_frq_!=SINGLE_FRQ) qc_.PriResCheck(OMCs_[OMCs_.size()-2],OMCs_.back(),R_vec_.back());

            if(f==0) vaild_sat_num_++;
            if(sig_sat_data_->stat_!=EX_SLIP) sig_sat_data_->stat_=USED;
        }
    }

    // ionosphere constraint
    Vec3 rec_pos(rr[0],rr[1],rr[2]);
    Vec3 blh=Xyz2Blh(rec_pos,WGS84);
    if(popt_.GNSS_opt_.ion_opt_==ION_CONST){
        for(int i=0;i<n;i++){
            sig_sat_data_=&epoch_rover_.sat_infos_[i];
            if(sig_sat_data_->stat_>USED&&sig_sat_data_->stat_!=EX_SLIP) continue;
            int sat=sig_sat_data_->sat_.sat_no_;
            int iI=para_.IdxIonPar(sat);
            ion_corr_->sat_data_=sig_sat_data_;
            ion_corr_->isat_=isat_;
            ion_corr_->IonCorr(blh);
            double omc_v=isat_[sat-1].sion-x[iI];
            OMCs_.push_back(omc_v);
            for(int k=0;k<full_nx_;k++) A_coe_.push_back(k==iI?1.0:0.0);
            R_vec_.push_back(100);
            num_L_++;
            char ion_info[MAXBUFF]={'\0'};
            sprintf(ion_info,"%5d %s I  omc=%10.4f var=%7.4f el=%3.1f",epoch_count_,sig_sat_data_->sat_.sat_id_.c_str(),omc_v,isat_[sat-1].var_ion*3.0,sig_sat_data_->azel_[1]*R2D);
            LOG(DEBUG)<<ion_info;
        }
    }

    if(post&&lar_res_.size()>0){
        qc_flag=qc_.PostResCheck(lar_res_,lar_res_sat_,lar_res_frq_,lar_res_type_,epoch_rover_);
        lar_res_.clear();lar_res_sat_.clear();lar_res_frq_.clear();lar_res_type_.clear();
    }
    return post?qc_flag:num_L_;
}

int PPPSolver::InitPPPSolver() {
    popt_=kPrcOpt;

    pre_eph_= new PrecEph;
    pre_eph_->InitSatEph(nav_);
    pre_eph_->ant_corr_=&ant_corr_;

    if(popt_.GNSS_opt_.trp_opt_>=TRP_EST_WET) trp_corr_=new EstTrp;

    switch(popt_.adj_opt_.mode_){
        case ADJ_KF: ppp_adj_=new KFAdj; break;
        case ADJ_HEL:ppp_adj_=new HelmertAdj;break;
    }

    switch(popt_.GNSS_opt_.ion_opt_){
        case ION_OFF:  ion_corr_= nullptr; break;
        case ION_BRDC: ion_corr_=new KlobIon; break;
        case ION_TEC:
        case ION_CONST:
            ion_corr_=new GIMIon;
            reader_gim_.InitReadGIM(kFileOpt.ion_);
            reader_gim_.ReadTec(ion_corr_->tecs_);
            break;
        case ION_IF:   ion_corr_=new IFIon; break;
    }

    para_.popt_=popt_;
    full_nx_=para_.NumRPar()+para_.NumAmbPar();
    full_x_.assign(full_nx_,0.0);
    full_Px_.assign(full_nx_*full_nx_,0.0);
    zip_x_.assign(full_x_.begin(),full_x_.end());
    zip_Px_.assign(full_Px_.begin(),full_Px_.end());

    if(epoch_count_==1){
        for(int i=0;i<MAXSAT;i++){
            isat_[i].last_time=epoch_rover_.sat_infos_[0].sig_recep_;
        }
    }
}

int PPPSolver::Algorithm() {
    if(epoch_count_==1) InitPPPSolver();
    return PPPAlgorithm();
}

MainSolver::MainSolver() {
    site_count_=0;
    solver_= nullptr;
}

MainSolver::~MainSolver() {

};

int MainSolver::InitMainSolver() {
    int spp,ppp;
    spp=kPrcOpt.mode_==SPP||kPrcOpt.mode_opt_==MDOPT_SPP;
    ppp=kPrcOpt.mode_==PPP||kPrcOpt.mode_opt_==MDOPT_PPP;


    if(solver_) delete solver_;
    if(spp) solver_=new SPPSolver;
    else if(ppp) solver_=new PPPSolver;

    if(!InitReader()){
        LOG(ERROR)<<"No exist rover observation";
        return 0;
    }

    InitSolver();

    if(out_sol_) delete out_sol_;
    out_sol_=new OutSol;
    out_sol_->InitOutSol(kFileOpt.sol_);
    if(kSolOpt.out_head_) out_sol_->WriteSolHead(solver_->epoch_rover_.sta_);
}

int MainSolver::InitReader() {
    int ppp,dgps,ppk;
    ppp=kPrcOpt.mode_==PPP||kPrcOpt.mode_opt_==MDOPT_PPP;
    dgps=kPrcOpt.mode_==DGPS||kPrcOpt.mode_opt_==MDOPT_DGPS;
    ppk=kPrcOpt.mode_==PPK||kPrcOpt.mode_opt_==MDOPT_PPK;

    if(reader_rov_) delete reader_rov_;
    reader_rov_=new ReadRnxObs;
    reader_rov_->InitReadObs(kFileOpt.rover_,ts_,te_,tint_);
    reader_rov_->ReadHead(&nav_,&solver_->epoch_rover_.sta_);
    if(dgps||ppk){
        reader_bas_->InitReadObs(kFileOpt.base_,ts_,te_,tint_);reader_bas_->ReadHead(&nav_,&solver_->epoch_base_.sta_);
        int stat=reader_bas_->TestOpen();
        if(stat==0) return 0;
    }
    int stat=reader_rov_->TestOpen();
    if(!stat) return 0;

    if(site_count_==1){
        // 读取导航星历
        reader_nav_.InitReadNav(kFileOpt.brdc_);reader_nav_.ReadHead(nav_);
        reader_nav_.ReadBody(nav_);reader_nav_.CloseFile();nav_.UpdateNav();

        // 读取码偏差
        if(kPrcOpt.GNSS_opt_.cb_prd_type_>CBC_OFF){
            if(kPrcOpt.GNSS_opt_.cb_prd_type_==CBC_DCB){
                reader_bia_.InitReadBias();reader_bia_.ReadBiasFile(bia_corr_.dcb_);
                bia_corr_.type_=reader_bia_.type_;
            }
        }

        // 读取erp/blq
        if(kPrcOpt.GNSS_opt_.tid_opt_>TIDE_OFF){
            reader_erp_.InitReadErp(kFileOpt.erp_);reader_erp_.ReadErpPara(tid_corr_.erp_);
            reader_blq_.InitReadBlq(kFileOpt.blq_);reader_blq_.ReadBlqPara(kSite,tid_corr_.ocean_par_[0]);
        }

        // 读取精密星历
        if(ppp){
            for(int i=0;i<3;i++){
                if(reader_pre_.InitReadPreEph(kFileOpt.sp3_[i])) reader_pre_.ReadPre(nav_);
                else continue;
            }
            reader_pre_.InitReadPreEph(kFileOpt.clk_);
            reader_pre_.ReadPre(nav_);
        }

        // 读取atx
        if(kPrcOpt.GNSS_opt_.rec_pcv_||kPrcOpt.GNSS_opt_.sat_pcv_){
            reader_ant_.InitReadAtx(kFileOpt.atx_);
            reader_ant_.ReadAtx(ant_corr_.ants_);
        }
    }

    // 设置天线信息
    if(kPrcOpt.GNSS_opt_.rec_pcv_||kPrcOpt.GNSS_opt_.sat_pcv_){
        ant_corr_.sta_[0]=solver_->epoch_rover_.sta_;
        ant_corr_.ant_time_=kPrcOpt.prc_date_;
        ant_corr_.SetAntPara();
    }
    return 1;
}

int MainSolver::InitSolver() {
    solver_->ant_corr_=ant_corr_;
    solver_->bia_corr_=bia_corr_;
    solver_->tid_corr_=tid_corr_;
    solver_->nav_=nav_;
}


int MainSolver::EpochSol() {
    reader_rov_->sys_mask_=kPrcOpt.GNSS_opt_.nav_sys_;
    while(reader_rov_->ReadBody(solver_->epoch_rover_))
    {
        if(!solver_->SortObs(solver_->epoch_rover_)){
            LOG(ERROR)<<"Sort rover observation error";
            return 0;
        }

        solver_->site_count_=site_count_;
        solver_->epoch_count_++;
        solver_->epoch_sol_.i_epoch_=solver_->epoch_count_;

        if(!solver_->Algorithm()){
            LOG(ERROR)<<solver_->epoch_count_<<"th epoch resolver fail";
        }

        out_sol_->isat_=solver_->isat_;
        out_sol_->WriteSols(solver_->epoch_sol_,solver_->epoch_rover_);
    }

    return 1;
}
