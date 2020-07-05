/**************************************************************************

Copyright(c): 2020, by Chao Chen, All rights reserved.
              China University of Mining and Technology, XuZhou, P.R. China

Author: Chao Chen(cchen@cumt.edu.cn)

Date: 2020-05-30

Description:

**************************************************************************/

#include "QuaCtrl.h"
#include "CmnFunc.h"

QuaCtrl::QuaCtrl() {
    epoch_=0;
    isat_= nullptr;
}

QuaCtrl::~QuaCtrl() {

}

double QuaCtrl::MWMeas(Isat_t &isat, ObsData data,double dcb[MAXSAT][MAXCBIASPAIR], int flag) {
    // 探测周跳只需修正C码
    int sat=data.sat_.sat_no_;
    int i=0,j=1;
    if(flag==1) j=2;
    double P1=data.P_[i],P2=data.P_[j];
    double L1=data.L_[i],L2=data.L_[j];
    double lami=data.lam_[i],lamj=data.lam_[j];
    double P1_C1=0.0,P2_C2=0.0;
    if(data.sat_.sat_sys_==GPS){
        P1_C1=dcb[data.sat_.sat_no_-1][GPS_C1CC1W];
        P2_C2=dcb[data.sat_.sat_no_-1][GPS_C2CC2W];
        P1-=P1_C1;P2-=P2_C2;
    }

    if(P1==0.0||P2==0.0||L1==0.0||L2==0.0) return 0.0;
    double mw=(L1-L2)-(lamj-lami)/(lami+lamj)*(P1/lami+P2/lamj);
    return mw;
}

double QuaCtrl::GFMeas(Isat_t &isat, ObsData data, int flag) {
    int i=0,j=1;
    if(flag==1) j=2;
    double L1=isat.corr_L[i],L2=isat.corr_L[j];
    if(L1==0.0||L2==0.0) return 0.0;
    return L1-L2;
}

void QuaCtrl::GFCycSlip(PrcOpt popt, Obss &obss, double dt) {
    int del_ep=1,sat;
    ObsData *data= nullptr;
    Isat_t *isat= nullptr;
    del_ep=Round(fabs(dt/popt.ti_));

    for(int i=0;i<obss.num_;i++){
        data=&obss.sat_infos_[i];
        sat=data->sat_.sat_no_;
        isat=&isat_[sat-1];
        double el=data->azel_[1]*R2D;

#if 0
        double g1=GFMeas(*isat,*data, 0);
        if(g1==0.0) continue;
        double thres=0.0;
        if(el<popt.GNSS_opt_.ele_min_) continue;
        if(el>15.0) thres=popt.QC_opt_.slip_opt_.thres_gf_;
        else{
            thres=-popt.QC_opt_.slip_opt_.thres_gf_/15.0*el+2.0*popt.QC_opt_.slip_opt_.thres_gf_;
        }

        double g0=isat->gf;
        if(g0==0.0) continue;

        double dgf=g1-g0;
        if(fabs(dgf)>MIN(thres*del_ep,1.5)){
            data->stat_=EX_SLIP;
            char buff[MAXBUFF]={'\0'};
            sprintf(buff,"%d %s %s dgf=%5.3f - %5.3f el=%3.1f thres=%4.1f",epoch_,data->sig_recep_.time_str_.c_str(),
                    data->sat_.sat_id_.c_str(),g1,g0,el,MIN(thres,1.5));
            LOG(WARNING)<<buff;
        }
#else
        int f=data->frq3_?2:1;
        for(int k=0;k<f;k++){
            double g1=GFMeas(*isat,*data, k);
            if(g1==0.0) continue;
            double thres=0.0;
            if(el<popt.GNSS_opt_.ele_min_) continue;
            if(el>15.0) thres=popt.QC_opt_.slip_opt_.thres_gf_;
            else{
                thres=-popt.QC_opt_.slip_opt_.thres_gf_/15.0*el+2.0*popt.QC_opt_.slip_opt_.thres_gf_;
            }

            double g0=isat->gf[k];
            if(g0==0.0) continue;

            double dgf=g1-g0;
            if(fabs(dgf)>MIN(thres*del_ep,1.5)){
                data->stat_=EX_SLIP;
                char buff[MAXBUFF]={'\0'};
                sprintf(buff,"%d %s %s dgf%d=%5.3f - %5.3f el=%3.1f thres=%4.1f",epoch_,data->sig_recep_.time_str_.c_str(),
                        data->sat_.sat_id_.c_str(),k+1,g1,g0,el,MIN(thres,1.5));
                LOG(WARNING)<<buff;
            }
        }

#endif

    }
}

void QuaCtrl::MWCycSlip(PrcOpt popt, Obss &obss,double dt,double dcb[MAXSAT][MAXCBIASPAIR]) {
#if 0
    int del_ep=1;
    double fact=1.0;
    ObsData *data= nullptr;
    Isat_t *isat= nullptr;
    if(popt.ti_>0.0) del_ep=Round(dt/popt.ti_);
    if(popt.ti_<=1.5){
        if(fabs(dt)<=10.0) del_ep=1;
        else if(fabs(dt)<=15.0) del_ep=2;
        else if(fabs(dt)<=22.0) del_ep=3;
    }
    if(popt.ti_>=29.5){
        if(del_ep<=2) fact=1.0;
        else if(del_ep<=4) fact=1.25;
        else if(del_ep<=6) fact=1.5;
        else fact=2.0;
    }
    for(int i=0;i<obss.num_;i++){
        data=&obss.sat_infos_[i];
        int sat=data->sat_.sat_no_;
        double el=data->azel_[1]*R2D;
        isat=&isat_[sat-1];

        if(data->sig_recep_.TimeDiff(isat->last_time)>popt.ti_){
            isat->smw=isat->mw_idx=0;
        }
        double w1=MWMeas(*isat,*data,dcb,0);
        double w0=isat->smw;
        if(w1==0.0||w0==0.0){
            isat->smw=isat->mw_idx=0;
            continue;
        }
        double thres=0.0;
        if(el>20.0) thres=popt.QC_opt_.slip_opt_.thres_mw_;
        else thres=-popt.QC_opt_.slip_opt_.thres_mw_*0.1*el+3.0*popt.QC_opt_.slip_opt_.thres_mw_;
        if(el<popt.GNSS_opt_.ele_min_) continue;
        double dmw=w1-w0;
        if(fabs(dmw)>MIN(thres*fact,6.0)){
            data->stat_=EX_SLIP;
            char buff[MAXBUFF]={'\0'};
            sprintf(buff,"%s %s dmw=%5.3f - %5.3f el=%3.1f thres=%4.1f",data->sig_recep_.time_str_.c_str(),
                    data->sat_.sat_id_.c_str(),w1,w0,el,MIN(thres,6.0));
            LOG(WARNING)<<buff;
        }
    }
#else
    int del_ep=1;
    double fact=1.0;
    ObsData *data= nullptr;
    Isat_t *isat= nullptr;
    if(popt.ti_>0.0) del_ep=Round(dt/popt.ti_);
    if(popt.ti_<=1.5){
        if(fabs(dt)<=10.0) del_ep=1;
        else if(fabs(dt)<=15.0) del_ep=2;
        else if(fabs(dt)<=22.0) del_ep=3;
    }
    if(popt.ti_>=29.5){
        if(del_ep<=2) fact=1.0;
        else if(del_ep<=4) fact=1.25;
        else if(del_ep<=6) fact=1.5;
        else fact=2.0;
    }

    for(int i=0;i<obss.num_;i++){
        data=&obss.sat_infos_[i];
        int sat=data->sat_.sat_no_;
        double el=data->azel_[1]*R2D;
        isat=&isat_[sat-1];

        // f=0 L1_L2 f=1 L1_L5
        int f=data->frq3_?2:1;
        for(int k=0;k<f;k++){
            if(data->sig_recep_.TimeDiff(isat->last_time)>popt.ti_){
                isat->smw[k]=isat->mw_idx[k]=0;
            }
            double w1=MWMeas(*isat,*data,dcb,k);
            double w0=isat->smw[k];
            if(w1==0.0||w0==0.0){
                isat->smw[k]=isat->mw_idx[k]=0;
                continue;
            }
            double thres=0.0;
            if(el>20.0) thres=popt.QC_opt_.slip_opt_.thres_mw_;
            else thres=-popt.QC_opt_.slip_opt_.thres_mw_*0.1*el+3.0*popt.QC_opt_.slip_opt_.thres_mw_;
            if(el<popt.GNSS_opt_.ele_min_) continue;
            double dmw=w1-w0;
            if(fabs(dmw)>MIN(thres*fact,6.0)){
                data->stat_=EX_SLIP;
                char buff[MAXBUFF]={'\0'};
                sprintf(buff,"%s %s dmw%d=%5.3f - %5.3f el=%3.1f thres=%4.1f",data->sig_recep_.time_str_.c_str(),
                        data->sat_.sat_id_.c_str(),k+1,w1,w0,el,MIN(thres,6.0));
                LOG(WARNING)<<buff;
            }
        }
    }

#endif
}

void QuaCtrl::PriResCheck(double omc_CP, double omc_PR, double& var_PR) {
    if(fabs(omc_CP-omc_PR)>=0.2){
        var_PR*=100;
    }
}

int QuaCtrl::PostResCheck(vector<double> res, vector<int> sats, vector<int> frqs, vector<int>types,Obss &obss) {
    double max_v=0.0;
    ObsData *data= nullptr;
    int max_sat=0,max_frq=0,type;
    max_v=res[0];max_sat=sats[0];max_frq=frqs[0];
    for(int i=1;i<res.size();i++){
        if(fabs(max_v)>=fabs(res[i])) continue;
        max_v=res[i];max_sat=sats[i];max_frq=frqs[i];type=types[i];
    }
    for(int i=0;i<obss.num_;i++){
        if(obss.sat_infos_[i].sat_.sat_no_==max_sat){
            data=&obss.sat_infos_[i];
            data->stat_=type?EX_CP_POST_RES:EX_PR_POST_RES;
            char buff[MAXBUFF]={'\0'};
            sprintf(buff,"%d %s %s reject by large post %s%d residual: v=%6.3f",epoch_,data->sig_recep_.time_str_.c_str(),
                    data->sat_.sat_id_.c_str(),(type?"L":"P"),max_frq+1,max_v);
            LOG(WARNING)<<buff;
            return 1;
        }else continue;
    }
    return 0;
}

int QuaCtrl::PostResQC(PrcOpt popt,Obss &obss,vector<double>v,vector<double>norm_v, int nv,int type) {
    ObsData *data= nullptr;
    int qc_flag=0;
    double thres_PR=3.0,thres_CP=0.03,thres_normv=2.0;
    vector<int>info;
    if(type==OBS_PR) info.assign(PRv_info_.begin(),PRv_info_.end());
    else info.assign(CPv_info_.begin(),CPv_info_.end());

    if(popt.GNSS_opt_.use_frq_==SINGLE_FRQ&&popt.GNSS_opt_.ion_opt_==ION_IF) thres_CP=thres_PR;

    for(int i=0;i<nv;i++){v[i]=fabs(v[i]);norm_v[i]=fabs(norm_v[i]);}
    auto max_v=max_element(begin(v),end(v));
    auto max_norm_v=max_element(begin(norm_v),end(norm_v));
    int idx_max_v=distance(begin(v),max_v);
    int idx_max_norm_v=distance(begin(norm_v),max_norm_v);
    int max_v_f=info[idx_max_v]&0x0F;
    int max_norm_v_f=info[idx_max_norm_v]&0x0F;
    int max_v_sat=(info[idx_max_v]>>8)&0xFF;
    int max_norm_v_sat=(info[idx_max_norm_v]>>8)&0xFF;

    char buff[MAXBUFF]={'\0'};
    double thres_v=type?thres_CP:thres_PR;
    for(int i=0;i<obss.num_;i++){
        data=&obss.sat_infos_[i];
        if(data->stat_>USED&&data->stat_!=EX_SLIP) continue;
        if(data->sat_.sat_no_==max_v_sat){
            if(data->sat_.sat_sys_==BD2&&(data->sat_.sat_prn_>=1&&data->sat_.sat_prn_<=5)){
                thres_v*=5.0;
            }
            if(data->sat_.sat_sys_==GLO){
                thres_v*=2.0;
            }
            double el=data->azel_[1];
            if(el*R2D>15.0) thres_v/=sin(el);

            if(*max_v>thres_v){
                sprintf(buff,"%4d %20s: %4s %s%d exclude by large post residual v=%7.3f, thres=%7.3f",epoch_,
                        data->sig_recep_.time_str_.c_str(),data->sat_.sat_id_.c_str(),(type?"L":"P"),max_v_f+1,*max_v,thres_v/(sin(el)));
                LOG(WARNING)<<buff;
                data->stat_=EX_PR_POST_RES;qc_flag=true;
                return qc_flag;
            }
        }
#if 0
        buff[0]={'\0'};
        if(data->sat_.sat_no_==max_norm_v_sat){
            if(*max_norm_v>thres_normv){
                sprintf(buff,"%4d %20s: %4s %s%d exclude by large norm pseudorange residual v=%7.3f, thres=%7.3f",epoch_,
                        data->sig_recep_.time_str_.c_str(),data->sat_.sat_id_.c_str(),(type?"L":"P"),max_norm_v_f+1,*max_norm_v,thres_normv);
                LOG(WARNING)<<buff;
                data->stat_=EX_PR_NOR_RES;qc_flag=true;
                return qc_flag;
            }
        }
#endif
    }
    return false;
}

#define MAX_VAR_EPH SQR(300.0)
void QuaCtrl::ExcludeSat(ObsData& obs) {
    int sat=obs.sat_.sat_no_;
    if(kPrcOpt.QC_opt_.exclude_sat_[sat-1]==EX_MAN){
        obs.stat_=EX_MAN;
    }
    if(obs.svh_<0){
        obs.stat_=EX_SVH;
    }
    if(obs.sat_var_>MAX_VAR_EPH){
        obs.stat_=EX_URA;
    }
}

