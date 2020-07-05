/**************************************************************************

Copyright(c): 2020, by Chao Chen, All rights reserved.
              China University of Mining and Technology, XuZhou, P.R. China

Author: Chao Chen(cchen@cumt.edu.cn)

Date: 2020-05-30

Description:

**************************************************************************/

#include "PPPLibGlo.h"
#include "ParSetting.h"
#include "TrpModel.h"
#include "AmbModel.h"
#include "IonModel.h"

Parameter::Parameter() {
    tt_=0;
    x_=Px_= nullptr;
    glo_fcn_= nullptr;
    epoch_=0;
}

Parameter::~Parameter()=default;

int Parameter::NumUseFrq() {
    if(popt_.GNSS_opt_.use_frq_==DUAL_FRQ){
        if(popt_.GNSS_opt_.ion_opt_==ION_IF) return 1;
        else return popt_.GNSS_opt_.use_frq_;
    }
    else if(popt_.GNSS_opt_.use_frq_==TRIPLE_FRQ){
        if(popt_.GNSS_opt_.ion_opt_==ION_IF){
            if(popt_.GNSS_opt_.IF3_format_==IF3_DUAL) return 2;
            // 三频单无电离层组合可能需要估计两种模糊度,0~MAXSAT: 双频无电离层模糊度f=0;MAXSAT~2*MAXSAT: 双频无电离层模糊度f=1
            else if(popt_.GNSS_opt_.IF3_format_==IF3_SINGLE) return 2;
        }
        else return popt_.GNSS_opt_.use_frq_;
    }
}

int Parameter::NumPosPar() {
    return popt_.GNSS_opt_.dynamic_modeling_?9:3;
}

int Parameter::NumClkPar() {
    int dgnss=0;
    dgnss=popt_.mode_==DGPS||popt_.mode_==PPK||popt_.mode_opt_==DGPS||popt_.mode_opt_==PPK;
    if(dgnss) return 0;
    else return NSYS;
}

int Parameter::NumSPPPar() {
    return NumPosPar()+NumClkPar();
}

// 电离层约束时才需要估计接收机端DCB参数
int Parameter::NumDcbPar() {
    if(popt_.GNSS_opt_.ion_opt_==ION_CONST){
        if(popt_.GNSS_opt_.use_frq_==SINGLE_FRQ) return 0;
        else return NSYS;
    }else return 0;
}

int Parameter::NumIfbPar() {
    return popt_.GNSS_opt_.use_frq_>=TRIPLE_FRQ?NSYS:0;
}

int Parameter::NumGloPar() {
    int spp=popt_.mode_==SPP||popt_.mode_opt_==MDOPT_SPP;
    if(spp) return 0;
    if(!(popt_.GNSS_opt_.nav_sys_&GLO)) return 0;
    if(!exist_sys_mark_[INDEXGLO]) return 0;
    if(popt_.adj_opt_.glo_icb_modeling_==GLOICB_OFF) return 0;
    else if(popt_.adj_opt_.glo_icb_modeling_==GLOICB_LNF) return 1;
    else if(popt_.adj_opt_.glo_icb_modeling_==GLOICB_QUAD) return 2;
    else return 0;
}

int Parameter::NumTrpPar() {
    return popt_.GNSS_opt_.trp_opt_<=TRP_SAAS?0:(popt_.GNSS_opt_.trp_opt_==TRP_EST_WET?1:3);
}

int Parameter::NumIonPar() {
    return popt_.GNSS_opt_.ion_opt_>ION_IF?MAXSAT:0;

}

int Parameter::NumRPar() {
    return NumPosPar()+NumClkPar()+NumDcbPar()+NumIfbPar()+NumGloPar()+NumTrpPar()+NumIonPar();
}

int Parameter::NumAmbPar() {
    if(popt_.mode_==SPP||popt_.mode_opt_==MDOPT_SPP) return 0;
    return NumUseFrq()*MAXSAT;
}

int Parameter::IdxClkPar(int i) {
    return NumPosPar()+i;
}

int Parameter::IdxDcbPar() {
    return NumPosPar()+NumClkPar();
}

int Parameter::IdxIfbPar(int i) {
    return NumPosPar()+NumClkPar()+NumDcbPar()+i;
}

int Parameter::IdxGloPar() {
    return NumPosPar()+NumClkPar()+NumDcbPar()+NumIfbPar();
}

int Parameter::IdxTrpPar() {
    int a= NumGloPar();
    return NumPosPar()+NumClkPar()+NumDcbPar()+NumIfbPar()+NumGloPar();
}

int Parameter::IdxIonPar(int sat) {
    return NumPosPar()+NumClkPar()+NumDcbPar()+NumIfbPar()+NumGloPar()+NumTrpPar()+sat-1;
}

int Parameter::IdxAmbPar(int sat,int f) {
    return NumRPar()+f*MAXSAT+sat-1;
}

void Parameter::CoePosPar(vector<double>& A_coe) {
    A_coe.push_back(-sat_info_.sig_vec_.i_);
    A_coe.push_back(-sat_info_.sig_vec_.j_);
    A_coe.push_back(-sat_info_.sig_vec_.k_);
}

void Parameter::CoeClkPar(vector<double>& A_coe) {
    if(!NumClkPar()) return;
    int sys=sat_info_.sat_.sat_sys_;
    int npos=NumPosPar();
    int idx_clock=0;
    for(int i=0;i<NSYS;i++){
        if(exist_sys_mark_[i]!=0){ idx_clock=i;break;}
    }
    for(int j=npos;j<nx_;j++) A_coe.push_back(j==idx_clock+npos?1.0:0.0);
    int idx=0;
    if(sys==BD2) idx=INDEXBD2;
    else if(sys==BD3) idx=INDEXBD3;
    else if(sys==GAL) idx=INDEXGAL;
    else if(sys==GLO) idx=INDEXGLO;
    else if(sys==QZS) idx=INDEXQZS;
    A_coe[idx+npos+num_L_*nx_]=1.0;
}

void Parameter::CoeDcbPar(vector<double> &A_coe) {
    if(!NumDcbPar()) return;

    int id=IdxDcbPar();

    A_coe[id+num_L_*nx_]=1.0;
}

void Parameter::CoeIfbPar(vector<double> &A_coe) {
    if(!NumIfbPar()) return;

    int sys=sat_info_.sat_.sat_sys_idx_;
    int ifb=IdxIfbPar(sys);
    A_coe[ifb+num_L_*nx_]=1.0;
}

void Parameter::CoeGloPar(vector<double> &A_coe) {
    if(!NumGloPar()) return;

    int idx=IdxGloPar();
    int prn=sat_info_.sat_.sat_prn_;
    int frq=glo_fcn_[prn-1];
    A_coe[idx+num_L_*nx_]=frq;
    if(popt_.adj_opt_.glo_icb_modeling_==GLOICB_QUAD){
        A_coe[idx+1+num_L_*nx_]=frq*frq;
    }
}

void Parameter::CoeTrpPar(vector<double> &A_coe) {
    if(!NumTrpPar()) return;
    int it=IdxTrpPar();
    int sat=sat_info_.sat_.sat_no_;

    A_coe[it+num_L_*nx_]=isat_[sat-1].map_trp_w;
}

void Parameter::CoeIonPar(vector<double> &A_coe,int use_f,int type) {
    if(!NumIonPar()) return;

    double C=SQR(sat_info_.lam_[use_f]/sat_info_.lam_[sat_info_.use_f1_])*(type?-1.0:1.0);
    int sat=sat_info_.sat_.sat_no_;
    int ii=IdxIonPar(sat);

    A_coe[ii+num_L_*nx_]=C;
}

void Parameter::CoeAmbPar(int f, vector<double> &A_coe) {
    if(!NumAmbPar()) return;
    int sat=sat_info_.sat_.sat_no_;
    int ia=IdxAmbPar(sat,f);

    if(popt_.GNSS_opt_.use_frq_==SINGLE_FRQ&&popt_.GNSS_opt_.ion_opt_==ION_IF) A_coe[ia+num_L_*nx_]=0.5;
    else A_coe[ia+num_L_*nx_]=1.0;
}

int Parameter::InitX(double xi, double var, int i) {
    x_[i]=xi;
    for(int j=0;j<nx_;j++){
        Px_[i+j*nx_]=Px_[j+i*nx_]=i==j?var:0.0;
    }
}

#define POS_VAR SQR(60.0)
#define VEL_VAR SQR(60.0)
void Parameter::PosParUpdate() {
    int i=0;
    if(popt_.mode_==FIX){
        PosOpt posopt=popt_.GNSS_opt_.pos_opt_;
        double rr[3]={posopt.rover_.i_,posopt.rover_.j_,posopt.rover_.k_};
        for(i=0;i<3;i++) InitX(rr[i],1E-8,i);
        return;
    }
    if(popt_.mode_opt_==MDOPT_KINE||popt_.mode_opt_==MDOPT_KINE_SIM){
        // 动态定位使用SPP解初始化位置信息
        double rr[3]={last_sol_.pos_.i_,last_sol_.pos_.j_,last_sol_.pos_.k_};
        double rv[3]={last_sol_.vel_.i_,last_sol_.vel_.j_,last_sol_.vel_.k_};
        for(i=0;i<3;i++) InitX(rr[i],POS_VAR,i);
        return;
    }

    if(popt_.mode_opt_==MDOPT_STATIC){
        if(epoch_==1){
            //静态定位的首历元使用SPP解进行初始化，并方差初始化
            double rr[3]={last_sol_.pos_.i_,last_sol_.pos_.j_,last_sol_.pos_.k_};
            for(i=0;i<NumPosPar();i++) InitX(rr[i],POS_VAR,i);
        }
        return;
    }
}

#define CLK_VAR SQR(60.0)
void Parameter::ClkParUpdate() {
    if(!NumClkPar()) return;
    // 钟差参数每个历元需要使用spp进行初始化,isb可以建模成随机游走或常数模型
    // 需要先判断哪个系统是基准钟差
    int npos=NumPosPar();
    double clk[NSYS]={last_sol_.clk_G_,last_sol_.B2_isb_,last_sol_.B3_isb_,
                      last_sol_.E_isb_,last_sol_.R_isb_,last_sol_.J_isb_};
    int idx_clock=0;

    if(epoch_==1){
        for(int i=0;i<NSYS;i++){
            InitX(clk[i],CLK_VAR,npos+idx_clock+i);
        }
    }
    else{
        for(int i=0;i<NSYS;i++){
            if(exist_sys_mark_[i]!=0){ idx_clock=i;break;}
        }
        // 钟差建模为白噪声每历元重新初始化
        InitX(clk[idx_clock],CLK_VAR,npos+idx_clock);

        for(int i=0;i<NSYS;i++){
            if(i==idx_clock) continue;
            if(!exist_sys_mark_[i]) continue;
            if(popt_.adj_opt_.gnss_isb_modeling_==STO_RW){
                int ic=i+npos;
                Px_[ic+ic*nx_]+=SQR(0.001)*tt_;
            }
            else if(popt_.adj_opt_.gnss_isb_modeling_==STO_WN){
                InitX(clk[idx_clock+i],CLK_VAR,npos+i);
            }
        }
    }
}

#define DCB_VAR  SQR(30.0)
void Parameter::DcbParUpdate() {
    if(!NumDcbPar()) return;

    int id=IdxDcbPar();

    if(epoch_==1){
        InitX(1E-6,DCB_VAR,id);
    }
    else{
        //接收机的硬件延迟建模
    }
}

#define IFB_VAR  SQR(30.0)
void Parameter::IfbParUpdate() {
    if(!NumIfbPar()) return;

    if(epoch_==1){
        for(int i=0;i<NSYS;i++){
            if(!exist_sys_mark_[i]) continue;
            if(i==INDEXGLO) continue;
            int ifb=IdxIfbPar(i);
            InitX(0.01,IFB_VAR,ifb);
        }
    }
    else{
        for(int i=0;i<NSYS;i++){
            if(!exist_sys_mark_[i]) continue;
            int ifb=IdxIfbPar(i);
            Px_[ifb+ifb*nx_]+=1.0e-8*fabs(tt_);
        }
    }
}

#define GLO_ICB_VAR  SQR(60.0)
void Parameter::GloParUpdate() {
    if(!NumGloPar()) return;

    int idx=IdxGloPar();

    if(epoch_==1){
        InitX(0.01,GLO_ICB_VAR,idx);
        if(popt_.adj_opt_.glo_icb_modeling_==GLOICB_QUAD){
            InitX(0.01,GLO_ICB_VAR,idx+1);
        }
    }else{
        Px_[idx+idx*nx_]+=1e-6*fabs(tt_);
        if(popt_.adj_opt_.glo_icb_modeling_==GLOICB_QUAD){
            Px_[idx+1+(idx+1)*nx_]+=1e-6*fabs(tt_);
        }
    }
}

#define TRP_GRA_VAR   SQR(0.01)
void Parameter::TrpParUpdate() {
    if(!NumTrpPar()) return;
    int it=IdxTrpPar();

    if(epoch_==1){
        double var=SQR(0.3);
        InitX(0.175,var,it);
        if(popt_.GNSS_opt_.trp_opt_==TRP_EST_GRAD){
            for(int j=it+1;j<it+3;j++) InitX(1E-6,TRP_GRA_VAR,j);
        }
    }
    else{
        // 对流层延迟建模成随机游走
        Px_[it+it*nx_]+=popt_.adj_opt_.gnss_psd_[2]*tt_;
        if(popt_.GNSS_opt_.trp_opt_==TRP_EST_GRAD){
            for(int j=it+1;j<it+3;j++) Px_[j+j*nx_]+=popt_.adj_opt_.gnss_psd_[2]*0.1*tt_;
        }
    }
}

#define GAP_RESION  120
#define VAR_ION  SQR(60.0)
void Parameter::IonParUpdate() {
    if(!NumIonPar()) return;
    int i,j;
    double ion=0.0;
    ObsData *data= nullptr;
    Isat_t *isat= nullptr;

    for(i=0;i<MAXSAT;i++){
        j=IdxIonPar(i+1);
        if(x_[j]!=0.0&&isat_[i].outc[0]>GAP_RESION) x_[j]=0.0;
    }
    for(i=0;i<obss_.num_;i++){
        int sat=obss_.sat_infos_[i].sat_.sat_no_;
        data=&obss_.sat_infos_[i];
        isat=&isat_[sat-1];
        j=IdxIonPar(sat);

        if(x_[j]==0.0||!last_sol_.ppp_flag_){
            if(isat->corr_P[data->use_f1_]==0.0||isat->corr_P[data->use_f2_]==0.0){
                data->stat_=EX_NOPOBS;
                continue;
            }

            Vec3 blh=Xyz2Blh(last_sol_.pos_,WGS84);
            if(popt_.GNSS_opt_.use_frq_==SINGLE_FRQ&&popt_.GNSS_opt_.ion_opt_==ION_UC){
                KlobIon klob;
                klob.sat_data_=data;
                klob.isat_=isat_;
                ion=klob.IonCorr(blh);
            }
            else if(popt_.GNSS_opt_.ion_opt_==ION_CONST){
                gim_model_.sat_data_=data;
                gim_model_.IonCorr(blh);
                ion=isat_[sat-1].sion;
            }
            else{
                ion=(isat->corr_P[data->use_f1_]-isat->corr_P[data->use_f2_])/(1.0-SQR(data->lam_[data->use_f2_]/data->lam_[data->use_f1_]));
//                ion=(data->P_[0]-data->P_[2])/(1.0-SQR(data->lam_[2]/data->lam_[0]));
                IonModel ion_model;
                double map=ion_model.IonMapFun(blh,data->azel_[1]);
                ion/=map;
            }
            InitX(ion,VAR_ION,j);
        }
        else{
            if(isat->stat==NO_USE){
                InitX(x_[j],VAR_ION,j);
            }else{
                double sinel=sin(MAX(data->azel_[1],5.0*D2R));
                Px_[j+j*nx_]+=popt_.adj_opt_.gnss_psd_[1]*tt_;
            }
        }
    }
}

#if 1
#define AMB_VAR SQR(60.0)
void Parameter::AmbParUpdate() {
    if(!NumAmbPar()) return;
    int i,f,k;
    int ia;
    int sat=0;
    vector<double>amb_bias;
    double bias=0.0;
    ObsData *data= nullptr;
    Isat_t *isat=nullptr;

    for(f=0;f<NumUseFrq();f++){
        amb_bias.assign(obss_.num_,0.0);
        k=f;
        for(i=0;i<MAXSAT;i++){
            ia=IdxAmbPar(i+1,f);
            //除非发生周跳，重置模糊度，因其他原因在当前历元没有使用该卫星，需要保持模糊度
            if(++isat_[i].outc[f]>popt_.QC_opt_.max_outage2reset_amb_||
               popt_.GNSS_opt_.ar_opt_.mode_==1){
                InitX(0.0,0.0,ia);
            }
        }

        for(i=0;i<obss_.num_;i++){
            data=&obss_.sat_infos_[i];
            sat=data->sat_.sat_no_;
            ia=IdxAmbPar(data->sat_.sat_no_,f);

            if(popt_.GNSS_opt_.use_frq_>SINGLE_FRQ&&popt_.GNSS_opt_.ion_opt_==ION_IF){
                if(popt_.GNSS_opt_.use_frq_==DUAL_FRQ){
                    // 只使用L1_L2/L1_L3 两种无电离层组合,L2_L3噪声太大不适合使用

                    int frq=f;
                    if(data->use_f2_==2&&f==0) frq=1;
                    double Pc=isat_[sat-1].IF_P[frq],Lc=isat_[sat-1].IF_L[frq];
                    if(Pc==0.0||Lc==0.0) continue;
                    else bias=Lc-Pc;
                    amb_bias[i]=isat_[sat-1].amb[frq]=bias;

                }
                if(popt_.GNSS_opt_.use_frq_==TRIPLE_FRQ){
                    if(popt_.GNSS_opt_.IF3_format_==IF3_DUAL){
                        // f==0: L1_L2  f==1 L1_L5
                        double Pc=isat_[sat-1].IF_P[f],Lc=isat_[sat-1].IF_L[f];
                        if(Pc==0.0||Lc==0.0) continue;
                        else bias=Lc-Pc;
                        amb_bias[i]=isat_[sat-1].amb[f]=bias;
                    }
                    else if(popt_.GNSS_opt_.IF3_format_==IF3_SINGLE){
                        // 先判断三频数据是否存在,如果存在三频观测，则这颗卫星只使用一个无电离层组合
                        // f=0 优先存三频无电离层组合
                        double Pc=isat_[sat-1].IF_P[3],Lc=isat_[sat-1].IF_L[3];
                        if(f==0&&Pc!=0.0&&Lc!=0.0){
                            bias=Lc-Pc;
                            amb_bias[i]=bias;
                        }
                        else if(f==1&&(Pc==0.0||Lc==0.0)){
                            Pc=isat_[sat-1].IF_P[0];Lc=isat_[sat-1].IF_L[0];
                            if(Pc==0.0||Lc==0.0) continue;
                            bias=Lc-Pc;
                            amb_bias[i]=bias;
                        }
                    }
                }
            }
            else {
                isat=&isat_[sat-1];
                double ion=0.0;
                double P1=0.0,P2=0.0,lam1=data->lam_[k],lam2=1.0;
                if(popt_.GNSS_opt_.use_frq_==SINGLE_FRQ){
                    if(popt_.GNSS_opt_.ion_opt_==ION_CONST) ion=isat->sion;
                    else if(popt_.GNSS_opt_.ion_opt_==ION_IF) ion=0.0;
                }
                else{
                    P1=isat->corr_P[data->use_f1_],lam1=data->lam_[data->use_f1_];
                    P2=isat->corr_P[data->use_f2_],lam2=data->lam_[data->use_f2_];

                    if(P1==0.0||P2==0.0) {ion=0.0;}
                    else ion=(P1-P2)/(1.0-SQR(lam2/lam1));
                    if(f==1&&data->use_f2_!=f){
                        k=data->use_f2_;
                    }
                    if(f==2&&data->use_f3_!=f){
                        k=data->use_f3_;
                    }
                    if(data->L_[k]==0.0||data->P_[k]==0.0){
                        amb_bias[i]=0.0;
                        continue;
                    }
                }

                bias=data->L_[k]*data->lam_[k]-isat->corr_P[k]+2.0*ion*SQR(data->lam_[k]/lam1);
                amb_bias[i]=isat_[sat-1].amb[k]=bias;
            }
            if(x_[ia]==0.0||(amb_bias.size()>0&&amb_bias[i]==0.0)) continue;
        }

        for(i=0;i<obss_.num_;i++){
            if(amb_bias[i]==0.0) continue;
            data=&obss_.sat_infos_[i];
            int sys=data->sat_.sat_sys_;
            int sat=data->sat_.sat_no_;
            int ia=IdxAmbPar(sat,f);

#if 0
            if(epoch_==1669&&f==1){
                InitX(amb_bias[i],AMB_VAR,ia);
                continue;
            }
#endif
            if((epoch_!=1&&isat_[sat-1].stat==NO_USE)||!last_sol_.ppp_flag_){
                double amb=x_[ia];
                if(amb==0.0) amb=amb_bias[i];
                InitX(amb,AMB_VAR,ia);
            }else{
                if(popt_.GNSS_opt_.use_frq_==TRIPLE_FRQ&&popt_.GNSS_opt_.ion_opt_==ION_UC&&f==2&&sys==GPS){
                    Px_[ia+ia*nx_]+=3e-7*fabs(tt_);
                }
                else if(popt_.GNSS_opt_.use_frq_==TRIPLE_FRQ&&popt_.GNSS_opt_.ion_opt_==ION_IF&&f==1&&sys==GPS){
                    Px_[ia+ia*nx_]+=3e-7*fabs(tt_);
                }
                else{
                    Px_[ia+ia*nx_]+=popt_.adj_opt_.gnss_psd_[0]*fabs(tt_);
                }
                if(amb_bias[i]==0.0||(x_[ia]!=0.0&&data->stat_!=EX_SLIP)) continue;
                InitX(amb_bias[i],AMB_VAR,ia);
                LOG(WARNING)<<epoch_<<" "<<data->sig_recep_.time_str_<<": "<<data->sat_.sat_id_<<" reset ambiguity x[ia]="<<x_[ia]<<" amb="<<amb_bias[i];
            }
        }

    }
}


#endif

