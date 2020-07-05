/**************************************************************************

Copyright(c): 2020, by Chao Chen, All rights reserved.
              China University of Mining and Technology, XuZhou, P.R. China

Author: Chao Chen(cchen@cumt.edu.cn)

Date: 2020-05-30

Description:

**************************************************************************/

#include "PPPLibGlo.h"
#include "BiaModel.h"

BiaModel::BiaModel() {
    type_=kPrcOpt.GNSS_opt_.cb_prd_type_;
    ac_=kPrcOpt.GNSS_opt_.cb_prd_ac_;
    cbias_=0.0;
    sat_data_= nullptr;
    for(int i=0;i<MAXSAT;i++){
        for(int j=0;j<MAXCBIASPAIR;j++){
            dcb_[i][j]=0.0;
        }
    }

}


BiaModel::~BiaModel() {

}

double BiaModel::BDMultiPathCorr() {
    const static double kBDIGSOCoef[3][10]={
            {-0.55,-0.40,-0.34,-0.23,-0.15,-0.04,0.09,0.19,0.27,0.35},	//B1
            {-0.71,-0.36,-0.33,-0.19,-0.14,-0.03,0.08,0.17,0.24,0.33},	//B2
            {-0.27,-0.23,-0.21,-0.15,-0.11,-0.04,0.05,0.14,0.19,0.32},	//B3
    };
    const static double kBDMEOCoef[3][10]={
            {-0.47,-0.38,-0.32,-0.23,-0.11,0.06,0.34,0.69,0.97,1.05},	//B1
            {-0.40,-0.31,-0.26,-0.18,-0.06,0.09,0.28,0.48,0.64,0.69},	//B2
            {-0.22,-0.15,-0.13,-0.10,-0.04,0.05,0.14,0.27,0.36,0.47},	//B3
    };

    double el=sat_data_->azel_[1]*R2D*0.1;
    int int_el=(int)el;
    int prn=sat_data_->sat_.sat_prn_;
    int BD2_IGSO,BD2_MEO;
    double mp[3]={0};
    cbias_=0.0;

    BD2_IGSO=std::binary_search(kBDS2_IGSO,kBDS2_IGSO+NUM_BDS2_IGSO,prn);
    BD2_MEO=std::binary_search(kBDS2_MEO,kBDS2_MEO+NUM_BDS2_MEO,prn);
    if(BD2_IGSO){
        if(el<0) for(int i=0;i<3;i++) mp[i]=kBDIGSOCoef[i][0];
        else if(el>=9) for(int i=0;i<3;i++) mp[i]=kBDIGSOCoef[i][9];
        else for(int i=0;i<3;i++) mp[i]=kBDIGSOCoef[i][int_el]*(1.0-el+int_el)+kBDIGSOCoef[i][int_el+1]*(el-int_el);
    }
    else if(BD2_MEO){
        if(el<0) for(int i=0;i<3;i++) mp[i]=kBDMEOCoef[i][0];
        else if(el>=9) for(int i=0;i<3;i++) mp[i]=kBDMEOCoef[i][9];
        else for(int i=0;i<3;i++) mp[i]=kBDMEOCoef[i][int_el]*(1.0-el+int_el)+kBDMEOCoef[i][int_el+1]*(el-int_el);
    }
    else return 0.0;
    for(int i=0;i<3;i++) if(sat_data_->P_[i]!=0.0) cbias_-=mp[i];
    return cbias_;
}

int BiaModel::MatchOsb2Sat() {
    int i,j,sys=sat_data_->sat_.sat_sys_;
    int sat=sat_data_->sat_.sat_no_;

    if(!sat_bias_.dt) return 0;
    if(sat_data_->sig_recep_.TimeDiff(sat_bias_.tmin)<0.0) return 0;
    if(sat_data_->sig_recep_.TimeDiff(sat_bias_.tmax)>0.0) return 0;
    i=(int)(sat_data_->sig_recep_.TimeDiff(sat_bias_.tmin)/sat_bias_.dt);
    for(j=0;j<NFREQ;j++){
        sat_data_->code_osb_[j]=sat_bias_.bias[i].code[sat-1][sat_data_->code_[j]];
    }

    if(sys==GPS){
        sat_data_->datum_brdc=sat_bias_.bias[i].code[sat-1][CODE_L1W]-sat_bias_.bias->code[sat-1][CODE_L2W];
    }
    else if(sys==BD2||sys==BD3){
        sat_data_->datum_brdc=sat_bias_.bias[i].code[sat-1][CODE_L6I];
    }
    else if(sys==GAL){
        sat_data_->datum_brdc=sat_bias_.bias[i].code[sat-1][CODE_L1X]-sat_bias_.bias->code[sat-1][CODE_L5X];
    }
    else if(sys==GLO){
        sat_data_->datum_brdc=sat_bias_.bias[i].code[sat-1][CODE_L1P]-sat_bias_.bias->code[sat-1][CODE_L2P];
    }
    else if(sys==QZS){
        sat_data_->datum_brdc=sat_bias_.bias[i].code[sat-1][CODE_L1C]-sat_bias_.bias->code[sat-1][CODE_L2L];
    }
}

int BiaModel::AlignOsb2SatBias(SatOsb_t sat_osb) {
    int i,j,num_epoch,epoch_s=0,epoch_e=0;
    double dt=0.0;
    Time t_min,t_max;
    SatMask sat;

    for(i=0;i<sat_osb.nb;i++){
        if(i==0){
            t_min=sat_osb.osbs[0].ts;
            t_max=sat_osb.osbs[0].te;
            dt=t_max.TimeDiff(t_min);
        }
        if(sat_osb.osbs[i].ts.TimeDiff(t_min)<0.0) t_min=sat_osb.osbs[i].ts;
        if(sat_osb.osbs[i].te.TimeDiff(t_max)>0.0) t_max=sat_osb.osbs[i].te;
        if(sat_osb.osbs[i].te.TimeDiff(sat_osb.osbs[i].ts)<dt) dt=sat_osb.osbs[i].te.TimeDiff(sat_osb.osbs[i].ts);
    }

    if(!dt){
        sat_bias_.dt=0.0;
        sat_bias_.bias=NULL;
        return 0;
    }
    sat_bias_.dt=dt;
    sat_bias_.tmin=t_min;
    sat_bias_.tmax=t_max;
    num_epoch=t_max.TimeDiff(t_min)/dt;
    sat_bias_.bias=(Bias_t *)calloc(num_epoch,sizeof(Bias_t));
    if(sat_bias_.bias==NULL){
        cout<<"satellite bias malloc error"<<endl;
        return 0;
    }
    int code_idx;
    for(i=0;i<sat_osb.nb;i++){
        sat=sat_osb.osbs[i].sat;
        code_idx=sat_osb.osbs[i].code;
        epoch_s=(int)(sat_osb.osbs[i].ts.TimeDiff(t_min)/dt);
        epoch_e=(int)(sat_osb.osbs[i].te.TimeDiff(t_min)/dt);
        for(j=epoch_s;j<epoch_e;j++){
            if(sat_osb.osbs[i].type==PHASE_BIAS){
                sat_bias_.bias[j].phase[sat.sat_no_-1][code_idx]=sat_osb.osbs[i].val*CLIGHT*1E-9;
            }
            else{
                sat_bias_.bias[j].code[sat.sat_no_-1][code_idx]=sat_osb.osbs[i].val*CLIGHT*1E-9;
            }
        }
    }
    return 1;
}

double BiaModel::GetTGD() {
//    int sys=sat_data_->sat_.sat_sys_;
//    int prn=sat_data_->sat_.sat_prn_;
//    tgd_=tgd1_=tgd2_=0.0;
//    if(sys==GPS){
//        return tgd_=CLIGHT*eph_[sat_data_->eph_idx_].tgd_[0];  //d_p1-d_p2
//    }
//    else if(sys==BD2){
//        tgd1_=CLIGHT*eph_[sat_data_->eph_idx_].tgd_[0]; //d_b1-d_b3;
//        tgd2_=CLIGHT*eph_[sat_data_->eph_idx_].tgd_[2]; //d_b2-d_b3;
//    }
//    else if(sys==BD3){
//        tgd1_=CLIGHT*eph_[sat_data_->eph_idx_].tgd_[0];
//    }
//    else if(sys==GAL){
//
//    }
//    else if(sys==GLO){
//
//    }
//    else if(sys==QZS){
//
//    }
}

double BiaModel::DCBCorr(PrcOpt popt,int f) {
    // TODO GPS P1C1没有修正
    int sys=sat_data_->sat_.sat_sys_,sat=sat_data_->sat_.sat_no_;
    int prn=sat_data_->sat_.sat_prn_;
    int bd3_flag=0;
    double *frq=sat_data_->frq_;
    double alpha=0.0,beta=0.0;
    double P1_P2=0.0,P1_P3=0.0,P2_P3=0.0;
    cbias_=0.0;

    if(sys==BD2&&prn>=19) bd3_flag=1;

    if(sys==GPS){
        P1_P2=dcb_[sat-1][GPS_C1WC2W];
        alpha=SQR(frq[GPS_L1])/(SQR(frq[GPS_L1])-SQR(frq[GPS_L2]));
        beta=-SQR(frq[GPS_L2])/(SQR(frq[GPS_L1])-SQR(frq[GPS_L2]));
        if(f==GPS_L1){
           cbias_=beta*P1_P2;
        }
        else if(f==GPS_L2){
            cbias_=-alpha*P1_P2;
        }
        else if(f==GPS_L5){
            if(sat_data_->code_[f]==CODE_L5Q) P1_P3=dcb_[sat-1][GPS_C1CC5Q]-dcb_[sat-1][GPS_C1CC1W];
            else if(sat_data_->code_[f]==CODE_L5X) P1_P3=dcb_[sat-1][GPS_C1CC5X]-dcb_[sat-1][GPS_C1CC1W];
            cbias_=beta*P1_P2-P1_P3;
        }
        else cbias_=0.0;
        if(popt.GNSS_opt_.use_frq_==DUAL_FRQ&&popt.GNSS_opt_.ion_opt_==ION_IF){
            cbias_=0.0;
        }
        if(sat_data_->code_[f]==CODE_L1C) cbias_+=dcb_[sat-1][GPS_C1CC1W];
        else if(sat_data_->code_[f]==CODE_L2C) cbias_+=dcb_[sat-1][GPS_C2CC2W];
        else if(sat_data_->code_[f]==CODE_L2S) cbias_-=dcb_[sat-1][GPS_C2WC2S];
        else if(sat_data_->code_[f]==CODE_L2L) cbias_-=dcb_[sat-1][GPS_C2WC2L];
        else if(sat_data_->code_[f]==CODE_L2X) cbias_-=dcb_[sat-1][GPS_C2WC2X];
        return cbias_;
    }
    else if(sys==BD2){
        if(popt.GNSS_opt_.sat_eph_==BRDC){
            if(bd3_flag){
                P1_P3=dcb_[sat-1][BD3_C2IC6I];
                if(f==BD3_B1I) cbias_=P1_P3;
                else if(f==BD3_B2a){
                    if(sat_data_->code_[f]==CODE_L5X) P2_P3=dcb_[sat-1][BD3_C1XC6I]-dcb_[sat-1][BD3_C1XC5X];
                    else if(sat_data_->code_[f]==CODE_L5P) P2_P3=dcb_[sat-1][BD3_C1PC6I]-dcb_[sat-1][BD3_C1PC5P];
                    else if(sat_data_->code_[f]==CODE_L5D) P2_P3=dcb_[sat-1][BD3_C1DC6I]-dcb_[sat-1][BD3_C1DC5D];
                    cbias_=P2_P3;
                }
                else if(f==BD3_B3I) cbias_=0.0;
                else if(f==BD3_B1C){
                    if(sat_data_->code_[f]==CODE_L1X) P1_P3=dcb_[sat-1][BD3_C1XC6I];
                    else if(sat_data_->code_[f]==CODE_L1P) P1_P3=dcb_[sat-1][BD3_C1PC6I];
                    else if(sat_data_->code_[f]==CODE_L1D) P1_P3=dcb_[sat-1][BD3_C1DC6I];
                    cbias_=P1_P3;
                }
                else cbias_=0.0;
                return cbias_;
            }
            else{
                P1_P3=dcb_[sat-1][BD2_C2IC6I];
                P1_P2=dcb_[sat-1][BD2_C2IC7I];
                if(f==BD2_B1I) cbias_=P1_P3;
                else if(f==BD2_B2I) cbias_=P1_P3-P1_P2;
                else cbias_=0.0;
                return cbias_;
            }
        }
        else{
            if(bd3_flag){
                P1_P3=dcb_[sat-1][BD3_C2IC6I];
                if(sat_data_->code_[f]==CODE_L5X) P2_P3=dcb_[sat-1][BD3_C1XC6I]-dcb_[sat-1][BD3_C1XC5X];
                else if(sat_data_->code_[f]==CODE_L5P) P2_P3=dcb_[sat-1][BD3_C1PC6I]-dcb_[sat-1][BD3_C1PC5P];
                else if(sat_data_->code_[f]==CODE_L5D) P2_P3=dcb_[sat-1][BD3_C1DC6I]-dcb_[sat-1][BD3_C1DC5D];
                P1_P2=P1_P3-P2_P3;
            }
            else{
                P1_P3=dcb_[sat-1][BD2_C2IC6I];
                P1_P2=dcb_[sat-1][BD2_C2IC7I];
            }
            if(popt.GNSS_opt_.prods_ac_==COD){
                // CODE的BD2是基于B1/B2频点
                alpha=SQR(frq[BD2_B1I])/(SQR(frq[BD2_B1I])-SQR(frq[BD2_B2I]));
                beta=-SQR(frq[BD2_B2I])/(SQR(frq[BD2_B1I])-SQR(frq[BD2_B2I]));
                if(f==BD2_B1I) cbias_=beta*P1_P2;
                else if(f==BD2_B2I) cbias_=-alpha*P1_P2;
                else if(f==BD2_B3I){
                    cbias_=beta*P1_P2-P1_P3;
                }
                return cbias_;
            }
            else if(popt.GNSS_opt_.prods_ac_==WUM||popt.GNSS_opt_.prods_ac_==GBM){
                // WUM的BD2是基于B1/B3频点
                if(bd3_flag){
                    alpha=SQR(frq[BD3_B1I])/(SQR(frq[BD3_B1I])-SQR(frq[BD3_B3I]));
                    beta=-SQR(frq[BD3_B3I])/(SQR(frq[BD3_B1I])-SQR(frq[BD3_B3I]));
                    if(f==BD3_B1I) cbias_=beta*P1_P3;
                    else if(f==BD3_B2a) cbias_=beta*P1_P3-P1_P2;
                    else if(f==BD3_B3I) cbias_=-alpha*P1_P3;
                }
                else{
                    alpha=SQR(frq[BD2_B1I])/(SQR(frq[BD2_B1I])-SQR(frq[BD2_B3I]));
                    beta=-SQR(frq[BD2_B3I])/(SQR(frq[BD2_B1I])-SQR(frq[BD2_B3I]));
                    if(f==BD2_B1I) cbias_=beta*P1_P3;
                    else if(f==BD2_B2I) cbias_=beta*P1_P3-P1_P2;
                    else if(f==BD2_B3I) cbias_=-alpha*P1_P3;
                }
            }
            return cbias_;
        }
    }
    else if(sys==BD3){
        P1_P3=dcb_[sat-1][BD3_C2IC6I];
        if(sat_data_->code_[f]==CODE_L5X) P2_P3=dcb_[sat-1][BD3_C1XC6I]-dcb_[sat-1][BD3_C1XC5X];
        else if(sat_data_->code_[f]==CODE_L5P) P2_P3=dcb_[sat-1][BD3_C1PC6I]-dcb_[sat-1][BD3_C1PC5P];
        else if(sat_data_->code_[f]==CODE_L5D) P2_P3=dcb_[sat-1][BD3_C1DC6I]-dcb_[sat-1][BD3_C1DC5D];
        P1_P2=P1_P3-P2_P3;
        if(popt.GNSS_opt_.sat_eph_==BRDC){
            // 广播星历是基于B3频点
            P1_P3=dcb_[sat-1][BD3_C2IC6I];
            if(f==BD3_B1I) cbias_=P1_P3;
            else if(f==BD3_B2a){
                if(sat_data_->code_[f]==CODE_L5X) P2_P3=dcb_[sat-1][BD3_C1XC6I]-dcb_[sat-1][BD3_C1XC5X];
                else if(sat_data_->code_[f]==CODE_L5P) P2_P3=dcb_[sat-1][BD3_C1PC6I]-dcb_[sat-1][BD3_C1PC5P];
                else if(sat_data_->code_[f]==CODE_L5D) P2_P3=dcb_[sat-1][BD3_C1DC6I]-dcb_[sat-1][BD3_C1DC5D];
                cbias_=P2_P3;
            }
            else if(f==BD3_B3I) cbias_=0.0;
            else if(f==BD3_B1C){
                if(sat_data_->code_[f]==CODE_L1X) P1_P3=dcb_[sat-1][BD3_C1XC6I];
                else if(sat_data_->code_[f]==CODE_L1P) P1_P3=dcb_[sat-1][BD3_C1PC6I];
                else if(sat_data_->code_[f]==CODE_L1D) P1_P3=dcb_[sat-1][BD3_C1DC6I];
                cbias_=P1_P3;
            }
            else cbias_=0.0;
            return cbias_;
        }
        else{
            // GFZ/WUM使用的是B1与B3频点生成BDS的钟差
            alpha=SQR(frq[BD3_B1I])/(SQR(frq[BD3_B1I])-SQR(frq[BD3_B3I]));
            beta=-SQR(frq[BD3_B3I])/(SQR(frq[BD3_B1I])-SQR(frq[BD3_B3I]));
            if(f==BD3_B1I) cbias_=beta*P1_P3;
            else if(f==BD3_B2a) cbias_=beta*P1_P3-P1_P2;
            else if(f==BD3_B3I) cbias_=-alpha*P1_P3;
            return cbias_;
        }
    }
    else if(sys==GAL){
        alpha=SQR(frq[GAL_E1])/(SQR(frq[GAL_E1])-SQR(frq[GAL_E5a]));
        beta=-SQR(frq[GAL_E5a])/(SQR(frq[GAL_E1])-SQR(frq[GAL_E5a]));
        if(sat_data_->code_[f]==CODE_L1C) P1_P2=dcb_[sat-1][GAL_C1CC5Q];
        else if(sat_data_->code_[f]==CODE_L1X) P1_P2=dcb_[sat-1][GAL_C1XC5X];
        if(f==GAL_E1) cbias_=beta*P1_P2;
        else if(f==GAL_E5a) cbias_=-alpha*P1_P2;
        else if(f==GAL_E5b){
            if(sat_data_->code_[f]==CODE_L7X) P1_P3=dcb_[sat-1][GAL_C1XC7X];
            else if(sat_data_->code_[f]==CODE_L7Q) P1_P3=dcb_[sat-1][GAL_C1CC7Q];
            cbias_=beta*P1_P2-P1_P3;
        }
        else cbias_=0.0;
        return cbias_;
    }
    else if(sys==GLO){
        alpha=SQR(frq[GLO_G1])/(SQR(frq[GLO_G1])-SQR(frq[GLO_G2]));
        beta=-SQR(frq[GLO_G2])/(SQR(frq[GLO_G1])-SQR(frq[GLO_G2]));
        P1_P2=dcb_[sat-1][GLO_C1PC2P];
        if(f==GLO_G1){
            cbias_=beta*P1_P2;
        }
        else if(f==GLO_G2){
            cbias_=-alpha*P1_P2;
        }
        else cbias_=0.0;
        if(sat_data_->code_[f]==CODE_L1C) cbias_+=dcb_[sat-1][GLO_C1CC1P];
        else if(sat_data_->code_[f]==CODE_L2C) cbias_+=dcb_[sat-1][GLO_C2CC2P];
    }
    else if(sys==QZS){
        alpha=SQR(frq[QZS_L1])/(SQR(frq[QZS_L1])-SQR(frq[QZS_L2]));
        beta=-SQR(frq[QZS_L2])/(SQR(frq[QZS_L1])-SQR(frq[QZS_L2]));
        P1_P2=dcb_[sat-1][QZS_C1CC2L];
        if(f==QZS_L1){
            cbias_=beta*P1_P2;
        }
        else if(f==QZS_L2){
            cbias_=-alpha*P1_P2;
        }
        else if(f==QZS_L5){
            if(sat_data_->code_[sat_data_->use_f1_]==CODE_L5Q) P1_P3=dcb_[sat-1][QZS_C1CC5Q];
            else if(sat_data_->code_[sat_data_->use_f1_]==CODE_L5X) P1_P3=dcb_[sat-1][QZS_C1CC5X];
            cbias_=beta*P1_P2-P1_P3;
        }
        else cbias_=0.0;
        if(sat_data_->code_[f]==CODE_L1X) cbias_-=dcb_[sat-1][QZS_C1CC1X];
        else if(sat_data_->code_[f]==CODE_L2X) cbias_-=(dcb_[sat-1][QZS_C1CC1X]+dcb_[sat-1][QZS_C1XC2X]);
        return cbias_;
    }
}

int BiaModel::OSBCorr() {
    cbias_=0.0;
    int sys=sat_data_->sat_.sat_sys_;
    double f1=0.0,f2=0.0,alpha=0.0,beta=0.0;
}

void BiaModel::InitBiaCorr(ObsData &obs) {
    sat_data_=&obs;
    type_=kPrcOpt.GNSS_opt_.cb_prd_type_;
}

double BiaModel::CBiaCorr(PrcOpt popt,int f,ObsData &obs) {
    InitBiaCorr(obs);
    if(type_==CBC_DCB||type_==CBC_TGD){
        return DCBCorr(popt,f);
    }
    else if(type_==CBC_OSB){
        return OSBCorr();
    }
    else return 0.0;
}
