/**************************************************************************

Copyright(c): 2020, by Chao Chen, All rights reserved.
              China University of Mining and Technology, XuZhou, P.R. China

Author: Chao Chen(cchen@cumt.edu.cn)

Date: 2020-05-30

Description:

**************************************************************************/
#include "PPPLibGlo.h"
#include "AntModel.h"

AntModel::AntModel() = default;

AntModel::~AntModel() = default;

Ant_t* AntModel::SearchAntPara(int sat, const string& type) {
    int i,j;
    Ant_t *ant=nullptr;
    if(sat){
        for(i=0;i<ants_.size();i++){
            ant=&ants_.at(i);
            if(ant->sat_.sat_no_!=sat) continue;
            if(ant->ts.time_!=0.0&&(ant->ts.TimeDiff(ant_time_)>0.0)) continue;
            int b=ant->ts.TimeDiff(ant_time_)<0.0;
            if(ant->te.time_!=0.0&&(ant->te.TimeDiff(ant_time_)<0.0)) continue;
            return ant;
        }
    }
    else{
        int n=0;
        char *p,*types[2];
        const string& buff=type;
        for(p=strtok((char*)buff.c_str()," ");p&&n<2;p=strtok(nullptr," ")) types[n++]=p;
        if(n<=0) return nullptr;

        for(i=0;i<ants_.size();i++){
            ant=&ants_.at(i);
            for(j=0;j<n;j++) if(!strstr(ant->ant_type.c_str(),types[j])) break;
            if(j>=n) return ant;
        }

        for(i=0;i<ants_.size();i++){
            ant=&ants_.at(i);
            if(strstr(ant->ant_type.c_str(),types[0])!=ant->ant_type) continue;
            return ant;
        }
    }
    return nullptr;
}

double AntModel::InterpPcv(double ang, const double *var) {
    double a=ang/5.0;
    int i=(int)a;

    if(i<0) return var[0];
    else if(i>=18) return var[18];

    return var[i]*(1.0-a+i)+var[i+1]*(a-i);
}

double AntModel::InterpAziPcv(const Ant_t& ant,double az, double ze, int f) {
    double p,q,rpcv=0.0;
    int ize,iaz;
    int i=(int)((ant.zen2-ant.zen1)/ant.dzen)+1;

    ize=(int)((ze-ant.zen1)/(ant.dzen));
    iaz=(int)(az/ant.dazi);

    p=ze/ant.dzen-ize;
    q=az/ant.dazi-iaz;

    rpcv=(1.0-p)*(1.0-q)*ant.pcv_[f][(iaz+0)*i+(ize+0)]
          +p*(1.0-q)*ant.pcv_[f][(iaz+0)*i+(ize+1)]
          +q*(1.0-p)*ant.pcv_[f][(iaz+1)*i+(ize+0)]
          +p*q*ant.pcv_[f][(iaz+1)*i+(ize+1)];
    return rpcv;
}

void AntModel::SatPcvModel(ObsData &data, double nadir, double *dant) {
    int i,ii,sys=data.sat_.sat_sys_;
    Ant_t ant=sat_ant_[data.sat_.sat_no_-1];

    for(i=0;i<NFREQ;i++){
        switch (sys) {
            case GPS:
                ii=i;
                if(i==2) ii=1;
                break;
            case BD2:
            case BD3:
                ii=i+NFREQ;
                if(i==2) {ii=1+NFREQ;}
                break;
            case GAL:
                ii=i+2*NFREQ;
                if(i==2) {ii=1+2*NFREQ;}
                break;
            case GLO:
                ii=i+3*NFREQ;
                if(i==2) {ii=1+3*NFREQ;}
                break;
            case QZS:
                ii=i+4*NFREQ;
                if(i==2) ii=1+4*NFREQ;
                break;
        }
        dant[i]=InterpPcv(nadir*R2D*5.0,ant.pcv_[ii]);
    }
}

void AntModel::SetAntPara() {
    int i;
    SatMask sat;
    Ant_t *ant;
    int dgnss=(kPrcOpt.mode_==DGPS||kPrcOpt.mode_==PPK)||
              (kPrcOpt.mode_opt_==MDOPT_DGPS||kPrcOpt.mode_opt_==MDOPT_PPK)||kPrcOpt.mode_==FIX;
    for(i=0;i<MAXSAT;i++){
        sat=SatMask(i+1);
        if(sat.sat_sys_|kPrcOpt.GNSS_opt_.nav_sys_){
            if(!(ant=SearchAntPara(sat.sat_no_,""))){
                LOG(WARNING)<<sat.sat_id_<<" antenna correction parameter no exist";
                continue;
            }
            sat_ant_[i]=*ant;
        }
    }

    for(i=0;i<(dgnss?2:1);i++){
        if(sta_[i].del_fmt_==1){
            if(Norm(sta_->pos_,3)>0.0){
                Vec3 xyz(sta_->pos_[0],sta_->pos_[1],sta_->pos_[2]);
                Vec3 blh=Xyz2Blh(xyz,WGS84);
                double pos[3]={blh.i_,blh.j_,blh.k_};
                double del[3]={0};
                Ecef2Enu(pos,sta_[i].del_,del);
                for(int j=0;j<3;j++) ant_del[i][j]=del[j];
            }
        }
        else for(int j=0;j<3;j++) ant_del[i][j]=sta_[i].del_[j];
    }

    if(!(ant=SearchAntPara(0,sta_[0].ant_desc_))){
        LOG(WARNING)<<"Rover "<<sta_[0].ant_desc_<<" antenna correction parameter no exist";
    }
    else rec_ant_[0]=*ant;

    if(dgnss){
        if(!(ant=SearchAntPara(0,sta_[1].ant_desc_))){
            LOG(WARNING)<<"Reference "<<sta_[0].ant_desc_<<" antenna correction parameter no exist";
        }else rec_ant_[1]=*ant;
    }
}

void AntModel::SatPcoCorr(ObsData &data, double *rs, double *dant) {
   double sun_pos[3]={0},erp_val[5]={0};
   const Ant_t *sat_ant=&sat_ant_[data.sat_.sat_no_-1];

   SunMoonPos(*data.sig_trans_.GPST2UTC(),erp_val,sun_pos,NULL);

   double r[3],ez[3],es[3],ey[3],ex[3];
   for(int i=0;i<3;i++) r[i]=-rs[i];
   if(!NormV3(r,ez)) return;
   for(int i=0;i<3;i++) r[i]=sun_pos[i]-rs[i];
   if(!NormV3(r,es)) return;
   CrossVec3(ez,es,r);
   if(!NormV3(r,ey)) return;
   CrossVec3(ey,ez,ex);

   //TODO 精密星历基于无电离层组合，要细分不同分析中心产品
   int sys=data.sat_.sat_sys_;
   int i=0,j=1;
   double gamma=SQR(data.lam_[j])/SQR(data.lam_[i]);
   double C1=gamma/(gamma-1.0),C2=-1.0/(gamma-1.0);

   if(sys==GPS){i=0;j=1;}
   else if(sys==BD2||sys==BD3){
       i=0+NFREQ;j=2+NFREQ;
       if(kPrcOpt.GNSS_opt_.prods_ac_==COD){
           i=0+NFREQ;j=1+NFREQ;
       }
   }
   else if(sys==GAL){i=0+2*NFREQ;j=1+2*NFREQ;}
   else if(sys==GLO){i=0+3*NFREQ;j=1+3*NFREQ;}
   else if(sys==QZS){i=0+4*NFREQ;j=1+4*NFREQ;}

   for(int k=0;k<3;k++){
       double dant1=sat_ant->pco_[i][0]*ex[k]+sat_ant->pco_[i][1]*ey[k]+sat_ant->pco_[i][2]*ez[k];
       double dant2=sat_ant->pco_[j][0]*ex[k]+sat_ant->pco_[j][1]*ey[k]+sat_ant->pco_[j][2]*ez[k];
       dant[k]=C1*dant1+C2*dant2;
   }
}

void AntModel::SatPcvCorr(ObsData &data, double *rr, double *pcv_dant) {
    double ru[3],rz[3],eu[3],ez[3],nadir,cosa;
    int i;

    double rs[3]={data.pos_.i_,data.pos_.j_,data.pos_.k_};
    for(i=0;i<3;i++){
        ru[i]=rr[i]-rs[i];
        rz[i]=-rs[i];
    }
    if(!NormV3(ru,eu)||!NormV3(rz,ez)) return;

    cosa=Dot(eu,ez,3);
    cosa=cosa<-1.0?-1.0:(cosa>1.0?1.0:cosa);
    nadir=acos(cosa);

    SatPcvModel(data,nadir,pcv_dant);
}

void AntModel::RecAntCorr(int rov_bas,ObsData& data,double *dant) {
    double e[3],off[3],cosel=cos(data.azel_[1]);
    int i,j,k=0;

    e[0]=sin(data.azel_[0])*cosel;
    e[1]=cos(data.azel_[0])*cosel;
    e[2]=sin(data.azel_[1]);

    int sys=data.sat_.sat_sys_,ii;
    for(i=0;i<NFREQ;i++){
        if(sys==GPS||sys==BD2||sys==BD3||sys==GAL||sys==QZS){
            //使用GPS的接收机天线PCV近似替代
            ii=i;
            if(i==2) ii=1; //L5-->L2
        }
        else if(sys==GLO){
            ii=i+3*NFREQ;
            if(ii==2) ii=1+3*NFREQ;
        }
        for(j=0;j<3;j++) off[j]=rec_ant_[rov_bas].pco_[ii][j]+ant_del[rov_bas][j];
        if(rec_ant_[rov_bas].dazi!=0.0){
            double a=-Dot(off,e,3);
            dant[i]=-Dot(off,e,3)+InterpAziPcv(rec_ant_[rov_bas],data.azel_[0]*R2D,90-data.azel_[1]*R2D,ii);
        }
        else{
            dant[i]=-Dot(off,e,3)+InterpPcv(90-data.azel_[1]*R2D,rec_ant_[rov_bas].pcv_[ii]);
        }
    }
}