/**************************************************************************

Copyright(c): 2020, by Chao Chen, All rights reserved.
              China University of Mining and Technology, XuZhou, P.R. China

Author: Chao Chen(cchen@cumt.edu.cn)

Date: 2020-05-30

Description:

**************************************************************************/

#include "AttModel.h"

AttModel::AttModel() {
    phw_=0.0;
    for(int i=0;i<3;i++) exs_[i]=eys_[i]=0.0;
}

AttModel::~AttModel() = default;

int AttModel::SatYaw(ObsData obs) {
    double sun_pos[3]={0},erp[5]={0},n[3]={0},p[3]={0};
    double es[3],en[3],ep[3],ex[3],esun[3];
    SunMoonPos(*obs.sig_recep_.GPST2UTC(),erp,sun_pos,NULL);

    double sat_pos[3]={0},sat_vel[3]={0};
    sat_pos[0]=obs.pos_.i_;sat_pos[1]=obs.pos_.j_;sat_pos[2]=obs.pos_.k_;
    sat_vel[0]=obs.vel_.i_;sat_vel[1]=obs.vel_.j_;sat_vel[2]=obs.vel_.k_;
    sat_vel[0]-=OMGE_GPS*sat_vel[1];sat_vel[1]+=OMGE_GPS*sat_vel[0];
    CrossVec3(sat_pos,sat_vel,n);
    CrossVec3(sun_pos,n,p);
    if(!NormV3(sat_pos,es)||!NormV3(sun_pos,esun)||!NormV3(n,en)||!NormV3(p,ep)) return 0;
    double beta=PI/2.0-acos(Dot(esun,en,3));
    double E=acos(Dot(es,ep,3));
    double mu=PI/2.0+(Dot(es,esun,3)<=0?-E:E);
    if(mu<-PI/2.0) mu+=2.0*PI;
    else if(mu>=PI/2.0) mu-=2.0*PI;

    double yaw=0.0;
    if(fabs(beta)<1E-12&&fabs(mu)<1E-12) yaw=PI;
    else yaw=atan2(-tan(beta),sin(mu))+PI;

    CrossVec3(en,es,ex);
    double cosy=cos(yaw);
    double siny=sin(yaw);
    for(int i=0;i<3;i++){
        exs_[i]=-siny*en[i]+cosy*ex[i];
        eys_[i]=-cosy*en[i]-siny*ex[i];
    }
    return 1;
}

double AttModel::SatPhw(ObsData &obs,double *rec_pos) {

    SatYaw(obs);

    double r[3]={0};
    r[0]=rec_pos[0]-obs.pos_.i_;
    r[1]=rec_pos[1]-obs.pos_.j_;
    r[2]=rec_pos[2]-obs.pos_.k_;
    double ek[3]={0};
    if(!NormV3(r,ek)) return 0;

    Vec3 rr(rec_pos[0],rec_pos[1],rec_pos[2]);
    Vec3 blh=Xyz2Blh(rr,WGS84);
    double E[9];
    double pos[3]={blh.i_,blh.j_,blh.k_};
    Xyz2Enu(pos,E);
    double exr[3],eyr[3];
    exr[0]=E[1];exr[1]=E[4];exr[2]=E[7];
    eyr[0]=-E[0];eyr[1]=-E[3];eyr[2]=-E[6];

    double eks[3]={0},ekr[3]={0};
    CrossVec3(ek,eys_,eks);
    CrossVec3(ek,eyr,ekr);
    double ds[3]={0},dr[3]={0};
    for(int i=0;i<3;i++){
        ds[i]=exs_[i]-ek[i]*Dot(ek,exs_,3)-eks[i];
        dr[i]=exr[i]-ek[i]*Dot(ek,exr,3)+ekr[i];
    }
    double cosp=Dot(ds,dr,3)/Norm(ds,3)/Norm(dr,3);
    if(cosp<-1.0) cosp=-1.0;
    else if(cosp>1.0) cosp=1.0;
    double ph=acos(cosp)/2.0/PI;
    double drs[3]={0};
    CrossVec3(ds,dr,drs);
    if(Dot(ek,drs,3)<0.0) ph=-ph;
    phw_=ph+floor(phw_-ph+0.5);
    return phw_;
}