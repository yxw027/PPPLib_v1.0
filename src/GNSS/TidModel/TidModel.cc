/**************************************************************************

Copyright(c): 2020, by Chao Chen, All rights reserved.
              China University of Mining and Technology, XuZhou, P.R. China

Author: Chao Chen(cchen@cumt.edu.cn)

Date: 2020-05-30

Description:

**************************************************************************/

#include "PPPLibGlo.h"
#include "TidModel.h"

TidModel::TidModel() {
    tut_=Time();
    for(int i=0;i<3;i++){
        E_[i*3]=E_[i*3+1]=E_[i*3+2]=0.0;
        tid_dr_[0][i]=tid_dr_[1][i]=tid_dr_[2][i]=0.0;
        denu_[0][i]=denu_[1][i]=0.0;
        sun_pos_[i]=moon_pos_[i]=0.0;
    }
    for(int i=0;i<66;i++){
        ocean_par_[0][i]=ocean_par_[1][i]=0.0;
    }
    gmst_=0.0;
    for(double & i : erp_val_) i=0.0;

    tid_opt_=TIDE_OFF;
}

TidModel::~TidModel() {
}

int TidModel::GetErpVal(Time time) {
    const double ep[]={2000,1,1,12,0,0};
    const Time t(ep);
    double mjd,day,a;
    int i,j,k;

    if(erp_.size()<=0) return 0;

    mjd=51544.5+time.GPST2UTC()->TimeDiff(t)/86400.0;

    if(mjd<=erp_[0].mjd){
        day=mjd-erp_[0].mjd;
        erp_val_[0]=erp_[0].xp+erp_[0].xpr*day;
        erp_val_[1]=erp_[0].yp+erp_[0].ypr*day;
        erp_val_[2]=erp_[0].ut1_utc-erp_[0].lod*day;
        erp_val_[3]=erp_[0].lod;
        return 1;
    }
    int n=erp_.size()-1;
    if(mjd>=erp_[n].mjd){
        day=mjd-erp_[0].mjd;
        erp_val_[0]=erp_[n].xp+erp_[n].xpr*day;
        erp_val_[1]=erp_[n].yp+erp_[n].ypr*day;
        erp_val_[2]=erp_[n].ut1_utc-erp_[n].lod*day;
        erp_val_[3]=erp_[n].lod;
        return 1;
    }
    for(j=0,k=n;j<k-1;){
        i=(j+k)/2;
        if(mjd<erp_[i].mjd) k=i;else j=i;
    }
    if(erp_[j].mjd==erp_[j+1].mjd) a=0.5;
    else{
        a=(mjd-erp_[j].mjd)/(erp_[j+1].mjd-erp_[j].mjd);
    }
    erp_val_[0]=(1.0-a)*erp_[j].xp+a*erp_[j+1].xp;
    erp_val_[1]=(1.0-a)*erp_[j].yp+a*erp_[j+1].yp;
    erp_val_[2]=(1.0-a)*erp_[j].ut1_utc+a*erp_[j+1].ut1_utc;
    erp_val_[3]=(1.0-a)*erp_[j].lod+a*erp_[j+1].lod;
    return 1;
}

#define GME         3.986004415E+14 /* earth gravitational constant */
void TidModel::TideSl(const double *eu, const double *rp, double GMp, Vec3 blh, double *dr) {
    const double H3=0.292,L3=0.015;
    double r,ep[3],latp,lonp,p,K2,K3,a,H2,L2,dp,du,cosp,sinl,cosl;
    int i;

    if ((r=Norm(rp,3))<=0.0) return;

    for (i=0; i<3; i++) ep[i]=rp[i]/r;

    K2=GMp/GME*SQR(RE_WGS84)*SQR(RE_WGS84)/(r*r*r);
    K3=K2*RE_WGS84/r;
    latp=asin(ep[2]); lonp=atan2(ep[1],ep[0]);
    cosp=cos(latp); sinl=sin(blh.i_); cosl=cos(blh.i_);

    /* step1 in phase (degree 2) */
    p=(3.0*sinl*sinl-1.0)/2.0;
    H2=0.6078-0.0006*p;
    L2=0.0847+0.0002*p;
    a=Dot(ep,eu,3);
    dp=K2*3.0*L2*a;
    du=K2*(H2*(1.5*a*a-0.5)-3.0*L2*a*a);

    /* step1 in phase (degree 3) */
    dp+=K3*L3*(7.5*a*a-1.5);
    du+=K3*(H3*(2.5*a*a*a-1.5*a)-L3*(7.5*a*a-1.5)*a);

    /* step1 out-of-phase (only radial) */
    du+=3.0/4.0*0.0025*K2*sin(2.0*latp)*sin(2.0*blh.i_)*sin(blh.j_-lonp);
    du+=3.0/4.0*0.0022*K2*cosp*cosp*cosl*cosl*sin(2.0*(blh.j_-lonp));

    dr[0]=dp*ep[0]+du*eu[0];
    dr[1]=dp*ep[1]+du*eu[1];
    dr[2]=dp*ep[2]+du*eu[2];
}
#define GMS         1.327124E+20    /* sun gravitational constant */
#define GMM         4.902801E+12    /* moon gravitational constant */
void TidModel::TidSolid() {
    double dr1[3],dr2[3],eu[3],du,dn,sinl,sin2l;

    /* step1: time domain */
    eu[0]=E_[2]; eu[1]=E_[5]; eu[2]=E_[8];

    /* step2: frequency domain, only K1 radial */
    sin2l=sin(2.0*blh_.i_);
    du=-0.012*sin2l*sin(gmst_+blh_.j_);
    TideSl(eu,sun_pos_,GMS,blh_,dr1);
    TideSl(eu,moon_pos_,GMM,blh_,dr2);

    tid_dr_[0][0]=dr1[0]+dr2[0]+du*E_[2];
    tid_dr_[0][1]=dr1[1]+dr2[1]+du*E_[5];
    tid_dr_[0][2]=dr1[2]+dr2[2]+du*E_[8];

    /* eliminate permanent deformation */
    if (tid_opt_&8) {
        sinl=sin(blh_.i_);
        du=0.1196*(1.5*sinl*sinl-0.5);
        dn=0.0247*sin2l;
        tid_dr_[0][0]+=du*E_[2]+dn*E_[1];
        tid_dr_[0][1]+=du*E_[5]+dn*E_[4];
        tid_dr_[0][2]+=du*E_[8]+dn*E_[7];
    }

}

void TidModel::TidOcean(int rov_bas) {
    if (!ocean_par_[rov_bas]) return ;
    const double args[][5]={
            { 1.40519E-4, 2.0,-2.0, 0.0, 0.00 },  /* M2 */
            { 1.45444E-4, 0.0, 0.0, 0.0, 0.00 },  /* S2 */
            { 1.37880E-4, 2.0,-3.0, 1.0, 0.00 },  /* N2 */
            { 1.45842E-4, 2.0, 0.0, 0.0, 0.00 },  /* K2 */
            { 0.72921E-4, 1.0, 0.0, 0.0, 0.25 },  /* K1 */
            { 0.67598E-4, 1.0,-2.0, 0.0,-0.25 },  /* O1 */
            { 0.72523E-4,-1.0, 0.0, 0.0,-0.25 },  /* P1 */
            { 0.64959E-4, 1.0,-3.0, 1.0,-0.25 },  /* Q1 */
            { 0.53234E-5, 0.0, 2.0, 0.0, 0.00 },  /* Mf */
            { 0.26392E-5, 0.0, 1.0,-1.0, 0.00 },  /* Mm */
            { 0.03982E-5, 2.0, 0.0, 0.0, 0.00 }   /* Ssa */
    };
    const double ep1975[]={ 1975,1,1,0,0,0 };
    double fday,days,t,t2,t3,a[5],ang,dp[3]={ 0 };
    int i,j;

    /* angular argument: see subroutine arg.f for reference [1] */
    Time time=tut_,t1975(ep1975);
    time.Time2Epoch();
    fday=time.epoch_[3]*3600.0+time.epoch_[4]*60.0+time.epoch_[5];
    time.epoch_[3]=time.epoch_[4]=time.epoch_[5]=0.0;
    days=(time.Epoch2Time(time.epoch_)->TimeDiff(t1975))/86400.0+1.0;
    t=(27392.500528+1.000000035*days)/36525.0;
    t2=t*t; t3=t2*t;

    a[0]=fday;
    a[1]=(279.69668+36000.768930485*t+3.03E-4*t2)*D2R; /* H0 */
    a[2]=(270.434358+481267.88314137*t-0.001133*t2+1.9E-6*t3)*D2R; /* S0 */
    a[3]=(334.329653+4069.0340329577*t-0.010325*t2-1.2E-5*t3)*D2R; /* P0 */
    a[4]=2.0*PI;

    /* displacements by 11 constituents */
    for (i=0; i<11; i++) {
        ang=0.0;
        for (j=0; j<5; j++) ang+=a[j]*args[i][j];
        for (j=0; j<3; j++) dp[j]+=ocean_par_[rov_bas][j+i*6]*cos(ang-ocean_par_[rov_bas][j+3+i*6]*D2R);
    }
    denu_[0][0]=-dp[1];
    denu_[0][1]=-dp[2];
    denu_[0][2]= dp[0];
}

void TidModel::IersMeanPole(double *xp, double *yp) {
    const double ep2000[]={ 2000,1,1,0,0,0 };
    double y,y2,y3;

    Time t2000(ep2000);
    y=tut_.TimeDiff(t2000)/86400.0/365.25;

    if (y<3653.0/365.25) { /* until 2010.0 */
        y2=y*y; y3=y2*y;
        *xp= 55.974+1.8243*y+0.18413*y2+0.007024*y3; /* (mas) */
        *yp=346.346+1.7896*y-0.10729*y2-0.000908*y3;
    }
    else { /* after 2010.0 */
        *xp= 23.513+7.6141*y; /* (mas) */
        *yp=358.891-0.6287*y;
    }
}

void TidModel::TidPole() {
    double xp_bar,yp_bar,m1,m2,cosl,sinl;

    /* iers mean pole (mas) */
    IersMeanPole(&xp_bar,&yp_bar);

    /* ref [7] eq.7.24 */
    m1= erp_val_[0]/AS2R-xp_bar*1E-3; /* (as) */
    m2=-erp_val_[1]/AS2R+yp_bar*1E-3;

    /* sin(2*theta) = sin(2*phi), cos(2*theta)=-cos(2*phi) */
    cosl=cos(blh_.j_);
    sinl=sin(blh_.j_);
    denu_[1][0]=  9E-3*cos(blh_.i_)    *(m1*sinl-m2*cosl); /* de= Slambda (m) */
    denu_[1][1]= -9E-3*cos(2.0*blh_.i_)*(m1*cosl+m2*sinl); /* dn=-Stheta  (m) */
    denu_[1][2]=-32E-3*sin(2.0*blh_.i_)*(m1*cosl+m2*sinl); /* du= Sr      (m) */
}

void TidModel::TidCorr(Time time, int rov_bas, double *xyz, double *dr) {
    tid_opt_=kPrcOpt.GNSS_opt_.tid_opt_;
    for(int i=0;i<3;i++) dr[i]=0.0;
    if(erp_.size()>0) GetErpVal(time);
    tut_=time;tut_.TimeAdd(erp_val_[2]);
    if((re_=Norm(xyz,3))<=0.0) return;

    blh_.i_=asin(xyz[2]/re_);
    blh_.j_=atan2(xyz[1],xyz[0]);
    blh_.k_=0.0;
    double pos[3]={0};
    pos[0]=blh_.i_;pos[1]=blh_.j_;
    Xyz2Enu(pos,E_);

    if(tid_opt_&TIDE_SOLID){
        gmst_=SunMoonPos(time,erp_val_,sun_pos_,moon_pos_);
        TidSolid();
        for(int i=0;i<3;i++) dr[i]+=tid_dr_[0][i];
    }
    if(tid_opt_&TIDE_OCEAN){
        TidOcean(rov_bas);
        MatMulPnt("TN",3,1,3,1.0,E_,denu_[0],0.0,tid_dr_[1]);
        for (int i=0; i<3; i++) dr[i]+=tid_dr_[1][i];
    }
    if(tid_opt_&TIDE_POLE){
        TidPole();
        MatMulPnt("TN",3,1,3,1.0,E_,denu_[1],0.0,tid_dr_[2]);
        for (int i=0; i<3; i++) dr[i]+=tid_dr_[2][i];
    }
}