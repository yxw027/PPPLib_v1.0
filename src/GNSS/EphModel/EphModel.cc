/**************************************************************************

Copyright(c): 2020, by Chao Chen, All rights reserved.
              China University of Mining and Technology, XuZhou, P.R. China

Author: Chao Chen(cchen@cumt.edu.cn)

Date: 2020-05-30

Description:

**************************************************************************/

#include "PPPLibGlo.h"
#include "EphModel.h"

SatEph::SatEph() {
    eph_.clear();glo_eph_.clear();
    peph_.clear();pclk_.clear();
}

SatEph::~SatEph() {

}

int SatEph::SatClk(ObsData& data) {
    return 0;
}

void SatEph::InitSatEph(const Nav &nav) {
    ne_=nav.n_;ng_=nav.ng_;
    eph_=nav.eph_;
    glo_eph_=nav.glo_eph_;
    np_=nav.np_;nc_=nav.nc_;
    peph_=nav.pre_eph_;
    pclk_=nav.pre_clk_;
}

int SatEph::SatPos(ObsData& data,int iode) {
    return 0;
}

int SatEph::BrdcPosClk(Obss &obs) {
    return 1;
}

int SatEph::PrecPosClk(Obss &obs) {
    return 1;
}

int SatEph::SatPosClk(Obss &obs) {
    int num=0;
    for(int i=0;i<obs.num_&&i<MAXOBS;i++){

        obs.sat_infos_[i].ReSetSat();

        // approximate signal translate time
        if(!obs.sat_infos_[i].SigTransTime()) continue;

        if(!SatClk(obs.sat_infos_[i])){
            cout<<"no satellite clock correction,skipped"<<endl;
            continue;
        }
        // signal translate time with clock correction
        obs.sat_infos_[i].ClkCorrection();

        if(!SatPos(obs.sat_infos_[i],999)){
            cout<<"no ephemeris,skipped"<<endl;
            continue;
        }

        if(obs.sat_infos_[i].svh_!=-1) num++;

        cout<<std::fixed<<setprecision(3)<<obs.sat_infos_[i].sig_recep_.time_str_<<" "<<obs.sat_infos_[i].sat_.sat_id_<<" ";
        cout<<obs.sat_infos_[i].pos_.i_<<" "<<obs.sat_infos_[i].pos_.j_<<" "<<obs.sat_infos_[i].pos_.k_<<" ";
        cout<<obs.sat_infos_[i].clk_err_[0]*1E9<<" "<<obs.sat_infos_[i].sat_var_<<endl;
    }

    return num;
}

BrdcEph::BrdcEph() {

}

BrdcEph::~BrdcEph() {

}

static const double kNormalUra[]={     /* ura nominal values */
        2.0,2.8,4.0,5.7,8.0,11.3,16.0,32.0,64.0,128.0,256.0,512.0,1024.0,
        2048.0,4096.0,8192.0
};
double BrdcEph::UraEph(int ura) {

    const double ura_value[]={
            2.4,3.4,4.85,6.85,9.65,13.65,24.0,48.0,96.0,192.0,384.0,768.0,1536.0,
            3072.0,6144.0
    };
    return ura<0||15<ura?SQR(6144.0):SQR(ura_value[ura]);
#if 0
    return ura<0||15<ura ? SQR(8192.0) : SQR(kNormalUra[ura]);
#endif
}

int BrdcEph::SelectEph(ObsData &data, int iode) {
    double t,tmax,tmin;
    int neph=-1;
    int sys=data.sat_.sat_sys_;
    switch (sys){
        case QZS: tmax=MAXDTOE_QZS+1.0; break;
        case GAL: tmax=MAXDTOE_GAL+1.0; break;
        case BD2: tmax=MAXDTOE_BDS+1.0; break;
        case BD3: tmax=MAXDTOE_BDS+1.0; break;
        default:  tmax=MAXDTOE_GPS+1.0; break;
    }
    tmin=tmax+1.0;

    for (int i=0; i<ne_; i++) {
        if (eph_[i].sat_!=data.sat_.sat_no_) continue;
        if (iode>=0&&eph_[i].iode_!=iode) continue;
        if ((t=fabs(eph_[i].toe_.TimeDiff(data.sig_recep_)))>tmax) continue;
        if (iode>=0) return i;
        if (t<=tmin) { neph=i; tmin=t; } /* toe closest to time */
    }
    if (iode>=0||neph<0) {
        return -1;
    }
    return neph;
}

int BrdcEph::SelGloEph(ObsData &data, int iode) {
    double t,tmax=MAXDTOE_GLO,tmin=tmax+1.0;
    int ngeph=-1;
    char buff[MAXBUFF]={'\0'};

    for (int i=0; i<ng_; i++) {
        if (glo_eph_[i].sat_!=data.sat_.sat_no_) continue;
        if (iode>=0&&glo_eph_[i].iode_!=iode) continue;
        if ((t=fabs(glo_eph_[i].toe_.TimeDiff(data.sig_recep_)))>tmax) continue;
        if (iode>=0) return i;
        if (t<=tmin) { ngeph=i; tmin=t; } /* toe closest to time */
    }
    if (iode>=0||ngeph<0) {
        sprintf(buff,"%4d %20s: %4s no match glo broadcast ephemeris",data.epoch_,data.sig_recep_.time_str_.c_str(),data.sat_.sat_id_.c_str());
        LOG(WARNING)<<buff;
        return -1;
    }
    return ngeph;
}

void BrdcEph::Eph2Clk(ObsData &data, const Eph eph) {
    double t;

    t=data.sig_trans_.TimeDiff(eph.toc_);
    for(int i=0;i<2;i++){
        t-=eph.f0_+eph.f1_*t+eph.f2_*t*t;
    }
    data.clk_err_[0]=eph.f0_+eph.f1_*t+eph.f2_*t*t;
}

void BrdcEph::GloEph2Clk(ObsData &data, const GloEph glo_eph) {
    double t;

    t=data.sig_trans_.TimeDiff(glo_eph.toe_);

    for (int i=0; i<2; i++) {
        t-=-glo_eph.taun_+glo_eph.gamn_*t;
    }
    data.clk_err_[0]=-glo_eph.taun_+glo_eph.gamn_*t;
}

int BrdcEph::SatBrdcClk(ObsData &data) {
    SatMask sat=data.sat_;
    int sys=sat.sat_sys_,idx_eph=0,idx_geph=0;
    if(sys==GPS||sys==BD2||sys==BD3||sys==GAL||sys==QZS){
        if((data.eph_idx_=idx_eph=SelectEph(data,-1))==-1) return 0;
        else{Eph2Clk(data,eph_[idx_eph]);}
    }
    else if(sys==GLO){
        if((data.eph_idx_=idx_geph=SelGloEph(data,-1))==-1) return 0;
        else{GloEph2Clk(data,glo_eph_[idx_geph]);}
    }

    return 1;
}

#define RTOL_KEPLER 1E-13         /* relative tolerance for Kepler equation */
#define MAX_ITER_KEPLER 30        /* max number of iteration of Kelpler */
#define SIN_5 -0.0871557427476582 /* sin(-5.0 deg) */
#define COS_5  0.9961946980917456 /* cos(-5.0 deg) */
void BrdcEph::Eph2Pos(ObsData &data, int iode, const Eph eph) {
    double tk,M,E,Ek,sinE,cosE,u,r,i,O,sin2u,cos2u,x,y,sinO,cosO,cosi,mu,omge;
    double xg,yg,zg,sino,coso;
    SatMask sat=data.sat_;
    int sys=sat.sat_sys_,n;

    if (eph.A_<=0.0) {
        data.ReSetSat();
        return;
    }
    tk=data.sig_trans_.TimeDiff(eph.toe_);

    switch (sys) {
        case GAL: mu=MU_GAL; omge=OMGE_GAL; break;
        case BD2:
        case BD3: mu=MU_BDS; omge=OMGE_BDS; break;
        default:  mu=MU_GPS; omge=OMGE_GPS; break;
    }
    M=eph.M0_+(sqrt(mu/(eph.A_*eph.A_*eph.A_))+eph.deln_)*tk;

    for (n=0,E=M,Ek=0.0; fabs(E-Ek)>RTOL_KEPLER&&n<MAX_ITER_KEPLER; n++) {
        Ek=E; E-=(E-eph.e_*sin(E)-M)/(1.0-eph.e_*cos(E));
    }
    if (n>=MAX_ITER_KEPLER) {
        return;
    }
    sinE=sin(E); cosE=cos(E);

    u=atan2(sqrt(1.0-eph.e_*eph.e_)*sinE,cosE-eph.e_)+eph.omg_;
    r=eph.A_*(1.0-eph.e_*cosE);
    i=eph.i0_+eph.idot_*tk;
    sin2u=sin(2.0*u); cos2u=cos(2.0*u);
    u+=eph.cus_*sin2u+eph.cuc_*cos2u;
    r+=eph.crs_*sin2u+eph.crc_*cos2u;
    i+=eph.cis_*sin2u+eph.cic_*cos2u;
    x=r*cos(u); y=r*sin(u); cosi=cos(i);

    /* beidou geo satellite (ref [9]) */
    if (sys==BD2&&sat.sat_prn_<=5) {
        O=eph.OMG0_+eph.OMGd_*tk-omge*eph.toes_;
        sinO=sin(O); cosO=cos(O);
        xg=x*cosO-y*cosi*sinO;
        yg=x*sinO+y*cosi*cosO;
        zg=y*sin(i);
        sino=sin(omge*tk); coso=cos(omge*tk);
        data.pos_.i_= xg*coso+yg*sino*COS_5+zg*sino*SIN_5;
        data.pos_.j_=-xg*sino+yg*coso*COS_5+zg*coso*SIN_5;
        data.pos_.k_=-yg*SIN_5+zg*COS_5;
    }
    else {
        O=eph.OMG0_+(eph.OMGd_-omge)*tk-omge*eph.toes_;
        sinO=sin(O); cosO=cos(O);
        data.pos_.i_=x*cosO-y*cosi*sinO;
        data.pos_.j_=x*sinO+y*cosi*cosO;
        data.pos_.k_=y*sin(i);
    }
    tk=data.sig_trans_.TimeDiff(eph.toc_);
    data.clk_err_[0]=eph.f0_+eph.f1_*tk+eph.f2_*tk*tk;

    // relativity correction
    // the onboad clock will be affected by a relativistic clock correction
    // caused by the motion of the satellite as well as the change in the
    // gravitational potential, the effect due to the orbit eccentricity can
    // be computed by
    data.clk_rel_=eph.RelClkCorr(mu,sinE);
    data.clk_err_[0]-=data.clk_rel_;

    /* position and clock error variance */
    data.sat_var_=UraEph(eph.sva_);
}

#define J2_GLO   1.0826257E-3     /* 2nd zonal harmonic of geopot   ref [2] */
#define RE_GLO   6378136.0        /* radius of earth (m)            ref [2] */
void BrdcEph::GloDefEq(const double *x, double *xdot, const double *acc) {
    double a,b,c,r2=Dot(x,x,3),r3=r2*sqrt(r2),omg2=SQR(OMGE_GLO);

    if (r2<=0.0) {
        xdot[0]=xdot[1]=xdot[2]=xdot[3]=xdot[4]=xdot[5]=0.0;
        return;
    }
    /* ref [2] A.3.1.2 with bug fix for xdot[4],xdot[5] */
    a=1.5*J2_GLO*MU_GLO*SQR(RE_GLO)/r2/r3; /* 3/2*J2*mu*Ae^2/r^5 */
    b=5.0*x[2]*x[2]/r2;                    /* 5*z^2/r^2 */
    c=-MU_GLO/r3-a*(1.0-b);                /* -mu/r^3-a(1-b) */
    xdot[0]=x[3]; xdot[1]=x[4]; xdot[2]=x[5];
    xdot[3]=(c+omg2)*x[0]+2.0*OMGE_GLO*x[4]+acc[0];
    xdot[4]=(c+omg2)*x[1]-2.0*OMGE_GLO*x[3]+acc[1];
    xdot[5]=(c-2.0*a)*x[2]+acc[2];
}

void BrdcEph::GloOrbit(double t, double *x, const double *acc) {
    double k1[6],k2[6],k3[6],k4[6],w[6];
    int i;

    GloDefEq(x,k1,acc); for (i=0; i<6; i++) w[i]=x[i]+k1[i]*t/2.0;
    GloDefEq(w,k2,acc); for (i=0; i<6; i++) w[i]=x[i]+k2[i]*t/2.0;
    GloDefEq(w,k3,acc); for (i=0; i<6; i++) w[i]=x[i]+k3[i]*t;
    GloDefEq(w,k4,acc);
    for (i=0; i<6; i++) x[i]+=(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i])*t/6.0;
}

#define TSTEP      60.0             /* integration step glonass ephemeris (s) */
#define ERREPH_GLO 5.0            /* error of glonass ephemeris (m) */
void BrdcEph::GloEph2Pos(ObsData &data, const GloEph glo_eph) {
    double t,tt,x[6];
    int i;

    t=data.sig_trans_.TimeDiff(glo_eph.toe_);

    data.clk_err_[0]=-glo_eph.taun_+glo_eph.gamn_*t;

    for (i=0; i<3; i++) {
        x[i]=glo_eph.pos_[i];
        x[i+3]=glo_eph.vel_[i];
    }
    for (tt=t<0.0 ? -TSTEP : TSTEP; fabs(t)>1E-9; t-=tt) {
        if (fabs(t)<TSTEP) tt=t;
        GloOrbit(tt,x,glo_eph.acc_);
    }
    //for (i=0; i<3; i++) data.pos_[i]=x[i];
    data.pos_=Vec3(x[0],x[1],x[2]);

    data.sat_var_=SQR(ERREPH_GLO);
}

int BrdcEph::SatBrdcPos(ObsData &data, int iode) {
    data.svh_=-1;
    SatMask sat=data.sat_;
    ObsData d1E_3;
    int sys=sat.sat_sys_,idx_eph=0,idx_geph=0;
    if(sys==GPS||sys==BD2||sys==BD3||sys==GAL||sys==QZS){
        if((idx_eph=SelectEph(data,-1))==-1) return 0;
        else{Eph2Pos(data,-1,eph_[idx_eph]);}
        if(data.sat_var_<129.0) data.svh_=eph_[idx_eph].svh_;

        d1E_3.sig_trans_=data.sig_trans_;d1E_3.sig_trans_.TimeAdd(1E-3);
        Eph2Pos(d1E_3,-1,eph_[idx_eph]);
    }
    else if(sys==GLO){
        if((idx_geph=SelGloEph(data,-1))==-1) return 0;
        else{GloEph2Pos(data,glo_eph_[idx_geph]);}

        d1E_3.sig_trans_=data.sig_trans_;d1E_3.sig_trans_.TimeAdd(1E-3);
        GloEph2Pos(d1E_3,glo_eph_[idx_geph]);
        data.svh_=glo_eph_[idx_eph].svh_;
    }
    else{
        cout<<"invalid satellite system"<<endl;
        return 0;
    }
    data.vel_=Vec3(data.pos_.i_-d1E_3.pos_.i_,data.pos_.j_-d1E_3.pos_.j_,data.pos_.k_-d1E_3.pos_.k_)/1E-3;
    data.clk_err_[1]=(data.clk_err_[0]-d1E_3.clk_err_[1])/1E-3;

    return 1;
}

int BrdcEph::BrdcPosClk(Obss& obs) {
    int num=0;
    for(int i=0;i<obs.num_&&i<MAXOBS;i++){

        obs.sat_infos_[i].ReSetSat();

        // approximate signal translate time
        if(!obs.sat_infos_[i].SigTransTime()) continue;

        if(!SatBrdcClk(obs.sat_infos_[i])){
            char buff[MAXBUFF]={'\0'};
            sprintf(buff,"%4d %20s: %4s no broadcast clock",obs.sat_infos_[i].epoch_,
                    obs.sat_infos_[i].sig_recep_.time_str_.c_str(),obs.sat_infos_[i].sat_.sat_id_.c_str());
            LOG(WARNING)<<buff;
            obs.sat_infos_[i].stat_=EX_NOPRED;
            continue;
        }
        // signal translate time with clock correction
        obs.sat_infos_[i].ClkCorrection();

        if(!SatBrdcPos(obs.sat_infos_[i],999)){
            LOG(INFO)<<obs.sat_infos_[i].sig_recep_.time_str_<<": "<<obs.sat_infos_[i].sat_.sat_id_<<" no broadcast orbit";
            obs.sat_infos_[i].stat_=EX_SVH;
            continue;
        }

        if(obs.sat_infos_[i].svh_!=-1) num++;

        LOG_IF(i==0,DEBUG)<<"Broadcast ephemeris";
        char sat_info[MAXBUFF]={'\0'};
        sprintf(sat_info,"%s %14.3f %14.3f %14.3f %12.3f",
                         obs.sat_infos_[i].sat_.sat_id_.c_str(),obs.sat_infos_[i].pos_.i_,obs.sat_infos_[i].pos_.j_,\
                         obs.sat_infos_[i].pos_.k_,obs.sat_infos_[i].clk_err_[0]*CLIGHT);
        LOG(DEBUG)<<sat_info;
    }

    return num;
}

PrecEph::PrecEph() {
    clk_var_=0.0;
}

PrecEph::~PrecEph() {

}

#define NMAX        10              /* order of polynomial interpolation */
#define MAXDTE      900.0           /* max time difference to ephem time (s) */
#define EXTERR_CLK  1E-3            /* extrapolation error for clock (m/s) */
#define EXTERR_EPH  5E-7            /* extrapolation error for ephem (m/s^2) */
static void EarthRotCorr(int k, double *pos, double p[3][NMAX+1], double dt) {
    double sinl,cosl;
    sinl=sin(OMGE_GPS*dt);
    cosl=cos(OMGE_GPS*dt);
    p[0][k]=cosl*pos[0]-sinl*pos[1];
    p[1][k]=sinl*pos[0]+cosl*pos[1];
    p[2][k]=pos[2];
}

double PrecEph::InterpoLagr(const double* dt, Time* ptime, const double* ppos, int n) {
    int i,j;
    double item=1.0,value=0.0;

    for (i=0;i<n;i++) {
        item=1.0;
        for (j=0; j<i; j++) item*=dt[j]/ptime[j].TimeDiff(ptime[i]);
        for (j=i+1; j<n; j++) item*=dt[j]/ptime[j].TimeDiff(ptime[i]);
        item*=ppos[i];
        value+=item;
    }
    return value;
}

double PrecEph::InterpolNev(const double* x, double *y,int n) {
    // Neville algorithm
    int i,j;
    for(j=1;j<n;j++){
        for(i=0;i<n-j;i++){
            y[i]=(x[i+j]*y[i]-x[i]*y[i+1])/(x[i+j]-x[i]);
        }
    }
    return y[0];
}

int PrecEph::PreEph2Pos(ObsData &data) {
    double t[NMAX+1],p[3][NMAX+1],c[2],std=0.0,s[3],sinl,cosl,*pos;
    int i,j,k,index;
    int id[NMAX+1],idx;

    data.pos_=Vec3(0.0);
    data.clk_err_[0]=0.0;
    if (np_<NMAX||
        data.sig_trans_.TimeDiff(peph_[0].time_)<-MAXDTE||
        data.sig_trans_.TimeDiff(peph_[np_-1].time_)>MAXDTE) {
        LOG(WARNING)<<data.sig_recep_.time_str_<<": "<<data.sat_.sat_id_<<" no precise ephemeris";
        return 0;
    }
    /* binary search */
    for (i=0,j=np_-1; i<j;) {
        k=(i+j)/2;
        if (peph_[k].time_.TimeDiff(data.sig_trans_)<0.0) i=k+1;
        else j=k;
    }
    index=i<=0?0:i-1;

    /* polynomial interpolation for orbit */
    i=index-(NMAX+1)/2;
    if (i<0) i=0; else if (i+NMAX>np_) i=np_-NMAX-1;

    for(j=k=0;j<NMAX*50;j++){
        if(index+j>=0&&index+j<np_&&k<=NMAX){
            id[k]=index+j;
            t[k]=peph_[id[k]].time_.TimeDiff(data.sig_trans_);
            pos=peph_[id[k]].pos_[data.sat_.sat_no_-1];
            if(Norm(pos,3)>0.0){
                EarthRotCorr(k,pos,p,t[k]);
                k++;
            }
        }
        if(k==NMAX+1) break;
        if(index-j>=0&&index-j<np_&&k<=NMAX&&j!=0){
            id[k]=index-j;
            t[k]=peph_[id[k]].time_.TimeDiff(data.sig_trans_);
            if(t[k]==10){
                int a=1;
            }
            pos=peph_[id[k]].pos_[data.sat_.sat_no_-1];
            if(Norm(pos,3)>0.0){
                EarthRotCorr(k,pos,p,t[k]);
                k++;
            }
        }
        if(k==NMAX+1) break;
    }
    if(k<=NMAX) return 0;

    for (i=0;i<NMAX+1;i++) {
        if(i==9){
            int a=1;
        }
        for(j=i+1;j<NMAX+1;j++){
            if(t[i]<=t[j]) continue;
            sinl=t[j];t[j]=t[i];t[i]=sinl;
            k=id[j];id[j]=id[i];id[i]=k;
            for(k=0;k<3;k++){
                sinl=p[k][j];
                p[k][j]=p[k][i];
                p[k][i]=sinl;
            }
        }
    }

    idx=0;
    for(i=0;i<=NMAX;i++){
        if(fabs(t[idx])<=fabs(t[i])) idx=i;
    }
    index=id[idx];
    if(t[0]>900.0||t[NMAX]<-900.0) return 0;

    data.pos_.i_=InterpolNev(t,p[0],NMAX+1);
    data.pos_.j_=InterpolNev(t,p[1],NMAX+1);
    data.pos_.k_=InterpolNev(t,p[2],NMAX+1);

    /* satellite position variance */
    for(i=0;i<3;i++) s[i]=peph_[index].std_[data.sat_.sat_no_-1][i];
    std=Norm(s,3);
    if(t[0]>0.0)         std+=EXTERR_EPH*SQR(t[0])/2.0;
    else if(t[NMAX]<0.0) std+=EXTERR_EPH*SQR(t[NMAX])/2.0;
    orb_var_=SQR(std);
    return 1;
}


int PrecEph::SatPrecClk(ObsData &data) {
    int i,j,k,idx;
    double t[2],c[2],std;

    if(nc_<2||
       data.sig_trans_.TimeDiff(pclk_[0].time_)<-MAXDTE||
       data.sig_trans_.TimeDiff(pclk_[nc_-1].time_)>MAXDTE){
        LOG(WARNING)<<data.sig_recep_.time_str_<<": "<<data.sat_.sat_id_<<" missing precise clock";
        return 0;
    }
    for(i=0,j=nc_-1;i<j;){
        k=(i+j)/2;
        if(pclk_[k].time_.TimeDiff(data.sig_trans_)<0.0) i=k+1;
        else j=k;
    }
    idx=i<=0?0:i-1;

    t[0]=data.sig_trans_.TimeDiff(pclk_[idx].time_);
    t[1]=data.sig_trans_.TimeDiff(pclk_[idx+1].time_);
    c[0]=pclk_[idx].clk_[data.sat_.sat_no_-1];
    c[1]=pclk_[idx+1].clk_[data.sat_.sat_no_-1];

    if(t[0]<=0.0){
        if((data.clk_err_[0]=c[0])==0.0) return 0;
        std=pclk_[idx].std_[data.sat_.sat_no_-1]*CLIGHT-EXTERR_CLK*t[0];
    }
    else if(t[1]>=0.0){
        if((data.clk_err_[0]=c[1])==0.0) return 0;
        std=pclk_[idx+1].std_[data.sat_.sat_no_-1]*CLIGHT+EXTERR_CLK*t[1];
    }
    else if(c[0]!=0.0&&c[1]!=0.0){
        data.clk_err_[0]=(c[1]*t[0]-c[0]*t[1])/(t[0]-t[1]);
        i=t[0]<-t[1]?0:1;
        std=pclk_[idx+i].std_[data.sat_.sat_no_-1];

        if(std*CLIGHT>0.05) std+=EXTERR_CLK*fabs(t[i]);
        else                std=std*CLIGHT+EXTERR_CLK*fabs(t[i]);
    }
    else{
        LOG(WARNING)<<data.sig_recep_.time_str_<<": "<<data.sat_.sat_id_<<" precise clock outage";
        return 0;
    }
    clk_var_=SQR(std);
    data.sat_var_=orb_var_+clk_var_;

    return 1;
}

int PrecEph::SatPrecPos(ObsData &data, int iode) {
    double tt=1E-3;
    int i;

    ObsData data1;
    data1.sig_trans_=data.sig_trans_;
    data1.sig_trans_.TimeAdd(tt);
    data1.sat_=data.sat_;

    if(data.sat_.sat_no_<=0||MAXSAT<data.sat_.sat_no_) return 0;

    if(!PreEph2Pos(data)||!SatPrecClk(data)) return 0;
    if(!PreEph2Pos(data1)||!SatPrecClk(data1)) return 0;

    // satellite antenna pco
    double pos[3]={data.pos_.i_,data.pos_.j_,data.pos_.k_};
    double dant[3]={0};
    ant_corr_->SatPcoCorr(data,pos,dant);
    data.pos_.i_+=dant[0];data.pos_.j_+=dant[1];data.pos_.k_+=dant[2];

    double vel[3]={0};
    vel[0]=data.vel_.i_=(data1.pos_.i_-pos[0])/tt;
    vel[1]=data.vel_.j_=(data1.pos_.j_-pos[1])/tt;
    vel[2]=data.vel_.k_=(data1.pos_.k_-pos[2])/tt;
    pos[0]=data.pos_.i_;pos[1]=data.pos_.j_;pos[2]=data.pos_.k_;

    // relativistic effect correction
    if(data.clk_err_[0]!=0.0){
        data.clk_err_[0]=data.clk_err_[0]-2.0*Dot(pos,vel,3)/CLIGHT/CLIGHT;
        data.clk_err_[1]=(data1.clk_err_[0]-data.clk_err_[0])/tt;
    }
    else{
        data.clk_err_[0]=data.clk_err_[1]=0.0;
    }
    return 1;
}

int PrecEph::PrecPosClk(Obss &obs) {
    int num=0;
    for(int i=0;i<obs.num_&&i<MAXOBS;i++){

        obs.sat_infos_[i].ReSetSat();

        // approximate signal translate time
        if(!obs.sat_infos_[i].SigTransTime()) continue;
        if(!SatBrdcClk(obs.sat_infos_[i])){
            LOG(WARNING)<<obs.sat_infos_[i].sig_recep_.time_str_<<": "<<obs.sat_infos_[i].sat_.sat_id_<<" no precise clock";
            continue;
        }
        // signal translate time with broadcast clock correction
        obs.sat_infos_[i].ClkCorrection();

        if(!SatPrecPos(obs.sat_infos_[i],999)){
            obs.sat_infos_[i].stat_=EX_NOPRED;
            LOG(DEBUG)<<obs.sat_infos_[i].epoch_<<" "<<obs.sat_infos_[i].sig_recep_.time_str_<<": "<<obs.sat_infos_[i].sat_.sat_id_<<" no precise orbit";
            continue;
        }

        if(obs.sat_infos_[i].svh_!=-1) num++;

        LOG_IF(i==0,DEBUG)<<"Precise ephemeris";
        char sat_info[MAXBUFF]={'\0'};
        sprintf(sat_info,"%s %14.3f %14.3f %14.3f %12.3f",
                obs.sat_infos_[i].sat_.sat_id_.c_str(),obs.sat_infos_[i].pos_.i_,obs.sat_infos_[i].pos_.j_,\
                         obs.sat_infos_[i].pos_.k_,obs.sat_infos_[i].clk_err_[0]*CLIGHT);
        LOG(DEBUG)<<sat_info;
    }

    return num;
}