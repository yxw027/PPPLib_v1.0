/**************************************************************************

Copyright(c): 2020, by Chao Chen, All rights reserved.
              China University of Mining and Technology, XuZhou, P.R. China

Author: Chao Chen(cchen@cumt.edu.cn)

Date: 2020-05-30

Description:

**************************************************************************/

#include "PPPLibGlo.h"
#include "IonModel.h"

Tec::Tec() {
    t_=Time();re_=0;
    for(int i=0;i<3;i++){
        lats_[i]=lons_[i]=hgts_[i]=ndata_[i]=0;
    }
}

Tec::~Tec() {
    data_.clear();rms_.clear();
}

IonModel::IonModel() {
    map_=ion_=ion_var_=0.0;
    for(int i=0;i<8;i++) ion_paras_[i]=0.0;
    sat_data_= nullptr;
}

IonModel::~IonModel() {
    isat_= nullptr;
    sat_data_= nullptr;
    tecs_.clear();
}

#define HION  350000.0   // ionosphere height(m)
double IonModel::IonMapFun(const Vec3 &blh, double el) {
    if(blh.k_>=HION) return 1.0;
    return 1.0/cos(asin((RE_WGS84+blh.k_)/(RE_WGS84+HION)*sin(PI/2.0-el)));
}

void IonModel::GetBrdIonPara(double *paras) {
    for(int i=0;i<8;i++) ion_paras_[i]=paras[i];
}

#define ERR_BRDCI   0.5					/* broadcast iono model error factor */
int IonModel::KlobModel(const Vec3 &blh) {
    int sat=sat_data_->sat_.sat_no_;
    const double ion_default[]={ /* 2004/1/1 */
            0.1118E-07,-0.7451E-08,-0.5961E-07, 0.1192E-06,
            0.1167E+06,-0.2294E+06,-0.1311E+06, 0.1049E+07
    };
    double tt,psi,phi,lam,amp,per,x;
    int week;
    Time obs_time=sat_data_->sig_recep_;
    Vec para(8,ion_paras_);

    if (blh.k_<-1E3||sat_data_->azel_[1]<=0) return 0;
    if (VecNorm(para)<=0.0) para=Vec(8,ion_default);

    /* earth centered angle (semi-circle) */
    psi=0.0137/(sat_data_->azel_[1]/PI+0.11)-0.022;

    /* subionospheric latitude/longitude (semi-circle) */
    phi=blh.i_/PI+psi*cos(sat_data_->azel_[0]);
    // latitude boundary protection
    if (phi> 0.416) phi= 0.416;
    else if (phi<-0.416) phi=-0.416;
    lam=blh.j_/PI+psi*sin(sat_data_->azel_[0])/cos(phi*PI);

    /* geomagnetic latitude (semi-circle) */
    phi+=0.064*cos((lam-1.617)*PI);

    /* local time (s) */
    tt=43200.0*lam+obs_time.Time2GPST(&week);
    tt-=floor(tt/86400.0)*86400.0; /* 0<=tt<86400 */

    /* slant factor */
    double map_ion=0.0;
    map_ion=isat_[sat-1].map_ion=1.0+16.0*pow(0.53-sat_data_->azel_[1]/PI,3.0);

    /* ionospheric delay */
    amp=para.dd_[0]+phi*(para.dd_[1]+phi*(para.dd_[2]+phi*para.dd_[3]));
    per=para.dd_[4]+phi*(para.dd_[5]+phi*(para.dd_[6]+phi*para.dd_[7]));
    amp=amp<0.0?0.0:amp;
    per=per<72000.0?72000.0:per;
    x=2.0*PI*(tt-50400.0)/per;

    double sion=CLIGHT*map_ion*(fabs(x)<1.57 ? 5E-9+amp*(1.0+x*x*(-0.5+x*x/24.0)):5E-9);
    double ion_var_=SQR(sion*ERR_BRDCI);
    double lam_G1=CLIGHT/FREQ_GPS_L1;
    double l=lam_G1;
    l=sat_data_->lam_[sat_data_->use_f1_];

    isat_[sat-1].sion=sion*SQR(l/lam_G1);
    isat_[sat-1].var_ion=ion_var_*SQR(l/lam_G1)*SQR(l/lam_G1);

    return 1;
}

int IonModel::IonCorr(const Vec3 blh) {
    ion_=map_=ion_var_=0.0;
    return 1;
}

KlobIon::KlobIon() = default;

KlobIon::~KlobIon() = default;

int KlobIon::IonCorr(const Vec3 blh) {
    return KlobModel(blh);
}

GIMIon::GIMIon() = default;

GIMIon::~GIMIon() = default;

int GIMIon::DataIndex(int i, int j, int k, const int *ndata) {
    if (i<0||ndata[0]<=i||j<0||ndata[1]<=j||k<0||ndata[2]<=k) return -1;
    return i+ndata[0]*(j+ndata[1]*k);
}

int GIMIon::InterpTec(Tec& tec, int k, const double *posp, double &val, double &rms) {
    double dlat,dlon,a,b,d[4]={0},r[4]={0};
    int i,j,n,index;

    val=rms=0.0;
    if (tec.lats_[2]==0.0||tec.lons_[2]==0.0) return 0;
    dlat=posp[0]*R2D-tec.lats_[0];
    dlon=posp[1]*R2D-tec.lons_[0];
    if (tec.lons_[2]>0.0) dlon-=floor(dlon/360)*360.0;  /*  0<=dlon<360 */
    else                  dlon+=floor(-dlon/360)*360.0; /* -360<dlon<=0 */

    a=dlat/tec.lats_[2];
    b=dlon/tec.lons_[2];
    i=(int)floor(a); a-=i;
    j=(int)floor(b); b-=j;

    /* get gridded tec data */
    for (n=0; n<4; n++) {
        if ((index=DataIndex(i+(n%2),j+(n<2 ? 0 : 1),k,tec.ndata_))<0) continue;
        d[n]=tec.data_[index];
        r[n]=tec.rms_[index];
    }
    if (d[0]>0.0&&d[1]>0.0&&d[2]>0.0&&d[3]>0.0) {
        /* bilinear interpolation (inside of grid) */
        val=(1.0-a)*(1.0-b)*d[0]+a*(1.0-b)*d[1]+(1.0-a)*b*d[2]+a*b*d[3];
        rms=(1.0-a)*(1.0-b)*r[0]+a*(1.0-b)*r[1]+(1.0-a)*b*r[2]+a*b*r[3];
    }
    /* nearest-neighbour extrapolation (outside of grid) */
    else if (a<=0.5&&b<=0.5&&d[0]>0.0) {val=d[0];rms=r[0];}
    else if (a> 0.5&&b<=0.5&&d[1]>0.0) {val=d[1];rms=r[1];}
    else if (a<=0.5&&b> 0.5&&d[2]>0.0) {val=d[2];rms=r[2];}
    else if (a> 0.5&&b> 0.5&&d[3]>0.0) {val=d[3];rms=r[3];}
    else {
        i=0;
        for (n=0; n<4; n++) if (d[n]>0.0) { i++;val+=d[n];rms+=r[n];}
        if (i==0) return 0;
        val/=i; rms/=i;
    }
    return 1;
}

double GIMIon::IonPpp(const Vec3& blh, const double re, const double hion, double *posp) {
    double cosaz,rp,ap,sinap,tanap;

    /* single layer mapping function(M-SLM) */
    rp=re/(re+hion)*sin(0.9782*(PI/2.0-sat_data_->azel_[1]));
    ap=PI/2.0-sat_data_->azel_[1]-asin(rp);
    sinap=sin(ap);
    tanap=tan(ap);
    cosaz=cos(sat_data_->azel_[0]);
    posp[0]=asin(sin(blh.i_)*cos(ap)+cos(blh.i_)*sinap*cosaz);

    if ((blh.i_> 70.0*D2R&& tanap*cosaz>tan(PI/2.0-blh.i_))||
        (blh.i_<-70.0*D2R&&-tanap*cosaz>tan(PI/2.0+blh.i_))) {
        posp[1]=blh.j_+PI-asin(sinap*sin(sat_data_->azel_[0])/cos(posp[0]));
    }
    else {
        posp[1]=blh.j_+asin(sinap*sin(sat_data_->azel_[0])/cos(posp[0]));
    }
    return 1.0/sqrt(1.0-rp*rp);
}

int GIMIon::GIMModel(const Vec3& blh,Tec& tec,double *delay,double *var) {
    const double fact=40.30E16/SQR(FREQ_GPS_L1);
    double vtec,rms,hion,posp[2]={0};
    int i,sat=sat_data_->sat_.sat_no_;

    *delay=*var=0.0;

    for(i=0;i<tec.ndata_[2];i++){
        hion=tec.hgts_[0]+tec.hgts_[2]*i;
        isat_[sat-1].map_ion=IonPpp(blh,tec.re_,hion,posp);
        if(!InterpTec(tec,i,posp,vtec,rms)) return 0;

        *delay+=fact*vtec*isat_[sat-1].map_ion;
        *var+=fact*fact*rms*rms*isat_[sat-1].map_ion;
    }
    return 1;
}

#define MIN_EL  0.0
#define MIN_HGT -1000.0
#define VAR_NOTEC SQR(30.0)
int GIMIon::IonCorr(Vec3 blh) {
    int i,stat[2];
    int sat=sat_data_->sat_.sat_no_;
    double tt,dels[2],vars[2];
    if(sat_data_->azel_[1]<MIN_EL||blh.k_<MIN_HGT){
        map_=ion_=0.0;
        isat_[sat-1].sion=isat_[sat-1].map_ion=0.0;
        isat_[sat-1].var_ion=VAR_NOTEC;
        return 1;
    }
    for(i=0;i<tecs_.size();i++){
        if(tecs_[i].t_.TimeDiff(sat_data_->sig_recep_)>0.0) break;
    }
    if(i==0||i>=tecs_.size()){
        cout<<"tec grid out of period"<<endl;
        return 0;
    }
    if((tt=tecs_[i].t_.TimeDiff(tecs_[i-1].t_))==0.0){
        cout<<"tec grid time interval error"<<endl;
        return 0;
    }
    stat[0]=GIMModel(blh,tecs_[i-1],dels,vars);
    stat[1]=GIMModel(blh,tecs_[i],dels+1,vars+1);

    if(!stat[0]&&!stat[1]){
        cout<<"tec grid out of area pos"<<endl;
        return 0;
    }
    if(stat[0]&&stat[1]){
        double a=sat_data_->sig_recep_.TimeDiff(tecs_[i-1].t_)/tt;
        isat_[sat-1].sion=dels[0]*(1.0-a)+dels[1]*a;
        isat_[sat-1].var_ion=vars[0]*(1.0-a)+vars[1]*a;
    }
    else if(stat[0]){
        isat_[sat-1].sion=dels[0];isat_[sat-1].var_ion=vars[0];
    }
    else{
        isat_[sat-1].sion=dels[1];isat_[sat-1].var_ion=vars[1];
    }
    return 1;
}

IFIon::IFIon() = default;

IFIon::~IFIon() = default;

int IFIon::IonCorr(Vec3 blh) {
    int sat=sat_data_->sat_.sat_no_;
    map_=ion_var_=0.0;
    isat_[sat-1].sion=0.0;
    isat_[sat-1].var_ion=SQR(0.01);
    return 1;
}

