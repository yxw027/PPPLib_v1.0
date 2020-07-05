/**************************************************************************

Copyright(c): 2020, by Chao Chen, All rights reserved.
              China University of Mining and Technology, XuZhou, P.R. China

Author: Chao Chen(cchen@cumt.edu.cn)

Date: 2020-05-30

Description:

**************************************************************************/

#include "CmnFunc.h"

// time class
#define MAXLEAPS    64                                  /* max number of leap seconds table */
static const double kGPSTimeStart[]={1980,1, 6,0,0,0};
static const double kBDSTimeStart[]={2006,1, 1,0,0,0};
static double kLeapSec[MAXLEAPS+1][7]={                 /* leap seconds (y,m,d,h,m,s,utc-gpst) */
        {2017,1,1,0,0,0,-18},
        {2015,7,1,0,0,0,-17},
        {2012,7,1,0,0,0,-16},
        {2009,1,1,0,0,0,-15},
        {2006,1,1,0,0,0,-14},
        {1999,1,1,0,0,0,-13},
        {1997,7,1,0,0,0,-12},
        {1996,1,1,0,0,0,-11},
        {1994,7,1,0,0,0,-10},
        {1993,7,1,0,0,0, -9},
        {1992,7,1,0,0,0, -8},
        {1991,1,1,0,0,0, -7},
        {1990,1,1,0,0,0, -6},
        {1988,1,1,0,0,0, -5},
        {1985,7,1,0,0,0, -4},
        {1983,7,1,0,0,0, -3},
        {1982,7,1,0,0,0, -2},
        {1981,7,1,0,0,0, -1},
        {0}
};

Time::Time() {
    time_=0; sec_=0.0;
    for(int i=0;i<6;i++) epoch_[i]=0.0;
    epoch_[0]=1980.0;epoch_[1]=12.0;epoch_[2]=6;
    time_str_='\0';
    doy_=0;
    time_sys_=0;
    mjd_={0};
}

Time::Time(const double *epoch) {
    Epoch2Time(epoch);
}

Time::~Time() {

}

Time* Time::Epoch2Time(const double * epoch) {
    const int doy[]={ 1,32,60,91,121,152,182,213,244,274,305,335 };

    int days, dsec, year=int(epoch[0]), mon=(int)epoch[1], day=(int)epoch[2];

    if (year<1970||2099<year||mon<1||12<mon) {
        time_=0;sec_=0; return this;
    }
    /* leap year if year%4==0 in 1901-2099 */
    days=(year-1970)*365+(year-1969)/4+doy[mon-1]+day-2+(year%4==0&&mon>=3 ? 1 : 0);
    dsec=(int)floor(epoch[5]);
    time_=(time_t)days*86400+(time_t)epoch[3]*3600+(time_t)epoch[4]*60+dsec;
    sec_=epoch[5]-dsec;

    return this;    // 在类的非静态成员函数中返回类对象本身的时候 | 静态成员函数和非静态成员函数的根本区别在于有无this指针
}

void Time::Time2Epoch() {
    const int mday[]={ /* # of days in a month */
            31,28,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31,
            31,29,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31
    };
    int days, dsec, mon, day;

    /* leap year if year%4==0 in 1901-2099 */
    days=(int)(time_/86400);
    dsec=(int)(time_-(time_t)days*86400);
    for (day=days%1461, mon=0; mon<48; mon++) {
        if (day>=mday[mon]) day-=mday[mon]; else break;
    }
    epoch_[0]=1970+days/1461*4+mon/12; epoch_[1]=mon%12+1; epoch_[2]=day+1;
    epoch_[3]=dsec/3600; epoch_[4]=dsec%3600/60; epoch_[5]=dsec%60+sec_;
}

int Time::Str2Time(std::__cxx11::string s) {
    if (sscanf(s.c_str(),"%lf %lf %lf %lf %lf %lf",epoch_,epoch_+1,epoch_+2,epoch_+3,epoch_+4,epoch_+5)<6)
        return -1;
    if (epoch_[0]<100) epoch_[0]+=2000;
    if (epoch_[0]<=1900||epoch_[1]==0||epoch_[2]==0) return -1;

    Epoch2Time(epoch_);

    time_str_=s;

    return 0;
}

string Time::Time2Str(int n) {
    if (n<0) n=0; else if (n>12) n=12;
    string str;
    if (1.0-sec_<0.5/pow(10.0,n)) { time_++; sec_=0.0; };
    Time2Epoch();
    time_str_=Int2Str(4,"0",(int)epoch_[0],str)+"/"+Int2Str(2,"0",(int)epoch_[1],str)+"/"+
              Int2Str(2,"0",(int)epoch_[2],str)+" "+Int2Str(2,"0",(int)epoch_[3],str)+":"+
              Int2Str(2,"0",(int)epoch_[4],str)+":"+Doul2Str(2+n+1,n,"0",epoch_[5],str);
    return time_str_;
}

Time* Time::TimeAdd(double sec) {
    double tt;
    sec_+=sec;tt=floor(sec_);time_+=(int)tt;sec_-=tt;
    return this;
}

double Time::TimeDiff(class Time t2) const {
    return difftime(time_,t2.time_)+sec_-t2.sec_;
}

double Time::Time2Sec(Time& day) {
    double sec;
    double epoch_0[6]={0};
    int i;
    Time2Epoch();

    sec=epoch_[3]*3600.0+epoch_[4]*60.0+epoch_[5];
    for (i=0;i<3;i++) epoch_0[i]=epoch_[i];
    day.Epoch2Time(epoch_0);
    return sec;
}

double Time::Time2Doy() {
    double epoch_0[6]={0};
    Time t0;

    Time2Epoch();
    epoch_0[0]=epoch_[0]; epoch_0[1]=epoch_0[2]=1.0; epoch_0[3]=epoch_0[4]=epoch_0[5]=0.0;
    return this->TimeDiff(*t0.Epoch2Time(epoch_0))/86400.0+1.0;
}

Time* Time::GPST2Time(int week, double sec) {
    Epoch2Time(kGPSTimeStart);

    if (sec<-1E9||1E9<sec) sec=0.0;
    time_+=(time_t)86400*7*week+(int)sec;
    sec_=sec-(int)sec;

    return this;
}

double Time::Time2GPST(int * week) {
    Time t0;
    time_t sec;

    t0.Epoch2Time(kGPSTimeStart);
    sec=time_-t0.time_;

    int w=(int)(sec/(86400*7));
    if (week) *week=w;

    return (double)(sec-(double)w*86400*7)+sec_;
}

Time* Time::UTC2GPST() {
    int i;
    Time t0;

    for (i=0;kLeapSec[i][0]>0;i++) {
        if (TimeDiff(*t0.Epoch2Time(kLeapSec[i]))>=0.0)
            return TimeAdd(-kLeapSec[i][6]);
    }
    return this;
}

Time* Time::GPST2UTC() {
    int i;
    Time tu,t0;

    for (i=0; kLeapSec[i][0]>0; i++) {
        tu=*this;
        tu.TimeAdd(kLeapSec[i][6]);
        if (tu.TimeDiff(*t0.Epoch2Time(kLeapSec[i]))>=0.0) { *this=tu; return this; }
    }
    return this;
}

Time* Time::GPST2BDT() {
    return TimeAdd(-14.0);
}

Time* Time::BDT2GPST() {
    return TimeAdd(14.0);
}

Time* Time::BDT2Time(int week, double wos) {
    Epoch2Time(kBDSTimeStart);

    if (wos<-1E9||1E9<wos) wos=0.0;
    time_+=(time_t)86400*7*week+(int)wos;
    sec_=wos-(int)wos;

    return this;
}

double Time::Time2BDT(int *week) {
    Time t0;
    time_t sss;

    t0.Epoch2Time(kBDSTimeStart);
    sss=time_-t0.time_;

    int w=(int)(sss/(86400*7));
    if (week) *week=w;

    return (double)(sss-(double)w*86400*7)+sec_;
}

int Time::Monitor(Time ts, Time te, double tint) {
    return (tint<=0.0||fmod(Time2GPST(NULL)+DTTOL,tint)<=DTTOL*2.0)&&
           (ts.time_==0||TimeDiff(ts)>=-DTTOL)&&
           (te.time_==0||TimeDiff(te)<  DTTOL);
}

Time* Time::AdjWeek(Time t) {
    double dt=TimeDiff(t);
    if (dt < -302400.0) return TimeAdd( 604800.0);
    if (dt >  302400.0) return TimeAdd(-604800.0);
    return this;
}

Time* Time::AdjDay(Time t) {
    double dt=TimeDiff(t);
    if (dt < -43200.0) return TimeAdd( 86400.0);
    if (dt >  43200.0) return TimeAdd(-86400.0);
    return this;
}

Time* Time::CopyTime(Time t) {
    time_=t.time_;
    sec_=t.sec_;
    time_str_=t.time_str_;
    time_sys_=t.time_sys_;
    for (int i=0; i<6; i++) epoch_[i]=t.epoch_[i];

    return this;
}

void Time::Time2Mjd() {
    this->Time2Epoch();
    int year=(int)Round(epoch_[0]);
    int mon=(int)Round(epoch_[1]);
    if(mon<=2){
        year=year-1;
        mon=mon+12;
    }

    int a=(int)(365.25*year);
    int b=(int)(30.6001*(mon+1));
    mjd_.day=a+b+Round(epoch_[2])-679019;
    mjd_.sod.sn=Round(epoch_[3])*3600+Round(epoch_[4])*60+Round(epoch_[5]);
    mjd_.sod.tos=epoch_[5]-Round(epoch_[5]);
}

double Time::UTC2Gmst(double ut1_utc) {
    const double ep2000[]={ 2000,1,1,12,0,0 };
    Time tut,tut0,t2000;
    double ut,t1,t2,t3,gmst0,gmst;

    tut=*this; tut.TimeAdd(ut1_utc);
    ut=tut.Time2Sec(tut0);
    t1=tut0.TimeDiff(*t2000.Epoch2Time(ep2000))/86400.0/36525.0;
    t2=t1*t1; t3=t2*t1;
    gmst0=24110.54841+8640184.812866*t1+0.093104*t2-6.2E-6*t3;
    gmst=gmst0+1.002737909350795*ut;

    return fmod(gmst,86400.0)*PI/43200.0; /* 0 <= gmst <= 2*PI */
}

// Core math related class
int Round(double d) {
    int i;
    if(d>=0)
        i=(int)(d+0.5);
    else
        i=(int)(d-0.5);
    return i;
}


static int QuickSortOnce(double *array,int low,int high) {
    double pivot=array[low];
    int i=low,j=high;

    while(i<j){
        while(array[j]>=pivot&&i<j) j--;
        array[i]=array[j];

        while(array[i]<=pivot&&i<j) i++;
        array[j]=array[i];
    }
    array[i]=pivot;
    return i;
}

static void QuickSort(double *array,int low,int high) {
    if(low>=high) return;

    int pivot=QuickSortOnce(array,low,high);

    QuickSort(array,low,pivot-1);
    QuickSort(array,pivot+1,high);
}

double GetMedian(const double *array,int n) {
    double a[n];
    for(int i=0;i<n;i++) a[i]=array[i];
    QuickSort(a,0,n-1);

    if(n%2!=0) return a[n/2];
    else return (a[n/2]+a[n/2-1])/2;
}

Vec3::Vec3(void) {
    i_=j_=k_=0.0;
}

Vec3::Vec3(double xyz) {
    i_=j_=k_=xyz;
}

Vec3::Vec3(double x, double y, double z) {
    i_=x;j_=y;k_=z;
}

Vec3::~Vec3() {

}

Vec3 Vec3::operator+(const Vec3 &v) const {
    return Vec3(this->i_+v.i_,this->j_+v.j_,this->k_+v.k_);
}

Vec3 Vec3::operator-(const Vec3 &v) const {
    return Vec3(this->i_-v.i_,this->j_-v.j_,this->k_-v.k_);
}

Vec3 Vec3::operator*(double f) const {
    return Vec3(this->i_*f,this->j_*f,this->k_*f);
}

Vec3 Vec3::operator*(const Vec3 &v) const {
    return Vec3(this->j_*v.k_-this->k_*v.j_,this->k_*v.i_-this->i_*v.k_,this->i_*v.j_-this->j_*v.i_);
}

Vec3 Vec3::operator*(const Mat3 &m) const {
    return Vec3(i_*m.e00_+j_*m.e10_+k_*m.e20_,i_*m.e01_+j_*m.e11_+k_*m.e21_,i_*m.e02_+j_*m.e12_+k_*m.e22_);
}

Vec3 Vec3::operator/(double f) const {
    return Vec3(i_/f,j_/f,k_/f);
}

Vec3 Vec3::operator/(const Vec3 & v) const {
    return Vec3(i_/v.i_,j_/v.j_,k_/v.k_);
}

Vec3& Vec3::operator+=(const Vec3 &v) {
    i_+=v.i_;j_+=v.j_;k_+=v.k_;
    return *this;
}

Vec3& Vec3::operator-=(const Vec3 &v) {
    i_-=v.i_;j_-=v.j_;k_-=v.k_;
    return *this;
}

Vec3& Vec3::operator*=(double f) {
    i_*=f;j_*=f;k_*=f;
    return *this;
}

Vec3& Vec3::operator/=(double f) {
    i_/=f;j_/=f;k_/=f;
    return *this;
}

Vec3 & Vec3::operator/=(const Vec3 & v) {
    i_/=v.i_;j_/=v.j_;k_/=v.k_;
    return *this;
}

double Vec3Norm(const Vec3 &v) {
    return sqrt(SQR(v.i_)+SQR(v.j_)+SQR(v.k_));
}

double Vec3NormXY(const class Vec3 & v) {
    return sqrt(SQR(v.i_)+SQR(v.j_));
}

double VecDot(const Vec3 &v1, const Vec3 &v2) {
    return (v1.i_*v2.i_+v1.j_*v2.j_+v1.k_*v2.k_);
}

double VecDot2(const Vec3& v1,const Vec3& v2) {
    return (v1.i_*v2.i_+v1.j_*v2.j_);
}

Vec3 VecPow(const Vec3 &v, int k) {
    Vec3 pp=v;
    for(int i=1;i<k;i++){
        pp.i_*=v.i_,pp.j_*=v.j_,pp.k_*=v.k_;
    }
    return pp;
}


Vec3 Xyz2Blh(const Vec3 &xyz,int coor) {
    double re,fe;
    Vec3 BLH;
    if(coor==CGCS2000){re=RE_CGCS2000,fe=EE_CGCS2000;}
    else if(coor==WGS84){re=RE_WGS84,fe=EE_WGS84;}
    double r2=VecDot2(xyz,xyz),e2=fe*(2.0-fe),z,zk,v=re,sinp;

    for(z=xyz.k_,zk=0.0;fabs(z-zk)>=1E-4;){
        zk=z;
        sinp=z/sqrt(r2+z*z);
        v=re/sqrt(1.0-e2*sinp*sinp);
        z=xyz.k_+v*e2*sinp;
    }
    BLH.i_=r2>1E-12?atan(z/sqrt(r2)):(xyz.k_>0.0?PI/2.0:-PI/2.0);;
    BLH.j_=r2>1E-12?atan2(xyz.j_,xyz.i_):0.0;
    BLH.k_=sqrt(r2+z*z)-v;
    return BLH;
}

Vec3 Blh2Xyz(const Vec3& blh,int coor) {

}

Mat3 Pos2Cen(const Vec3& pos) {
    double si=sin(pos.i_),ci=cos(pos.i_),sj=sin(pos.j_),cj=cos(pos.j_);
    return Mat3(-sj,-si*cj,ci*cj,
                 cj,-si*sj,ci*sj,
                 0.0,  ci,       si);	//C_n^e
//    return Mat3(-sj,cj,0,
//                -sj*cj,-si*sj,ci,
//                ci*cj,  ci*sj,       si   );	//C_n^e
}

Vec3 Xyz2Enu(const Vec3& xyz,const Vec3& pos) {
    // (C_e^n*V^e)=(V^e)^T*(C_e^n)^T=(V^e)^T*C_n^e;
    return xyz*Pos2Cen(pos);
}

const Vec3 I31(1.0),O31(0.0);


Mat3::Mat3() {
    e00_=e01_=e02_=0.0;
    e10_=e11_=e12_=0.0;
    e20_=e21_=e22_=0.0;
}

Mat3::Mat3(double xx, double xy, double xz, double yx, double yy, double yz, double zx, double zy, double zz) {
    e00_=xx;e01_=xy;e02_=xz;
    e10_=yx;e11_=yy;e12_=yz;
    e20_=zx;e21_=zy;e22_=zz;
}

Mat3::Mat3(const Vec3 &v0, const Vec3 &v1, const Vec3 &v2) {
    e00_=v0.i_;e01_=v0.j_;e02_=v0.k_;
    e10_=v1.i_;e11_=v1.j_;e12_=v1.k_;
    e20_=v2.i_;e21_=v2.j_;e22_=v2.k_;
}

Mat3::~Mat3() {

}

Mat3 Mat3::operator+(const Mat3 & m) const {
    Mat3 m_tmp;

    m_tmp.e00_=e00_+m.e00_;m_tmp.e01_=e01_+m.e01_;m_tmp.e02_=e02_+m.e02_;
    m_tmp.e10_=e10_+m.e10_;m_tmp.e11_=e11_+m.e11_;m_tmp.e12_=e12_+m.e12_;
    m_tmp.e20_=e20_+m.e20_;m_tmp.e21_=e21_+m.e21_;m_tmp.e22_=e22_+m.e22_;
    return m_tmp;
}

Mat3 Mat3::operator-(const class Mat3 & m) const {
    Mat3 m_tmp;

    m_tmp.e00_=e00_-m.e00_;m_tmp.e01_=e01_-m.e01_;m_tmp.e02_=e02_-m.e02_;
    m_tmp.e10_=e10_-m.e10_;m_tmp.e11_=e11_-m.e11_;m_tmp.e12_=e12_-m.e12_;
    m_tmp.e20_=e20_-m.e20_;m_tmp.e21_=e21_-m.e21_;m_tmp.e22_=e22_-m.e22_;
    return m_tmp;
}

Mat3 Mat3 ::operator*(const Mat3 & m) const {
    Mat3 m_tmp;
    m_tmp.e00_ = e00_*m.e00_ + e01_*m.e10_ + e02_*m.e20_;
    m_tmp.e00_ = e00_*m.e00_ + e01_*m.e10_ + e02_*m.e20_;
    m_tmp.e01_ = e00_*m.e01_ + e01_*m.e11_ + e02_*m.e21_;
    m_tmp.e02_ = e00_*m.e02_ + e01_*m.e12_ + e02_*m.e22_;
    m_tmp.e10_ = e10_*m.e00_ + e11_*m.e10_ + e12_*m.e20_;
    m_tmp.e11_ = e10_*m.e01_ + e11_*m.e11_ + e12_*m.e21_;
    m_tmp.e12_ = e10_*m.e02_ + e11_*m.e12_ + e12_*m.e22_;
    m_tmp.e20_ = e20_*m.e00_ + e21_*m.e10_ + e22_*m.e20_;
    m_tmp.e21_ = e20_*m.e01_ + e21_*m.e11_ + e22_*m.e21_;
    m_tmp.e22_ = e20_*m.e02_ + e21_*m.e12_ + e22_*m.e22_;
    return m_tmp;
}

Mat3 Mat3::operator*(double f) const {
    return Mat3(e00_*f,e01_*f,e02_*f, e10_*f,e11_*f,e12_*f, e20_*f,e21_*f,e22_*f);
}

Vec3 Mat3::operator*(const Vec3 & v) const {
    return Vec3(e00_*v.i_+e01_*v.j_+e02_*v.k_,e10_*v.i_+e11_*v.j_+e12_*v.k_,e20_*v.i_+e21_*v.j_+e22_*v.k_);
}

Mat3& Mat3::operator+=(const Mat3 &m) {
    e00_+=m.e00_;e01_+=m.e01_;e02_+=m.e02_;
    e10_+=m.e10_;e11_+=m.e11_;e12_+=m.e12_;
    e20_+=m.e20_;e21_+=m.e21_;e22_+=m.e22_;
    return *this;
}

Mat3& Mat3::operator-=(const Mat3 &m) {
    e00_-=m.e00_;e01_-=m.e01_;e02_-=m.e02_;
    e10_-=m.e10_;e11_-=m.e11_;e12_-=m.e12_;
    e20_-=m.e20_;e21_-=m.e21_;e22_-=m.e22_;
    return *this;
}

Mat3 operator-(const Mat3& m) {
    return Mat3(-m.e00_,-m.e01_,-m.e02_,
                -m.e10_,-m.e11_,-m.e12_,
                -m.e20_,-m.e21_,-m.e22_);
}

Mat3 operator~(const Mat3& m) {
    return Mat3(m.e00_,m.e10_,m.e20_,
                m.e01_,m.e11_,m.e21_,
                m.e02_,m.e12_,m.e22_);
}

Mat3 operator*(double f,const Mat3& m) {
    return Mat3(m.e00_*f,m.e01_*f,m.e02_*f,
                m.e10_*f,m.e11_*f,m.e12_*f,
                m.e20_*f,m.e21_*f,m.e22_*f);
}

Mat3 Mat3Pow(const Mat3& m, int k) {
    Mat3 mat=m;
    for(int i=1;i<k;i++) mat=mat*m;
    return mat;
}

double Mat3Trace(const Mat3& m) {
    return (m.e00_+m.e11_+m.e22_);
}

double Mat3Det(const Mat3& m) {
    return (m.e00_*(m.e11_*m.e22_-m.e12_*m.e21_)-m.e01_*(m.e10_*m.e22_-m.e12_*m.e20_)+m.e02_*(m.e10_*m.e21_-m.e11_*m.e20_));
}

Vec3 Mat3DigV(const Mat3& m) {
    return Vec3(m.e00_,m.e11_,m.e22_);
}

Mat3 Mat3DigM(const Mat3& m) {
    return Mat3(m.e00_,0.0,0.0,
                 0.0,m.e11_,0.0,
                 0.0,0.0,m.e22_);
}

Mat3 Mat3Adj(const Mat3& m) {
    Mat3 m_tmp;
    m_tmp.e00_= (m.e11_*m.e22_-m.e12_*m.e21_);
    m_tmp.e10_=-(m.e10_*m.e22_-m.e12_*m.e20_);
    m_tmp.e20_= (m.e10_*m.e21_-m.e11_*m.e20_);
    m_tmp.e01_=-(m.e01_*m.e22_-m.e02_*m.e21_);
    m_tmp.e11_= (m.e00_*m.e22_-m.e02_*m.e20_);
    m_tmp.e21_=-(m.e00_*m.e21_-m.e01_*m.e20_);
    m_tmp.e02_= (m.e01_*m.e12_-m.e02_*m.e11_);
    m_tmp.e12_=-(m.e00_*m.e12_-m.e02_*m.e10_);
    m_tmp.e22_= (m.e00_*m.e11_-m.e01_*m.e10_);
    return m_tmp;
}

Mat3 Mat3Inv(const Mat3& m) {
    Mat3 adj_m=Mat3Adj(m);
    double det_m=m.e00_*adj_m.e00_+m.e01_*adj_m.e10_+m.e02_*adj_m.e20_;
    return adj_m*(1.0/det_m);
}

Vec::Vec() {

}

Vec::Vec(int row, int clm) {
    if(clm==1) {row_=row;clm_=1;}
    else       {row_=1;  clm_=clm_;}
    rc_=row_*clm_;
}

Vec::Vec(int row, double f) {
    row_=row,clm_=1;rc_=row_*clm_;
    for(int i=0;i<row_;i++) dd_[i]=f;
}

Vec::Vec(int row, const double *pf) {
    row_=row;clm_=1;rc_=row_*clm_;
    memcpy(dd_,pf,row* sizeof(double));
}

Vec::Vec(const Vec3 &v) {
    row_=3;clm_=1;rc_=row_*clm_;
    dd_[0]=v.i_;dd_[1]=v.j_;dd_[2]=v.k_;
}

Vec::Vec(const Vec3 &v1, const Vec3 &v2) {
    row_=6;clm_=1;rc_=row_*clm_;
    dd_[0]=v1.i_;dd_[1]=v1.j_;dd_[2]=v1.k_;
    dd_[3]=v2.i_;dd_[4]=v2.j_;dd_[5]=v2.k_;
}

Vec::Vec(int row, vector<double> v1) {
    row_=row;clm_=1;rc_=row_*clm_;

    double *tp=dd_;
    for(vector<double>::iterator iter=v1.begin();iter!=v1.end();++iter,++tp){
        *tp=*iter;
    }
    tp=dd_;
}

Vec::~Vec() {

}

Vec Vec::operator+(const Vec &v) const {
    assert(row_==v.row_&&clm_==v.clm_);
    const double *p1=dd_, *p2=v.dd_,*p1_e=&dd_[rc_];
    Vec v_tmp(row_,clm_);
    for(double *p=v_tmp.dd_;p1<p1_e;p++,p1++,p2++){
        *p=*p1+*p2;
    }
    return v_tmp;
}

Vec Vec::operator-(const Vec &v) const {
    assert(row_!=v.row_&&clm_!=v.clm_);
    const double *p1=dd_, *p2=v.dd_,*p1_e=&dd_[rc_];
    Vec v_tmp(row_,clm_);
    for(double *p=v_tmp.dd_;p1<p1_e;p++,p1++,p2++){
        *p=*p1-*p2;
    }
    return v_tmp;
}

Vec Vec::operator*(double f) const {
    Vec v_tmp(row_,clm_);
    const double *p1=dd_,*p1_e=&dd_[rc_];
    for(double *p=v_tmp.dd_;p1<p1_e;p++,p1++){
        *p=*p1*f;
    }
    return v_tmp;
}

Vec& Vec::operator=(double f) {
    for(double *p=dd_,*p_e=&dd_[rc_];p<p_e;p++){*p=f;}

    return *this;
}

Vec& Vec::operator=(const double * pf) {
    for(double *p=dd_,*p_e=&dd_[rc_];p<p_e,pf++;){*p=*pf;};
    return *this;
}

Vec& Vec::operator+=(const class Vec & v) {
    assert(row_==v.row_&&clm_==v.clm_);
    const double *p1=v.dd_;
    for(double *p=dd_,*p_e=&dd_[rc_];p<p_e;p++,p1++){*p+=*p1;}
    return *this;
}

Vec& Vec::operator-=(const class Vec & v) {
    assert(row_==v.row_&&clm_==v.clm_);
    const double *p1=v.dd_;
    for(double *p=dd_,*p_e=&dd_[rc_];p<p_e;p++,p1++){*p-=*p1;}
    return *this;
}

Vec& Vec::operator*=(double f) {
    for(double *p=dd_,*p_e=&dd_[rc_];p<p_e; p++)  {*p*=f;}
    return *this;
}

Vec Vec::VecPow(const class Vec & v, int k) {
    Vec pp=v;
    double *p,*p_e=&pp.dd_[pp.rc_];
    for(int i=1;i<k;i++){
        p=pp.dd_;
        for(const double *p1=v.dd_;p<p_e;p++,p1++){
            *p*=*p1;
        }
    }
    return pp;
}

Vec Vec::VecAbs(const Vec &v) {
    Vec res(v.row_,v.clm_);
    const double *p=v.dd_,*p_e=&v.dd_[v.rc_];
    for(double *p1=res.dd_;p<p_e;p++,p1++){
        *p1=*p>0?*p:-*p;
    }
    return res;
}

void Vec::TraceVec(int p,int q) {
    int i,j;
    for(i=0;i<row_;i++){
        for(j=0;j<clm_;j++) fprintf(stderr," %*.*f",p,q,dd_[j+i*clm_]);
        fprintf(stderr,"\n");
    }
}

double VecNorm(const Vec & v) {
    const double *p=v.dd_,*p_e=&v.dd_[v.rc_];
    double f=0.0;
    for(;p<p_e;p++){
        f+=(*p)*(*p);
    }
    return sqrt(f);
}

Mat::Mat() {

}

Mat::Mat(int row0, int clm0) {
    row_=row0;clm_=clm0;
    rc_=row0*clm0;
}

Mat::Mat(int row0, int clm0, double f) {
    row_=row0;clm_=clm0;
    rc_=row0*clm0;
    for(double *p_s=dd_,*p_e=&dd_[rc_];p_s<p_e;p_s++) *p_s=f;
}

Mat::Mat(int row0, int clm0, double *pf) {
    row_=row0;clm_=clm0;rc_=row_*clm_;
    memcpy(dd_,pf,rc_* sizeof(double));
}

Mat::Mat(int row0, int clm0, vector<double> v) {
    row_=row0;clm_=clm0;rc_=row_*clm_;

    double *tp=dd_;
    for(vector<double>::iterator iter=v.begin();iter!=v.end();++iter,++tp){
        *tp=*iter;
    }
    tp=dd_;
}

Mat::~Mat() {

}

void Mat::MatClear() {
    for(double *p_s=dd_,*p_e=&dd_[rc_];p_s<p_e;p_s++) *p_s=0.0;
}

Mat Mat::operator+(const Mat &m) const {
    assert(row_==m.row_&&clm_==m.clm_);
    Mat m_tmp(row_,clm_);
    double *p=m_tmp.dd_,*p_e=&m_tmp.dd_[m_tmp.rc_];
    const double *p1=this->dd_,*p2=m.dd_;
    while(p<p_e){
        *p++=(*p1++)+(*p2++);
    }
    return m_tmp;
}

Mat Mat::operator-(const Mat &m) const {
    assert(row_==m.row_&&clm_==m.clm_);
    Mat m_tmp(row_,clm_);
    double *p=m_tmp.dd_,*p_e=&m_tmp.dd_[m_tmp.rc_];
    const double *p1=this->dd_,*p2=m.dd_;
    while(p<p_e){
        *p++=(*p1++)-(*p2++);
    }
    return m_tmp;
}

Mat Mat::operator*(double f) const {
    Mat m_tmp(row_,clm_);
    double *p=m_tmp.dd_,*p_e=&m_tmp.dd_[m_tmp.rc_];
    const double *p1=this->dd_;
    while(p<p_e){
        *p++=(*p1++)*f;
    }
    return m_tmp;
}

Mat Mat::operator*(Mat & mat) const {
    assert(this->clm_==mat.row_);
    Mat m_tmp(this->row_,mat.clm_);
    int m=this->row_,k=this->clm_,n=mat.clm_;
    double *p=m_tmp.dd_;
    const double *p1_i=this->dd_,*p2=mat.dd_;
    for(int i=0;i<m;i++,p1_i+=k){
        for(int j=0;j<n;j++){
            double f=0.0;
            const double *p1_ij=p1_i,*p1_ij_e=&p1_i[k],*p2_ji=&p2[j];
            for(;p1_ij<p1_ij_e;p1_ij++,p2_ji+=n){
                f+=(*p1_ij)*(*p2_ji);
            }
            *p++=f;
        }
    }
    return m_tmp;
}

void Mat::TraceMat(int p,int q) {
    int i,j;
    for(i=0;i<row_;i++){
        for(j=0;j<clm_;j++) fprintf(stdout," %*.*f",p,q,dd_[j+i*clm_]);
        fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
}

// math lib from rtklib and HPRTK
void MatMulPnt(const char *TraFlag,int SizeA,int SizeB,int SizeAB,double CoeAB,
                const double *MatA,const double *MatB,double CoeC,double *MatC){
    double LineAB;
    int i,j,x,f=TraFlag[0]=='N' ? (TraFlag[1]=='N' ? 1 : 2) : (TraFlag[1]=='N' ? 3 : 4);

    for (i=0; i<SizeA; i++) for (j=0; j<SizeB; j++) {
            LineAB=0.0;
            switch (f) {
                case 1: for (x=0; x<SizeAB; x++) LineAB+= MatA[i+x*SizeA]  *MatB[x+j*SizeAB]; break;
                case 2: for (x=0; x<SizeAB; x++) LineAB+= MatA[i+x*SizeA]  *MatB[j+x*SizeB]; break;
                case 3: for (x=0; x<SizeAB; x++) LineAB+= MatA[x+i*SizeAB] *MatB[x+j*SizeAB]; break;
                case 4: for (x=0; x<SizeAB; x++) LineAB+= MatA[x+i*SizeAB] *MatB[j+x*SizeB]; break;
            }
            if (CoeC==0.0) MatC[i+j*SizeA] = CoeAB*LineAB;
            else MatC[i+j*SizeA] = CoeAB*LineAB + CoeC*MatC[i+j*SizeA];
        }
}
void MatMulVec(const char *TraFlag,int SizeA,int SizeB,int SizeAB,double CoeAB,
                const vector<double> &MatA,const vector<double> &MatB,double CoeC,vector<double> &MatC){
    double LineAB;
    int i,j,x,f=TraFlag[0]=='N' ? (TraFlag[1]=='N' ? 1 : 2) : (TraFlag[1]=='N' ? 3 : 4);

    for (i=0; i<SizeA; i++) for (j=0; j<SizeB; j++) {
            LineAB=0.0;
            switch (f) {
                case 1: for (x=0; x<SizeAB; x++) LineAB+= MatA[i+x*SizeA]  *MatB[x+j*SizeAB]; break;
                case 2: for (x=0; x<SizeAB; x++) LineAB+= MatA[i+x*SizeA]  *MatB[j+x*SizeB]; break;
                case 3: for (x=0; x<SizeAB; x++) LineAB+= MatA[x+i*SizeAB] *MatB[x+j*SizeAB]; break;
                case 4: for (x=0; x<SizeAB; x++) LineAB+= MatA[x+i*SizeAB] *MatB[j+x*SizeB]; break;
            }
            if (CoeC==0.0) MatC[i+j*SizeA] = CoeAB*LineAB;
            else MatC[i+j*SizeA] = CoeAB*LineAB + CoeC*MatC[i+j*SizeA];
        }
}
void MatMulVec(const char *TraFlag,int SizeA,int SizeB,int SizeAB,double CoeAB,
                const vector<double> &MatA,const vector<double> &MatB,double CoeC,vector<double> &MatC,
                const int StartA, const int StartB, const int StartC){
    double LineAB;
    int i,j,x,f=TraFlag[0]=='N' ? (TraFlag[1]=='N' ? 1 : 2) : (TraFlag[1]=='N' ? 3 : 4);

    for (i=0; i<SizeA; i++) for (j=0; j<SizeB; j++) {
            LineAB=0.0;
            switch (f) {
                case 1: for (x=0; x<SizeAB; x++)
                        LineAB+= MatA[StartA+i+x*SizeA]  *MatB[StartB+x+j*SizeAB]; break;
                case 2: for (x=0; x<SizeAB; x++)
                        LineAB+= MatA[StartA+i+x*SizeA]  *MatB[StartB+j+x*SizeB]; break;
                case 3: for (x=0; x<SizeAB; x++)
                        LineAB+= MatA[StartA+x+i*SizeAB] *MatB[StartB+x+j*SizeAB]; break;
                case 4: for (x=0; x<SizeAB; x++)
                        LineAB+= MatA[StartA+x+i*SizeAB] *MatB[StartB+j+x*SizeB]; break;
            }
            if (CoeC==0.0) MatC[StartC+i+j*SizeA] = CoeAB*LineAB;
            else MatC[StartC+i+j*SizeA] = CoeAB*LineAB + CoeC*MatC[StartC+i+j*SizeA];
        }
}

int LUDcmp(vector<double> &MatSrc,int Size,vector<int> &index){
    double big,s,tmp;
    int i,imax=0,j,k;
    vector<double> vv(Size,0.0);

    for (i=0; i<Size; i++) {
        big=0.0;
        for (j=0; j<Size; j++) if ((tmp=fabs(MatSrc[i+j*Size]))>big) big=tmp;
        if (big>0.0) vv[i]=1.0/big;
        else { vv.clear(); return -1; }
    }
    for (j=0; j<Size; j++) {
        for (i=0; i<j; i++) {
            s=MatSrc[i+j*Size];
            for (k=0; k<i; k++) s-=MatSrc[i+k*Size]*MatSrc[k+j*Size];
            MatSrc[i+j*Size]=s;
        }
        big=0.0;
        for (i=j; i<Size; i++) {
            s=MatSrc[i+j*Size];
            for (k=0; k<j; k++) s-=MatSrc[i+k*Size]*MatSrc[k+j*Size];
            MatSrc[i+j*Size]=s;
            if ((tmp=vv[i]*fabs(s))>=big) { big=tmp; imax=i; }
        }
        if (j!=imax) {
            for (k=0; k<Size; k++) {
                tmp=MatSrc[imax+k*Size];
                MatSrc[imax+k*Size]=MatSrc[j+k*Size]; MatSrc[j+k*Size]=tmp;
            }
            vv[imax]=vv[j];
        }
        index[j]=imax;
        if (MatSrc[j+j*Size]==0.0) { vv.clear(); return -1; }
        if (j!=Size-1) {
            tmp=1.0/MatSrc[j+j*Size];
            for (i=j+1; i<Size; i++) MatSrc[i+j*Size]*=tmp;
        }
    }
    vv.clear();
    return 0;
}

void LUBkSb(vector<double> &MatB,int Size,vector<int> &index,vector<double> &MatA,int startA){
    double s;
    int i,ii=-1,ip,j;

    for (i=0; i<Size; i++) {
        ip=index[i]; s=MatA[startA+ip]; MatA[startA+ip]=MatA[startA+i];
        if (ii>=0) for (j=ii; j<i; j++) s-=MatB[i+j*Size]*MatA[startA+j];
        else if (s) ii=i;
        MatA[startA+i]=s;
    }
    for (i=Size-1; i>=0; i--) {
        s=MatA[startA+i];
        for (j=i+1; j<Size; j++) s-=MatB[i+j*Size]*MatA[startA+j];
        MatA[startA+i]=s/MatB[i+i*Size];
    }
}

int MatInv(vector<double> &MatSrc,int Size){
    vector<double> matBuf;
    vector<int> index(Size,1);

    /* initialize matBuf use MatSrc */
    matBuf.assign(MatSrc.begin(),MatSrc.end());

    if (LUDcmp(matBuf,Size,index)==-1) { matBuf.clear(); index.clear(); return -1; }
    for (int j=0; j<Size; j++){
        for (int i=0; i<Size; i++) MatSrc[i+j*Size]=0.0;
        MatSrc[j+j*Size]=1.0;
        LUBkSb(matBuf,Size,index,MatSrc,j*Size);
    }

    matBuf.clear(); index.clear();
    return 0;
}


int SolverLineEq(const string trans_flag,const vector<double>& A,const vector<double>& y,vector<double>& x,int nx,int nl){
    vector<double>A_(A.begin(),A.end());

    string trans=trans_flag[0]=='N'?"NN":"TN";
    if(MatInv(A_,nx)==0)
        MatMulVec(trans.c_str(),nx,nl,nx,1.0,A_,y,0.0,x);
    else return -1;

    A_.clear();
    return 0;
}


// Core str related function
string Doul2Str(int str_len, int dec_len, const string str_filler, const double src_num, string &dst_str){
    dst_str=to_string(src_num); /* with 6 decimal digit */
    if(dec_len>0)dst_str=dst_str.substr(0,dst_str.length()-6+dec_len); /* decimal digits */
    else dst_str=dst_str.substr(0,dst_str.length()-7);
    while (str_len>dst_str.length())
        dst_str=str_filler+dst_str;
    return dst_str;
}

string Int2Str(int str_len, const string str_filler, const int src_num, string &dst_str){
    dst_str=to_string(src_num);
    while(str_len>dst_str.length())
        dst_str=str_filler+dst_str;
    return dst_str;
}

int Str2Double(string src_str,double &dst_num){
    int i,fnum;

    if (src_str.length()<=0) return 0;

    for (i=0,fnum=1; i<src_str.length(); i++){
        if (fnum && src_str[i]<='9' && src_str[i]>='0') fnum=0;
        if (fnum==0 && src_str[i]=='D') { src_str[i]='E'; break; }
    }
    if (fnum) return 0;

    dst_num=stod(src_str);

    return 1;
}

int Str2Int(const string src_str,int &dst_num){
    int i;
    for (i=0; i<src_str.length(); i++){
        if (src_str[i]<='9'&&src_str[i]>='0') break;
    }
    if (i>=src_str.length()) return 0;

    dst_num=stoi(src_str);

    return 1;
}

void SplitString(const string& s, vector<string>& v,string c)
{
    string::size_type pos1,pos2;
    pos2=s.find(c);
    pos1=0;
    while(string::npos!=pos2){
        v.push_back(s.substr(pos1,pos2-pos1));
        pos1=pos2+c.size();
        pos2=s.find(c,pos1);
    }
    if(pos1!=s.length())
        v.push_back(s.substr(pos1));
}

vector<string> MultiSplitStr(const string &s, const string &seperator) {
    vector<string> result;
    typedef string::size_type string_size;
    string_size i = 0;

    while(i!=s.size()){
        int flag = 0;
        while(i!=s.size()&&flag==0){
            flag = 1;
            for(string_size x =0;x<seperator.size();++x){
                if(s[i]==seperator[x]){
                    ++i;flag=0;break;
                }
            }
        }
        flag=0;string_size j=i;
        while(j!=s.size()&&flag==0){
            for(string_size x=0;x<seperator.size();++x){
                if(s[j]==seperator[x]){
                    flag=1;break;
                }
            }
            if(flag==0) ++j;
        }
        if(i!=j){
            result.push_back(s.substr(i,j-i));i=j;
        }
    }
    return result;
}

vector<string> TextSplit(const string &in, const string &delim){
    std::vector<std::string> ret;
    try
    {
        std::regex re{delim};
        return std::vector<std::string>{std::sregex_token_iterator(in.begin(), in.end(), re, -1), std::sregex_token_iterator()};
    }
    catch (const std::exception &e)
    {
        std::cout << "error:" << e.what() << std::endl;
    }
    return ret;
}

string &StrTrim(string &s){
    if(s.empty()) return s;
    s.erase(0,s.find_first_not_of(" "));
    s.erase(s.find_last_not_of(" ")+1);
    return s;
}

extern void CreateDir(const char *path)
{
    char buff[1024],*p;

    strcpy(buff,path);
    if (!(p=strrchr(buff,FILEPATHSEP))) return;
    *p='\0';

#ifdef WIN32
    CreateDirectory(buff,NULL);
#else
    mkdir(buff,0777);
#endif
}

extern int ExPath(char *path,char *paths[],int nmax) {
    int i,j,n=0;
    char tmp[1024];
#ifdef WIN32
    WIN32_FIND_DATA file;
    HANDLE h;
    char dir[1024]="",*p;

    if ((p=strrchr(path,'\\'))) {
        strncpy(dir,path,p-path+1); dir[p-path+1]='\0';
    }
    if ((h=FindFirstFile((LPCTSTR)path,&file))==INVALID_HANDLE_VALUE) {
        strcpy(paths[0],path);
        return 1;
    }
    sprintf(paths[n++],"%s%s",dir,file.cFileName);
    while (FindNextFile(h,&file)&&n<nmax) {
        if (file.dwFileAttributes&FILE_ATTRIBUTE_DIRECTORY) continue;
        sprintf(paths[n++],"%s%s",dir,file.cFileName);
    }
    FindClose(h);
#else
    struct dirent *d;
    DIR *dp;
    const char *file=path;
    char dir[1024]="",s1[1024],s2[1024],*p,*q,*r;

    if ((p=strrchr(path,'/'))||(p=strrchr(path,'\\'))) {
        file=p+1; strncpy(dir,path,p-path+1); dir[p-path+1]='\0';
    }
    if (!(dp=opendir(*dir?dir:"."))) return 0;
    while ((d=readdir(dp))) {
        if (*(d->d_name)=='.') continue;
        sprintf(s1,"^%s$",d->d_name);
        sprintf(s2,"^%s$",file);
        for (p=s1;*p;p++) *p=(char)tolower((int)*p);
        for (p=s2;*p;p++) *p=(char)tolower((int)*p);

        for (p=s1,q=strtok_r(s2,"*",&r);q;q=strtok_r(NULL,"*",&r)) {
            if ((p=strstr(p,q))) p+=strlen(q); else break;
        }
        if (p&&n<nmax) sprintf(paths[n++],"%s%s",dir,d->d_name);
    }
    closedir(dp);
#endif
    /* sort paths in alphabetical order */
    for (i=0;i<n-1;i++) {
        for (j=i+1;j<n;j++) {
            if (strcmp(paths[i],paths[j])>0) {
                strcpy(tmp,paths[i]);
                strcpy(paths[i],paths[j]);
                strcpy(paths[j],tmp);
            }
        }
    }
    return n;
}

void AstArgs(double t, double *f) {
    static const double fc[][5]={ /* coefficients for iau 1980 nutation */
            { 134.96340251, 1717915923.2178,  31.8792,  0.051635, -0.00024470 },
            { 357.52910918,  129596581.0481,  -0.5532,  0.000136, -0.00001149 },
            { 93.27209062, 1739527262.8478, -12.7512, -0.001037,  0.00000417 },
            { 297.85019547, 1602961601.2090,  -6.3706,  0.006593, -0.00003169 },
            { 125.04455501,   -6962890.2665,   7.4722,  0.007702, -0.00005939 }
    };
    double tt[4];
    int i,j;

    for (tt[0]=t,i=1; i<4; i++) tt[i]=tt[i-1]*t;
    for (i=0; i<5; i++) {
        f[i]=fc[i][0]*3600.0;
        for (j=0; j<4; j++) f[i]+=fc[i][j+1]*tt[j];
        f[i]=fmod(f[i]*AS2R,2.0*PI);
    }
}

static void nut_iau1980(double t,const double *f,double *dpsi,double *deps)
{
    static const double nut[106][10]={
            { 0,   0,   0,   0,   1, -6798.4, -171996, -174.2, 92025,   8.9 },
            { 0,   0,   2,  -2,   2,   182.6,  -13187,   -1.6,  5736,  -3.1 },
            { 0,   0,   2,   0,   2,    13.7,   -2274,   -0.2,   977,  -0.5 },
            { 0,   0,   0,   0,   2, -3399.2,    2062,    0.2,  -895,   0.5 },
            { 0,  -1,   0,   0,   0,  -365.3,   -1426,    3.4,    54,  -0.1 },
            { 1,   0,   0,   0,   0,    27.6,     712,    0.1,    -7,   0.0 },
            { 0,   1,   2,  -2,   2,   121.7,    -517,    1.2,   224,  -0.6 },
            { 0,   0,   2,   0,   1,    13.6,    -386,   -0.4,   200,   0.0 },
            { 1,   0,   2,   0,   2,     9.1,    -301,    0.0,   129,  -0.1 },
            { 0,  -1,   2,  -2,   2,   365.2,     217,   -0.5,   -95,   0.3 },
            { -1,   0,   0,   2,   0,    31.8,     158,    0.0,    -1,   0.0 },
            { 0,   0,   2,  -2,   1,   177.8,     129,    0.1,   -70,   0.0 },
            { -1,   0,   2,   0,   2,    27.1,     123,    0.0,   -53,   0.0 },
            { 1,   0,   0,   0,   1,    27.7,      63,    0.1,   -33,   0.0 },
            { 0,   0,   0,   2,   0,    14.8,      63,    0.0,    -2,   0.0 },
            { -1,   0,   2,   2,   2,     9.6,     -59,    0.0,    26,   0.0 },
            { -1,   0,   0,   0,   1,   -27.4,     -58,   -0.1,    32,   0.0 },
            { 1,   0,   2,   0,   1,     9.1,     -51,    0.0,    27,   0.0 },
            { -2,   0,   0,   2,   0,  -205.9,     -48,    0.0,     1,   0.0 },
            { -2,   0,   2,   0,   1,  1305.5,      46,    0.0,   -24,   0.0 },
            { 0,   0,   2,   2,   2,     7.1,     -38,    0.0,    16,   0.0 },
            { 2,   0,   2,   0,   2,     6.9,     -31,    0.0,    13,   0.0 },
            { 2,   0,   0,   0,   0,    13.8,      29,    0.0,    -1,   0.0 },
            { 1,   0,   2,  -2,   2,    23.9,      29,    0.0,   -12,   0.0 },
            { 0,   0,   2,   0,   0,    13.6,      26,    0.0,    -1,   0.0 },
            { 0,   0,   2,  -2,   0,   173.3,     -22,    0.0,     0,   0.0 },
            { -1,   0,   2,   0,   1,    27.0,      21,    0.0,   -10,   0.0 },
            { 0,   2,   0,   0,   0,   182.6,      17,   -0.1,     0,   0.0 },
            { 0,   2,   2,  -2,   2,    91.3,     -16,    0.1,     7,   0.0 },
            { -1,   0,   0,   2,   1,    32.0,      16,    0.0,    -8,   0.0 },
            { 0,   1,   0,   0,   1,   386.0,     -15,    0.0,     9,   0.0 },
            { 1,   0,   0,  -2,   1,   -31.7,     -13,    0.0,     7,   0.0 },
            { 0,  -1,   0,   0,   1,  -346.6,     -12,    0.0,     6,   0.0 },
            { 2,   0,  -2,   0,   0, -1095.2,      11,    0.0,     0,   0.0 },
            { -1,   0,   2,   2,   1,     9.5,     -10,    0.0,     5,   0.0 },
            { 1,   0,   2,   2,   2,     5.6,      -8,    0.0,     3,   0.0 },
            { 0,  -1,   2,   0,   2,    14.2,      -7,    0.0,     3,   0.0 },
            { 0,   0,   2,   2,   1,     7.1,      -7,    0.0,     3,   0.0 },
            { 1,   1,   0,  -2,   0,   -34.8,      -7,    0.0,     0,   0.0 },
            { 0,   1,   2,   0,   2,    13.2,       7,    0.0,    -3,   0.0 },
            { -2,   0,   0,   2,   1,  -199.8,      -6,    0.0,     3,   0.0 },
            { 0,   0,   0,   2,   1,    14.8,      -6,    0.0,     3,   0.0 },
            { 2,   0,   2,  -2,   2,    12.8,       6,    0.0,    -3,   0.0 },
            { 1,   0,   0,   2,   0,     9.6,       6,    0.0,     0,   0.0 },
            { 1,   0,   2,  -2,   1,    23.9,       6,    0.0,    -3,   0.0 },
            { 0,   0,   0,  -2,   1,   -14.7,      -5,    0.0,     3,   0.0 },
            { 0,  -1,   2,  -2,   1,   346.6,      -5,    0.0,     3,   0.0 },
            { 2,   0,   2,   0,   1,     6.9,      -5,    0.0,     3,   0.0 },
            { 1,  -1,   0,   0,   0,    29.8,       5,    0.0,     0,   0.0 },
            { 1,   0,   0,  -1,   0,   411.8,      -4,    0.0,     0,   0.0 },
            { 0,   0,   0,   1,   0,    29.5,      -4,    0.0,     0,   0.0 },
            { 0,   1,   0,  -2,   0,   -15.4,      -4,    0.0,     0,   0.0 },
            { 1,   0,  -2,   0,   0,   -26.9,       4,    0.0,     0,   0.0 },
            { 2,   0,   0,  -2,   1,   212.3,       4,    0.0,    -2,   0.0 },
            { 0,   1,   2,  -2,   1,   119.6,       4,    0.0,    -2,   0.0 },
            { 1,   1,   0,   0,   0,    25.6,      -3,    0.0,     0,   0.0 },
            { 1,  -1,   0,  -1,   0, -3232.9,      -3,    0.0,     0,   0.0 },
            { -1,  -1,   2,   2,   2,     9.8,      -3,    0.0,     1,   0.0 },
            { 0,  -1,   2,   2,   2,     7.2,      -3,    0.0,     1,   0.0 },
            { 1,  -1,   2,   0,   2,     9.4,      -3,    0.0,     1,   0.0 },
            { 3,   0,   2,   0,   2,     5.5,      -3,    0.0,     1,   0.0 },
            { -2,   0,   2,   0,   2,  1615.7,      -3,    0.0,     1,   0.0 },
            { 1,   0,   2,   0,   0,     9.1,       3,    0.0,     0,   0.0 },
            { -1,   0,   2,   4,   2,     5.8,      -2,    0.0,     1,   0.0 },
            { 1,   0,   0,   0,   2,    27.8,      -2,    0.0,     1,   0.0 },
            { -1,   0,   2,  -2,   1,   -32.6,      -2,    0.0,     1,   0.0 },
            { 0,  -2,   2,  -2,   1,  6786.3,      -2,    0.0,     1,   0.0 },
            { -2,   0,   0,   0,   1,   -13.7,      -2,    0.0,     1,   0.0 },
            { 2,   0,   0,   0,   1,    13.8,       2,    0.0,    -1,   0.0 },
            { 3,   0,   0,   0,   0,     9.2,       2,    0.0,     0,   0.0 },
            { 1,   1,   2,   0,   2,     8.9,       2,    0.0,    -1,   0.0 },
            { 0,   0,   2,   1,   2,     9.3,       2,    0.0,    -1,   0.0 },
            { 1,   0,   0,   2,   1,     9.6,      -1,    0.0,     0,   0.0 },
            { 1,   0,   2,   2,   1,     5.6,      -1,    0.0,     1,   0.0 },
            { 1,   1,   0,  -2,   1,   -34.7,      -1,    0.0,     0,   0.0 },
            { 0,   1,   0,   2,   0,    14.2,      -1,    0.0,     0,   0.0 },
            { 0,   1,   2,  -2,   0,   117.5,      -1,    0.0,     0,   0.0 },
            { 0,   1,  -2,   2,   0,  -329.8,      -1,    0.0,     0,   0.0 },
            { 1,   0,  -2,   2,   0,    23.8,      -1,    0.0,     0,   0.0 },
            { 1,   0,  -2,  -2,   0,    -9.5,      -1,    0.0,     0,   0.0 },
            { 1,   0,   2,  -2,   0,    32.8,      -1,    0.0,     0,   0.0 },
            { 1,   0,   0,  -4,   0,   -10.1,      -1,    0.0,     0,   0.0 },
            { 2,   0,   0,  -4,   0,   -15.9,      -1,    0.0,     0,   0.0 },
            { 0,   0,   2,   4,   2,     4.8,      -1,    0.0,     0,   0.0 },
            { 0,   0,   2,  -1,   2,    25.4,      -1,    0.0,     0,   0.0 },
            { -2,   0,   2,   4,   2,     7.3,      -1,    0.0,     1,   0.0 },
            { 2,   0,   2,   2,   2,     4.7,      -1,    0.0,     0,   0.0 },
            { 0,  -1,   2,   0,   1,    14.2,      -1,    0.0,     0,   0.0 },
            { 0,   0,  -2,   0,   1,   -13.6,      -1,    0.0,     0,   0.0 },
            { 0,   0,   4,  -2,   2,    12.7,       1,    0.0,     0,   0.0 },
            { 0,   1,   0,   0,   2,   409.2,       1,    0.0,     0,   0.0 },
            { 1,   1,   2,  -2,   2,    22.5,       1,    0.0,    -1,   0.0 },
            { 3,   0,   2,  -2,   2,     8.7,       1,    0.0,     0,   0.0 },
            { -2,   0,   2,   2,   2,    14.6,       1,    0.0,    -1,   0.0 },
            { -1,   0,   0,   0,   2,   -27.3,       1,    0.0,    -1,   0.0 },
            { 0,   0,  -2,   2,   1,  -169.0,       1,    0.0,     0,   0.0 },
            { 0,   1,   2,   0,   1,    13.1,       1,    0.0,     0,   0.0 },
            { -1,   0,   4,   0,   2,     9.1,       1,    0.0,     0,   0.0 },
            { 2,   1,   0,  -2,   0,   131.7,       1,    0.0,     0,   0.0 },
            { 2,   0,   0,   2,   0,     7.1,       1,    0.0,     0,   0.0 },
            { 2,   0,   2,  -2,   1,    12.8,       1,    0.0,    -1,   0.0 },
            { 2,   0,  -2,   0,   1,  -943.2,       1,    0.0,     0,   0.0 },
            { 1,  -1,   0,  -2,   0,   -29.3,       1,    0.0,     0,   0.0 },
            { -1,   0,   0,   1,   1,  -388.3,       1,    0.0,     0,   0.0 },
            { -1,  -1,   0,   2,   1,    35.0,       1,    0.0,     0,   0.0 },
            { 0,   1,   0,   1,   0,    27.3,       1,    0.0,     0,   0.0 }
    };
    double ang;
    int i,j;

    *dpsi=*deps=0.0;

    for (i=0; i<106; i++) {
        ang=0.0;
        for (j=0; j<5; j++) ang+=nut[i][j]*f[j];
        *dpsi+=(nut[i][6]+nut[i][7]*t)*sin(ang);
        *deps+=(nut[i][8]+nut[i][9]*t)*cos(ang);
    }
    *dpsi*=1E-4*AS2R; /* 0.1 mas -> rad */
    *deps*=1E-4*AS2R;
}
#define Rx(t,X) do { \
    (X)[0]=1.0; (X)[1]=(X)[2]=(X)[3]=(X)[6]=0.0; \
    (X)[4]=(X)[8]=cos(t); (X)[7]=sin(t); (X)[5]=-(X)[7]; \
} while (0)

#define Ry(t,X) do { \
    (X)[4]=1.0; (X)[1]=(X)[3]=(X)[5]=(X)[7]=0.0; \
    (X)[0]=(X)[8]=cos(t); (X)[2]=sin(t); (X)[6]=-(X)[2]; \
} while (0)

#define Rz(t,X) do { \
    (X)[8]=1.0; (X)[2]=(X)[5]=(X)[6]=(X)[7]=0.0; \
    (X)[0]=(X)[4]=cos(t); (X)[3]=sin(t); (X)[1]=-(X)[3]; \
} while (0)

void Eci2Ecef(Time tutc,const double *erpv,double *U,double *gmst) {
    const double ep2000[]={ 2000,1,1,12,0,0 };
    static Time tutc_;
    static double U_[9],gmst_;
    Time tgps,baseTime;
    double eps,ze,th,z,t,t2,t3,dpsi,deps,gast,f[5];
    double R1[9],R2[9],R3[9],R[9],W[9],N[9],P[9],NP[9];
    int i;

    if (fabs(tutc.TimeDiff(tutc_))<0.01) { /* read cache */
        for (i=0; i<9; i++) U[i]=U_[i];
        if (gmst) *gmst=gmst_;
        return;
    }
    tutc_=tutc;

    /* terrestrial time */
    tgps=tutc_;
    tgps.UTC2GPST();
    baseTime.Epoch2Time(ep2000);
    t=(tgps.TimeDiff(baseTime)+19.0+32.184)/86400.0/36525.0;
    t2=t*t; t3=t2*t;

    /* astronomical arguments */
    AstArgs(t,f);

    /* iau 1976 precession */
    ze=(2306.2181*t+0.30188*t2+0.017998*t3)*AS2R;
    th=(2004.3109*t-0.42665*t2-0.041833*t3)*AS2R;
    z =(2306.2181*t+1.09468*t2+0.018203*t3)*AS2R;
    eps=(84381.448-46.8150*t-0.00059*t2+0.001813*t3)*AS2R;
    Rz(-z,R1); Ry(th,R2); Rz(-ze,R3);
    MatMulPnt("NN",3,3,3,1.0,R1,R2,0.0,R);
    MatMulPnt("NN",3,3,3,1.0,R,R3,0.0,P); /* P=Rz(-z)*Ry(th)*Rz(-ze) */

    /* iau 1980 nutation */
    nut_iau1980(t,f,&dpsi,&deps);
    Rx(-eps-deps,R1); Rz(-dpsi,R2); Rx(eps,R3);
    MatMulPnt("NN",3,3,3,1.0,R1,R2,0.0,R);
    MatMulPnt("NN",3,3,3,1.0,R,R3,0.0,N); /* N=Rx(-eps)*Rz(-dspi)*Rx(eps) */

    /* greenwich aparent sidereal time (rad) */
    gmst_=tutc_.UTC2Gmst(erpv[2]);
    gast=gmst_+dpsi*cos(eps);
    gast+=(0.00264*sin(f[4])+0.000063*sin(2.0*f[4]))*AS2R;

    /* eci to ecef transformation matrix */
    Ry(-erpv[0],R1); Rx(-erpv[1],R2); Rz(gast,R3);
    MatMulPnt("NN",3,3,3,1.0,R1,R2,0.0,W);
    MatMulPnt("NN",3,3,3,1.0,W,R3,0.0,R); /* W=Ry(-xp)*Rx(-yp) */
    MatMulPnt("NN",3,3,3,1.0,N,P,0.0,NP);
    MatMulPnt("NN",3,3,3,1.0,R,NP,0.0,U_); /* U=W*Rz(gast)*N*P */

    for (i=0; i<9; i++) U[i]=U_[i];
    if (gmst) *gmst=gmst_;
}

static void SunMoonPosEci(Time tut, double *rsun, double *rmoon) {
    const double ep2000[]={ 2000,1,1,12,0,0 };
    double t,f[5],eps,Ms,ls,rs,lm,pm,rm,sine,cose,sinp,cosp,sinl,cosl;
    Time baseTime;

    baseTime.Epoch2Time(ep2000);
    t=tut.TimeDiff(baseTime)/86400.0/36525.0;

    /* astronomical arguments */
    AstArgs(t,f);

    /* obliquity of the ecliptic */
    eps=23.439291-0.0130042*t;
    sine=sin(eps*D2R); cose=cos(eps*D2R);

    /* sun position in eci */
    if (rsun) {
        Ms=357.5277233+35999.05034*t;
        ls=280.460+36000.770*t+1.914666471*sin(Ms*D2R)+0.019994643*sin(2.0*Ms*D2R);
        rs=AU*(1.000140612-0.016708617*cos(Ms*D2R)-0.000139589*cos(2.0*Ms*D2R));
        sinl=sin(ls*D2R); cosl=cos(ls*D2R);
        rsun[0]=rs*cosl;
        rsun[1]=rs*cose*sinl;
        rsun[2]=rs*sine*sinl;
    }
    /* moon position in eci */
    if (rmoon) {
        lm=218.32+481267.883*t+6.29*sin(f[0])-1.27*sin(f[0]-2.0*f[3])+
           0.66*sin(2.0*f[3])+0.21*sin(2.0*f[0])-0.19*sin(f[1])-0.11*sin(2.0*f[2]);
        pm=5.13*sin(f[2])+0.28*sin(f[0]+f[2])-0.28*sin(f[2]-f[0])-
           0.17*sin(f[2]-2.0*f[3]);
        rm=RE_WGS84/sin((0.9508+0.0518*cos(f[0])+0.0095*cos(f[0]-2.0*f[3])+
                         0.0078*cos(2.0*f[3])+0.0028*cos(2.0*f[0]))*D2R);
        sinl=sin(lm*D2R); cosl=cos(lm*D2R);
        sinp=sin(pm*D2R); cosp=cos(pm*D2R);
        rmoon[0]=rm*cosp*cosl;
        rmoon[1]=rm*(cose*cosp*sinl-sine*sinp);
        rmoon[2]=rm*(sine*cosp*sinl+cose*sinp);
    }
}

double SunMoonPos(Time ut1t,const double *erp_val,double *sun_pos,double *moon_pos) {
    double rs[3],rm[3],U[9],gmst0;
    ut1t.TimeAdd(erp_val[2]); //utc->ut1

    SunMoonPosEci(ut1t,sun_pos?rs:NULL,moon_pos?rm:NULL);
    Eci2Ecef(ut1t,erp_val,U,&gmst0);

    if(sun_pos) MatMulPnt("NN",3,1,3,1.0,U,rs,0.0,sun_pos);
    if(moon_pos) MatMulPnt("NN",3,1,3,1.0,U,rm,0.0,moon_pos);
    return gmst0;
}


