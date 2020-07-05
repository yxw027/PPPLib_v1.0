//
// Created by cc on 4/9/20.
//

#ifndef MSNAVS_CMNFUNC_H
#define MSNAVS_CMNFUNC_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <ctime>
#include <cstdarg>
#include <cstdlib>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <vector>
#include <cctype>
#include <assert.h>
#include <regex>
#ifdef WIN32
#include <winsock2.h>
#include <windows.h>
#include <unistd.h>
#else
#include <dirent.h>
#include <pthread.h>
#include <termio.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/types.h>
#endif

#ifdef WIN32
#define thread_t    HANDLE
#define lock_t      CRITICAL_SECTION
#define initlock(f) InitializeCriticalSection(f)
#define lock(f)     EnterCriticalSection(f)
#define unlock(f)   LeaveCriticalSection(f)
#define FILEPATHSEP '\\'
#define DEFPRC_DIR  "E:\\PhdWorks\\DataProcess\\data"
#define LOGINI_PATH "..\\conf\\log.ini"
#else
#define thread_t    pthread_t
#define lock_t      pthread_mutex_t
#define initlock(f) pthread_mutex_init(f,NULL)
#define lock(f)     pthread_mutex_lock(f)
#define unlock(f)   pthread_mutex_unlock(f)
#define FILEPATHSEP '/'
//#define DEFPRC_DIR  "/mnt/hgfs/DataProcess/data"
#define DEFPRC_DIR "/home/cc/test_data"
#define LOGINI_PATH "/home/cc/Codings/PPPLib/conf/log.ini"
#endif

// log settings
//#define  ELPP_DISABLE_LOGS
//#define ELPP_DISABLE_WARNING_LOGS
//#define ELPP_NO_DEFAULT_LOG_FILE
#define  ELPP_STL_LOGGING
#define  ELPP_FEATURE_PERFORMANCE_TRACKING
#include "easylogging++.h"

// constant
#define PI		3.14159265358979
#define PI_2	(PI/2.0)
#define PI_4	(PI/4.0)
#define _2PI	(2.0*PI)
#define CLIGHT  299792458.0         /* speed of light (m/s) */
#define D2R     (PI/180.0)          /* deg to rad */
#define R2D     (180.0/PI)          /* rad to deg */
#define AS2R    (D2R/3600.0)        /* arc sec to radian */
#define AU      149597870691.0      /* 1 AU (m) */
#define MIN(x,y) ((x)<=(y)?(x):(y))
#define MAX(x,y) ((x)>=(y)?(x):(y))
#define SQR(x)  ((x)*(x))
#define SQRT(x) ((x)<0.0||(x)!=(x)?0.0:sqrt(x))
#define DTTOL       0.005               /* tolerance of time difference (s) */
#define MAXSTRPATH  1024
#define MAXBUFF     4096
using std::cin;
using std::cout;
using std::endl;
using std::setw;
using std::setprecision;
using std::istringstream;
using std::ostringstream;
using std::ifstream;
using std::fstream;
using std::ios;
using std::string;
using std::to_string;
using std::vector;
using std::ofstream;

typedef struct { // seconds of day
    long sn;
    double tos;
}Sod_t;

typedef struct { //Julian day
    long day;
    Sod_t sod;
}Mjd_t;

// MANAVS base time class
class Time {
public:
    Time();
    Time(const double *epoch);
    ~Time();

public:
    Time* Epoch2Time(const double *ep);
    void Time2Epoch();
    int Str2Time(string s);
    string Time2Str(int n);
    Time* TimeAdd(double sec);
    double TimeDiff(Time t2) const;
    double Time2Sec(Time& day);
    double  Time2Doy();
    Time* GPST2Time(int week,double sss);
    double Time2GPST(int *week);
    Time* UTC2GPST();
    Time* GPST2UTC();
    Time* GPST2BDT();
    Time* BDT2GPST();
    Time* BDT2Time(int week,double wos);
    double Time2BDT(int *week);
    int Monitor(Time ts,Time te,double tint);
    Time* AdjWeek(Time t);
    Time* AdjDay(Time t);
    Time* CopyTime(Time t);
    void Time2Mjd();
    double UTC2Gmst(double ut1_utc);

public:
    time_t time_;
    double sec_;
    double epoch_[6];
    string  time_str_;
    int doy_;
    int time_sys_;
    Mjd_t mjd_;
};

// Core math related class
int Round(double d);
double GetMedian(const double *array,int n);

class Vec3;
class Mat3;
class Vec;
class Mat;

#define MAT_MAX_DIM   200
#define MAT_MAX_DIM2  (MAT_MAX_DIM*MAT_MAX_DIM)
#define RE_WGS84    6378137.0           /* earth semimajor axis (WGS84) (m) */
#define EE_WGS84    (1.0/298.257223563) /* earth eccentricity (WGS84) */
#define RE_CGCS2000 6378137.0			/* earth semimajor axis (CGCS2000) (m) */
#define EE_CGCS2000 (1.0/298.257222101) /* earth eccentricity (CGCS2000) */
enum COOR {CGCS2000,WGS84};
class Vec3 {
public:
    Vec3(void);
    ~Vec3();
    Vec3(double xyz);
    Vec3(double x,double y,double z);

public:
    Vec3 operator+(const Vec3& v) const;       // vector addition              // const 表示成员函数隐含传入的this指针为const指针，决定了在该成员函数中，任意修改它所在的类的成员的操作都是不允许的
    Vec3 operator-(const Vec3& v) const;       // vector subtraction
    Vec3 operator*(const Vec3& v) const;       // vector cross multiplication
    Vec3 operator*(double f) const;            // vector multiply scale
    Vec3 operator*(const Mat3& m) const;
    Vec3 operator/(double f) const;
    Vec3 operator/(const Vec3& v) const;
    Vec3& operator+=(const Vec3& v);
    Vec3& operator-=(const Vec3& v);
    Vec3& operator*=(double f);
    Vec3& operator/=(double f);
    Vec3& operator/=(const Vec3& v);

    friend double Vec3Norm(const Vec3& v);
    friend double Vec3NormXY(const Vec3& v);
    friend double VecDot(const Vec3& v1,const Vec3& v2);
    friend double VecDot2(const Vec3& v1,const Vec3& v2);
    friend Vec3 VecPow(const Vec3& v,int k=2);

    friend Vec3 Xyz2Blh(const Vec3& xyz,int coor);
    friend Vec3 Blh2Xyz(const Vec3& blh,int coor);
    friend Mat3 Pos2Cen(const Vec3& pos);
    friend Vec3 Xyz2Enu(const Vec3& xyz,const Vec3& pos);

public:
    double i_,j_,k_;
};

extern const Vec3 I31, O31;


class Mat3 {
public:
    Mat3(void);
    ~Mat3();
    Mat3(double xx,double xy,double xz,
          double yx,double yy,double yz,
          double zx,double zy,double zz);
    Mat3(const Vec3 &v0,const Vec3 &v1,const Vec3 &v2);

public:
    Mat3  operator+(const Mat3 &m) const;     // matrix addition
    Mat3  operator-(const Mat3 &m) const;     // matrix subtraction
    Mat3  operator*(const Mat3 &m) const;     // matrix multiplication
    Mat3  operator*(double f) const;          // matrix multiply scale
    Vec3  operator*(const Vec3 &v) const;     // matrix multiply vector
    Mat3& operator+=(const Mat3 &m);
    Mat3& operator-=(const Mat3 &m);

    friend Mat3 operator-(const Mat3& m);
    friend Mat3 operator~(const Mat3& m);
    friend Mat3 operator*(double f,const Mat3& m);
    friend Mat3 Mat3Pow(const Mat3& m, int k);
    friend double Mat3Trace(const Mat3& m);
    friend double Mat3Det(const Mat3& m);
    friend Vec3 Mat3DigV(const Mat3& m);
    friend Mat3 Mat3DigM(const Mat3& m);
    friend Mat3 Mat3Adj(const Mat3& m);
    friend Mat3 Mat3Inv(const Mat3& m);

public:
    double e00_,e01_,e02_;
    double e10_,e11_,e12_;
    double e20_,e21_,e22_;
};

class Vec {
public:
    Vec();
    ~Vec();
    Vec(int row,int clm=1);
    Vec(int row,double f);
    Vec(int row,const double *pf);
    Vec(const Vec3& v);
    Vec(const Vec3& v1,const Vec3& v2);
    Vec(int row,vector<double>v1);

public:
    Vec  operator+(const Vec &v) const;
    Vec  operator-(const Vec &v) const;
    Vec  operator*(double f) const;
    Vec& operator=(double f);
    Vec& operator=(const double *pf);
    Vec& operator+=(const Vec& v);
    Vec& operator-=(const Vec& v);
    Vec& operator*=(double f);
    Vec& operator*(const Mat& m) const;
    Mat& operator*(const Vec &v) const;
    double& operator()(int r);
//    Vec operator~(const Vec& v);
    Vec VecAbs(const Vec& v);
    friend double VecNorm(const Vec &v);
    Vec VecPow(const Vec& v,int k=2);

public:
    void TraceVec(int p,int q);

public:
    int row_,clm_,rc_;
    double dd_[MAT_MAX_DIM];
};

class Mat {
public:
    Mat(void);
    ~Mat();
    Mat(int row0,int clm0);
    Mat(int row0,int clm0,double f);
    Mat(int row0,int clm0,double *pf);
    Mat(int row0,int clm0,vector<double> m);

    void MatClear(void);
    void SetDiag(double f,...);
    Mat operator+(const Mat &m) const;
    Mat operator-(const Mat &m) const;
    Mat operator*(double f) const;
    Mat operator*(Mat &mat) const;

public:
    void TraceMat(int p,int q);
public:
    int row_,clm_,rc_;
    double dd_[MAT_MAX_DIM2];
};

// math lib from rtklib and HPRTK
template <typename Iter1,typename Iter2>
double Dot(const Iter1 VecA,const Iter2 VecB,int SizeVec){
    double dInn=0.0;

    while (--SizeVec>=0){
        dInn+=VecA[SizeVec]*VecB[SizeVec];
    }
    return dInn;
}

template <typename Iter>
double Norm(const Iter VecA,int SizeVec){
    return sqrt(Dot(VecA,VecA,SizeVec));
}

template <typename Iter1,typename Iter2>
int NormV3(const Iter1 vec1,Iter2 vec2){
   double r;
   if((r=Norm(vec1,3))<=0.0) return 0;
   vec2[0]=vec1[0]/r;
   vec2[1]=vec1[1]/r;
   vec2[2]=vec1[2]/r;
   return 1;
}

template <typename Iter1,typename Iter2, typename Iter3>
void CrossVec3(const Iter1 vec1,const Iter2 vec2,Iter3 vec3){
    vec3[0]=vec1[1]*vec2[2]-vec1[2]*vec2[1];
    vec3[1]=vec1[2]*vec2[0]-vec1[0]*vec2[2];
    vec3[2]=vec1[0]*vec2[1]-vec1[1]*vec2[0];
}

template <typename Iter>
void EyeMat(Iter A,const int n){
    for(int i=0;i<n;i++){     //row
        for(int j=0;j<n;j++){ //clm
            A[j+i*n]=i==j?1.0:0.0;
        }
    }
}


void MatMulPnt(const char *Tra_Flag,int sizeA,int SizeB,int SizeAB,double CoeAB,
               const double *MatA,const double *MatB,double CoeC,double *MatC);
void MatMulVec(const char *TraFlag,int SizeA,int SizeB,int SizeAB,double CoeAB,
               const vector<double> &MatA,const vector<double> &MatB,double CoeC,vector<double> &MatC);
void MatMulVec(const char *TraFlag,int SizeA,int SizeB,int SizeAB,double CoeAB,
               const vector<double> &MatA,const vector<double> &MatB,double CoeC,vector<double> &MatC,
               const int StartA, const int StartB, const int StartC);
template <typename Iter1,typename Iter2,typename Iter3>
void MatMul(const char *TraFlag,int SizeA,int SizeB,int SizeAB,double CoeAB,
            const Iter1 MatA,const Iter2 MatB,double CoeC,Iter3 MatC) {
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
int LUDcmp(vector<double> &MatSrc,int Size,vector<int> &index);
void LUBkSb(vector<double> &MatB,int Size,vector<int> &index,vector<double> &MatA,int startA);
int MatInv(vector<double> &MatSrc,int Size);

template <typename Iter>
double MatTrace(const Iter A,const int size){
    double tra=0.0;
    for(int i=0;i<size;i++)
        tra+=A[i+i*size];

    return tra;
}
int SolverLineEq(const string trans_flag,const vector<double>& A,const vector<double>& y,vector<double>& x,int nx,int nl);

// MSNAV str related function
string Doul2Str(int str_len, int dec_len, const string str_filler, const double src_num, string &dst_str);
string Int2Str(int str_len, const string str_filler, const int src_num, string &dst_str);
int Str2Double(string src_str,double &dst_num);
int Str2Int(const string src_str,int &dst_num);
void SplitString(const string& s, vector<string>& v,string c);
vector<string> MultiSplitStr(const string &s, const string &seperator);
vector<string> TextSplit(const string &in, const string &delim);
string &StrTrim(string &s);

extern void CreateDir(const char *path);
extern int ExPath(char *path,char *paths[],int nmax);

void AstArgs(double t, double *f);
void Eci2Ecef(Time tutc,const double *erpv,double *U,double *gmst);
double SunMoonPos(Time ut1t,const double *erp_val,double *sun_pos,double *moon_pos);
template <typename Iter1,typename Iter2>
void Xyz2Enu(const Iter1 BlhPos,Iter2 TranMat){
    double sinLat=sin(BlhPos[0]),cosLat=cos(BlhPos[0]),sinLon=sin(BlhPos[1]),cosLon=cos(BlhPos[1]);

    TranMat[0]=-sinLon;			TranMat[3]=cosLon;			TranMat[6]=0.0;
    TranMat[1]=-sinLat*cosLon;	TranMat[4]=-sinLat*sinLon;	TranMat[7]=cosLat;
    TranMat[2]=cosLat*cosLon;	TranMat[5]=cosLat*sinLon;	TranMat[8]=sinLat;
}

template <typename Iter1,typename Iter2,typename Iter3>
void Ecef2Enu(const Iter1 BlhPos,const Iter2 SightVec,Iter3 EnuPos){
    double E[9];

    Xyz2Enu(BlhPos,E);
    MatMul("NN",3,1,3,1.0,E,SightVec,0.0,EnuPos);
}


#endif //MSNAVS_CMNFUNC_H
