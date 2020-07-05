//
// Created by cc on 4/9/20.
//

#ifndef MSNAVS_GNSS_H
#define MSNAVS_GNSS_H

#include "CmnFunc.h"
#include "TidModel.h"
// Core option related definition

#define ENABD2
#define ENABD3
#define ENAGAL
#define ENAGLO
#define ENAQZS
// GNSS related definition
enum {MINPRNGPS=1,MAXPRNGPS=32,NSATGPS=(MAXPRNGPS-MINPRNGPS+1),NSYSGPS=1};
#ifdef ENABD2
enum {MINPRNBD2=1,MAXPRNBD2=18,NSATBD2=(MAXPRNBD2-MINPRNBD2+1),NSYSBD2=1};
#else
enum {MINPRNBD2=0,MAXPRNBD2=0,NSATBD2=(MAXPRNBD2-MINPRNBD2+1),NSYSBD2=0};
#endif
#ifdef ENABD3
enum {MINPRNBD3=1,MAXPRNBD3=27,NSATBD3=(MAXPRNBD3-MINPRNBD3+1),NSYSBD3=1};
#else
enum {MINPRNBD3=0,MAXPRNBD3=0,NSATBD3=(MAXPRNBD3-MINPRNBD3+1),NSYSBD3=1};
#endif
#ifdef ENAGAL
enum {MINPRNGAL=1,MAXPRNGAL=36,NSATGAL=(MAXPRNGAL-MINPRNGAL+1),NSYSGAL=1};
#else
enum {MINPRNGAL=0,MAXPRNGAL=0,NSATGAL=(MAXPRNGAL-MINPRNGAL+1),NSYSGAL=0};
#endif
#ifdef ENAGLO
enum {MINPRNGLO=1,MAXPRNGLO=28,NSATGLO=(MAXPRNGLO-MINPRNGLO+1),NSYSGLO=1};
#else
enum {MINPRNGLO=0,MAXPRNGLO=0,NSATGLO=(MAXPRNGLO-MINPRNGLO+1),NSYSGLO=0};
#endif
#ifdef ENAQZS
enum {MINPRNQZS=193,MAXPRNQZS=201,NSATQZS=(MAXPRNQZS-MINPRNQZS+1),NSYSQZS=1};
#else
enum {MINPRNQZS=0,MAXPRNQZS=0,NSATQZS=(MAXPRNQZS-MINPRNQZS+1),NSYSQZS=0};
#endif
enum {NSYS=NSYSGPS+NSYSBD2+NSYSBD3+NSYSGAL+NSYSGLO+NSYSQZS,MAXSAT=NSATGPS+NSATBD2+NSATBD3+NSATGAL+NSATGLO+NSATQZS};
enum {INDEXGPS,INDEXBD2,INDEXBD3,INDEXGAL,INDEXGLO,INDEXQZS};
enum {TIMESYSGPS,TIMESYSBDS,TIMESYSGAL,TIMESYSUTC,TIMESYSQZS};

#define FREQ_NONE    0.0
#define FREQ_GPS_L1  1.575420E9       /*  L1 */
#define FREQ_GPS_L2  1.227600E9       /*  L2 */
#define FREQ_GPS_L5  1.176450E9       /*  L5 */
#define FREQ_BD2_B1I 1.561098E9
#define FREQ_BD2_B2I 1.207140E9
#define FREQ_BD2_B3I 1.268520E9
#define FREQ_BD3_B1I 1.561098E9
#define FREQ_BD3_B2a 1.176450E9
#define FREQ_BD3_B3I 1.268520E9
#define FREQ_BD3_B1C 1.575420E9
#define FREQ_BD3_B2b 1.207140E9
#define FREQ_GAL_E1  1.575420E9       /*  L1 */
#define FREQ_GAL_E5a 1.176450E9       /*  L5 */
#define FREQ_GAL_E5b 1.207140E9       /*  L7 */
#define FREQ_GAL_E5  1.191795E9       /*  L8 */
#define FREQ_GAL_E6  1.278750E9       /*  L6 */
#define FREQ_GLO_G1  1.602000E9       /*  L1 */
#define FREQ_GLO_D1  0.562500E6       /*  L1 */
#define FREQ_GLO_G1a 1.600995E9       /*  L4 */
#define FREQ_GLO_G2  1.246000E9       /*  L2 */
#define FREQ_GLO_D2  0.437500E6       /*  L1 */
#define FREQ_GLO_G2a 1.248060E9       /*  L6 */
#define FREQ_GLO_G3  1.202025E9       /*  L3 */
#define FREQ_QZS_L1  1.575420E9       /*  L1 */
#define FREQ_QZS_L2  1.227600E9       /*  L2 */
#define FREQ_QZS_L5  1.176450E9       /*  L5 */
#define FREQ_QZS_L6  1.278750E9       /*  L6 */
#define MINFREQ_GLO -7                /* min frequency number glonass */
#define MAXFREQ_GLO 13                /* max frequency number glonass */

#define CODE_NONE    0
#define CODE_L1C     1               /* GPS L1 GLO G1 GAL-E1 QZS-L1 */
#define CODE_L1S     2               /* GPS L1 */
#define CODE_L1L     3
#define CODE_L1X     4
#define CODE_L1P     5
#define CODE_L1W     6
#define CODE_L1Y     7
#define CODE_L1M     8
#define CODE_L1A     9
#define CODE_L1B     10
#define CODE_L1Z     11
#define CODE_L1D     12
#define CODE_L2C     13
#define CODE_L2D     14
#define CODE_L2S     15
#define CODE_L2L     16
#define CODE_L2X     17
#define CODE_L2P     18
#define CODE_L2W     19
#define CODE_L2Y     20
#define CODE_L2M     21
#define CODE_L2I     22
#define CODE_L2Q     23
#define CODE_L3I     24
#define CODE_L3Q     25
#define CODE_L3X     26
#define CODE_L4A     27
#define CODE_L4B     28
#define CODE_L4X     29
#define CODE_L5I     30
#define CODE_L5Q     31
#define CODE_L5X     32
#define CODE_L5D     33
#define CODE_L5P     34
#define CODE_L5Z     35
#define CODE_L6A     36
#define CODE_L6B     37
#define CODE_L6X     38
#define CODE_L6C     39
#define CODE_L6Z     40
#define CODE_L6S     41
#define CODE_L6L     42
#define CODE_L6E     43
#define CODE_L6I     44
#define CODE_L6Q     45
#define CODE_L7I     46
#define CODE_L7Q     47
#define CODE_L7X     48
#define CODE_L7D     49
#define CODE_L7P     50
#define CODE_L7Z     51
#define CODE_L8I     52
#define CODE_L8Q     53
#define CODE_L8X     54
#define CODE_L8D     55
#define CODE_L8P     56
#define CODE_L8Z     57
#define MAXCODE      57
#define MAXOBSTYPE   64

#define GPS_C1CC2W   0
#define GPS_C1CC5Q   1
#define GPS_C1CC5X   2
#define GPS_C1WC2W   3
#define GPS_C1CC1W   4
#define GPS_C2CC2W   5
#define GPS_C2WC2S   6
#define GPS_C2WC2L   7
#define GPS_C2WC2X   8
#define GPS_MAXCPAIR 8

#define BD2_C2IC7I   0
#define BD2_C2IC6I   1
#define BD2_MAXCPAIR 1

#define BD3_C1XC5X   0
#define BD3_C1PC5P   1
#define BD3_C1DC5D   2
#define BD3_C1XC6I   3
#define BD3_C1PC6I   4
#define BD3_C1DC6I   5
#define BD3_C2IC6I   6
#define BD3_C1XC7Z   7
#define BD3_C1XC8X   8
#define BD3_MAXCPAIR 8

#define GAL_C1CC5Q   0
#define GAL_C1CC6C   1
#define GAL_C1CC7Q   2
#define GAL_C1CC8Q   3
#define GAL_C1XC5X   4
#define GAL_C1XC7X   5
#define GAL_C1XC8X   6
#define GAL_MAXCPAIR 6

#define GLO_C1CC2C   0
#define GLO_C1CC2P   1
#define GLO_C1PC2P   2
#define GLO_C1CC1P   3
#define GLO_C2CC2P   4
#define GLO_MAXCPAIR 4

#define QZS_C1CC2L   0
#define QZS_C1CC5X   1
#define QZS_C1CC5Q   2
#define QZS_C1XC2X   3
#define QZS_C1XC5X   4
#define QZS_C1CC1X   5
#define QZS_MAXCPAIR 5


#define MAXDTOE_GPS  7200.0              /* max time difference to GPS Toe (s) */
#define MAXDTOE_BDS  21600.0             /* max time difference to BeiDou Toe (s) */
#define MAXDTOE_GAL  10800.0             /* max time difference to Galileo Toe 14400 (s) */
#define MAXDTOE_GLO  1800.0              /* max time difference to GLONASS Toe (s) */
#define MAXDTOE_QZS  7200.0              /* max time difference to QZSS Toe (s) */
#define MU_GPS   3.9860050E14            /* gravitational constant         ref [1] */
#define MU_BDS   3.986004418E14          /* earth gravitational constant   ref [9] */
#define MU_GAL   3.986004418E14          /* earth gravitational constant   ref [7] */
#define MU_GLO   3.9860044E14            /* gravitational constant         ref [2] */
#define OMGE_GPS 7.2921151467E-5     /* earth angular velocity (IS-GPS) (rad/s) */
#define OMGE_BDS 7.292115E-5             /* earth angular velocity (rad/s) ref [9] */
#define OMGE_GAL 7.2921151467E-5         /* earth angular velocity (rad/s) ref [7] */
#define OMGE_GLO 7.292115E-5             /* earth angular velocity (rad/s) ref [2] */

const int kGPS_II_R[]={2,11,13,14,16,19,20,21,22,28};   //(10) L1 L2
#define NUM_GPS_IIR  10
const int kGPS_IIR_M[]={5,7,12,15,17,29,31};            //(7)
#define NUM_GPS_IIRM  7
const int kGPS_II_F[]={1,3,6,8,9,10,24,25,26,27,30,32}; //(12)
#define NUM_GPS_IIF  12
const int kGPS_III_A[]={4,18};                          //(2) L1 L2 L5
#define NUM_GPS_IIIA  2
const int kBDS2_GEO[]={1,2,3,4,5};
#define NUM_BDS2_GEO  5
const int kBDS2_IGSO[]={6,7,8,9,10,13,16};
#define NUM_BDS2_IGSO 7
const int kBDS2_MEO[]={11,12,14};   //B1-2 B2b B3
#define NUM_BDS2_MEO  3
const int kBDS3_GEO[]={59,60};
#define NUM_BDS3_GEO  2
const int kBDS3_IGSO[]={38,39,40};
#define NUM_BDS3_IGSO 3
const int kBDS3_MEO[]={19,20,21,22,23,24,25,26,
                       27,28,29,30,32,33,34,35,
                       36,37,41,42,43,44,45,46}; //B1-2 B3 B1 B2a B2b B2
#define NUM_BDS3_MEO  24

const int kGAL_IOV[]={11,12,19,20};
#define NUM_GAL_IOV   4
const int kGAL_FOV[]={ 1, 2, 3, 4, 5, 7, 8, 9,
                      13,14,15,18,21,22,24,
                      25,26,27,30,31,33,36}; //E1 E5a E5b E5
#define NUM_GAL_FOV  22
const int kGLO_M[]={ 1, 2, 3, 4, 5, 6, 7, 8,
                     10,11,12,13,14,15,16,
                     17,18,19,22,23,24};
const int kGLO_M_PLUS[]={21};
const int kGLO_K1[]={9,20};  // L1 L2
const int kQZS_GEO[]={7};
#define NUM_QZS_GEO  1
const int kQZS_IGSO[]={1,2,3};   // L1 L2 L5
#define NUM_QZS_IGSO 3

enum GPSFrq{GPS_L1,GPS_L2,GPS_L5};
enum BD2Frq{BD2_B1I,BD2_B2I,BD2_B3I};
enum BD3Frq{BD3_B1I,BD3_B2a,BD3_B3I,BD3_B1C};
enum GALFrq{GAL_E1,GAL_E5a,GAL_E5b};
enum GLOFrq{GLO_G1,GLO_G2,GLO_G3};
enum QZSFrq{QZS_L1,QZS_L2,QZS_L5};

enum ObsType{OBS_PR,OBS_CP,OBS_DP};

// global consant
#define MAXOBS   64
#define MAXFREQ  5
#define MAXCBIASPAIR 10
extern const double kGNSSFreqs[NSYS][MAXFREQ];
extern string kSignalCodes[];
extern string kCodeInterBiasPair[NSYS+1][MAXCBIASPAIR];

class Nav;
// signal transition related class
class Signal {
public:
    Signal();
    ~Signal();

public:
    unsigned char Signal2Code(string signal,int *frq,int sys);
    string Code2Signal(unsigned char code, int *frq,int sys);
    int GetSignalPri(int sys, unsigned char code);

public:
    int sat_sys_;
    int n_;
    int frq_[MAXOBSTYPE];
    int pos_[MAXOBSTYPE];
    unsigned char pri_[MAXOBSTYPE];
    unsigned char type_[MAXOBSTYPE];
    unsigned char code_[MAXOBSTYPE];
    double shift_[MAXOBSTYPE];
};

// sat num/prn/sys class
class SatMask {
public:
    SatMask();
    SatMask(int sat_no);
    SatMask(int sat_sys,int sat_prn);
    SatMask(string sat_id);
    ~SatMask();

protected:
    void SatSysPrn2No();
    void SatNo2SysPrn();
    void SatID2No();
    void SatNo2ID();

public:
    int sat_no_,sat_sys_,sat_prn_;
    string sat_id_;
    string str_sat_sys_;
    int sat_sys_idx_;
    int bd3_flag_;
};

// gnss obs class
#define NFREQ  4
#define NEXOBS 0

typedef struct {
    Time last_time;
    int last_ep;
    int stat;
    double corr_P[NFREQ]; // corrected by code bias/satellite pcv/receiver ant model
    double corr_L[NFREQ]; // corrected by satellite pcv/receiver ant model/phase wind-up
    double IF_P[4];       // L1_L2 L1_L3 L2_L3 L1_L2_L3
    double IF_L[4];

    double mw[2];
    double smw[2];
    int mw_idx[2];
    double mw_var[2];
    double gf[2];

    int outc[NFREQ];
    int lock[NFREQ];
    int rejc[NFREQ];
    int slip[NFREQ];

    double dist;
    double clk_err[2];

    double strp_h;
    double strp_w;
    double strp;
    double map_trp_h;
    double map_trp_w;
    double map_grad_e;
    double map_grad_n;
    double var_trp;

    double sion;
    double map_ion;
    double var_ion;

    double tide;
    double sagnac;
    double shapiro;
    double phw;
    double dants[NFREQ];
    double dantr[NFREQ];

    double amb[NFREQ];
    double cbias[NFREQ];
    double pbias[NFREQ];

    double post_res[2][NFREQ];
    double pri_res[2][NFREQ];

}Isat_t;

class ObsData {
public:
    ObsData();
    ~ObsData();

public:
    void ReSet();
    void ReSetSat();
    int SigTransTime();
    void ClkCorrection();
    double ShapiroCorr(Vec3 rec_pos);
    double SagnacCorr(Vec3 rec_pos);

public:
    int epoch_;
    unsigned char rcv_;
    int stat_;
    SatMask sat_;
    Time sig_recep_,sig_trans_;
    unsigned char code_[NFREQ+NEXOBS];
    double P_[NFREQ+NEXOBS];
    double L_[NFREQ+NEXOBS];
    float  D_[NFREQ+NEXOBS];
    unsigned char SNR_[NFREQ+NEXOBS];
    unsigned char LLI_[NFREQ+NEXOBS];
    double frq_[NFREQ];
    double lam_[NFREQ];
    int use_f1_,use_f2_,use_f3_,use_f4_;
    int eph_idx_;
    int svh_;
    int frq3_;
    Vec3 pos_;
    Vec3 vel_;
    double clk_err_[2];
    double clk_rel_;
    double dcb_[5];
    double code_osb_[NFREQ];
    double datum_brdc;        // GPS L1L2(C1WC2W)  BD3 L3(B3I)  GAL E1E5a(C1XC7X) GLO L1L2(C1PC2P) QZS
    double phase_osb_[NFREQ];
    Vec3 ant_corr_;

    double azel_[2];
    Vec3 sig_vec_;
    double pri_res_[2][NFREQ+NEXOBS];
    double post_res_[2][NFREQ+NEXOBS];

    double sat_var_;
    double meas_var_[NFREQ+NEXOBS];
};

class Sta {
public:
    Sta();
    ~Sta();

public:
    string name_;
    string marker_;
    string ant_desc_;
    string ant_seri_;
    string rec_type_;
    string firm_ver_;
    int ant_setup_;
    int itrf_;
    int del_fmt_;        // 0:enu,1:xyz
    double pos_[3];
    double del_[3];
    double hgt_;
};

class Obss {
public:
    Obss();
    ~Obss();

public:
    void ReSet();

public:
    int num_;
    unsigned char rcv_;
    Sta sta_;
    vector<ObsData> sat_infos_;
};

//
typedef struct {
    double utc[4];
}UtcPara_t;

typedef struct {
    double ion[8];
}IonPara_t;

class PreClk {
public:
    PreClk();
    ~PreClk();

public:
    Time time_;             /* time (GPST) */
    double clk_[MAXSAT]; /* satellite clock (s) */
    float  std_[MAXSAT];  /* satellite clock std (s) */
};

class PreEph {
public:
    PreEph();
    ~PreEph();

public:
    Time time_;             /* time (GPST) */
    double pos_[MAXSAT][4]; /* satellite position/clock (ecef) (m|s) */
    float  std_[MAXSAT][4]; /* satellite position/clock std (m|s) */
    double vel_[MAXSAT][4]; /* satellite velocity/clk-rate (m/s|s/s) */
    float  vst_[MAXSAT][4]; /* satellite velocity/clk-rate std (m/s|s/s) */
    float  cov_[MAXSAT][3]; /* satellite position covariance (m^2) */
    float  vco_[MAXSAT][3]; /* satellite velocity covariance (m^2) */
};

class GloEph {
public:
    GloEph();
    ~GloEph();
public:
    int sat_;            /* satellite number */
    int iode_;           /* IODE (0-6 bit of tb field) */
    int frq_;            /* satellite frequency number */
    int svh_,sva_,age_;    /* satellite health, accuracy, age of operation */
    Time toe_;         /* epoch of epherides (gpst) */
    Time tof_;         /* message frame time (gpst) */
    double pos_[3];      /* satellite position (ecef) (m) */
    double vel_[3];      /* satellite velocity (ecef) (m/s) */
    double acc_[3];      /* satellite acceleration (ecef) (m/s^2) */
    double taun_,gamn_;   /* SV clock bias (s)/relative freq bias */
    double dtaun_;       /* delay between L1 and L2 (s) */
};

class Eph {
public:
    Eph();
    ~Eph();

public:
    double RelClkCorr(double mu,double sinE) const;

public:
    int sat_;            /* satellite number */
    int iode_,iodc_;      /* IODE,IODC */
    int sva_;            /* SV accuracy (URA index) */
    int svh_;            /* SV health (0:ok) */
    int week_;           /* GPS/QZS: gps week, GAL: galileo week */
    int code_;           /* GPS/QZS: code on L2 */
                         /* GAL: data source defined as rinex 3.03 */
                         /* BDS: data source (0:unknown,1:B1I,2:B1Q,3:B2I,4:B2Q,5:B3I,6:B3Q) */
    int flag_;           /* GPS/QZS: L2 P data flag */
                         /* BDS: nav type (0:unknown,1:IGSO/MEO,2:GEO) */
    Time toe_,toc_,ttr_; /* Toe,Toc,T_trans */
                         /* SV orbit parameters */
    double A_,e_,i0_,OMG0_,omg_,M0_,deln_,OMGd_,idot_;
    double crc_,crs_,cuc_,cus_,cic_,cis_;
    double toes_;          /* Toe (s) in week */
    double fit_;           /* fit interval (h) */
    double f0_,f1_,f2_;    /* SV clock parameters (af0,af1,af2) */
    double tgd_[4];        /* group delay parameters */
                           /* GPS/QZS:tgd[0]=TGD */
                           /* GAL    :tgd[0]=BGD E5a/E1,tgd[1]=BGD E5b/E1 */
                           /* CMP    :tgd[0]=BGD1,tgd[1]=BGD2 */
    double Adot_,ndot_;    /* Adot,ndot for CNAV */
};

class Nav {
public:
    using Ptr_=std::shared_ptr<Nav>;

public:
    Nav();
    ~Nav();

public:
    double SatWaveLen(int sat_no,int frq);
    void UpdateNav();
    static Ptr_  GetInstance();

private:
    static Ptr_ nav_;

public:
    int n_;         /* number of broadcast ephemeris */
    int ng_;       /* number of glonass ephemeris */
    int np_;       /* number of precise orbit */
    int nc_;       /* number of precise clock */
    vector<Eph> eph_;
    vector<GloEph> glo_eph_;
    vector<PreEph> pre_eph_;
    vector<PreClk> pre_clk_;
    int glo_frq_num_[MAXPRNGLO+1];
    double glo_cp_bias_[4];
    int leaps_;
    UtcPara_t utc_para_[NSYS];
    IonPara_t ion_para_[NSYS];
    double lam_[MAXSAT][NFREQ];
    double dcb_[MAXSAT][MAXCBIASPAIR];
};

double GeoDist(const Vec3 SatPos,const Vec3 RecPos,Vec3& SigVec);
double SatAzEl(const Vec3 blh,const Vec3 SigVec,double *azel);

#endif //MSNAVS_GNSS_H
