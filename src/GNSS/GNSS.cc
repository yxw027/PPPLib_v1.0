/**************************************************************************

Copyright(c): 2020, by Chao Chen, All rights reserved.
              China University of Mining and Technology, XuZhou, P.R. China

Author: Chao Chen(cchen@cumt.edu.cn)

Date: 2020-05-30

Description:

**************************************************************************/

#include "CmnFunc.h"
#include "PPPLibGlo.h"
#include "GNSS.h"

extern const double kGNSSFreqs[NSYS][MAXFREQ]={
        {FREQ_GPS_L1,  FREQ_GPS_L2, FREQ_GPS_L5,    FREQ_NONE,      FREQ_NONE},
        {FREQ_BD2_B1I, FREQ_BD2_B2I, FREQ_BD2_B3I,  FREQ_NONE,      FREQ_NONE},
        {FREQ_BD3_B1I, FREQ_BD3_B2a, FREQ_BD3_B3I, FREQ_BD3_B1C, FREQ_BD3_B2b},
        {FREQ_GAL_E1,  FREQ_GAL_E5a,FREQ_GAL_E5b,   FREQ_GAL_E5,    FREQ_GAL_E6},
        {FREQ_GLO_G1,  FREQ_GLO_G2,  FREQ_GLO_G3,  FREQ_GLO_G1a,  FREQ_GLO_G2a},
        {FREQ_QZS_L1,  FREQ_QZS_L2,  FREQ_QZS_L5,   FREQ_QZS_L6,      FREQ_NONE},
};

extern string kSignalCodes[]={
        "",  "1C","1S","1L","1X","1P","1W","1Y","1M","1A",
        "1B","1Z","1D","2C","2D","2S","2L","2X","2P","2W",
        "2Y","2M","2I","2Q","3I","3Q","3X","4A","4B","4X",
        "5I","5Q","5X","5D","5P","5Z","6A","6B","6X","6C",
        "6Z","6S","6L","6E","6I","6Q","7I","7Q","7X","7D",
        "7P","7Z","8I","8Q","8X","8D","8P", "8Z" , "" , "" ,
};

static string kCodePriors[NSYS][MAXFREQ]={
        /* GPS L1           L2      L5                        */
        {"CWPYMSLX", "CWPYMDSLX",  "IQX",     "",     ""},
        /* BD2 B1I        B2I      B3I        */
        {"IQX",           "IQX",  "IQXA",     "",     ""},
        /* BD3 B1I        B2a      B3I      B1C    B2b*/
        {"IQX",           "DPX",   "IQX",  "DPXA","DPZ" },
        /* GAL E1          E5a      E5b     E5      E6        */
        {"CABXZ",          "IQX",  "IQX",  "IQX","CABXZ"},
        /* GLO G1           G2       G3      G1a    G2a       */
        {   "PC",           "PC",  "IQX",  "ABX", "ABXP"},
        /* QZS L1           L2       L5       L6              */
        {"CSLXZ",          "SLX", "IQXDPZ", "SLXEZ",  ""},
};

extern string kCodeInterBiasPair[NSYS+1][MAXCBIASPAIR]= {
        //GPS
        {"C1C-C2W", "C1C-C5Q", "C1C-C5X", "C1W-C2W", "C1C-C1W","C2C-C2W","C2W-C2S","C2W-C2L","C2W-C2X",""},
        //BDS-2
        {"C2I-C7I","C2I-C6I","","","","","","","",""},
        //BDS-3
        {"C1X-C5X", "C1P-C5P", "C1D-C5D", "C1X-C6I", "C1P-C6I", "C1D-C6I", "C2I-C6I", "C1X-C7Z", "C1X-C8X", ""},
        //GAL
        {"C1C-C5Q", "C1C-C6C", "C1C-C7Q", "C1C-8Q", "C1X-C5X", "C1X-C7X", "C1X-C8X", "", "", ""},
        //GLO
        {"C1C-C2C","C1C-C2P","C1P-C2P","C1C-C1P","C2C-C2P","","","","",""},
        //QZS
        {"C1C-C2L","C1C-C5X","C1C-C5Q","C1X-C2X","C1X-C5X","C1C-C1X","","","",""},
};


extern const unsigned char kGPSFreqBand[]={
        0, 1, 1, 1, 1, 1, 1, 1, 1, 0,
        0, 0, 0, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 0, 0, 0, 0, 0, 0, 0, 0,
        3, 3, 3, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};

extern const unsigned  char kBD2FreqBand[]={
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
        0, 0, 1, 1, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 3, 0, 3, 0,
        0, 0, 0, 0, 3, 3, 2, 2, 2, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};

extern const unsigned  char kBD3FreqBand[]={
        0, 0, 0, 0, 4, 4, 0, 0, 0, 4,
        0, 0, 4, 0, 0, 0, 0, 1, 0, 0,
        0, 0, 1, 1, 0, 0, 0, 0, 0, 0,
        0, 0, 2, 2, 2, 0, 3, 0, 3, 0,
        0, 0, 0, 0, 3, 3, 0, 0, 0, 5,
        5, 5, 0, 0, 0, 0, 0, 0, 0, 0,
};

extern const unsigned  char kGALFreqBand[]={
        0, 1, 0, 0, 1, 0, 0, 0, 0, 1,
        1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        2, 2, 2, 0, 0, 0, 5, 5, 5, 5,
        5, 0, 0, 0, 0, 0, 3, 3, 3, 0,
        0, 0, 4, 4, 4, 0, 0, 0, 0, 0,
};
extern const unsigned  char kGLOFreqBand[]={
        0, 1, 0, 0, 0, 1, 0, 0, 0, 0,
        0, 0, 0, 2, 0, 0, 0, 0, 2, 0,
        0, 0, 0, 0, 3, 3, 3, 4, 4, 4,
        0, 0, 0, 0, 0, 0, 5, 5, 5, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};
extern const unsigned  char kQZSFreqBand[]={
        0, 1, 1, 1, 1, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 2, 2, 2, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        3, 3, 3, 3, 3, 0, 0, 0, 4, 0,
        4, 4, 4, 4, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};



// signal transition related class
Signal::Signal() {
    sat_sys_=GPS;
    n_=0;
    for(int &i:frq_) i=0;
    for(int &i:pos_) i=-1;
    for(unsigned char &i:pri_) i=0;
    for(unsigned char &i:type_) i=0;
    for(unsigned char &i:code_) i=CODE_NONE;
    for(double &i:shift_) i=0.0;
}

Signal::~Signal() {

}

unsigned char Signal::Signal2Code(string signal, int *frq, int sys) {
    int i;
    if(frq) *frq=0;
    for(i=1;i<MAXCODE;i++){
        if(kSignalCodes[i].compare(signal)!=0) continue;
        if(frq) {
            switch(sys) {
                case GPS: *frq=kGPSFreqBand[i];break;
                case BD2: *frq=kBD2FreqBand[i];break;
                case BD3: *frq=kBD3FreqBand[i];break;
                case GAL: *frq=kGALFreqBand[i];break;
                case GLO: *frq=kGLOFreqBand[i];break;
                case QZS: *frq=kQZSFreqBand[i];break;
                default:
                    *frq=0;break;
            }
        }
        return (unsigned char) i;
    }
    return CODE_NONE;
}

std::__cxx11::string Signal::Code2Signal(unsigned char code, int * frq, int sys) {
    if(frq) *frq=0;
    if(code<=CODE_NONE||MAXCODE<code) return "";
    if(frq) {
        switch(sys) {
            case GPS: *frq=kGPSFreqBand[code];break;
            case BD2: *frq=kBD2FreqBand[code];break;
            case BD3: *frq=kBD3FreqBand[code];break;
            case GAL: *frq=kGALFreqBand[code];break;
            case GLO: *frq=kGLOFreqBand[code];break;
            case QZS: *frq=kQZSFreqBand[code];break;
            default:
                *frq=0;break;
        }
    }
    return kSignalCodes[code];
}

int Signal::GetSignalPri(int sys, unsigned char code) {
    if(code==CODE_NONE) return 0;
    size_t str_pos;
    string str_signal,str_buff;
    int i,j;

    str_signal=Code2Signal(code,&j,sys);
    switch(sys){
        case GPS: i=INDEXGPS;break;
        case BD2: i=INDEXBD2;break;
        case BD3: i=INDEXBD3;break;
        case GAL: i=INDEXGAL;break;
        case GLO: i=INDEXGLO;break;
        case QZS: i=INDEXQZS;break;
        default:i=0;
    }

    return ((str_pos=kCodePriors[i][j-1].find(str_signal[1]))!=string::npos)?14-(int)str_pos:0;
}

// satellite mask class
SatMask::SatMask() {
    sat_no_=sat_prn_=1;
    sat_sys_=GPS;
    str_sat_sys_="G";
    sat_id_="G01";
    sat_sys_idx_=0;
    bd3_flag_=0;
}

SatMask::SatMask(int sat_no) {
    if(sat_no<=1) sat_no=1;
    else if(sat_no>=MAXSAT) sat_no=MAXSAT-1;
    sat_no_=sat_no;
    SatNo2SysPrn();
    SatNo2ID();
    SatID2No();
}

SatMask::SatMask(int sat_sys, int sat_prn) {
    sat_sys_=sat_sys;
    sat_prn_=sat_prn;
    SatSysPrn2No();
    SatNo2ID();
    SatID2No();
}

SatMask::SatMask(string sat_id) {
    sat_no_=0;
    sat_id_=sat_id;
    SatID2No();
    SatNo2SysPrn();
}

SatMask::~SatMask() {

}

void SatMask::SatSysPrn2No() {
    if (sat_prn_<=0) return;
    switch (sat_sys_) {
        case GPS:
            if (sat_prn_<MINPRNGPS||MAXPRNGPS<sat_prn_) {sat_no_=0;return;}
            sat_no_=sat_prn_-MINPRNGPS+1;break;
        case BD2:
            if(sat_prn_<MINPRNBD2||MAXPRNBD2<sat_prn_){sat_no_=0;return;}
            sat_no_=NSATGPS+sat_prn_-MINPRNBD2+1;break;
        case BD3:
            if (sat_prn_<MINPRNBD3||MAXPRNBD3+MAXPRNBD2<sat_prn_) {sat_no_=0;return;}
            sat_no_=NSATGPS+sat_prn_-MINPRNBD3+1;break;
        case GAL:
            if (sat_prn_<MINPRNGAL||MAXPRNGAL<sat_prn_) {sat_no_=0;return;}
            sat_no_=NSATGPS+NSATBD2+NSATBD3+sat_prn_-MINPRNGAL+1;break;
        case GLO:
            if (sat_prn_<MINPRNGLO||MAXPRNGLO<sat_prn_) {sat_no_=0;return;}
            sat_no_=NSATGPS+NSATBD2+NSATBD3+NSATGAL+sat_prn_-MINPRNGLO+1;break;
        case QZS:
            if (sat_prn_<MINPRNQZS||MAXPRNQZS<sat_prn_) {sat_no_=0;return;}
            sat_no_= NSATGPS+NSATBD2+NSATBD3+NSATGLO+NSATGAL+sat_prn_-MINPRNQZS+1;break;
        case IRN:
            sat_no_=0;break;
        case SBS:
            sat_no_=0;break;
    }
    return;
}

void SatMask::SatNo2SysPrn() {
    int sat_num=sat_no_;
    sat_sys_=NONE;
    if (sat_num<=0||MAXSAT<sat_num) sat_num=0;
    else if (sat_num<=NSATGPS) {
        sat_sys_=GPS; sat_num+=MINPRNGPS-1;
    }
    else if ((sat_num-=NSATGPS)<=NSATBD2) {
        sat_sys_=BD2; sat_num+=MINPRNBD2-1;
    }
    else if((sat_num-=NSATBD2)<=NSATBD3){
        sat_sys_=BD3;sat_num+=MINPRNBD3-1+NSATBD2;
    }
    else if ((sat_num-=NSATBD3)<=NSATGAL) {
        sat_sys_=GAL; sat_num+=MINPRNGAL-1;
    }
    else if ((sat_num-=NSATGAL)<=NSATGLO) {
        sat_sys_=GLO; sat_num+=MINPRNGLO-1;
    }
    else if ((sat_num-=NSATGLO)<=NSATQZS) {
        sat_sys_=QZS; sat_num+=MINPRNQZS-1;
    }
    else sat_num=0;
    sat_prn_=sat_num;
}

void SatMask::SatNo2ID() {
    string buf;
    switch (sat_sys_){
        case GPS: sat_id_="G"+Int2Str(2,"0",sat_prn_-MINPRNGPS+1,buf); return;
        case BD2: sat_id_="C"+Int2Str(2,"0",sat_prn_-MINPRNBD2+1,buf); return;
        case BD3: sat_id_="C"+Int2Str(2,"0",sat_prn_-MINPRNBD3+1,buf); return;
        case GAL: sat_id_="E"+Int2Str(2,"0",sat_prn_-MINPRNGAL+1,buf); return;
        case GLO: sat_id_="R"+Int2Str(2,"0",sat_prn_-MINPRNGLO+1,buf); return;
        case QZS: sat_id_="J"+Int2Str(2,"0",sat_prn_,buf); return;
    }
    sat_id_="";
}

void SatMask::SatID2No() {
    char code;
    string aaa=" 0123456789";
    sat_sys_=NONE;

    if (aaa.find(sat_id_[0])!=string::npos) {
        if (Str2Int((sat_id_.substr(1,2)),sat_prn_)==0) return;
        if      (MINPRNGPS<=sat_prn_&&sat_prn_<=MAXPRNGPS) sat_sys_=GPS;
        else if (MINPRNQZS<=sat_prn_&&sat_prn_<=MAXPRNQZS) sat_sys_=QZS;
        else return;
        SatSysPrn2No();
    }

    if(Str2Int((sat_id_.substr(1,2)),sat_prn_)==0) return;
    code=sat_id_[0];

    switch (code) {
        case 'G': sat_sys_=GPS;str_sat_sys_='G';sat_sys_idx_=INDEXGPS;sat_prn_+=MINPRNGPS-1; break;
        case 'C':
            if(sat_prn_<=18){
                sat_sys_=BD2;str_sat_sys_='C';sat_sys_idx_=INDEXBD2;sat_prn_+=MINPRNBD2-1; break;
            }else{
                sat_sys_=BD3;str_sat_sys_='C';sat_sys_idx_=INDEXBD3;sat_prn_+=MINPRNBD3-1; break;
            }
        case 'E': sat_sys_=GAL;str_sat_sys_='E';sat_sys_idx_=INDEXGAL;sat_prn_+=MINPRNGAL-1; break;
        case 'R': sat_sys_=GLO;str_sat_sys_='R';sat_sys_idx_=INDEXGLO;sat_prn_+=MINPRNGLO-1; break;
        case 'J': sat_sys_=QZS;str_sat_sys_='J';sat_sys_idx_=INDEXQZS;sat_prn_+=MINPRNQZS-1; break;
        case 'I': sat_sys_=IRN;str_sat_sys_='I';sat_prn_=0;sat_no_=0;break;
        case 'S': sat_sys_=SBS;str_sat_sys_='S';sat_prn_=0;sat_no_=0;break;
        default: return;
    }
    SatSysPrn2No();
}

// obs class
ObsData::ObsData() {
    rcv_=0;
    stat_=USED;
    sat_=SatMask();
    sig_recep_=sig_trans_=Time();
    for(int i=0;i<NFREQ+NEXOBS;i++){
        code_[i]=0;
        P_[i]=L_[i]=D_[i]=SNR_[i]=LLI_[i]=meas_var_[i]=0.0;
        post_res_[0][i]=pri_res_[0][i]=0.0;
        post_res_[1][i]=pri_res_[1][i]=0.0;
        dcb_[i]=code_osb_[i]=phase_osb_[i]=0.0;
    }
    use_f1_=use_f2_=use_f3_=use_f4_=-1;
    eph_idx_=-1;
    svh_=0;
    frq3_=0;
    sig_vec_=ant_corr_=Vec3(0.0);
    pos_=vel_=Vec3(0.0);
    for(int i=0;i<2;i++){
        clk_err_[i]=azel_[i]=0.0;
    }
    clk_rel_=0.0;
    for(auto &i:lam_) i=-1.0;
    for(auto &i:frq_) i=-1.0;
    sat_var_=0.0;
}

ObsData::~ObsData() {

}

void ObsData::ReSet() {
    int i;
    sat_=SatMask();

    for(i=0;i<NFREQ+NEXOBS;i++){
        SNR_[i]=LLI_[i]=code_[i]='\0';
        L_[i]=P_[i]=D_[i]=0.0;
    }
}

void ObsData::ReSetSat() {
    svh_=0;
    sat_var_=0.0;
    pos_=Vec3(0.0);
    vel_=Vec3(0.0);
//    azel_[0]=azel_[1]=0.0;
    for(auto &i:clk_err_) i=0.0;
}

int ObsData::SigTransTime() {
    double pr;
    int i;

    sig_trans_=sig_recep_;
    for(i=0,pr=0.0;i<NFREQ;i++) if((pr=P_[i])!=0.0) break;
    if(i<NFREQ){
        sig_trans_.TimeAdd(-pr/CLIGHT);
        return 1;
    }
    else{
        char buff[MAXBUFF]={'\0'};
        sprintf(buff,"%4d %20s: %4s no pseudorange observation for signal transmit time compute",epoch_,sig_recep_.time_str_.c_str(),sat_.sat_id_.c_str());
        LOG(WARNING)<<buff;
        stat_=EX_NOPOBS;
        return 0;
    }
}

void ObsData::ClkCorrection() {
    sig_trans_.TimeAdd(-clk_err_[0]);
}

double ObsData::ShapiroCorr(Vec3 rec_pos) {
    if(rec_pos.i_*rec_pos.j_*rec_pos.k_==0.0) rec_pos=Vec3(100.0);
    double rr=Vec3Norm(rec_pos);
    double rs=Vec3Norm(pos_);
    double l=Vec3Norm(pos_-rec_pos);
    double mu=0.0;
    int sys=sat_.sat_sys_;

    switch(sys){
        case GPS: mu=MU_GPS;break;
        case BD2:
        case BD3: mu=MU_BDS;break;
        case GAL: mu=MU_GAL;break;
        case GLO: mu=MU_GLO;break;
        case QZS: mu=MU_GPS;break;
    }
    return -2.0*mu/SQR(CLIGHT)*log((rr+rs+l)/(rr+rs-l));
}

double ObsData::SagnacCorr(Vec3 rec_pos) {
    return OMGE_GPS*(pos_.i_*rec_pos.j_-pos_.j_*rec_pos.i_)/CLIGHT;
}

Sta::Sta() {
    ant_setup_=itrf_=del_fmt_=0;
    hgt_=0.0;
    for(double &i:pos_) i=0.0;
    for(double &i:del_) i=0.0;
}

Sta::~Sta(){

}

Obss::Obss() {
    num_=0;
    rcv_=0;
    sta_=Sta();
}

Obss::~Obss() {
    sat_infos_.clear();
}

void Obss::ReSet() {
    num_=0;
    sat_infos_.clear();
}

PreClk::PreClk() {
    time_=Time();
    for(double &i:clk_) i=0.0;
    for(float  &j:std_) j=0.0;
}

PreClk::~PreClk() {

}

PreEph::PreEph() {
    time_=Time();
    for(int i=0;i<MAXSAT;i++){
        for(int j=0;j<4;j++){
            pos_[i][j]=std_[i][j]=vel_[i][j]=vst_[i][j]=0.0;
        }
        for(int j=0;j<3;j++){
            cov_[i][j]=vco_[i][j]=0.0;
        }
    }
}

PreEph::~PreEph() {

}

GloEph::GloEph() {
    sat_=1;
    iode_=frq_=svh_=sva_=age_=0;
    toe_=tof_=Time();
    for(int i=0;i<3;i++) pos_[i]=vel_[i]=acc_[i]=0.0;
    taun_=gamn_=dtaun_=0.0;
}

GloEph::~GloEph() {

}

Eph::Eph() {
    sat_=1;
    iode_=iodc_=sva_=svh_=week_=code_=flag_=0;
    toc_=toc_=ttr_=Time();
    A_=e_=i0_=OMG0_=omg_=M0_=deln_=OMGd_=idot_=0.0;
    crc_=crs_=cuc_=cus_=cic_=cis_=0.0;
    toes_=fit_=f0_=f1_=f2_=0.0;
    for(int i=0;i<4;i++) tgd_[i]=0.0;
}

Eph::~Eph() {

}

double Eph::RelClkCorr(double mu,double sinE) const {
    return 2.0*sqrt(mu*A_)*e_*sinE/SQR(CLIGHT);
}

Nav::Ptr_ Nav::nav_(new Nav());

Nav::Nav() {
    n_=0;ng_=0;
    np_=0; nc_=0;
    for(int i=0;i<MAXPRNGLO+1;i++) glo_frq_num_[i]=0;
    for(double &i:glo_cp_bias_) i=0.0;
    leaps_=0;

    for(int i=0;i<NSYS;i++){
        for(int j=0;j<8;j++){
            ion_para_[i].ion[j]=0.0;
        }
    }
    for(int j=0;j<MAXSAT;j++) {
        for(int k=0;k<MAXCBIASPAIR;k++) dcb_[j][k]=0.0;
    }
}

Nav::~Nav() {
    eph_.clear();
    glo_eph_.clear();
    pre_eph_.clear();pre_clk_.clear();
}

double Nav::SatWaveLen(int sat_no, int frq) {
    SatMask sat(sat_no);
    int sys=sat.sat_sys_,idx_sys=sat.sat_sys_idx_;
    const double dfreq[]={FREQ_GLO_D1,FREQ_GLO_D2};

    if(kGNSSFreqs[idx_sys][frq]==0.0){
        return 0.0;
    }
    if(sys==GLO){
        if(frq==2){
            return CLIGHT/kGNSSFreqs[idx_sys][frq];
        }
        else{
            for(int i=0;i<ng_;i++){
                if(glo_eph_[i].sat_!=sat.sat_no_) continue;
                glo_frq_num_[sat.sat_prn_-1]=glo_eph_[i].frq_;
                return CLIGHT/(kGNSSFreqs[idx_sys][frq]+dfreq[frq]*glo_eph_[i].frq_);
            }
        }
    }

    return CLIGHT/kGNSSFreqs[idx_sys][frq];
}

void Nav::UpdateNav() {
    int i,j;
    SatMask sat;

    for(i=0;i<MAXSAT;i++){
        sat=SatMask(i+1);
        for(j=0;j<NFREQ;j++){
            lam_[i][j]=SatWaveLen(i+1,j);
        }
        if(sat.sat_sys_==GPS){
            LOG_N_TIMES(1,DEBUG)<<"GPS wavelength: "<<"L1   "<<setprecision(6)<<lam_[i][0]<<" L2  "<<lam_[i][1]\
                               <<" L5  "<<lam_[i][2];
        }
        else if(sat.sat_sys_==BD2){
            LOG_N_TIMES(1,DEBUG)<<"BD2 wavelength: "<<"B1I  "<<setprecision(6)<<lam_[i][0]<<" B2I "<<lam_[i][1]\
                               <<" B3I "<<lam_[i][2];
        }
        else if(sat.sat_sys_==BD3){
            LOG_N_TIMES(1,DEBUG)<<"BD3 wavelength: "<<"B1I  "<<setprecision(6)<<lam_[i][0]<<" B2a "<<lam_[i][1]\
                               <<" B3I "<<lam_[i][2];
        }
        else if(sat.sat_sys_==GAL){
            LOG_N_TIMES(1,DEBUG)<<"GAL wavelength: "<<"E1   "<<setprecision(6)<<lam_[i][0]<<" E5a "<<lam_[i][1]\
                               <<" E5b "<<lam_[i][2];
        }
        else if(sat.sat_sys_==GLO){
            LOG_N_TIMES(1,DEBUG)<<"GLO wavelength: "<<"G1   "<<setprecision(6)<<lam_[i][0]<<" G2  "<<lam_[i][1]\
                               <<" G3  "<<lam_[i][2];
        }
        else if(sat.sat_sys_==QZS){
            LOG_N_TIMES(1,DEBUG)<<"QZS wavelength: "<<"L1   "<<setprecision(6)<<lam_[i][0]<<" L2  "<<lam_[i][1]\
                               <<" L5  "<<lam_[i][2];
        }

    }
}

Nav::Ptr_ Nav::GetInstance() {
    return nav_;
}

double GeoDist(const Vec3 sat_pos,const Vec3 rec_pos,Vec3& sig_vec) {
    double r;
    Vec3 vec;

    if(Vec3Norm(sat_pos)<RE_WGS84) return -1.0;
    vec=sat_pos-rec_pos;
    r=Vec3Norm(vec);
    sig_vec=vec/r;

    return r;
}

double SatAzEl(const Vec3 blh,const Vec3 sig_vec,double *az_el) {
    double az=0.0,el=PI/2.0;
    Vec3 enu;

    if(blh.k_>-RE_WGS84){
        enu=Xyz2Enu(sig_vec,blh);
        az=VecDot2(enu,enu)<1E-12?0.0:atan2(enu.i_,enu.j_);
        if(az<0.0) az+=2*PI;
        el=asin(enu.k_);
    }
    if(az_el) {az_el[0]=az;az_el[1]=el;}
    return el;
}