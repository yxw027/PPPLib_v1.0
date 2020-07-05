/**************************************************************************

Copyright(c): 2020, by Chao Chen, All rights reserved.
              China University of Mining and Technology, XuZhou, P.R. China

Author: Chao Chen(cchen@cumt.edu.cn)

Date: 2020-05-30

Description:

**************************************************************************/

#include "PPPLibGlo.h"
#include "BiaModel.h"
#include "ReadFiles.h"

MatchFile::MatchFile() {
    Time t;

    dir_=kPrcOpt.dir_;
    sep_=(char)FILEPATHSEP;
    t=kPrcOpt.prc_date_;
    year_=t.epoch_[0];
    yy_=year_>=2000.0?Round(year_-2000.0):Round(year_-1900.0);;
    doy_=(int)t.Time2Doy();
    wod_=(int)(t.Time2GPST(&week_)/86400);
}

MatchFile::~MatchFile() {

}

int MatchFile::MatchProd() {
    string prods_dir;
    char f[MAXSTRPATH]={'\0'};

    if(!kFileOpt.customize_file_){
        sprintf(f,"%s%c%d%c%d%c%03d%c%s",dir_.c_str(),sep_,year_,sep_,week_,sep_,doy_,sep_,"prods");
        prods_dir=f;f[0]='\0';
    }else prods_dir=kFileOpt.customize_dir_;

    LOG(DEBUG)<<"Prods dir: "<<prods_dir;

    sprintf(f,"%s%cbrdm%03d0.%02dp",prods_dir.c_str(),sep_,doy_,yy_);
    if((access(f,0))==-1){
        cout<<"[ERROR]-"<<__FUNCTION__<<": "<<"no found "<<f<<endl;
        defaultLogger->error("Broadcast ephemeris no exist(%v)",f);
        return 0;
    }
    kFileOpt.brdc_=f;
    LOG(DEBUG)<<"Broadcast ephemeris file: "<<kFileOpt.brdc_;
    f[0]='\0';

    int ion=kPrcOpt.GNSS_opt_.ion_opt_==ION_TEC||kPrcOpt.GNSS_opt_.ion_opt_==ION_CONST;
    if(ion){
        sprintf(f,"%s%cCODG%03d0.%02dI",prods_dir.c_str(),sep_,doy_,yy_);
        if((access(f,0))==-1){
            defaultLogger->error("Ionsphere rinex no exist(%v)",f);
            return 0;
        }
        kFileOpt.ion_=f;f[0]='\0';
        LOG(DEBUG)<<"Ionsphere rinex file: "<<kFileOpt.ion_;
    }

    if(kPrcOpt.GNSS_opt_.cb_prd_type_==CBC_DCB&&kPrcOpt.GNSS_opt_.cb_prd_ac_==BIA_CAS){
        sprintf(f,"%s%cCAS0MGXRAP_%d%03d0000_01D_01D_DCB.BSX",prods_dir.c_str(),sep_,year_,doy_);
        if(access(f,0)==-1){
            defaultLogger->error("CAS Multi-GNSS differential code biases no exist(%v)",f);
            return 0;
        }
        kFileOpt.bia_=f;f[0]='\0';
        LOG(DEBUG)<<"CAS Multi-GNSS differential code biases file: "<<kFileOpt.bia_;
    }
    else if(kPrcOpt.GNSS_opt_.cb_prd_type_==CBC_OSB){
        if(kPrcOpt.GNSS_opt_.cb_prd_ac_==BIA_CODE){
            sprintf(f,"%s%cCOD0MGXFIN_%d%03d0000_01D_01D_OSB.BIA",prods_dir.c_str(),sep_,year_,doy_);
        }
        else if(kPrcOpt.GNSS_opt_.cb_prd_ac_==BIA_CAS){
            sprintf(f,"%s%cCAS0MGXRAP_%d%03d0000_01D_01D_OSB.BIA",prods_dir.c_str(),sep_,year_,doy_);
        }
        else if(kPrcOpt.GNSS_opt_.cb_prd_ac_==BIA_CNES){
            sprintf(f,"%s%ccnt%d%d",prods_dir.c_str(),sep_,week_,wod_);
        }
        if(access(f,0)==-1){
            defaultLogger->error("CAS Multi-GNSS Observation-Specific code biases no exist(%v)",f);
            return 0;
        }
        kFileOpt.bia_=f;f[0]='\0';
        LOG(DEBUG)<<"Multi-GNSS Observation-Specific code biases file: "<<kFileOpt.ion_;
    }
    if(kPrcOpt.mode_opt_==MDOPT_KINE){
        if(kSolOpt.err_fmt_==true){
            sprintf(f,"%s%c%s%03d0.%s",prods_dir.c_str(),sep_,kSite.c_str(),doy_,"ref");
            if((access(f,0))==-1){
                cout<<"[ERROR]-"<<__FUNCTION__<<": "<<"no found "<<f<<endl;
                defaultLogger->error("Reference solution file(%v) no exist, change to normal solution format",f);
            }
            kFileOpt.ref_=f;
            LOG(DEBUG)<<"Reference solution file: "<<f;
        }
    }
    return 1;
}

int MatchFile::MatchPrec() {
    string pres_dir;
    int wod1,wk1;
    int wod2,wk2;
    int ppp;
    char f[MAXSTRPATH]={'\0'};

    if(!kFileOpt.customize_file_){
        sprintf(f,"%s%c%d%c%d%c%s",dir_.c_str(),sep_,year_,sep_,week_,sep_,kListACsStr[kPrcOpt.GNSS_opt_.prods_ac_].c_str());
        pres_dir=f;f[0]='\0';
    }else pres_dir=kFileOpt.customize_dir_;

    LOG(DEBUG)<<"Precise products dir: "<<pres_dir;
    if(kPrcOpt.GNSS_opt_.sat_eph_==PRECISE){
        sprintf(f,"%s%c%s%d%d%s",pres_dir.c_str(),sep_,kListACsStr[kPrcOpt.GNSS_opt_.prods_ac_].c_str(),week_,wod_,".sp3");
        if((access(f,0))==-1){
            defaultLogger->error("On the day precise orbit file(%v) no exist",f);
            return 0;
        }
        kFileOpt.sp3_[1]=f;f[0]='\0';
        LOG(DEBUG)<<"On the day precise orbit file: "<<kFileOpt.sp3_[1];

        if(wod_==0) {
            wod1=6;wk1=week_-1;
        }
        else {
            wod1=wod_-1;wk1=week_;
        }
        sprintf(f,"%s%c%s%d%d%s",pres_dir.c_str(),sep_,kListACsStr[kPrcOpt.GNSS_opt_.prods_ac_].c_str(),wk1,wod1,".sp3");
        if((access(f,0))==-1){
            LOG(DEBUG)<<"Before the day precise orbit file no exist: "<<f;
        }
        else{
            kFileOpt.sp3_[0]=f;f[0]='\0';
            LOG(DEBUG)<<"Before the day precise orbit file: "<<kFileOpt.sp3_[0];
        }

        if(wod_==6) {
            wod2=0;wk2=week_+1;
        }
        else {
            wod2=wod_+1;wk2=week_;
        }
        sprintf(f,"%s%c%s%d%d%s",pres_dir.c_str(),sep_,kListACsStr[kPrcOpt.GNSS_opt_.prods_ac_].c_str(),wk2,wod2,".sp3");
        if((access(f,0))==-1){
            LOG(DEBUG)<<"After the day precise orbit file no exist: "<<f;
        }
        else{
            kFileOpt.sp3_[2]=f;f[0]='\0';
            LOG(DEBUG)<<"After the day precise orbit file: "<<kFileOpt.sp3_[2];
        }

        sprintf(f,"%s%c%s%d%d%s",pres_dir.c_str(),sep_,kListACsStr[kPrcOpt.GNSS_opt_.prods_ac_].c_str(),week_,wod_,".clk");
        if((access(f,0))==-1){
            LOG(ERROR)<<"On the day precise clock file no exist: "<<f;
            return 0;
        }
        kFileOpt.clk_=f;f[0]='\0';
        LOG(DEBUG)<<"On the day precise clock file: "<<kFileOpt.clk_;
    }
    if(kPrcOpt.GNSS_opt_.tid_opt_>TIDE_OFF){
        sprintf(f,"%s%c%s%d%d%s",pres_dir.c_str(),sep_,kListACsStr[kPrcOpt.GNSS_opt_.prods_ac_].c_str(),week_,wod_,".erp");
        if((access(f,0))!=-1){
            kFileOpt.erp_=f;f[0]='\0';
            LOG(DEBUG)<<"Erp file: "<<kFileOpt.erp_;
        }
        else{
            sprintf(f,"%s%c%d%c%d%igs%2dP%d.erp",dir_.c_str(),sep_,year_,sep_,week_,yy_,week_);
            if((access(f,0))==-1){
                LOG(DEBUG)<<"Erp file no exist: "<<f;
            }
            kFileOpt.erp_=f;
            LOG(DEBUG)<<"Erp file: "<<kFileOpt.erp_;
        }
    }
    return 1;
}

int MatchFile::MatchComn() {
    string cmn_prods;
    char f[MAXSTRPATH]={'\0'};

    if(!kFileOpt.customize_file_){
        sprintf(f,"%s%ccmnprods",dir_.c_str(),sep_);
        cmn_prods=f;f[0]='\0';
    }else cmn_prods=kFileOpt.customize_dir_;

    LOG(DEBUG)<<"cmnprods dir: "<<cmn_prods;
    if(kPrcOpt.GNSS_opt_.sat_pcv_||kPrcOpt.GNSS_opt_.rec_pcv_)
    {
        sprintf(f,"%s%cigs14_2097.atx",cmn_prods.c_str(),sep_);
        if((access(f,0))==-1){
            LOG(DEBUG)<<"Antex file no exist: "<<f;
        }
        kFileOpt.atx_=f;f[0]='\0';
        LOG(DEBUG)<<"Antex file: "<<kFileOpt.atx_;
    }
    if(kPrcOpt.GNSS_opt_.tid_opt_){
        sprintf(f,"%s%cocnload.blq",cmn_prods.c_str(),sep_);
        if((access(f,0))==-1){
            LOG(DEBUG)<<"Blq file no exist: "<<f;
        }
        kFileOpt.blq_=f;
        LOG(DEBUG)<<"Blq file: "<<kFileOpt.blq_;
    }
    return 1;
}

int MatchFile::MatchDCB() {
    int dgnss=(kPrcOpt.mode_==DGPS||kPrcOpt.mode_==PPK)||
              (kPrcOpt.mode_opt_==MDOPT_DGPS||kPrcOpt.mode_opt_==MDOPT_PPK)||kPrcOpt.mode_==FIX;
    string dcb_dir;
    char f[MAXSTRPATH]={'\0'};

    if(!dgnss){
        sprintf(f,"%s%c%d%cdcb%c*.DCB",dir_.c_str(),sep_,year_,sep_,sep_);
        kFileOpt.bia_=f;
        LOG(DEBUG)<<"Normal DCb file: "<<kFileOpt.bia_;
    }
}

int MatchFile::MatchPos() {
    int dgnss=(kPrcOpt.mode_==DGPS||kPrcOpt.mode_==PPK)||
              (kPrcOpt.mode_opt_==MDOPT_DGPS||kPrcOpt.mode_opt_==MDOPT_PPK)||kPrcOpt.mode_==FIX;
    int flag;
    string crd_file;
    char f[MAXSTRPATH]={'\0'};

    if(dgnss){
        if(kFileOpt.customize_file_){
            sprintf(f,"%s%c%d%c%d%c%03d%c%s%csite.sta",dir_.c_str(),sep_,year_,sep_,week_,sep_,doy_,sep_,"prods",sep_);
        }else sprintf(f,"%s%c%d%c%d%c%03d%c%s%csite.sta",dir_.c_str(),sep_,year_,sep_,week_,sep_,doy_,sep_,"prods",sep_);
        crd_file=f;f[0]='\0';
        if((access(crd_file.c_str(),0))!=-1){
            if(!GetPosFrmSta(crd_file)){
                defaultLogger->error("Reference station position no setting");
                return 0;
            }
        }
    }
    else{
        if(kFileOpt.customize_file_){
            sprintf(f,"%s%cigs%2dP%d.snx",dir_.c_str(),sep_,yy_,week_);
        }else sprintf(f,"%s%c%d%c%d%cigs%2dP%d.snx",dir_.c_str(),sep_,year_,sep_,week_,sep_,yy_,week_);
        crd_file=f;f[0]='\0';
        if((access(crd_file.c_str(),0))!=-1){
            flag=GetPosFrmSnx(crd_file);
            if(flag) return 1;
        }
        else flag=0;

        if(kFileOpt.customize_file_){
            sprintf(f,"%s%csite.crd",dir_.c_str(),sep_);
        }else sprintf(f,"%s%ccmnprods%csite.crd",dir_.c_str(),sep_,sep_);
        crd_file.clear();
        crd_file=f;
        if((access(crd_file.c_str(),0))!=-1){
            flag=GetPosFrmCrd(crd_file);
        }
        else flag=0;

        if(!flag && kSolOpt.err_fmt_){
            LOG(WARNING)<<"Rover station position no setting";
            kSolOpt.err_fmt_= false;
        }
    }
    return 1;
}

int MatchFile::GetPosFrmSnx(string file) {
    ifstream inf;
    int flag=0;
    Vec3 *pos;
    string buff,site;

    inf.open(file);
    if(inf.is_open()){
        if(kPrcOpt.mode_==FIX){
            pos=&kPrcOpt.GNSS_opt_.pos_opt_.rover_;
        }
        else{
            pos=&kSolOpt.ref_sol_.pos;
        }
        while(getline(inf,buff)&&!inf.eof()){
            if(buff.find("+SOLUTION/ESTIMATE")!=string::npos) flag=1;
            if(buff.find("-SOLUTION/ESTIMATE")!=string::npos) {flag=0;break;}
            if(flag){
                site=buff.substr(14,4);
                int a=buff.find("STAX");
                if((strcasecmp(site.c_str(),kSite.c_str())==0)&&(buff.find("STAX")!=-1)){
                    Str2Double(buff.substr(47,21),pos->i_);
                    continue;
                }
                if((strcasecmp(site.c_str(),kSite.c_str())==0)&&(buff.find("STAY")!=-1)){
                    Str2Double(buff.substr(47,21),pos->j_);
                    continue;
                }
                if((strcasecmp(site.c_str(),kSite.c_str())==0)&&(buff.find("STAZ")!=-1)){
                    Str2Double(buff.substr(47,21),pos->k_);
                    break;
                }
            }
        }
    }
    inf.close();
    char ref_coor[MAXBUFF]={'\0'};
    sprintf(ref_coor,"%s %12.3f %12.3f %12.3f",kSite.c_str(),pos->i_,pos->j_,pos->k_);
    LOG(INFO)<<ref_coor;
    return 1;
}

int MatchFile::GetPosFrmCrd(string file) {
    ifstream inf;
    Vec3 *pos;
    string buff,site;

    inf.open(file);
    if(inf.is_open()) {
        if (kPrcOpt.mode_ == FIX) {
            pos = &kPrcOpt.GNSS_opt_.pos_opt_.rover_;
        } else {
            pos = &kSolOpt.ref_sol_.pos;
        }
        while(getline(inf,buff)&&!inf.eof()){
            site=buff.substr(0,4);
            if((strcasecmp(site.c_str(),kSite.c_str()))==0){
                sscanf(buff.c_str(),"%s %lf %lf %lf",(char *)site.c_str(),&pos->i_,&pos->j_,&pos->k_);
                break;
            }
        }
        if(pos->i_*pos->j_*pos->k_==0.0) return 0;
    }
    inf.close();
    char ref_coor[MAXBUFF]={'\0'};
    sprintf(ref_coor,"%s %12.3f %12.3f %12.3f",kSite.c_str(),pos->i_,pos->j_,pos->k_);
    LOG(INFO)<<ref_coor;
    return 1;
}

int MatchFile::GetPosFrmSta(string file) {
    int dgnss=(kPrcOpt.mode_==DGPS||kPrcOpt.mode_==PPK)||(kPrcOpt.mode_opt_==MDOPT_DGPS||kPrcOpt.mode_opt_==MDOPT_PPK);
    ifstream inf;
    int fmt;
    Vec3 pos;
    string buff,site;

    inf.open(file);
    if(inf.is_open()) {
        while(getline(inf,buff)&&!inf.eof()){
            if(buff.substr(0,1)=="#") continue;
            if(dgnss){
                if(buff.substr(0,4)=="base"){
                    sscanf((char*)buff.c_str(),"%s %d %lf %lf %lf",(char*)site.c_str(),&fmt,&pos.i_,&pos.j_,&pos.k_);
                }
                if(pos.i_*pos.j_*pos.k_==0.0) return 0;
                kPrcOpt.GNSS_opt_.pos_opt_.base_=pos;
            };

            if(kPrcOpt.mode_==FIX||kPrcOpt.mode_opt_==MDOPT_PPP){
                if(buff.substr(0,4)==kSite) {
                    sscanf(buff.c_str(),"%s %d %lf %lf %lf",(char*)site.c_str(),&fmt,&pos.i_,&pos.j_,&pos.k_);
                }
                if(pos.i_*pos.j_*pos.k_==0.0) return 0;
                kPrcOpt.GNSS_opt_.pos_opt_.rover_=pos;
            }
        }
    }
    char ref_coor[MAXBUFF]={'\0'};
    sprintf(ref_coor,"%s %12.3f %12.3f %12.3f",kSite.c_str(),pos.i_,pos.j_,pos.k_);
    LOG(INFO)<<ref_coor;
    inf.close();
    return 1;
}

int MatchFile::MatchOut() {
    string sys_str;
    string out_dir,mode_dir,par_dir,name;
    char f[MAXSTRPATH]={'\0'};

    if(kPrcOpt.GNSS_opt_.nav_sys_&GPS){
        sys_str+="G";
    }
    if(kPrcOpt.adj_opt_.bds_isb_){
        if(kPrcOpt.GNSS_opt_.nav_sys_&BD2){
            sys_str+="B2";
        }
        if(kPrcOpt.GNSS_opt_.nav_sys_&BD3){
            sys_str+="B3";
        }
    }
    else{
        if(kPrcOpt.GNSS_opt_.nav_sys_&BD2||kPrcOpt.GNSS_opt_.nav_sys_&BD3){
            sys_str+="B";
        }
    }
    if(kPrcOpt.GNSS_opt_.nav_sys_&GAL){
        sys_str+="E";
    }
    if(kPrcOpt.GNSS_opt_.nav_sys_&GLO){
        sys_str+="R";
    }
    if(kPrcOpt.GNSS_opt_.nav_sys_&QZS){
        sys_str+="J";
    }
    if(kPrcOpt.GNSS_opt_.sat_eph_==PRECISE){
        if(kFileOpt.customize_file_){
            sprintf(f,"%s%cresult_%s%c",dir_.c_str(),sep_,kListACsStr[kPrcOpt.GNSS_opt_.prods_ac_].c_str(),sep_);
        }else sprintf(f,"%s%c%d%c%d%c%03d%cresult_%s%c",dir_.c_str(),sep_,year_,sep_,week_,sep_,doy_,sep_,kListACsStr[kPrcOpt.GNSS_opt_.prods_ac_].c_str(),sep_);
    }
    else{
        if(kFileOpt.customize_file_){
            sprintf(f,"%s%cresult_brd%c",dir_.c_str(),sep_);
        }else sprintf(f,"%s%c%d%c%d%c%03d%cresult_brd%c",dir_.c_str(),sep_,year_,sep_,week_,sep_,doy_,sep_,sep_);
    }
    out_dir=f;
    CreateDir(out_dir.c_str());
    sprintf(f,"%s%c%s_%s%c",out_dir.c_str(),sep_,kListModesStr[kPrcOpt.mode_].c_str(),kListModOptsStr[kPrcOpt.mode_opt_].c_str(),sep_);
    mode_dir=f;f[0]='\0';
    CreateDir(mode_dir.c_str());
    sprintf(f,"%s%c%s_%s_%s%c",mode_dir.c_str(),sep_,sys_str.c_str(),kListFrqsStr[kPrcOpt.GNSS_opt_.use_frq_].c_str(),kListIonsStr[kPrcOpt.GNSS_opt_.ion_opt_].c_str(),sep_);
    par_dir=f;f[0]='\0';
    CreateDir(par_dir.c_str());
    sprintf(f,"%s%s.pos",par_dir.c_str(),kSite.c_str());
    kFileOpt.sol_=f;f[0]='\0';
    LOG(DEBUG)<<"Solution file path: "<<kFileOpt.sol_;
}

int MatchFile::MatchFileAuto() {
    if(!MatchProd()) return 0;
    if(!MatchPrec()) return 0;
    if(!MatchComn()) return 0;
    if(kPrcOpt.GNSS_opt_.cb_prd_type_==CBC_DCB&&kPrcOpt.GNSS_opt_.cb_prd_ac_==BIA_CODE)
        if(!MatchDCB()) return 0;
    if(!MatchPos())  return 0;
    MatchOut();
}

static const string kSysCodes = "GCERJ";	/* satellite system codes */
static const string kObsCodes = "CLDS";    /* obs type codes */
static const string kFrqCodes = "1256789"; /* frequency codes */
ReadRinex::ReadRinex(void) {
    ver_=3.0;
    ts_=Time();
    te_=Time();
    tint_=0.0;
    sys_mask_=kPrcOpt.GNSS_opt_.nav_sys_;
}

ReadRinex::ReadRinex(string& file) {
    inf_.open(file,ios::in);
    ver_=3.0;
    ts_=Time();
    te_=Time();
    tint_=0.0;
    SetSysMark();
}

ReadRinex::~ReadRinex() {
    inf_.close();
}

int ReadRinex::TestOpen() {
    return inf_.is_open();
}

void ReadRinex::CloseFile() {
    if(inf_.is_open()) inf_.close();
}

void ReadRinex::SetSysMark() {
    GNSSOpt gopt=kPrcOpt.GNSS_opt_;

    if(gopt.nav_sys_&GPS) {sys_mask_|=GPS;LOG(DEBUG)<<"Reading GPS Broadcast ephemeris...";}
    if(gopt.nav_sys_&BD2) {sys_mask_|=BD2;LOG(DEBUG)<<"Reading BD2 Broadcast ephemeris...";}
    if(gopt.nav_sys_&BD3) {sys_mask_|=BD3;LOG(DEBUG)<<"Reading BD3 Broadcast ephemeris...";}
    if(gopt.nav_sys_&GAL) {sys_mask_|=GAL;LOG(DEBUG)<<"Reading GAL Broadcast ephemeris...";}
    if(gopt.nav_sys_&GLO) {sys_mask_|=GLO;LOG(DEBUG)<<"Reading GLO Broadcast ephemeris...";}
    if(gopt.nav_sys_&QZS) {sys_mask_|=QZS;LOG(DEBUG)<<"Reading QZS Broadcast ephemeris...";}
}

int ReadRinex::ReadRnxHead() {
    while(getline(inf_,buff_)&&!inf_.eof()){

        if(buff_.find("RINEX VERSION / TYPE")!=string::npos){
            Str2Double(buff_.substr(0,9),ver_);
            type_=buff_.substr(20,1);
            switch(buff_[40]){
                case ' ':
                case 'G': sat_sys_=GPS; time_sys_=TIMESYSGPS;break;
                case 'C': sat_sys_=BD2|BD3; time_sys_=TIMESYSBDS;break;
                case 'E': sat_sys_=GAL; time_sys_=TIMESYSGAL;break;
                case 'R': sat_sys_=GLO; time_sys_=TIMESYSUTC;break;
                case 'J': sat_sys_=QZS; time_sys_=TIMESYSQZS;break;
                case 'M': sat_sys_=NONE;time_sys_=TIMESYSGPS;break;
                default: break;
            }
            break;
        }
    }
    return 1;
}

int ReadRinex::ReadHead() {
    return 0;
}

int ReadRinex::ReadBody() {
    return 0;
}

ReadRnxObs::ReadRnxObs() {

}

ReadRnxObs::ReadRnxObs(string &file):
    ReadRinex(file) {

}

ReadRnxObs::ReadRnxObs(string &file, int rcv):
    ReadRinex(file),rcv_(rcv){
}


ReadRnxObs::~ReadRnxObs() {
    inf_.close();
}

#define OBSTYPE_CODE    0
#define OBSTYPE_PHASE   1
#define OBSTYPE_DOPPLER 2
#define OBSTYPE_SNR     3

void ReadRnxObs::SaveSlips(unsigned char slips[][NFREQ],ObsData& data) {
    int i;
    for(i=0;i<NFREQ;i++) {
        if(data.LLI_[i]&1) slips[data.sat_.sat_no_-1][i]|=1;
    }
}

void ReadRnxObs::ReStoreSlips(unsigned char slips[][NFREQ], ObsData &data){
    int i;
    for(i=0;i<NFREQ;i++) {
        if(slips[data.sat_.sat_no_-1][i]&1) data.LLI_[i]|=1;
        slips[data.sat_.sat_no_-1][i]=0;
    }
}

void ReadRnxObs::SetSigIndex(int sys, string type_obs[MAXOBSTYPE], Signal& sig_idx) {
    size_t p;
    string str;
    double shift;
    int i,j,num_sig;
    int k_code=-1,k_phase=-1,k_doppler=-1,k_snr=-1;

    for(i=num_sig=0;type_obs[i][0];i++,num_sig++){
        sig_idx.code_[i]=sig_idx.Signal2Code(type_obs[i].substr(1),sig_idx.frq_+i,sys);
        sig_idx.type_[i]=((p=kObsCodes.find(type_obs[i][0]))!=string::npos)?(int)p:0;
        sig_idx.pri_[i]=sig_idx.GetSignalPri(sys,sig_idx.code_[i]);
    }

    for(i=0;i<NFREQ;i++){
        for(j=0;j<num_sig;j++){
            if(sig_idx.type_[j]==OBSTYPE_CODE){
                if(sig_idx.frq_[j]==i+1&&sig_idx.pri_[j]&&
                   (k_code<0||sig_idx.pri_[j]>=sig_idx.pri_[k_code])){
                    k_code=j;
                }
            }
            else if(sig_idx.type_[j]==OBSTYPE_PHASE){
                if(sig_idx.frq_[j]==i+1&&sig_idx.pri_[j]&&
                   (k_phase<0||sig_idx.pri_[j]>=sig_idx.pri_[k_phase])){
                    k_phase=j;
                }
            }
            else if(sig_idx.type_[j]==OBSTYPE_DOPPLER){
                if(sig_idx.frq_[j]==i+1&&sig_idx.pri_[j]&&
                   (k_doppler<0||sig_idx.pri_[j]>=sig_idx.pri_[k_doppler])){
                    k_doppler=j;
                }
            }
            else if(sig_idx.type_[j]==OBSTYPE_SNR){
                if(sig_idx.frq_[j]==i+1&&sig_idx.pri_[j]&&
                   (k_snr<0||sig_idx.pri_[j]>=sig_idx.pri_[k_snr])){
                    k_snr=j;
                }
            }
        }

        if(k_code>=0)    sig_idx.pos_[k_code]=i;
        if(k_phase>=0)   sig_idx.pos_[k_phase]=i;
        if(k_doppler>=0) sig_idx.pos_[k_doppler]=i;
        if(k_snr>=0)     sig_idx.pos_[k_snr]=i;
        k_code=k_phase=k_doppler=k_snr=-1;
    }
    sig_idx.n_=num_sig;
}


int ReadRnxObs::ReadObsBody(int &flag, vector<ObsData> &data) {
    ObsData obs;
    Time t;
    int i=0,n=0,num_sat=0;
    vector<SatMask> sats(MAXOBS,SatMask());

    // decode epoch obs line by line
    while(getline(inf_,buff_)&&!inf_.eof()){
        if(i==0){
            if((num_sat=DecodeObsEpoch(t,flag,sats))<=0){
                continue;
            }
        }
        else if(flag<=2||flag==6){
            obs.ReSet();
            obs.sig_recep_=t;

            // decode obs data
            if(DecodeObsData(obs)&&data.size()<MAXOBS){
                data.push_back(obs);n++;
            }
        }
        if(++i>num_sat) return n;
    }
    return 0;
}

int ReadRnxObs::DecodeObsEpoch(Time &t, int &flag, vector<SatMask> &sats) {
    int i,j,n=0;
    string sat_id;

    if (buff_.length()<32) return 0;

    if (ver_<=2.99) { /* ver.2 */
        Str2Int(buff_.substr(29,3),n);
        if (n<=0) return 0;

        /* epoch flag: 3:new site,4:header info,5:external event */
        Str2Int(buff_.substr(28,1),flag);
        if (3<=flag&&flag<=5) return n;

        if (t.Str2Time(buff_.substr(0,26))!=0) return 0;

        for (i=0,j=32;i<n;i++,j+=3) {
            if (j>=68) {
                if (!getline(inf_,buff_)) break; /* read next line */
                j=32;
            }
            if (i<MAXOBS) {
                sat_id=buff_.substr(j,3);
                sats[i]=SatMask(sat_id);
            }
        }
    }
    else { /* ver.3 */
        Str2Int(buff_.substr(32,3),n);
        if (n<=0) return 0;

        Str2Int(buff_.substr(31,1),flag);
        if (3<flag&&flag<=5) return n;

        if (buff_[0]!='>'||t.Str2Time(buff_.substr(1,28))!=0)
            return 0;
    }
    return n;
}

int ReadRnxObs::DecodeObsData(ObsData &obs) {
    Signal *ind;
    double val[MAXOBSTYPE] = { 0.0 };
    unsigned char lli[MAXOBSTYPE] = { 0 };
    string sat_id;
    int i,j,n,m,num,stat=1,p[MAXOBSTYPE],k[16],l[16],BD3_flag=0;

    if (ver_>2.99){ /* ver.3 */
        sat_id=buff_.substr(0,3);
        if(!sat_id.compare(1,1," ")) sat_id[1]='0';
        obs.sat_=SatMask(sat_id);
    }
    if (!obs.sat_.sat_no_) stat=0;
    else if (!(obs.sat_.sat_sys_&sys_mask_)) stat=0;

    switch(obs.sat_.sat_sys_){
        case BD2: ind=sig_index_+INDEXBD2; break;
        case BD3: ind=sig_index_+INDEXBD3; break;
        case GAL: ind=sig_index_+INDEXGAL; break;
        case GLO: ind=sig_index_+INDEXGLO; break;
        case QZS: ind=sig_index_+INDEXQZS; break;
        default:  ind=sig_index_;   break;
    }
    for (i=0,j=ver_<=2.99?0:3;i<ind->n_&&j+15<buff_.length();i++,j+=16){
        if (stat){
            Str2Double(buff_.substr(j,14),val[i]);
            val[i] += ind->shift_[i];
            Str2Int(buff_.substr(j+14,1),num);
            lli[i]=(unsigned char)num&3;
        }
        /* if the last data, read the next line */
        if (ver_<=2.99&&j+17>80&&i+1<ind->n_) { /* ver.2 */
            if (!getline(inf_,buff_)) break; /* read next line */
            j=-16;
        }
    }
    if (!stat) return 0;

    for (i=n=m=0;i<ind->n_;i++){
        p[i]=ver_<=2.11?ind->frq_[i]-1:ind->pos_[i];

        if (ind->type_[i]==0&&p[i]==0) k[n++]=i; /* C1? index */
        if (ind->type_[i]==0&&p[i]==1) l[m++]=i; /* C2? index */
    }
    if (ver_<=2.11){
        /* if multiple codes (C1/P1,C2/P2), select higher priority */
        if (n>=2) {
            if (val[k[0]]==0.0&&val[k[1]]==0.0) {
                p[k[0]]=-1; p[k[1]]=-1;
            }
            else if (val[k[0]]!=0.0&&val[k[1]]==0.0) {
                p[k[0]]=0; p[k[1]]=-1;
            }
            else if (val[k[0]]==0.0&&val[k[1]]!=0.0) {
                p[k[0]]=-1; p[k[1]]=0;
            }
            else if (ind->pri_[k[1]]>ind->pri_[k[0]]) {
                p[k[1]]=0; p[k[0]]=NEXOBS<1?-1:NFREQ;
            }
            else {
                p[k[0]]=0; p[k[1]]=NEXOBS<1?-1:NFREQ;
            }
        }
        if (m>=2) {
            if (val[l[0]]==0.0&&val[l[1]]==0.0) {
                p[l[0]]=-1; p[l[1]]=-1;
            }
            else if (val[l[0]]!=0.0&&val[l[1]]==0.0) {
                p[l[0]]=1; p[l[1]]=-1;
            }
            else if (val[l[0]]==0.0&&val[l[1]]!=0.0) {
                p[l[0]]=-1; p[l[1]]=1;
            }
            else if (ind->pri_[l[1]]>ind->pri_[l[0]]) {
                p[l[1]]=1; p[l[0]]=NEXOBS<2?-1:NFREQ+1;
            }
            else {
                p[l[0]]=1; p[l[1]]=NEXOBS<2?-1:NFREQ+1;
            }
        }
    }
    /* save obs data */
    for (i=0;i<ind->n_;i++){
        if(p[i]<0||val[i]==0.0) continue;
        switch(ind->type_[i]){
            case 0: obs.P_[p[i]]=val[i]; obs.code_[p[i]]=ind->code_[i]; break;
            case 1: obs.L_[p[i]]=val[i]; obs.LLI_ [p[i]]=lli[i];       break;
            case 2: obs.D_[p[i]]=(float)val[i];                        break;
            case 3: obs.SNR_[p[i]]=(unsigned char)(val[i]*4.0+0.5);    break;
        }
    }

    return 1;
}

void ReadRnxObs::InitReadObs(string file, Time ts, Time te, int tint) {
    inf_.open(file,ios::in);
    ts_=ts;te_=te;tint_=tint;
}

#define MAXPOSHEAD  1024                /* max head line position */
int ReadRnxObs::ReadHead(Nav* nav,Sta* sta) {
    int i,j,k;
    int num_sigs,idx_sig,prn,fcn,num_lines=0;
    if(!inf_.is_open()){
        return 0;
    }

    const string kDefCodes[]={
            "CWX    ",  /* GPS: L125____ */
            "X  XX  ",  /* BDS: L1__67__ */
            "X XXXX ",  /* GAL: L1_5678_ */
            "CC     ",  /* GLO: L12_____ */
            "CXXX   ",  /* QZS: L1256___ */
    };

    if(!ReadRnxHead()){return 0;}

    while(getline(inf_,buff_)&&!inf_.eof()){
        if (buff_.find("MARKER NAME")!=string::npos && sta)
            sta->name_=buff_.substr(0,10);
        else if (buff_.find("MARKER NUMBER")!=string::npos && sta)
            sta->marker_=buff_.substr(0,20);
        else if (buff_.find("MARKER TYPE")!=string::npos) continue;
        else if (buff_.find("OBSERVER / AGENCY")!=string::npos) continue;
        else if (buff_.find("REC # / TYPE / VERS")!=string::npos && sta){
            sta->ant_seri_=buff_.substr(0,20);
            sta->rec_type_=buff_.substr(20,20);
            sta->firm_ver_=buff_.substr(40,20);
        }
        else if (buff_.find("ANT # / TYPE")!=string::npos && sta){
            sta->ant_seri_=buff_.substr(0,20);
            sta->ant_desc_=buff_.substr(20,20);
        }
        else if (buff_.find("APPROX POSITION XYZ")!=string::npos && sta)
            for (i=0;i<3;i++){
                Str2Double(buff_.substr(i*14,14),sta->pos_[i]);
            }
        else if (buff_.find("ANTENNA: DELTA H/E/N")!=string::npos && sta){
            Str2Double(buff_.substr(0,14),sta->del_[2]);  /* h */
            Str2Double(buff_.substr(14,14),sta->del_[0]); /* e */
            Str2Double(buff_.substr(28,14),sta->del_[1]); /* n */
        }
        else if (buff_.find("ANTENNA: DELTA X/Y/Z")!=string::npos) continue; /* opt ver.3 */
        else if (buff_.find("ANTENNA: PHASECENTER")!=string::npos) continue; /* opt ver.3 */
        else if (buff_.find("ANTENNA: B.SIGHT XYZ")!=string::npos) continue; /* opt ver.3 */
        else if (buff_.find("ANTENNA: ZERODIR AZI")!=string::npos) continue; /* opt ver.3 */
        else if (buff_.find("ANTENNA: ZERODIR XYZ")!=string::npos) continue; /* opt ver.3 */
        else if (buff_.find("CENTER OF MASS: XYZ" )!=string::npos) continue; /* opt ver.3 */
        else if (buff_.find("SYS / # / OBS TYPES" )!=string::npos) { /* ver.3 */
            if((i=kSysCodes.find(buff_[0]))==string::npos){
                LOG(WARNING)<<"Invalid satellite system: "<<buff_[0];
                continue;
            }
            Str2Int(buff_.substr(3,3),num_sigs);
            for(j=idx_sig=0,k=7;j<num_sigs;j++,k+=4){
                if(k>58){
                    if(!getline(inf_,buff_)) break;
                    k=7;
                }
                if(idx_sig<MAXOBSTYPE-1) type_obs_[i][idx_sig++]=buff_.substr(k,3);
            }
            if(i==1){
                if(ver_<3.04){
                    for(j=0;j<num_sigs;j++){
                        if(type_obs_[i][j][1]=='1'){
                            type_obs_[i][j][1]='2';
                            LOG(INFO)<<"BD2 change C1x to C2x";
                        }
                    }
                }
            }
        }
        else if (buff_.find("WAVELENGTH FACT L1/2")!=string::npos) continue; /* opt ver.2 */
        else if (buff_.find("# / TYPES OF OBSERV")!=string::npos) { /* ver.2 */

        }
        else if (buff_.find("SIGNAL STRENGTH UNIT")!=string::npos) continue; /* opt ver.3 */
        else if (buff_.find("INTERVAL"			 )!=string::npos) continue; /* opt */
        else if (buff_.find("TIME OF FIRST OBS"   )!=string::npos) {
            if      (!buff_.compare(48,3,"GPS")) time_sys_=TIMESYSGPS;
            else if (!buff_.compare(48,3,"BDT")) time_sys_=TIMESYSBDS; /* ver.3.02 */
            else if (!buff_.compare(48,3,"GAL")) time_sys_=TIMESYSGAL;
            else if (!buff_.compare(48,3,"GLO")) time_sys_=TIMESYSUTC;
            else if (!buff_.compare(48,3,"QZS")) time_sys_=TIMESYSQZS; /* ver.3.02 */
        }
        else if (buff_.find("TIME OF LAST OBS"    )!=string::npos) continue; /* opt */
        else if (buff_.find("RCV CLOCK OFFS APPL" )!=string::npos) continue; /* opt */
        else if (buff_.find("SYS / DCBS APPLIED"  )!=string::npos) continue; /* opt ver.3 */
        else if (buff_.find("SYS / PCVS APPLIED"  )!=string::npos) continue; /* opt ver.3 */
        else if (buff_.find("SYS / SCALE FACTOR"  )!=string::npos) continue; /* opt ver.3 */
        else if (buff_.find("SYS / PHASE SHIFTS"  )!=string::npos) continue; /* ver.3.01 */
        else if (buff_.find("GLONASS SLOT / FRQ #")!=string::npos && nav) { /* ver.3.02 */
            for (i=0;i<8;i++) {
                if (buff_.compare(8*i+4,1,"R")!=0||!buff_.compare(8*i+8,2,"  ")) continue;
                Str2Int(buff_.substr(8*i+5,2),prn);
                Str2Int(buff_.substr(8*i+8,2),fcn);
                if (1<=prn&&prn<=MAXPRNGLO) nav->glo_frq_num_[prn-1]=fcn+8;
            }
        }
        else if (buff_.find("GLONASS COD/PHS/BIS" )!=string::npos && nav) {  /* ver.3.02 */
            for (i=0;i<4;i++) {
                if      (buff_.compare(13*i+1,3,"C1C"))
                    Str2Double(buff_.substr(13*i+5,8),nav->glo_cp_bias_[0]);
                else if (buff_.compare(13*i+1,3,"C1P"))
                    Str2Double(buff_.substr(13*i+5,8),nav->glo_cp_bias_[1]);
                else if (buff_.compare(13*i+1,3,"C2C"))
                    Str2Double(buff_.substr(13*i+5,8),nav->glo_cp_bias_[2]);
                else if (buff_.compare(13*i+1,3,"C2P"))
                    Str2Double(buff_.substr(13*i+5,8),nav->glo_cp_bias_[3]);
            }
        }
        else if (buff_.find("LEAP SECONDS")!=string::npos && nav) {/* opt */
            Str2Int(buff_.substr(0,6),nav->leaps_);
        }
        else if (buff_.find("# OF SALTELLITES")!=string::npos) continue;/* opt */
        else if (buff_.find("PRN / # OF OBS"  )!=string::npos) continue;/* opt */
        else if (buff_.find("PGM / RUN BY / DATE")!=string::npos) continue;
        else if (buff_.find("COMMENT" )!=string::npos) continue;
        if (buff_.find("END OF HEADER")!=string::npos)
            break;
        if (++num_lines>=MAXPOSHEAD && type_.compare(" ")==0) return 0; /* no rinex file */
    }

    SetSigIndex(GPS,type_obs_[0],sig_index_[INDEXGPS]);
    SetSigIndex(BD2,type_obs_[1],sig_index_[INDEXBD2]);
    SetSigIndex(BD3,type_obs_[1],sig_index_[INDEXBD3]);
    SetSigIndex(GAL,type_obs_[2],sig_index_[INDEXGAL]);
    SetSigIndex(GLO,type_obs_[3],sig_index_[INDEXGLO]);
    SetSigIndex(QZS,type_obs_[4],sig_index_[INDEXQZS]);
}

int ReadRnxObs::ReadBody(Obss& obss) {
    vector<ObsData>data;
    unsigned char slips[MAXSAT][NFREQ]={{0}};
    int i,n,flag=0;

    if(!&obss||!inf_.is_open()) return 0;
    obss.ReSet();
    // read rinex obs data body
    while((n=ReadObsBody(flag,data))>=0&&!inf_.eof()) {

        for(i=0;i<n;i++){

            // utc->gpst
            if(time_sys_==TIMESYSUTC) data[i].sig_recep_.UTC2GPST();
            data[i].sig_recep_.Time2Str(1);

            // save cycle slips
            SaveSlips(slips,data[i]);
        }
        if(n>0&&!data[0].sig_recep_.Monitor(ts_,te_,tint_)){
            obss.sat_infos_.clear();
            continue;
        }

        for(i=0;i<n;i++){
            // restore cycle-slip
            ReStoreSlips(slips,data[i]);

            obss.rcv_=(unsigned char)rcv_;

            // save obs data
            obss.sat_infos_.push_back(data[i]);
        }
        data.clear();

        if(n>0){
            obss.num_=obss.sat_infos_.size();
            break;
        }
    }
    return n;
}



ReadRnxNav::ReadRnxNav() {

}

ReadRnxNav::ReadRnxNav(string file) {
    inf_.open(file,ios::in);
    ver_=3.0;
}

ReadRnxNav::~ReadRnxNav() {

}

static const double kUraEph[]={
        2.4,3.4,4.85,6.85,9.65,13.65,24.0,48.0,96.0,192.0,384.0,768.0,1536.0,
        3072.0,6144.0,0.0
};

int ReadRnxNav::UraIndex(double val) {
    int i;
    for (i=0;i<15;i++) if (kUraEph[i]>=val) break;
    return i;
}

int CmpEph(const void* p1,const void* p2) {
    Eph *q1=(Eph *)p1,*q2=(Eph *)p2;
    return q1->ttr_.time_!=q2->ttr_.time_?(int)(q1->ttr_.time_-q2->ttr_.time_):
           (q1->toe_.time_!=q2->toe_.time_?(int)(q1->toe_.time_-q2->toe_.time_):
           (q1->sat_!= q2->sat_?(int)(q1->sat_-q2->sat_):(int)(q1->code_-q2->code_)));
}

int CmpGloEph(const void *p1,const void *p2){
    GloEph *q1=(GloEph *)p1,*q2=(GloEph *)p2;
    return q1->tof_.time_!=q2->tof_.time_?(int)(q1->tof_.time_-q2->tof_.time_):
           (q1->toe_.time_!=q2->toe_.time_?(int)(q1->toe_.time_-q2->toe_.time_):
           q1->sat_-q2->sat_);
}


void ReadRnxNav::UniqeGloEph(Nav &navs) {
    if(navs.ng_<=0) return;
    qsort(&navs.glo_eph_[0],navs.ng_, sizeof(GloEph),CmpGloEph);
    LOG(DEBUG)<<"Broadcast Glonass ephemeris number: "<<navs.ng_;
}

void ReadRnxNav::UniqeEph(Nav &navs) {
    int n=navs.n_;

    if(n<=0) return;
    qsort(&navs.eph_[0], navs.n_, sizeof(Eph), CmpEph);
    LOG(DEBUG)<<"Broadcast ephemeris number: "<<n;
}

int ReadRnxNav::ReadNavBody(Nav &nav,int &eph_type) {
    Time toc,last_toc;
    SatMask sat;
    int i=0,j,sp=3,flag=1,prn=0,stat=1,last_sat,last_iode;
    string id;

    while (getline(inf_,buff_)&&!inf_.eof()) {
        if (buff_.compare(0, 3, "   ")!=0) i=0;
        /* first line */
        if (i==0) {
            /* decode satellite field */
            if (ver_>=3.0||sat_sys_==GAL||sat_sys_==QZS||sat_sys_==NONE) { /* ver.3 or GAL/QZS */
                id=buff_.substr(0,3);
                sat=SatMask(id);
                sp=4;
                if (ver_>=3.0) sat_sys_=sat.sat_sys_;
            }
            else {
                Str2Int(buff_.substr(0,2),prn);

                if (sat_sys_==GLO) sat=SatMask(GLO,prn);

                else sat=SatMask(GPS, prn);
            }
            /* decode toc field */
            if (toc.Str2Time(buff_.substr(sp,19))) flag=0;
            else flag=1;
            /* decode data fields */
            for (j=0;j<3;j++) {
                Str2Double(buff_.substr(sp+19*(j+1),19),data_[i++]);
            }
        }
        /* next line */
        else if (flag==1){
            /* decode data fields */
            for (j=0;j<4;j++) {
                if (sp+19*(j+1)<=buff_.size()){
                    Str2Double(buff_.substr(sp+19*j,19),data_[i++]);
                }
                else data_[i++]=0.0;
            }
            /* decode ephemeris */
            if (sat_sys_==GLO&&i>=15) {
                if (!(sys_mask_&sat_sys_)) continue;

                stat=DecodeGloEph(toc,sat,glo_eph0_);
                nav.glo_eph_.push_back(glo_eph0_);
                nav.ng_++;eph_type=1;
            }
            else if (i>=31) {
                if (!(sys_mask_&sat_sys_)) continue;
                stat=DecodeEph(toc,sat,eph0_);
                nav.eph_.push_back(eph0_);
                nav.n_++;eph_type=0;
            }
        }
        buff_.clear();
    }
    UniqeEph(nav);
    if(nav.ng_>0) UniqeGloEph(nav);
    return 1;
}

int ReadRnxNav::DecodeEph(Time toc, SatMask sat, Eph &eph) {
    int sys=sat_sys_;

    if (!(sys&(GPS|GAL|QZS|BD2|BD3))) {
        return 0;
    }

    eph.sat_=sat.sat_no_;
    eph.toc_=toc;

    eph.f0_=data_[0];
    eph.f1_=data_[1];
    eph.f2_=data_[2];

    eph.A_=SQR(data_[10]); eph.e_=data_[ 8]; eph.i0_  =data_[15]; eph.OMG0_=data_[13];
    eph.omg_ =data_[17]; eph.M0_ =data_[ 6]; eph.deln_=data_[ 5]; eph.OMGd_=data_[18];
    eph.idot_=data_[19]; eph.crc_=data_[16]; eph.crs_ =data_[ 4]; eph.cuc_ =data_[ 7];
    eph.cus_ =data_[ 9]; eph.cic_=data_[12]; eph.cis_ =data_[14];

    if (sys==GPS||sys==QZS) {
        eph.iode_=(int)data_[ 3];      /* IODE */
        eph.iodc_=(int)data_[26];      /* IODC */
        eph.toes_=     data_[11];      /* toe (s) in gps week */
        eph.week_=(int)data_[21];      /* gps week */
        eph.toe_.GPST2Time(eph.week_,data_[11])->AdjWeek(toc);
        eph.ttr_.GPST2Time(eph.week_,data_[27])->AdjWeek(toc);

        eph.code_=(int)data_[20];      /* GPS: codes on L2 ch */
        eph.svh_ =(int)data_[24];      /* sv health */
        eph.sva_=UraIndex(data_[23]);  /* ura (m->index) */
        eph.flag_=(int)data_[22];      /* GPS: L2 P data flag */

        eph.tgd_[0]=   data_[25];      /* TGD */
        if (sys==GPS) {
            eph.fit_=data_[28];        /* fit interval (h) */
        }
        else {
            eph.fit_=data_[28]==0.0?1.0:2.0; /* fit interval (0:1h,1:>2h) */
        }
    }
    else if (sys==GAL) { /* GAL ver.3 */
        eph.iode_=(int)data_[ 3];      /* IODnav */
        eph.toes_=     data_[11];      /* toe (s) in galileo week */
        eph.week_=(int)data_[21];      /* gal week = gps week */
        eph.toe_.GPST2Time(eph.week_,data_[11])->AdjWeek(toc);
        eph.ttr_.GPST2Time(eph.week_,data_[27])->AdjWeek(toc);

        eph.code_=(int)data_[20];      /* data sources */
                                       /* bit 0 set: I/NAV E1-B */
                                       /* bit 1 set: F/NAV E5a-I */
                                       /* bit 2 set: F/NAV E5b-I */
                                       /* bit 8 set: af0-af2 toc are for E5a.E1 */
                                       /* bit 9 set: af0-af2 toc are for E5b.E1 */
        eph.svh_ =(int)data_[24];      /* sv health */
                                       /* bit     0: E1B DVS */
                                       /* bit   1-2: E1B HS */
                                       /* bit     3: E5a DVS */
                                       /* bit   4-5: E5a HS */
                                       /* bit     6: E5b DVS */
                                       /* bit   7-8: E5b HS */
        eph.sva_ =UraIndex(data_[23]); /* ura (m->index) */

        eph.tgd_[0]=   data_[25];      /* BGD E5a/E1 */
        eph.tgd_[1]=   data_[26];      /* BGD E5b/E1 */
    }
    else if (sys==BD2||sys==BD3) { /* BeiDou v.3.02 */
        eph.toc_.CopyTime(eph.toc_)->BDT2GPST();         /* bdt -> gpst */
        eph.iode_=(int)data_[ 3];                           /* AODE */
        eph.iodc_=(int)data_[28];                           /* AODC */
        eph.toes_=     data_[11];                           /* toe (s) in bdt week */
        eph.week_=(int)data_[21];                           /* bdt week */
        eph.toe_.BDT2Time(eph.week_,data_[11])->BDT2GPST(); /* bdt -> gpst */
        eph.ttr_.BDT2Time(eph.week_,data_[27])->BDT2GPST(); /* bdt -> gpst */
        eph.toe_.AdjWeek(toc);
        eph.ttr_.AdjWeek(toc);

        eph.svh_ =(int)data_[24];      /* satH1 */
        eph.sva_=UraIndex(data_[23]);  /* ura (m->index) */

        eph.tgd_[0]=   data_[25];      /* TGD1 B1/B3 */
        eph.tgd_[1]=   data_[26];      /* TGD2 B2/B3 */
    }
//    else if (sys==SYS_IRN) { /* IRNSS v.3.03 */
//        eph->iode=(int)data[ 3];      /* IODEC */
//        eph->toes=     data[11];      /* toe (s) in irnss week */
//        eph->week=(int)data[21];      /* irnss week */
//        eph->toe.gpst2time(eph->week,data[11])->adjweek(toc);
//        eph->ttr.gpst2time(eph->week,data[27])->adjweek(toc);
//        eph->svh =(int)data[24];      /* sv health */
//        eph->sva=uraindex(data[23]);  /* ura (m->index) */
//        eph->tgd[0]=   data[25];      /* TGD */
//    }

    if (eph.iode_<0||1023<eph.iode_) {
        eph.svh_=-1;
    }
    if (eph.iodc_<0||1023<eph.iodc_) {
        eph.svh_=-1;
    }
    eph.toc_.Time2Str(1);
    eph.toe_.Time2Str(1);

    return 1;
}

int ReadRnxNav::DecodeGloEph(Time toc, SatMask sat, GloEph &glo_eph) {
    Time tof;
    double tow,tod;
    int week,dow;

    if (sat.sat_sys_!=GLO) {
        return 0;
    }

    glo_eph.sat_=sat.sat_no_;

    /* toc rounded by 15 min in utc */
    tow=toc.Time2GPST(&week);
    toc.GPST2Time(week,floor((tow+450.0)/900.0)*900);
    dow=(int)floor(tow/86400.0);

    /* time of frame in utc */
    tod=ver_<=2.99?data_[2]:fmod(data_[2],86400.0); /* tod (v.2), tow (v.3) in utc */
    tof.GPST2Time(week,tod+dow*86400.0);
    tof.AdjDay(toc);

    glo_eph.toe_.CopyTime(toc)->UTC2GPST();   /* toc (gpst) */
    glo_eph.tof_.CopyTime(toc)->UTC2GPST();   /* tof (gpst) */

    /* iode = tb (7bit), tb =index of UTC+3H within current day */
    glo_eph.iode_=(int)(fmod(tow+10800.0,86400.0)/900.0+0.5);

    glo_eph.taun_=-data_[0];       /* -taun */
    glo_eph.gamn_= data_[1];       /* +gamman */

    glo_eph.pos_[0]=data_[3]*1E3; glo_eph.pos_[1]=data_[7]*1E3; glo_eph.pos_[2]=data_[11]*1E3;
    glo_eph.vel_[0]=data_[4]*1E3; glo_eph.vel_[1]=data_[8]*1E3; glo_eph.vel_[2]=data_[12]*1E3;
    glo_eph.acc_[0]=data_[5]*1E3; glo_eph.acc_[1]=data_[9]*1E3; glo_eph.acc_[2]=data_[13]*1E3;

    glo_eph.svh_=(int)data_[ 6];
    glo_eph.frq_=(int)data_[10];
    glo_eph.age_=(int)data_[14];

    /* some receiver output >128 for minus frequency number */
    if (glo_eph.frq_>128) glo_eph.frq_-=256;

    if (glo_eph.frq_<MINFREQ_GLO||MAXFREQ_GLO<glo_eph.frq_) {
        glo_eph.svh_=-1;
    }

    return 1;
}

void ReadRnxNav::InitReadNav(string file) {
    inf_.open(file,ios::in);
}

int ReadRnxNav::ReadHead(Nav& nav) {
    /* chekc inf */
    if (!inf_.is_open()){
        defaultLogger->error("Broadcast ephemeris file open error");
        return 0;
    }

    int nline=0,i,j;
    string str;

    /* base rinex  */
    if (ReadRnxHead()!=1) return 0;

    while (getline(inf_,buff_)&&!inf_.eof()) {
        if (buff_.find("ION ALPHA",60)!=string::npos) { /* opt ver.2 */
            if (&nav) {
                for (i=0,j=2;i<4;i++,j+=12) Str2Double(buff_.substr(j,12),nav.ion_para_[0].ion[i]);
            }
        }
        else if (buff_.find("ION BETA",60)!=string::npos) { /* opt ver.2 */
            if (&nav) {
                for (i=0,j=2;i<4;i++,j+=12) Str2Double(buff_.substr(j,12),nav.ion_para_[0].ion[i+4]);
            }
        }
        else if (buff_.find("DELTA-UTC: A0,A1,T,W",60)!=string::npos) { /* opt ver.2 */
            if (&nav) {
                for (i=0,j=3;i<2;i++,j+=19) Str2Double(buff_.substr(j,12),nav.utc_para_[0].utc[i]);
                for (;i<4;i++,j+=9) Str2Double(buff_.substr(j,9),nav.utc_para_[0].utc[i]);
            }
        }
        else if (buff_.find("IONOSPHERIC CORR",60)!=string::npos) { /* opt ver.3 */
            if (&nav) {
                if (buff_.compare(0,4,"GPSA")==0) {
                    for (i=0,j=5;i<4;i++,j+=12) Str2Double(buff_.substr(j,12),nav.ion_para_[0].ion[i]);
                }
                else if (buff_.compare(0,4,"GPSB")==0) {
                    for (i=0,j=5;i<4;i++,j+=12) Str2Double(buff_.substr(j,12),nav.ion_para_[0].ion[i+4]);
                }
                else if (buff_.compare(0,3,"GAL")==0) {
                    for (i=0,j=5;i<4;i++,j+=12) Str2Double(buff_.substr(j,12),nav.ion_para_[2].ion[i]);
                }
                else if (buff_.compare(0,4,"QZSA")==0) { /* v.3.02 */
                    for (i=0,j=5;i<4;i++,j+=12) Str2Double(buff_.substr(j,12),nav.ion_para_[4].ion[i]);
                }
                else if (buff_.compare(0,4,"QZSB")==0) { /* v.3.02 */
                    for (i=0,j=5;i<4;i++,j+=12) Str2Double(buff_.substr(j,12),nav.ion_para_[i+4].ion[i]);
                }
                else if (buff_.compare(0,4,"BDSA")==0) { /* v.3.02 */
                    for (i=0,j=5;i<4;i++,j+=12) Str2Double(buff_.substr(j,12),nav.ion_para_[1].ion[i]);
                }
                else if (buff_.compare(0,4,"BDSB")==0) { /* v.3.02 */
                    for (i=0,j=5;i<4;i++,j+=12) Str2Double(buff_.substr(j,12),nav.ion_para_[1].ion[i+4]);
                }
//                else if (buff_.compare(0,4,"IRNA")==0) { /* v.3.03 */
//                    for (i=0,j=5;i<4;i++,j+=12) str2double(buff_.substr(j,12),nav.ion_irn[i]);
//                }
//                else if (buff_.compare(0,4,"IRNB")==0) { /* v.3.03 */
//                    for (i=0,j=5;i<4;i++,j+=12) str2double(buff_.substr(j,12),nav.ion_irn[i+4]);
//                }
            }
        }
        else if (buff_.find("TIME SYSTEM CORR",60)!=string::npos) { /* opt ver.3 */
            if (&nav) {
                if (buff_.compare(0,4,"GPUT")==0) {
                    Str2Double(buff_.substr( 5,17),nav.utc_para_[0].utc[0]);
                    Str2Double(buff_.substr(22,16),nav.utc_para_[0].utc[1]);
                    Str2Double(buff_.substr(38, 7),nav.utc_para_[0].utc[2]);
                    Str2Double(buff_.substr(45, 5),nav.utc_para_[0].utc[3]);
                }
                else if (buff_.compare(0,4,"GLUT")==0) {
                    Str2Double(buff_.substr( 5,17),nav.utc_para_[3].utc[0]);
                    Str2Double(buff_.substr(22,16),nav.utc_para_[3].utc[1]);
                    Str2Double(buff_.substr(38, 7),nav.utc_para_[3].utc[2]);
                    Str2Double(buff_.substr(45, 5),nav.utc_para_[3].utc[3]);
                }
                else if (buff_.compare(0,4,"GAUT")==0) { /* v.3.02 */
                    Str2Double(buff_.substr( 5,17),nav.utc_para_[2].utc[0]);
                    Str2Double(buff_.substr(22,16),nav.utc_para_[2].utc[1]);
                    Str2Double(buff_.substr(38, 7),nav.utc_para_[2].utc[2]);
                    Str2Double(buff_.substr(45, 5),nav.utc_para_[2].utc[3]);
                }
                else if (buff_.compare(0,4,"QZUT")==0) { /* v.3.02 */
                    Str2Double(buff_.substr( 5,17),nav.utc_para_[4].utc[0]);
                    Str2Double(buff_.substr(22,16),nav.utc_para_[4].utc[1]);
                    Str2Double(buff_.substr(38, 7),nav.utc_para_[4].utc[2]);
                    Str2Double(buff_.substr(45, 5),nav.utc_para_[4].utc[3]);
                }
                else if (buff_.compare(0,4,"BDUT")==0) { /* v.3.02 */
                    Str2Double(buff_.substr( 5,17),nav.utc_para_[1].utc[2]);
                    Str2Double(buff_.substr(22,16),nav.utc_para_[1].utc[2]);
                    Str2Double(buff_.substr(38, 7),nav.utc_para_[1].utc[2]);
                    Str2Double(buff_.substr(45, 5),nav.utc_para_[1].utc[3]);
                }
//                else if (buff_.compare(0,4,"SBUT")==0) { /* v.3.02 */
//                    Str2Double(buff_.substr( 5,17),nav.utc_para_[0]);
//                    Str2Double(buff_.substr(22,16),nav.utc_para_[1]);
//                    Str2Double(buff_.substr(38, 7),nav.utc_para_[2]);
//                    Str2Double(buff_.substr(45, 5),nav.utc_para_[3]);
//                }
//                else if (buff_.compare(0,4,"IRUT")==0) { /* v.3.03 */
//                    Str2Double(buff_.substr( 5,17),nav.utc_irn[0]);
//                    Str2Double(buff_.substr(22,16),nav.utc_irn[1]);
//                    Str2Double(buff_.substr(38, 7),nav.utc_irn[2]);
//                    Str2Double(buff_.substr(45, 5),nav.utc_irn[3]);
//                }
            }
        }
        else if (buff_.find("LEAP SECONDS",60)!=string::npos) { /* opt */
            if (&nav) Str2Int(buff_.substr(0,6),nav.leaps_);
        }
        if (buff_.find("END OF HEADER")!=string::npos) return 1;
        if (++nline>=MAXPOSHEAD && type_.compare(" ")==0) break; /* no rinex file */
    }
    return 0;
}

int ReadRnxNav::ReadBody(Nav &nav) {
    int stat,eph_type;
    if (!&nav) return 0;

    SetSysMark();

    ReadNavBody(nav,eph_type);
    CloseFile();
    return nav.n_>0||nav.ng_>0;
}

ReadPreEph::ReadPreEph() {
    file_=type_="";
    ns_=0;
}

ReadPreEph::~ReadPreEph() {
    buff_.clear();inf_.close();
}

void ReadPreEph::SetSysMark() {
    GNSSOpt gopt=kPrcOpt.GNSS_opt_;

    if(gopt.nav_sys_&GPS) sys_mask_|=GPS;
    if(gopt.nav_sys_&BD2) sys_mask_|=BD2;
    if(gopt.nav_sys_&BD3) sys_mask_|=BD3;
    if(gopt.nav_sys_&GAL) sys_mask_|=GAL;
    if(gopt.nav_sys_&GLO) sys_mask_|=GLO;
    if(gopt.nav_sys_&QZS) sys_mask_|=QZS;
}

bool CmpPEph(PreEph p1,PreEph p2){
    double tt=p1.time_.TimeDiff(p2.time_);
    return tt<-1E-9||tt>1E-9;
}

void ReadPreEph::UniqePreEph(vector<PreEph>& peph) {
    sort(peph.begin(),peph.end(),CmpPEph);
    peph.erase(unique(peph.begin(),peph.end(),CmpPEph));
}

int ReadPreEph::ReadHead() {
    if(!inf_.is_open()) return 0;
    int clm=0,nsat=0;
    string last_block;

    if(type_=="pre_orb"){
        while(getline(inf_,buff_)&&!inf_.eof()){
            if(!buff_.substr(0,1).compare("*")) break;
            if(clm==0){
                time_.Str2Time(buff_.substr(3,28));
            }
            else if(!buff_.substr(0,2).compare("+ ")){
                if(!last_block.compare("##")||clm==2) Str2Int(buff_.substr(3,3),ns_);
            }
            else if(!buff_.substr(0,2).compare("%c")&&!last_block.compare("++")){

            }
            else if(!buff_.substr(0,2).compare("%f")&&!last_block.compare("%c")){

            }

            last_block=buff_.substr(0,2);
            clm++;
        }
        return 1;
    }
    else if(type_=="pre_clk"){
        while(getline(inf_,buff_)&&!inf_.eof()){
            if(!buff_.substr(60,13).compare("END OF HEADER")) break;
        }
    }
}

int ReadPreEph::ReadBody(Nav &nav) {
    if(!inf_.is_open()) return 0;
    Time t;
    SatMask sat;

    if(type_=="pre_orb"){
        while(!inf_.eof()){
            if(!buff_.substr(0,3).compare("EOF")) break;
            if(buff_[0]!='*'||t.Str2Time(buff_.substr(3,28))) continue;
            nav.pre_eph_.push_back(peph0_);
            nav.pre_eph_.back().time_=t;
            for(int i=0;i<ns_&&getline(inf_,buff_);i++){
                if(buff_.length()<4||buff_[0]!='P') continue;
                sat=SatMask(buff_.substr(1,3));
                if(!(sat.sat_no_)) continue;
                if(!(sat.sat_sys_&sys_mask_)) continue;
                for(int j=0;j<4;j++){
                    double std=0.0,val=0.0;
                    Str2Double(buff_.substr(4+j*14,14),val);
                    if(buff_.length()>=80) Str2Double(buff_.substr(61+j*3,j<3?2:3),std);
                    if(buff_[0]=='P'){
                        if(val!=0.0&&fabs(val-999999.999999)>=1E-6){
                            nav.pre_eph_.back().pos_[sat.sat_no_-1][j]=val*(j<3?1000.0:1E-6);
                        }
                    }
                }
            }
            getline(inf_,buff_);
        }
        nav.pre_eph_.pop_back();
        nav.np_=nav.pre_eph_.size();
//        UniqePreEph(nav.pre_eph_);
        inf_.close();
        return nav.np_;
    }
    else if(type_=="pre_clk"){
        int k=1;
        double data[2]={0};
        while(getline(inf_,buff_)&&!inf_.eof()){
            if(buff_.substr(0,2).compare("AS")||t.Str2Time(buff_.substr(8,26))) continue;
            if(nav.pre_clk_.size()<=0 || static_cast<bool>(fabs(t.TimeDiff(nav.pre_clk_.back().time_) > 1E-9))){
                nav.pre_clk_.push_back(pclk0_);
            }
            nav.pre_clk_.back().time_=t;
            sat=SatMask(buff_.substr(3,3));
            if(!(sat.sat_sys_&sys_mask_)) continue;
            if(!(sat.sat_no_)) continue;
            if(buff_.length()>60) k=2;
            for(int i=0,j=40;i<k;i++,j+=20) Str2Double(buff_.substr(j,19),data[i]);
            nav.pre_clk_.back().clk_[sat.sat_no_-1]=data[0];
            nav.pre_clk_.back().std_[sat.sat_no_-1]=(float)data[1];
        }
        nav.nc_=nav.pre_clk_.size();
        inf_.close();
        return nav.nc_;
    }
}

int ReadPreEph::InitReadPreEph(string file) {
    file_=file;
    sys_mask_=NONE;
    if(file.empty()) return 0;
    if(!file.substr(file.length()-3,file.length()).compare("sp3")) type_="pre_orb";
    else if(!file.substr(file.length()-3,file.length()).compare("clk")) type_="pre_clk";
    inf_.open(file_,ios::in);
    return 1;
}

int ReadPreEph::ReadPre(Nav &nav) {
    SetSysMark();
    ReadHead();
    ReadBody(nav);
    return 1;
}


ReadBias::ReadBias() {

}

ReadBias::~ReadBias() {

}

int ReadBias::InitReadBias() {
    type_=kPrcOpt.GNSS_opt_.cb_prd_type_;
    ac_=kPrcOpt.GNSS_opt_.cb_prd_ac_;
    return 1;
}

int ReadBias::DecodeDCB(string file,double dcb[MAXSAT][MAXCBIASPAIR]) {
    inf_.open(file);
    if(!inf_.is_open()){
        return 0;
    }
    int i,type=0,flag=1;
    double cbias=0.0;
    char str1[5]={0},str2[10]={0};
    SatMask sat;
    string buff;
    string code_i,code_j;
    string code_pair;

    while(getline(inf_,buff)&&!inf_.eof()){
        if(kPrcOpt.GNSS_opt_.cb_prd_ac_==BIA_CODE){
            if      (!buff.find("DIFFERENTIAL (P1-P2) CODE BIASES")) type=1;
            else if (!buff.find("DIFFERENTIAL (P1-C1) CODE BIASES")) type=2;
            else if (!buff.find("DIFFERENTIAL (P2-C2) CODE BIASES")) type=3;
            if (!type||sscanf(buff.c_str(),"%s %s",str1,str2)<1) continue;
            Str2Double(buff.substr(26,9),cbias);
            if(cbias==0) continue;
            if(buff.compare(6,4,"    ")) continue;
            sat=SatMask(str1);
            if(sat.sat_no_!=0){
                if(type==1){ // P1_P2
                    code_pair="C1W-C2W";
                    if(sat.sat_sys_==GPS){
                        for(i=0;i<MAXCBIASPAIR;i++){
                            if(code_pair==kCodeInterBiasPair[INDEXGPS][i]) break;
                        }
                        dcb[sat.sat_no_-1][i]=cbias*1E-9*CLIGHT;
                    }
                    else if(sat.sat_sys_==GLO){
                        code_pair="C1P-C2P";
                        for(i=0;i<MAXCBIASPAIR;i++){
                            if(code_pair==kCodeInterBiasPair[INDEXGLO][i]) break;
                        }
                        dcb[sat.sat_no_-1][i]=cbias*1E-9*CLIGHT;
                    }
                }
                else if(type==2){ // P1_C1
                    code_pair="C1C-C1W";  //转换与CAS一致
                    if(sat.sat_sys_==GPS){
                        for(i=0;i<MAXCBIASPAIR;i++){
                            if(code_pair==kCodeInterBiasPair[INDEXGPS][i]) break;
                        }
                        dcb[sat.sat_no_-1][i]=-cbias*1E-9*CLIGHT;
                    }
                    else if(sat.sat_sys_==GLO){
                        code_pair="C1C-C1P"; //转换与CAS一致
                        for(i=0;i<MAXCBIASPAIR;i++){
                            if(code_pair==kCodeInterBiasPair[INDEXGLO][i]) break;
                        }
                        dcb[sat.sat_no_-1][i]=-cbias*1E-9*CLIGHT;
                    }
                }
                else if(type==3){ // P2_C2
                    code_pair="C2W-C2L";
                    if(sat.sat_sys_==GPS){
                        for(i=0;i<MAXCBIASPAIR;i++){
                            if(code_pair==kCodeInterBiasPair[INDEXGPS][i]) break;
                        }
                        dcb[sat.sat_no_-1][i]=cbias*1E-9*CLIGHT;
                    }
                    else if(sat.sat_sys_==GLO){
                        code_pair="C2C-C2P";
                        for(i=0;i<MAXCBIASPAIR;i++){
                            if(code_pair==kCodeInterBiasPair[INDEXGLO][i]) break;
                        }
                        dcb[sat.sat_no_-1][i]=cbias*1E-9*CLIGHT;
                    }
                }
            }
        }
        else if(kPrcOpt.GNSS_opt_.cb_prd_ac_==BIA_CAS){
            if((!buff.compare(1,3,"DSB"))&&(!buff.compare(65,2,"ns"))&&
               (!buff.compare(15,4,"    "))){
                sat=SatMask(buff.substr(11,3));
                code_pair=buff.substr(25,3)+"-"+buff.substr(30,3);
                Str2Double(buff.substr(80,10),cbias);
                if(sat.sat_sys_==GPS){
                    for(i=0;i<MAXCBIASPAIR;i++){
                        if(code_pair==kCodeInterBiasPair[INDEXGPS][i]) break;
                    }

                    dcb[sat.sat_no_-1][i]=cbias*1E-9*CLIGHT;
                }
                else if(sat.sat_sys_==BD2){
                    for(i=0;i<MAXCBIASPAIR;i++){
                        if(code_pair==kCodeInterBiasPair[INDEXBD2][i]) break;
                    }
                    dcb[sat.sat_no_-1][i]=cbias*1E-9*CLIGHT;
                }
                else if(sat.sat_sys_==BD3){
                    for(i=0;i<MAXCBIASPAIR;i++){
                        if(code_pair==kCodeInterBiasPair[INDEXBD3][i]) break;
                    }
                    dcb[sat.sat_no_-1][i]=cbias*1E-9*CLIGHT;
                }
                else if(sat.sat_sys_==GAL){
                    for(i=0;i<MAXCBIASPAIR;i++){
                        if(code_pair==kCodeInterBiasPair[INDEXGAL][i]) break;
                    }
                    dcb[sat.sat_no_-1][i]=cbias*1E-9*CLIGHT;
                }
                else if(sat.sat_sys_==GLO){
                    for(i=0;i<MAXCBIASPAIR;i++){
                        if(code_pair==kCodeInterBiasPair[INDEXGLO][i]) break;
                    }
                    dcb[sat.sat_no_-1][i]=cbias*1E-9*CLIGHT;
                }
                else if(sat.sat_sys_==QZS){
                    for(i=0;i<MAXCBIASPAIR;i++){
                        if(code_pair==kCodeInterBiasPair[INDEXQZS][i]) break;
                    }
                    dcb[sat.sat_no_-1][i]=cbias*1E-9*CLIGHT;
                }
            }
            else continue;
        }
    }
    inf_.close();
}

int ReadBias::AddOsb(const Osb_t* osb) {
    Osb_t *osb_tmp= nullptr;
    if(sat_osbs_.nb>=sat_osbs_.nb_max){
        sat_osbs_.nb_max+=2048;
        if(!(osb_tmp=(Osb_t *)realloc(sat_osbs_.osbs, sizeof(Osb_t)*sat_osbs_.nb_max))){
            return 0;
        }
        sat_osbs_.osbs=osb_tmp;
    }
    sat_osbs_.osbs[sat_osbs_.nb++]=*osb;
    return 1;
}


int ReadBias::OsbTimeStr2Time(const char *s, int i, int n, Time& t) {
    double ep[6]={0};ep[1]=1.0;ep[2]=1.0;
    double day=0.0,sec=0.0;
    char str[256]={'\0'},*p=str;

    if(i<0||(int)strlen(s)<i) return -1;
    for(s+=i;*s&&--n>=0;) *p++=*s++;*p='\0';
    if(sscanf(str,"%lf:%lf:%lf",ep,&day,&sec)<3) return -1;
    if(ep[0]<100.0) ep[0]+=ep[0]<80.0?2000.0:1900.0;
    t.Epoch2Time(ep);t.TimeAdd(86400.0*(day-1.0)+sec);
    return 0;
}

int ReadBias::ReadCODDCB(double dcb[MAXSAT][MAXCBIASPAIR]) {
    int i,j,n,mon=kPrcOpt.prc_date_.epoch_[1],m;
    char *ex_files[36]={nullptr};
    string file_name,sep;
#if WIN32
    sep="\\.";
#else
    sep="/.";
#endif

    for(i=0;i<36;i++){
        if(!(ex_files[i]=(char *)malloc(1024))){
            for(i--;i>=0;i--) free(ex_files[i]);
            return 0;
        }
    }
    n=ExPath((char *)kFileOpt.bia_.c_str(),ex_files,36);

    for(i=0;i<n;i++){
        vector<string> splits=MultiSplitStr(ex_files[i],sep);
        m=atoi(((splits.end()-2)->substr(6,2)).c_str());
        if(mon!=m) continue;
        if(!DecodeDCB(ex_files[i],dcb)) return 0;
    }

    for(i=0;i<36;i++) free(ex_files[i]);
    return 1;
}

int ReadBias::ReadCASDCB(string file,double dcb[MAXSAT][MAXCBIASPAIR]) {
    DecodeDCB(file,dcb);
    return 1;
}

int ReadBias::ReadOSB(string file) {
    inf_.open(file);
    if(!inf_.is_open()){
        return 0;
    }

    Osb_t osb;
    SatMask sat;
    int i=0;
    string buff;

    while(getline(inf_,buff)&&!inf_.eof()){
        if((!buff.compare(1,3,"OSB"))&&(!buff.compare(65,2,"ns"))&&
          (!buff.compare(15,4,"    "))){
//          //读取第一个数据
            sat=SatMask(buff.substr(11,3));
            osb.sat=sat.sat_no_;
            osb.type=((buff.compare(25,1,"L"))?CODE_BIAS:PHASE_BIAS);
            for(i=0;kSignalCodes;i++){
                if(!buff.compare(26,2,kSignalCodes[i])){
                    osb.code=i;
                    break;
                }
            }
            if(OsbTimeStr2Time(buff.c_str(),35,14,osb.ts)) continue;
            if(OsbTimeStr2Time(buff.c_str(),50,14,osb.te)) continue;
            Str2Double(buff.substr(70,21),osb.val);
            if(!AddOsb( &osb)) return 0;
        }
    }
    inf_.close();

    return 1;
}

int ReadBias::ReadBiasFile(double dcb[MAXSAT][MAXCBIASPAIR]) {
    if(type_==CBC_DCB){
        if(ac_==BIA_CODE){
            ReadCODDCB(dcb);
        }
        else if(ac_==BIA_CAS){
            ReadCASDCB(kFileOpt.bia_,dcb);
        }
    }
    else if(type_==CBC_OSB){
        ReadOSB(kFileOpt.bia_);
    }
    return 1;
}

ReadGIM::ReadGIM() {
    hgts_[0]=hgts_[1]=450;hgts_[2]=0.0;
    lats_[0]=87.5;lats_[1]=-87.5;lats_[2]=-2.5;
    lons_[0]=-180.0;lons_[1]=180.0;lons_[2]=5.0;
    re_=6371.0;factor_=-1.0;
}

ReadGIM::~ReadGIM() = default;

int ReadGIM::DataIndex(int i, int j, int k, const int *ndata) {
    if (i<0||ndata[0]<=i||j<0||ndata[1]<=j||k<0||ndata[2]<=k) return -1;
    return i+ndata[0]*(j+ndata[1]*k);
}

int ReadGIM::GetIndex(double val, const double *range) {
    if(range[2]==0.0) return 0;
    if(range[2]>0.0&&(val<range[0]||range[1]<val)) return -1;
    if(range[2]<0.0&&(val<range[1]||range[0]<val)) return -1;
    return (int)floor((val-range[0])/range[2]+0.5);
}

int ReadGIM::GetNumItems(const double *range) {
    return GetIndex(range[1],range)+1;
}

Tec* ReadGIM::AddTec(vector<Tec>& tecs) {
    int num_data[3];

    num_data[0]=GetNumItems(lats_);
    num_data[1]=GetNumItems(lons_);
    num_data[2]=GetNumItems(hgts_);
    if(num_data[0]<=1||num_data[1]<=1||num_data[2]<=0) return NULL;

    tecs.push_back(Tec());
    tecs.back().re_=re_;
    for(int i=0;i<3;i++){
        tecs.back().ndata_[i]=num_data[i];
        tecs.back().lats_[i]=lats_[i];
        tecs.back().lons_[i]=lons_[i];
        tecs.back().hgts_[i]=hgts_[i];
    }
    int n=num_data[0]*num_data[1]*num_data[2];

    tecs.back().data_.assign(n,0.0);
    tecs.back().rms_.assign(n,0.0);

    for(int i=0;i<n;i++){
        tecs.back().data_[i]=0.0;
        tecs.back().rms_[i]=0.0f;
    }
    return &tecs.back();
}

int ReadGIM::ReadHead() {
    double ver=0;
    if(!inf_.is_open()) return 0;
    while(getline(inf_,buff_)&&!inf_.eof()){
        if (buff_.length()<60) continue;

        if (buff_.find("IONEX VERSION / TYPE")!=string::npos){
            if (buff_[20]=='I') Str2Double(buff_.substr(0,8),ver);
        }
        else if (buff_.find("BASE RADIUS")!=string::npos) {
            Str2Double(buff_.substr(0,8),re_);
        }
        else if (buff_.find("HGT1 / HGT2 / DHGT")!=string::npos){
            Str2Double(buff_.substr(2, 6),hgts_[0]);
            Str2Double(buff_.substr(8, 6),hgts_[1]);
            Str2Double(buff_.substr(14,6),hgts_[2]);
        }
        else if (buff_.find("LAT1 / LAT2 / DLAT")!=string::npos){
            Str2Double(buff_.substr(2, 6),lats_[0]);
            Str2Double(buff_.substr(8, 6),lats_[1]);
            Str2Double(buff_.substr(14,6),lats_[2]);
        }
        else if (buff_.find("LON1 / LON2 / DLON")!=string::npos) {
            Str2Double(buff_.substr(2, 6),lons_[0]);
            Str2Double(buff_.substr(8, 6),lons_[1]);
            Str2Double(buff_.substr(14,6),lons_[2]);
        }
        else if (buff_.find("EXPONENT")!=string::npos){
            Str2Double(buff_.substr(0,6),factor_);
        }
        else if (buff_.find("START OF AUX DATA")!=string::npos&&
                 buff_.find("DIFFERENTIAL CODE BIASES")!=string::npos){
            continue;
        }
        else if (buff_.find("PRN / BIAS / RMS")!=string::npos){
            continue;
        }
        else if (buff_.find("END OF HEADER")!=string::npos){
            return ver;
        }
    }
}

int ReadGIM::ReadBody(vector<Tec>& tecs) {
    if(!inf_.is_open()) return 0;

    Tec *tec=nullptr;
    int flag=0;

    while(getline(inf_,buff_)&&!inf_.eof()){
        if(buff_.length()<60) continue;
        if(buff_.find("START OF TEC MAP")!=string::npos){
            tec = AddTec(tecs);
            if(tec) flag=1;
        }
        else if(buff_.find("END OF TEC MAP")!=string::npos){
            flag=0;tec=nullptr;
        }
        else if(buff_.find("START OF RMS MAP")!=string::npos){
            flag=2;tec=nullptr;
        }
        else if(buff_.find("END OF RMS MAP")!=string::npos){
            flag=0;tec=nullptr;
        }
        else if(buff_.find("EPOCH OF CURRENT MAP")!=string::npos){
            if(ion_time_.Str2Time(buff_.substr(0,36))) continue;
            if(flag==2){
                for(int i=tecs.size()-1;i>=0;i--){
                    if(fabs(ion_time_.TimeDiff(tecs[i].t_))>=1.0) continue;
                    tec=&tecs[i];
                    break;
                }
            }
            else if(tec) tec->t_=ion_time_;
        }
        else if(buff_.find("LAT/LON1/LON2/DLON/H")!=string::npos&&tec){
            double lat,lon[3],hgt,x;
            Str2Double(buff_.substr(2,6), lat);
            Str2Double(buff_.substr(8,6), lon[0]);
            Str2Double(buff_.substr(14,6),lon[1]);
            Str2Double(buff_.substr(20,6),lon[2]);
            Str2Double(buff_.substr(26,6),hgt);
            int i=GetIndex(lat,tec->lats_);
            int k=GetIndex(hgt,tec->hgts_);
            int n=GetNumItems(lon);

            for (int m=0;m<n;m++) {
                int index;
                if (m%16==0&&!getline(inf_,buff_)) break;
                int j=GetIndex(lon[0]+lon[2]*m,tec->lons_);
                if ((index=DataIndex(i,j,k,tec->ndata_))<0) continue;
                Str2Double(buff_.substr(m%16*5,5),x);
                if (x==9999.0) continue;
                if (flag==1) tec->data_[index]=x*pow(10.0,factor_);
                else tec->rms_[index]=(float)(x*pow(10.0,factor_));
            }
        }
    }
}

void ReadGIM::InitReadGIM(string file) {
    hgts_[0]=hgts_[1]=450;hgts_[2]=0.0;
    lats_[0]=87.5;lats_[1]=-87.5;lats_[2]=-2.5;
    lons_[0]=-180.0;lons_[1]=180.0;lons_[2]=5.0;
    re_=6371.0;factor_=-1.0;
    ion_file_=file;
    inf_.open(ion_file_,ios::in);
}

int ReadGIM::ReadTec(vector<Tec> &tecs) {
    if(!ReadHead()) return 0;
    if(!ReadBody(tecs)) return 0;
    inf_.close();
    return 1;
}

ReadErp::ReadErp() = default;

ReadErp::~ReadErp() = default;

void ReadErp::InitReadErp(string file) {
    inf_.open(file,ios::in);
}

int ReadErp::ReadErpPara(vector<Erpd_t> &erps) {
    double v[14]={0};
    Erpd_t erp0={0};

    if(!inf_.is_open()) return 0;
    erps.clear();

    while(getline(inf_,buff_)&&!inf_.eof()){
        if (sscanf(buff_.c_str(),"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                   v,v+1,v+2,v+3,v+4,v+5,v+6,v+7,v+8,v+9,v+10,v+11,v+12,v+13)<5)
            continue;
        erp0.mjd=v[0];
        erp0.xp=v[1]*1E-6*AS2R;
        erp0.yp=v[2]*1E-6*AS2R;
        erp0.ut1_utc=v[3]*1E-7;
        erp0.lod=v[4]*1E-7;
        erp0.xpr=v[12]*1E-6*AS2R;
        erp0.ypr=v[13]*1E-6*AS2R;
        erps.push_back(erp0);
    }
    inf_.close();
    return 1;
}

ReadBlq::ReadBlq() = default;

ReadBlq::~ReadBlq() = default;

void ReadBlq::InitReadBlq(const string& file) {
    inf_.open(file,ios::in);
}

int ReadBlq::ReadBlqPara(string site, double *ocean_par) {
    string name="    ";
    if(!inf_.is_open()) return 0;

    transform(site.begin(),site.begin()+4,name.begin(),::toupper);

    while(inf_.is_open()&&getline(inf_,buff_)&&!inf_.eof()){
        if(buff_.size()<2||buff_.compare(0,2,"$$")==0) continue;
        if(buff_.compare(2,4,name)==0){
            double v[11]={0};
            int n=0;
            while(getline(inf_,buff_)&&!inf_.eof()){
                if(buff_.compare(0,2,"$$")==0) continue;
                if (sscanf(buff_.c_str(),"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                           v,v+1,v+2,v+3,v+4,v+5,v+6,v+7,v+8,v+9,v+10)<11) continue;
                for(int i=0;i<11;i++) ocean_par[n+i*6]=v[i];
                if(++n==6){inf_.close();return 1;}
            }
            inf_.close();return 0;
        }
    }
    inf_.close();
    return 1;
}


ReadAnt::ReadAnt() {

}

ReadAnt::~ReadAnt() {

}

int ReadAnt::DecodeAntPar(char *p, int n, double *v) {
    int i;
    for(i=0;i<n;i++) v[i]=0.0;
    for(i=0,p=strtok(p," ");p&&i<n;p=strtok(NULL," ")){
        v[i++]=atof(p)*1E-3;
    }
    return i;
}

void ReadAnt::InitReadAtx(const string &file) {
    inf_.open(file,ios::in);
}

int ReadAnt::ReadAtx(vector<Ant_t>& ants) {
    int stat=0,frq=0,f,i;
    Ant_t ant0={0};
    string sys;

    while(inf_.is_open()&&getline(inf_,buff_)&&!inf_.eof()){
        if(buff_.length()<60||buff_.find("COMMENT",60)!=string::npos) continue;
        if(buff_.find("START OF ANTENNA",60)!=string::npos){
            ants.push_back(ant0);stat=1;
        }
        if(buff_.find("END OF ANTENNA",60)!=string::npos) stat=0;
        if(!stat) continue;

        if(buff_.find("TYPE / SERIAL NO",60)!=string::npos){
            if(!buff_.find("JAVRINGANT_G5T")){
                int a=1;
            }
            ants.back().ant_type=buff_.substr(0,20);
            ants.back().ser_code=buff_.substr(20,20);
            if(ants.back().ser_code.compare(3,8,"        ")==0){
                ants.back().sat_=SatMask(ants.back().ser_code.substr(0,3));
            }
        }
        else if(buff_.find("VALID FROM",60)!=string::npos){
            if(!ants.back().ts.Str2Time(buff_.substr(0,43))) continue;
        }
        else if(buff_.find("VALID UNTIL",60)!=string::npos){
            if(!ants.back().te.Str2Time(buff_.substr(0,43))) continue;
        }
        else if(buff_.find("DAZI",60)!=string::npos){
            Str2Double(buff_.substr(2,6),ants.back().dazi); continue;
        }
        else if(buff_.find("ZEN1 / ZEN2 / DZEN")!=string::npos) {
            Str2Double(buff_.substr(2,6),ants.back().zen1);
            Str2Double(buff_.substr(8,6),ants.back().zen2);
            Str2Double(buff_.substr(14,6),ants.back().dzen);
            continue;
        }
        else if(buff_.find("START OF FREQUENCY",60)!=string::npos){
            if(sscanf(buff_.c_str()+4,"%d",&f)<1) continue;
            sys=buff_[3];
            if(sys=="G") frq=f;
            else if(sys=="C"){
                if(f==1) frq=f+1*NFREQ;
                else if(f==7) frq=2+1*NFREQ;
                else if(f==6) frq=3+1*NFREQ;
//                frq=f;
                //TODO 不同版本的atx文件对BDS定义的频率不一样
            }
            else if(sys=="E"){
                if(f==1) frq=f+2*NFREQ;
                else if(f==5) frq=2+2*NFREQ;
                else if(f==6) frq=3+2*NFREQ;
            }
            else if(sys=="R"){
                frq=f+3*NFREQ;
            }
            else if(sys=="J"){
                if(f==1) frq=f+4*NFREQ;
                else if(f==2) frq=2+4*NFREQ;
                else if(f==5) frq=3+4*NFREQ;
            }
            else frq=0;
        }
        else if(buff_.find("END OF FREQUENCY",60)!=string::npos){
            frq=0;
        }
        else if(buff_.find("NORTH / EAST / UP",60)!=string::npos){
            double neu[3]={0};
            if(frq<1) continue;
            if(sscanf(buff_.c_str(),"%lf %lf %lf",neu,neu+1,neu+2)<3) continue;
            ants.back().pco_[frq-1][0]=1E-3*neu[ants.back().sat_.sat_no_?0:1];
            ants.back().pco_[frq-1][1]=1E-3*neu[ants.back().sat_.sat_no_?1:0];
            ants.back().pco_[frq-1][2]=1E-3*neu[2];
        }
        else if(buff_.find("NOAZI")!=string::npos){
            if(frq<1) continue;
            double dd=(ants.back().zen2-ants.back().zen1)/ants.back().dzen+1;
            if(dd!=Round(dd)||dd<=1){
                LOG_N_TIMES(1,WARNING)<<"Number of PCV NOAZI parameter decode error";
                continue;
            }
            if(ants.back().dazi==0.0){
                i=DecodeAntPar(const_cast<char*>(buff_.c_str()+8),(int)dd,ants.back().pcv_[frq-1]);
                if(i<=0) {
                    LOG_N_TIMES(1,WARNING)<<"Number of PCV NOAZI parameter decode error";
                    continue;
                }
                else if(i!=(int)dd){
                    LOG_N_TIMES(1,WARNING)<<"Number of PCV NOAZI parameter decode error";
                    continue;
                }
            }
            else{
                int id=(int)(360-0)/ants.back().dazi+1;
                for(i=0;i<id;i++){
                    getline(inf_,buff_);
                    int j=DecodeAntPar(const_cast<char*>(buff_.c_str()+8),(int)dd,&ants.back().pcv_[frq-1][i*(int)dd]);
                    if(j<=0){
                        LOG_N_TIMES(1,WARNING)<<"Number of PCV AZI parameter decode error";
                        continue;
                    }
                    else if(j!=(int)dd){
                        LOG_N_TIMES(1,WARNING)<<"Number of PCV AZI parameter decode error";
                        continue;
                    }
                }
            }
        }
    }
    inf_.close();
    return 1;
}
