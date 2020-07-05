/**************************************************************************

Copyright(c): 2020, by Chao Chen, All rights reserved.
              China University of Mining and Technology, XuZhou, P.R. China

Author: Chao Chen(cchen@cumt.edu.cn)

Date: 2020-05-30

Description:

**************************************************************************/

#include "GNSS.h"
#include "ParSetting.h"
#include "OutSol.h"
SolStat::SolStat() {
    dt_=0.0;
    i_epoch_=0;
    stat_=SOL_NONE;
    ppp_flag_=0;
    clk_G_=0.0;
    B2_isb_=B3_isb_=E_isb_=R_isb_=J_isb_=0.0;
}

SolStat::~SolStat() {

}


OutSol::OutSol() {
    isat_= nullptr;
}

OutSol::OutSol(string file) {
    fout_.open(file,ios::out);
    sopt_=kSolOpt;
}

OutSol::~OutSol() {
    fout_.close();
}

void OutSol::OutSat(Obss obs) {
    int sat=0;
    ObsData data;
    Isat_t isat;
    Parameter para;

    for(int i=0;i<obs.num_;i++){
        data=obs.sat_infos_[i];
        sat=data.sat_.sat_no_;
        isat=isat_[sat-1];
        if(data.stat_>USED&&data.stat_!=EX_SLIP){
            isat.post_res[OBS_PR][0]=isat.post_res[OBS_PR][1]=isat.post_res[OBS_PR][2]=0.0;
            isat.post_res[OBS_CP][0]=isat.post_res[OBS_CP][1]=isat.post_res[OBS_CP][2]=0.0;
        }
        char info[MAXBUFF]={'\0'};
        char res[MAXBUFF]={'\0'};
        char amb[MAXBUFF]={'\0'};
        char trp[MAXBUFF]={'\0'};
        char ion[MAXBUFF]={'\0'};
        char loc[MAXBUFF]={'\0'};
        // info: id, stat, az, el;
        sprintf(info,"$SAT     %4s %10s %6.1f %6.1f",data.sat_.sat_id_.c_str(),kListSatEXOptStr[data.stat_].c_str(),data.azel_[1]*R2D,data.azel_[0]*R2D);

        // res: P1(IF_P),P2,P3(m), L1(IF_L),L2,L3(mm)
        sprintf(res,"%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f",isat.post_res[OBS_PR][0],isat.post_res[OBS_PR][1],
                isat.post_res[OBS_PR][2],isat.post_res[OBS_CP][0],isat.post_res[OBS_CP][1],isat.post_res[OBS_CP][2]);

        // amb: L1(IF_L), L2, L3, MW, SMW
        sprintf(amb,"%8.3f %8.3f %8.3f %8.3f %8.3f",isat.amb[0],isat.amb[1],isat.amb[2],isat.mw[0],isat.smw[0]);

        // trp: strp_h,strp_w,map_h,map_w
        sprintf(trp,"%6.3f %6.3f %6.3f %6.3f",isat.strp_h,isat.strp_w,isat.map_trp_h,isat.map_trp_w);

        // ion: sion
        sprintf(ion,"%6.3f",isat.sion);

        // lock and outc
        sprintf(loc,"%5d %5d",isat.lock[0],isat.outc[0]);

        fout_<<info<<" "<<res<<" "<<amb<<" "<<trp<<" "<<ion<<" "<<loc<<endl;
    }
}

void OutSol::InitOutSol(std::__cxx11::string file) {
   fout_.open(file,ios::out);
}

int OutSol::TestOpen() {
    return fout_.is_open();
}

void OutSol::WriteSolHead(Sta sta) {
    string comm="%";
    string use_sys,clk;

    if(kPrcOpt.GNSS_opt_.nav_sys_&GPS){
        use_sys+="GPS ";
    }
    if(kPrcOpt.adj_opt_.bds_isb_){
        if(kPrcOpt.GNSS_opt_.nav_sys_&BD2){
            use_sys+="BD2 ";
        }
        if(kPrcOpt.GNSS_opt_.nav_sys_&BD3){
            use_sys+="BD3 ";
        }
    }
    else{
        if(kPrcOpt.GNSS_opt_.nav_sys_&BD2||kPrcOpt.GNSS_opt_.nav_sys_&BD3){
            use_sys+="BDS ";
        }
    }
    if(kPrcOpt.GNSS_opt_.nav_sys_&GAL){
        use_sys+="GAL ";
    }
    if(kPrcOpt.GNSS_opt_.nav_sys_&GLO){
        use_sys+="GLO ";
    }
    if(kPrcOpt.GNSS_opt_.nav_sys_&QZS){
        use_sys+="QZS ";
    }

    fout_<<"+ PPPLib HEADER"<<endl;
    fout_<<"  Station : "<<sta.name_<<endl;
    fout_<<"  Rec Type: "<<sta.rec_type_<<endl;
    fout_<<"  Ant Type: "<<sta.ant_desc_<<endl;
    if(kPrcOpt.mode_opt_!=MDOPT_KINE){
        fout_<<"  Sta  Pos: "<<std::fixed<<setprecision(3)<<sta.pos_[0]<<"  "<<sta.pos_[1]<<"  "<<sta.pos_[2]<<endl;
    }
    fout_<<"  Nav  Sys: "<<use_sys<<endl;
    fout_<<"  Trp Mode: "<<kListTropStr[kPrcOpt.GNSS_opt_.trp_opt_]<<endl;
    fout_<<"  Ion Mode: "<<kListIonsStr[kPrcOpt.GNSS_opt_.ion_opt_]<<endl;
    fout_<<"  Sat  Eph: "<<kListEphStr[kPrcOpt.GNSS_opt_.sat_eph_]<<endl;
    fout_<<"  $EPOCH  : "<<"Y/M/D H:M:S | WEEK | WOS | Epoch | SolStat | ObsNum | ValidSat | PDOP | SIGMA0 "<<endl;
    fout_<<"  $POS    : "<<kListPosFmtStr[sopt_.pos_fmt_]<<" (m)"<<endl;
    fout_<<"  $CLK    : "<<"Clock | ISB"<<endl;
    fout_<<"  $TRP    : "<<"Hydrostatic | Wet "<<endl;
    fout_<<"  $SAT    : "<<"ID | Stat | Obs | Az | El | Res_P1 | Res_P2 | Res_P3 | Res_L1 | Res_L2 | Res_L3 | Amb1 | Amb2 | Amb3 | MwAmb | sMwAmb | sTrp_h | sTrp_w | map_h | map_w | sIon | Lock | Outc";
    if(kPrcOpt.GNSS_opt_.ion_opt_!=ION_IF) fout_<<"| Ion "<<endl;
    else fout_<<endl;
    fout_<<"- PPPLib HEADER"<<endl<<endl;
}

void OutSol::WriteSols(SolStat sols,Obss obs) {
    int week,nf=1;
    double wos;
    ObsData data;
    Vec3 base_blh(0.0),out_pos(0.0),enu[3],dr(0.0);
    int dgnss=(kPrcOpt.mode_==DGPS||kPrcOpt.mode_==PPK)||
              (kPrcOpt.mode_opt_==MDOPT_DGPS||kPrcOpt.mode_opt_==MDOPT_PPK)||kPrcOpt.mode_==FIX;

    if(kSolOpt.pos_fmt_==POSFMT_ENU){
        if(dgnss){
            base_blh=Xyz2Blh(kPrcOpt.GNSS_opt_.pos_opt_.base_,WGS84);
            dr=sols.pos_-kPrcOpt.GNSS_opt_.pos_opt_.base_;
        }
        else{
            base_blh=Xyz2Blh(kSolOpt.ref_sol_.pos,WGS84);
            dr=sols.pos_-kSolOpt.ref_sol_.pos;
        }
        out_pos=Xyz2Enu(dr,base_blh);
    }
    else{
        if(kSolOpt.err_fmt_) out_pos=sols.pos_-kSolOpt.ref_sol_.pos;
        else out_pos=sols.pos_;
    }

    wos=sols.sol_time_.Time2GPST(&week);
    char epoch[MAXBUFF]={'\0'};
    sprintf(epoch,  "$EPOCH    %12s %4d %10.2f %5d %6s %3d %3d %5.1f %5.1f",sols.sol_time_.time_str_.c_str(),week,wos,
            sols.i_epoch_,kListSolStatStr[sols.stat_].c_str(),obs.num_,sols.num_used_sat_,sols.pdop_,sols.sigma0_);
    fout_<<epoch<<endl;
    if(sols.stat_>SOL_NONE){
        char pos[MAXBUFF]={'\0'};
        sprintf(pos,"$POS   %12.3f %12.3f %12.3f\n",out_pos.i_,out_pos.j_,out_pos.k_);
        char clk[MAXBUFF]={'\0'};
        sprintf(clk,"$CLK   %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f\n",sols.clk_G_,sols.B2_isb_,sols.B3_isb_,sols.E_isb_,sols.R_isb_,sols.J_isb_);
        char trp[MAXBUFF]={'\0'};
        sprintf(trp,"$TRP   %12.3f %12.3f\n",sols.ztrp_dry_,sols.ztrp_wet_);
        fout_<<pos;
        fout_<<clk;
        fout_<<trp;
        if(kSolOpt.out_sat_) OutSat(obs);
        fout_<<endl;
    }
    else{
        //无解也可以输出sat信息
        fout_.setf(std::ios::right);
        if(kSolOpt.out_sat_) OutSat(obs);
        fout_<<endl;
    }
}
