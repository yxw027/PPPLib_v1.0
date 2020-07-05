//
// Created by cc on 4/13/20.
//

#ifndef MSNAVS_OUTSOL_H
#define MSNAVS_OUTSOL_H
#include "PPPLibGlo.h"


class SolStat {
public:
    SolStat();
    ~SolStat();
public:
    double dt_;
    Time sol_time_;
    int i_epoch_;
    int stat_;
    int sol_fmt_;
    int ppp_flag_;
    int num_used_sat_;
    double pdop_;
    double sigma0_;

    Vec3 pos_,vel_;
    Vec3 q_pos_,q_vel_;
    double clk_G_,B2_isb_,B3_isb_,E_isb_,R_isb_,J_isb_;
    double ztrp_dry_,ztrp_wet_;
};

class OutSol {
public:
    OutSol();
    OutSol(string file);
    ~OutSol();

protected:
    void OutSat(Obss obs);

public:
    void InitOutSol(string file);
    int TestOpen();
    void WriteSolHead(Sta sta);
    void WriteSols(SolStat sols,Obss obs);

public:
    Isat_t *isat_;
    ofstream fout_;
    string buff_;
    SolOpt sopt_;
};


#endif //MSNAVS_OUTSOL_H
