//
// Created by cc on 4/9/20.
//

#ifndef MSNAVS_ANTMODEL_H
#define MSNAVS_ANTMODEL_H

#include "GNSS.h"

typedef struct {
    SatMask sat_;
    string ant_type;
    string ser_code;
    Time ts,te;
//    double pco[NSYS*NFREQ][3];
//    double pcv[NSYS*NFREQ][80*50];
    double pco_[NFREQ*(NSYS-1)][3],pcv_[NFREQ*(NSYS-1)][80*30];
    double dazi;
    double zen1,zen2,dzen;
}Ant_t;

class AntModel {
public:
    AntModel();
    ~AntModel();

protected:
    Ant_t* SearchAntPara(int sat,const string& type);
    double InterpPcv(double ang,const double *var);
    double InterpAziPcv(const Ant_t& ant,double az,double ze,int f);
    void SatPcvModel(ObsData& data,double nadir,double *dant);

public:
    void SetAntPara();
    void SatPcoCorr(ObsData& data,double *rs,double *dants);
    void SatPcvCorr(ObsData& data,double *rr,double *pcv_dants);
    void RecAntCorr(int rov_bas,ObsData& data,double *dantr);

public:
    vector<Ant_t>ants_;
    Time ant_time_;
    Sta sta_[2];
    Ant_t sat_ant_[MAXSAT];
    Ant_t rec_ant_[2];
    double ant_del[2][3];
};


#endif //MSNAVS_ANTMODEL_H
