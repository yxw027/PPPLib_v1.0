//
// Created by cc on 4/9/20.
//

#ifndef MSNAVS_EPHMODEL_H
#define MSNAVS_EPHMODEL_H

#include "GNSS.h"
#include "AntModel.h"

class SatEph {
public:
    SatEph();
    virtual  ~SatEph();

protected:
    virtual int SatClk(ObsData& data);

public:
    void InitSatEph(const Nav& nav);
    virtual int SatPos(ObsData& data,int iode);
    virtual int BrdcPosClk(Obss& obs);
    virtual int PrecPosClk(Obss& obs);
    int SatPosClk(Obss& obs);

public:
    AntModel *ant_corr_;

public:
    int ne_,ng_;
    vector<Eph> eph_;
    vector<GloEph> glo_eph_;
    int np_,nc_;
    vector<PreEph> peph_;
    vector<PreClk> pclk_;
};

class BrdcEph: public SatEph {
public:
    BrdcEph();
    virtual ~BrdcEph();

protected:
    double UraEph(int ura);
    int SelectEph(ObsData& data,int iode);
    int SelGloEph(ObsData& data,int iode);
    void Eph2Clk(ObsData& data,const Eph eph);
    void GloEph2Clk(ObsData& data,const GloEph glo_eph);
    void Eph2Pos(ObsData &data, int iode, const Eph eph);
    void GloDefEq(const double *x,double *xdot,const double *acc);
    void GloOrbit(double t,double *x,const double *acc);
    void GloEph2Pos(ObsData& data,const GloEph glo_eph);
    int SatBrdcClk(ObsData& data);
    int SatBrdcPos(ObsData& data,int iode);

public:
    int BrdcPosClk(Obss& obs);
};

class PrecEph: public BrdcEph {
public:
    PrecEph();
    ~PrecEph();

protected:
    double InterpoLagr(const double* dt,Time* ptime,const double* ppos,int n);
    double InterpolNev(const double* x,double *y,int n);
    int PreEph2Pos(ObsData& data);
    int SatPrecClk(ObsData& data);
    int SatPrecPos(ObsData& data,int iode);
    
public:
    int PrecPosClk(Obss& obs);


protected:
    double orb_var_;
    double clk_var_;
};

#endif //MSNAVS_EPHMODEL_H
