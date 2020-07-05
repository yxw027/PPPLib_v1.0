//
// Created by cc on 4/9/20.
//

#ifndef MSNAVS_READFILES_H
#define MSNAVS_READFILES_H

#include "GNSS.h"
#include "BiaModel.h"
#include "IonModel.h"
#include "TidModel.h"
#include "AntModel.h"

class MatchFile {
public:
    MatchFile();
    ~MatchFile();

protected:
    int GetPosFrmSnx(string file);
    int GetPosFrmCrd(string file);
    int GetPosFrmSta(string file);
    int MatchProd();
    int MatchPrec();
    int MatchComn();
    int MatchDCB();
    int MatchPos();
    int MatchOut();

public:
    int MatchFileAuto();

private:
    string dir_;
    char sep_;
    int year_,yy_;
    int week_,doy_,wod_;
};

class ReadRinex {
public:
    ReadRinex();
    ReadRinex(string& file);
    virtual ~ReadRinex();

public:
    int TestOpen();
    void CloseFile();
    void SetSysMark();
    int ReadRnxHead();
    virtual int ReadHead();
    virtual int ReadBody();

public:
    ifstream inf_;
    double ver_;
    string type_;
    int sat_sys_;
    int time_sys_;
    Time ts_;
    Time te_;
    int sys_mask_;
    double tint_;
    string buff_;
};

class ReadRnxObs:public ReadRinex {
public:
    ReadRnxObs();
    ReadRnxObs(string& file);
    ReadRnxObs(string& file,int rcv);
    ~ReadRnxObs();

protected:
    static void SaveSlips(unsigned char slips[][NFREQ],ObsData& data);
    static void ReStoreSlips(unsigned char slips[][NFREQ],ObsData& data);
    static void SetSigIndex(int sys,string type_obs[MAXOBSTYPE],Signal& sig_idx);
    int ReadObsBody(int& flag,vector<ObsData>& data);
    int DecodeObsEpoch(Time& t, int& flag,vector<SatMask> &sats);
    int DecodeObsData(ObsData &obs);

public:
    void InitReadObs(string file,Time ts,Time te,int tint);
    virtual int ReadHead(Nav *nav,Sta* sta);
    virtual int ReadBody(Obss& obss);

public:
    int rcv_;
    string type_obs_[NSYS][MAXOBSTYPE];
    Signal sig_index_[NSYS];
};

class ReadRnxNav:public ReadRinex {
public:
    ReadRnxNav();
    ReadRnxNav(string file);
    ~ReadRnxNav();

protected:
    int UraIndex(double val);
    void UniqeGloEph(Nav& navs);
    void UniqeEph(Nav& navs);
    int ReadNavBody(Nav &navs,int& eph_type);
    int DecodeEph(Time toc,SatMask sat,Eph& eph);
    int DecodeGloEph(Time toc,SatMask sat,GloEph& glo_eph);

public:
    void InitReadNav(string file);
    virtual int ReadHead(Nav& navs);
    virtual int ReadBody(Nav& navs);

public:
    double data_[64];
    Eph eph0_;
    GloEph glo_eph0_;
};

class ReadPreEph {
public:
    ReadPreEph();
    ~ReadPreEph();

protected:
    void SetSysMark();
    void UniqePreEph(vector<PreEph>& peph);
    int ReadHead();
    int ReadBody(Nav &nav);

public:
    int InitReadPreEph(string file);
    int ReadPre(Nav &nav);

protected:
    PreEph peph0_;
    PreClk pclk0_;

public:
    Time time_;
    int sys_mask_;
    string file_;
    ifstream inf_;
    string buff_;
    string type_;
    int ns_;
};

// read bias products
// DCB: CODE,CAS
// OSB: CODE,CAS,CNES
class ReadBias {
public:
    ReadBias();
    ~ReadBias();

public:
    int AddOsb(const Osb_t* osb);
    int OsbTimeStr2Time(const char *s,int i,int n,Time& t);
    int DecodeDCB(string file,double dcb[MAXSAT][MAXCBIASPAIR]);
    int ReadCODDCB(double dcb[MAXSAT][MAXCBIASPAIR]);
    int ReadCASDCB(string file,double dcb[MAXSAT][MAXCBIASPAIR]);
    int ReadOSB(string file);

public:
    int InitReadBias();
    int ReadBiasFile(double dcb[MAXSAT][MAXCBIASPAIR]);

public:
    int type_;          // DCB OSB
    int ac_;            // CODE CAS CNES

    SatOsb_t sat_osbs_;
    ifstream inf_;
};

class ReadGIM {
public:
    ReadGIM();
    ~ReadGIM();

protected:
    int DataIndex(int i,int j,int k,const int *ndata);
    int GetIndex(double val,const double *range);
    int GetNumItems(const double *range);
    Tec* AddTec(vector<Tec>& tecs);
    int ReadHead();
    int ReadBody(vector<Tec>& tecs);

public:
    void InitReadGIM(string file);
    int ReadTec(vector<Tec>& tecs);

protected:
    Time ion_time_;
    double factor_;
    string buff_;

public:
    string ion_file_;
    ifstream inf_;
    double lats_[3],lons_[3],hgts_[3];
    double re_;
};

class ReadErp {
public:
    ReadErp();
    ~ReadErp();

public:
    void InitReadErp(string file);
    int ReadErpPara(vector<Erpd_t>& erps);

public:
    ifstream inf_;
    string buff_;
};

class ReadBlq {
public:
    ReadBlq();
    ~ReadBlq();

public:
    void InitReadBlq(const string& file);
    int ReadBlqPara(string site,double *ocean_par);

public:
    ifstream inf_;
    string buff_;
};

class ReadAnt {
public:
    ReadAnt();
    ~ReadAnt();

protected:
    int DecodeAntPar(char *p,int n,double *v);

public:
    void InitReadAtx(const string& file);
    int ReadAtx(vector<Ant_t>& ants);

public:
    ifstream inf_;
    string buff_;
};

#endif //MSNAVS_READFILES_H
