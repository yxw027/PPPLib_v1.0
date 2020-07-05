/**************************************************************************

Copyright(c): 2020, by Chao Chen, All rights reserved.
              China University of Mining and Technology, XuZhou, P.R. China

Author: Chao Chen(cchen@cumt.edu.cn)

Date: 2020-05-30

Description:

**************************************************************************/

#include "PPPLibGlo.h"
#include "Solver.h"
#include "dirent.h"

INITIALIZE_EASYLOGGINGPP
MainSolver *main_solver= new MainSolver;


static int PreTreatment(string obs_file) {
    char f[MAXSTRPATH]={'\0'};
    int dgnss=(kPrcOpt.mode_==DGPS||kPrcOpt.mode_==PPK)||(kPrcOpt.mode_opt_==MDOPT_DGPS||kPrcOpt.mode_opt_==MDOPT_PPK);
    if(dgnss) {
        vector<string> splits=MultiSplitStr(obs_file,".");
        sprintf(f,"%s_base.%s",splits[0].c_str(),splits[1].c_str());
        if((access(f,0))==-1){
            defaultLogger->fatal("Reference observation file(%v) no exist",f);
            return 0;
        }
        kFileOpt.base_=f;
        LOG(INFO)<<"Reference observation file :"<<f;
    }

    MatchFile match;
    match.MatchFileAuto();
}

static void Processor() {
    DIR *dir;
    struct dirent *file;
    int year,week,doy,i=0;
    char *ext;
    string sep;
#if WIN32
    sep="\\";
#else
    sep="/";
#endif

    kPrcOpt.prc_date_.Time2Epoch();
    year=kPrcOpt.prc_date_.epoch_[0];
    kPrcOpt.prc_date_.Time2GPST(&week);
    doy=(int)kPrcOpt.prc_date_.Time2Doy();

    string obs_dir;
    char f[MAXSTRPATH]={'\0'};

    if(!kFileOpt.customize_file_){
        sprintf(f,"%s%c%d%c%d%c%03d%cobs",kPrcOpt.dir_.c_str(),FILEPATHSEP,year,FILEPATHSEP,week,FILEPATHSEP,doy,FILEPATHSEP);
        obs_dir=f;
    }else obs_dir=kFileOpt.customize_dir_;

    LOG(DEBUG)<<"Observation data dir: "<< obs_dir;
    if(!(dir=opendir(obs_dir.c_str()))){
        LOG(FATAL)<<"Observation path not exist: "<<obs_dir;
        return;
    }

    while((file=readdir(dir))!=NULL){
        if(strncmp(file->d_name,".",1)==0)           continue;
        else if(strstr(file->d_name,"base"))  continue;
        else if(!(ext=strrchr(file->d_name,'.'))) continue;
        else if(!strstr(ext+3,"o"))   continue;
        f[0]='\0';
        sprintf(f,"%s%c%s",obs_dir.c_str(),FILEPATHSEP,file->d_name);
        kFileOpt.rover_=f;

        LOG(DEBUG)<<"Rover observation file: "<<f;
        vector<string> splits=MultiSplitStr(f,sep);
        kSite=(splits.end()-1)->substr(0,4);
        LOG(INFO)<<"------ Start processing  "<<"["<<kSite<<"]"<<" ------";

        long t1,t2;
        t1=clock();

        if(!PreTreatment(kFileOpt.rover_)){
            LOG(ERROR)<<"Data preparation fail";
            return;
        }

        main_solver->site_count_=++i;
        if(!main_solver->InitMainSolver()){
            LOG(INFO)<<"Initialize solver fail";
            return;
        }
        if(main_solver->EpochSol()){
            t2=clock();
            double t=(double)(t2-t1)/CLOCKS_PER_SEC;
            LOG(INFO)<<"Time consuming: "<<setprecision(2)<<t<<"s";
            LOG(INFO)<<"------ Processing complete  ^_^ ------";
        }else{
            LOG(INFO)<<"------ Processing fail      T_T ------";
            continue;
        }
    }

    closedir(dir);
}

static int ParsePara(int arc, char *arv[],string& conf_file) {
    double epoch[6]={0};
    const char *p;
    string pos_mode;
    int mask=NONE;

    for(int i=0;i<arc;i++) {
        if(!strcmp(arv[i],"-cf")&&i+1<arc)  conf_file=arv[++i];
        else if(!strcmp(arv[i],"-do")&&i+1<arc) kPrcOpt.use_def_opt_=atoi(arv[++i]);
        else if(!strcmp(arv[i],"-md")&&i+1<arc){
            pos_mode=arv[++i];
        }
        else if(!strcmp(arv[i],"-pd")&&i+1<arc){
            if(sscanf(arv[++i],"%lf/%lf/%lf",epoch,epoch+1,epoch+2)<3){
                fprintf(stderr,"Process date format error: yyyy/mm/dd \n");
                return 0;
            }
        }
        else if(!strcmp(arv[i],"-sys")&&i+1<arc){
            p=arv[++i];
            for(;*p&&*p!=' ';p++){
                switch(*p){
                    case 'G': mask|=GPS; break;
                    case 'B': mask|=(BD2|BD3);break;
                    case 'E': mask|=GAL;break;
                    case 'R': mask|=GLO;break;
                    case 'J': mask|=QZS;break;
                }
            }
            kPrcOpt.GNSS_opt_.nav_sys_=mask;
        }
        else if(!strcmp(arv[i],"-level")&&i+1<arc){
            kPrcOpt.debug_level_=atoi(arv[++i]);
        }
        else if(!strcmp(arv[i],"-frq")&&i+1<arc){
            kPrcOpt.GNSS_opt_.use_frq_=atoi(arv[++i]);
        }
        else if(!strcmp(arv[i],"-ion")&&i+1<arc){
            kPrcOpt.GNSS_opt_.ion_opt_=atoi(arv[++i]);
        }
    }
    kPrcOpt.prc_date_.Epoch2Time(epoch);

    vector<string>mode_opt;
    SplitString(pos_mode, mode_opt,  "-");
    for(int i=0;kListModesStr;i++){
        if(mode_opt[0]==kListModesStr[i]){
            kPrcOpt.mode_=i;break;
        }
    }
    for(int i=0;kListModOptsStr;i++){
        if(mode_opt[1]==kListModOptsStr[i]){
            kPrcOpt.mode_opt_=i;break;
        }
    }

    if(!kPrcOpt.use_def_opt_&&conf_file.empty()) {
        fprintf(stderr,"Custom configure file is missing, use default settings \n");
        kPrcOpt.use_def_opt_=true;
    }
    return 1;
}


void InitLog(int argc,char *argv[]){
    START_EASYLOGGINGPP(argc,argv);
    el::Configurations conf(LOGINI_PATH);
    el::Loggers::reconfigureAllLoggers(conf);
    el::Loggers::addFlag(el::LoggingFlag::StrictLogFileSizeCheck);
    el::Loggers::addFlag(el::LoggingFlag::ColoredTerminalOutput);
    el::Loggers::addFlag(el::LoggingFlag::HierarchicalLogging);
    el::Loggers::setLoggingLevel(static_cast<el::Level>(kPrcOpt.debug_level_));
    defaultLogger=el::Loggers::getLogger("default");
}

void PrintHelp() {
    fprintf(stderr,"PPPLib Usage: \n"
                   "-do       whether used default options(int 0:no 1:use, default use)\n"
                   "[-cf]     customize configure file path(string)\n"
                   "-md       PPPLib processing mode(string SPP-KINE, PPP-STATIC, PPP-KINE, ...)\n"
                   "-pd       process date(yyyy/mm/dd 2019/12/01)\n"
                   "[-sys]    use GNSS option(string GBERJ, default use GPS)\n"
                   "[-level]  output debug information level(int 4: Debug, 8: Fatal, 16: Error, 32: Warning, 128: Info, default Info)\n");
}

int main(int argc,char *argv[]){
    string conf_file;

    if(!ParsePara(argc,argv,conf_file)){
        PrintHelp();
        return 0;
    }

    InitLog(argc,argv);

    if(kPrcOpt.use_def_opt_){
        if(kPrcOpt.mode_==PPP) kPrcOpt.SetDefaultPPPOpt();
        else if(kPrcOpt.mode_==PPK) kPrcOpt.SetDefaultPPKOpt();
    }
    else{
        Config::Ptr_ config=Config::GetInstance();
        config->Open(conf_file);
        kPrcOpt.SetPrcOpt();
        kFileOpt.SetFileOpt();
        if(kFileOpt.customize_file_) kPrcOpt.dir_=kFileOpt.customize_dir_;
    }

    Processor();

    return 1;
}