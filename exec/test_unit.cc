#include "CmnFunc.h"
#include "Config.h"
#include "PPPLibGlo.h"
#include "ReadFiles.h"
INITIALIZE_EASYLOGGINGPP

void test_nav(){
    string path="../example/brdm3350.19p";
    ReadRnxNav read_nav;
    read_nav.InitReadNav(path);
//    Nav::Ptr_  nav=Nav::GetInstance();
//    read_nav.ReadHead(nav);
}

void test_config(){
    string path="../conf/PPPLib.ini";
    double ts,te;
    Config::Ptr_ config=Config::GetInstance();
    config->Open(path);
    kPrcOpt.SetPrcOpt();
    cout<<kPrcOpt.dir_<<endl;
    cout<<kPrcOpt.prc_date_.epoch_[0]<<" "<<kPrcOpt.prc_date_.epoch_[1]<<" "<<kPrcOpt.prc_date_.epoch_[2]<<endl;
    cout<<kPrcOpt.ts_.Time2Str(3)<<endl;
    cout<<kListModesStr[kPrcOpt.mode_]<<endl;
    cout<<kListModOptsStr[kPrcOpt.mode_opt_]<<endl;
    Isat_t isat {};
    int a=1+1;
    cout<<1;
}

int main(int argc,char const *argv[]){
    test_nav();
}