/**************************************************************************

Copyright(c): 2020, by Chao Chen, All rights reserved.
              China University of Mining and Technology, XuZhou, P.R. China

Author: Chao Chen(cchen@cumt.edu.cn)

Date: 2020-05-30

Description:

Note: refer to  https://github.com/2013fangwentao/Multi_Sensor_Fusion
**************************************************************************/

#include "Config.h"

Config::Ptr_ Config::config_info_(new Config());

bool Config::Open(std::string config_file) {
    ifstream inf(config_file,ios::in);
    if(!inf){
        cout<<"Configure file path error";
        return false;
    }
    while(!inf.eof()){
        string line,key;
        getline(inf,line);
        if(line.size()==0) continue;
        else if(line.substr(0,1)=="#"||line.substr(0,1)=="*"||line.substr(0,1)=="!") continue;
        auto data=TextSplit(line,";");
        key=StrTrim(data[0]);
        transform(key.begin(),key.end(),key.begin(),::tolower);
        if(data.size()>1){
            storage_[key]=StrTrim(data[1]);
        }
    }
    return true;
}

Config::Ptr_ Config::GetInstance() {
    return config_info_;
}