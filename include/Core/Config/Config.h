//
// Created by cc on 5/30/20.
//

#ifndef PPPLIB_CONFIG_H
#define PPPLIB_CONFIG_H

#include "CmnFunc.h"
#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <vector>
#include <exception>


class Config {
public:
    using Ptr_=std::shared_ptr<Config>;

private:
    Config() = default;
    Config(Config &&) = delete;
    Config(const Config &)= delete;
    Config &operator=(Config &&)= delete;
    Config &operator=(const Config &)= delete;
    static Ptr_ config_info_;
    std::map<std::string, std::string> storage_;

public:
    ~Config() = default;

public:
    static Ptr_ GetInstance();
    bool Open(std::string config_file);
    template <typename T>
    T Get(std::string key){
        transform(key.begin(),key.end(),key.begin(),::tolower);
        if(storage_.count(key)>0){
            try{
                double value=stod(storage_[key]);
                return static_cast<T>(value);
            }
            catch (const std::exception &e){
                std::cerr<<e.what()<<'\n';
            }
        }
        else{
//            LOG(ERROR)<<"The key of "<<key<<" does not exist";
//            getchar();
            return T(0x0);
        }
    }
    template <typename  T>
    std::vector<T> GetArray(std::string key){
        std::vector<T> data;
        transform(key.begin(),key.end(),key.begin(),::tolower);
        if(storage_.count(key)>0){
            try{
                auto text=TextSplit(storage_[key],",");
                for(auto index:text){
                    double value=stod(index);
                    data.emplace_back(static_cast<T>(value));
                }
            }
            catch(const std::exception &e){
                std::cerr<<e.what()<<'\n';
            }
        }
        else{
//            LOG(ERROR)<<"The key of "<<key<<" does not exist";
//            getchar();
        }
        return data;
    }
};

template <>
inline std::string Config::Get<std::string>(std::string key){
    transform(key.begin(),key.end(),key.begin(),::tolower);
    if(storage_.count(key)>0){
        try{
            return std::string(storage_[key]);
        }
        catch (const std::exception &e){
            std::cerr<<e.what()<<'\n';
        }
    }
    else{
//        LOG(ERROR)<<"The key of "<<key<<" does not exist"<<endl;
        getchar();
        return "";
    }
}

template <>
inline std::vector<std::string> Config::GetArray<std::string>(std::string key){
    std::vector<std::string> data;
    transform(key.begin(),key.end(),key.begin(),::tolower);
    if(storage_.count(key)>0){
        try{
            data=TextSplit(storage_[key],",");
        }
        catch(const std::exception &e){
            std::cerr<<e.what()<<'\n';
        }
    }
    else{
//        LOG(ERROR)<<"The key of "<<key<<" does not exist"<<endl;
        getchar();
    }
    return data;
}


#endif //PPPLIB_CONFIG_H
