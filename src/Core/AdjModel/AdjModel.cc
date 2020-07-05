/**************************************************************************

Copyright(c): 2020, by Chao Chen, All rights reserved.
              China University of Mining and Technology, XuZhou, P.R. China

Author: Chao Chen(cchen@cumt.edu.cn)

Date: 2020-05-30

Description:

**************************************************************************/

#include <GNSS/GNSS.h>
#include "AdjModel.h"

AdjModel::AdjModel() {
    epoch_=0;
    v_.clear();dx_.clear();
    v_CP_.clear();v_PR_.clear();
    norm_CPv_.clear();norm_PRv_.clear();
    Qx_.clear();Qvv_.clear();
    sigma0_=0.0;
    hel_iter_=0;
    use_sys_=0;
    for(int i=0;i<NSYS;i++){
        for(int j=0;j<2*NFREQ;j++){
            num_obs_[i][j]=0;
        }
    }
}

AdjModel::~AdjModel() {
    num_CP_=num_PR_=0;
    Qx_.clear();Qvv_.clear();
    dx_.clear();v_.clear();norm_PRv_.clear();norm_CPv_.clear();
    v_PR_.clear();v_CP_.clear();
}

int AdjModel::Lsq(const vector<double> &L, const vector<double> &A, int nl, int nx, const vector<double> &R,
                  vector<double> &X, vector<double> &Px) {
    v_.clear();
    v_PR_.clear();v_CP_.clear();norm_PRv_.clear();norm_CPv_.clear();
    dx_.clear();Qvv_.clear();Qx_.clear();CPv_info_.clear();PRv_info_.clear();

    int flag=0;
    vector<double> P,AP(nl*nx,0.0),APL(nx,0.0),Q(nx*nx,0.0);
    P.assign(R.begin(),R.end());
    dx_.assign(nx,0.0);
    if(MatInv(P,nl)==-1) return -1;

    // compute AP,APA^T 这里是列优先
    MatMulVec("NN",nx,nl,nl,1.0,A,P,0.0,AP); //AP
    MatMulVec("NT",nx,1,nl,1.0,AP,L,0.0,APL);  //APL
    MatMulVec("NT",nx,nx,nl,1.0,AP,A,0.0,Q);
    if (MatInv(Q,nx)==-1) return -1;
    /* compute dX */
    MatMulVec("NN",nx,1,nx,1.0,Q,APL,0.0,dx_);  //dX

    for(int i=0;i<nx;i++) X[i]+=dx_[i];

    vector<double> LP(nl,1),VPV(1,0.0);
    sigma0_=0.0;
    MatMulVec("NN",1,nl,nl,1.0,L,P,0.0,LP); //LP
    MatMulVec("NT",1,1,nl,1.0,LP,L,0.0,VPV); //LPL to VPV
    MatMulVec("NT",1,1,nx,-1.0,APL,dx_,1.0,VPV); //VPV = LPL-(APL)X

    // computer v/Qvv
    vector<double>AQ(nl*nx,0.0);
    Qvv_.assign(nl*nl,0.0);Qvv_.assign(R.begin(),R.end());
    v_.assign(nl,0.0);v_.assign(L.begin(),L.end());
    norm_PRv_.assign(nl,0.0);
    MatMulVec("TN",nl,1,nx,1.0,A,dx_,-1.0,v_);
    v_PR_.assign(v_.begin(),v_.end());
    MatMulVec("TN",nl,nx,nx,1.0,A,Q,0.0,AQ);
    MatMulVec("NN",nl,nl,nx,-1.0,AQ,A,1.0,Qvv_);
    for(int i=0;i<nl;i++) norm_PRv_[i]=v_[i]/sqrt(fabs(Qvv_[i+i*nl]));

    // sigma
    sigma0_=VPV[0]/(nl>nx ? nl-nx : 1);
    Px.assign(Q.begin(),Q.end());
    for(size_t i=0;i<Qx_.size();i++) Px[i]=sigma0_*Px[i];

    flag=Norm(dx_.begin(),nx)<1E-4?1:0;

    P.clear();AP.clear();APL.clear();Q.clear();
    LP.clear();VPV.clear();AQ.clear();
    return flag;
}

int AdjModel::Adjustment(const vector<double>& L,const vector<double>& A,int nl,int nx,
                         const vector<double>& R,vector<double>& X,vector<double>& Px) {

}

LSQAdj::LSQAdj() {

}

LSQAdj::~LSQAdj() {

}

int LSQAdj::Adjustment(const vector<double>& L,const vector<double>& A,int nl,int nx,
                       const vector<double>& R,vector<double>& X,vector<double>& Px) {
    return Lsq(L,A,nl,nx,R,X,Px);
}

KFAdj::KFAdj() {
    sigma0_=0.0;
}

KFAdj::~KFAdj() {

}

int KFAdj::Adjustment(const vector<double>& L,const vector<double>& A,int nl,int nx,
                      const vector<double>& R,vector<double>& X,vector<double>& Px) {

    v_.clear();
    v_PR_.clear();v_CP_.clear();norm_PRv_.clear();norm_CPv_.clear();
    dx_.clear();Qvv_.clear();Qx_.clear();CPv_info_.clear();PRv_info_.clear();

    vector<double> PxA(nl*nx,0.0),Q(R.begin(),R.end()),K(nl*nx,0.0);
    vector<double> E(nx*nx,0.0),Px0;

    EyeMat(E.begin(),nx);Px0.assign(Px.begin(),Px.end());
    dx_.assign(nx,0.0);
    MatMulVec("NN",nx,nl,nx,1.0,Px,A,0.0,PxA); //P*A'
    MatMulVec("TN",nl,nl,nx,1.0,A,PxA,1.0,Q);  //Q=A*P*A'
    // A 在程序里是列优先 nx*nl,实际应是nl*nx,两者互为倒置

    if(MatInv(Q,nl)==0){
        MatMulVec("NN",nx,nl,nl,1.0,PxA,Q,0.0,K); // K=P*A'*Q^-1;
        MatMulVec("NN",nx,1,nl,1.0,K,L,0.0,dx_); //dx=K*L
        for(int i=0;i<nx;i++) X[i]+=dx_[i];
        MatMulVec("NT",nx,nx,nl,-1.0,K,A,1.0,E);    //E=E-K*A
        MatMulVec("NN",nx,nx,nx,1.0,E,Px0,0.0,Px);  //Px=E*Px
        // v=Adx-l;
        v_.assign(L.begin(),L.end());
        MatMulVec("TN",nl,1,nx,1.0,A,dx_,-1.0,v_);
    }
    else return -1;

#if 0
    // Qvv=R-AN^(-1)A^T
    vector<double>P,AP,N,AN;
    P.assign(R.begin(),R.end());
    AP.assign(nx*nl,0.0);
    N.assign(nx*nx,0.0);
    AN.assign(nl*nx,0.0);
    Qvv_.assign(R.begin(),R.end());

    if(MatInv(P,nl)==0){
        MatMulVec("NN",nx,nl,nl,1.0,A,P,0.0,AP);
        MatMulVec("NT",nx,nx,nl,1.0,AP,A,0.0,N);
        if(MatInv(N,nx)==0){
            MatMulVec("TN",nl,nx,nx,1.0,A,N,0.0,AN);
            MatMulVec("NN",nl,nl,nx,1.0,AN,A,-1.0,Qvv_);
        }
    }
    else return -1;
    P.clear();VPV.clear();N.clear();AN.clear();
#endif

    // sigma
    vector<double>P, VPV ,VP;
    VP.assign(nl,0.0);
    VPV.assign(1,0.0);
    P.assign(R.begin(),R.end());
    MatInv(P,nl);
    MatMulVec("NN",1,nl,nl,1.0,v_,P,0.0,VP);
    MatMulVec("NT",1,1,nl,1.0,VP,v_,0.0,VPV);
    sigma0_=VPV[0]/(nl>nx?nl-nx:1);
    VP.clear();P.clear();

    int type=0;
    for(int i=0;i<nl;i++){
        type=(v_info_[i]>>4)&0x0F;
        if(type==OBS_CP){
            v_CP_.push_back(v_[i]);
            norm_CPv_.push_back(v_[i]/sqrt((R[i+i*nl])));
            CPv_info_.push_back(v_info_[i]);
        }
        else if(type==OBS_PR){
            v_PR_.push_back(v_[i]);
            norm_PRv_.push_back(v_[i]/sqrt(fabs(R[i+i*nl])));
            PRv_info_.push_back(v_info_[i]);
        }
    }

    PxA.clear();Q.clear();K.clear();E.clear();Px0.clear();
    return 0;
}

HelmertAdj::HelmertAdj() {

}

HelmertAdj::~HelmertAdj() {

}

int HelmertAdj::HelmertEst(const vector<double> &A, const vector<double> &Px, int nl, int nx) {

    vector<double>new_R(hel_R_.begin(),hel_R_.end()),new_Px(Px.begin(),Px.end());
    vector<double>QQ(nx*nx,0.0),Qi[NSYS],Ai,Pi,VPi,APi,W,S,sigma;

    if(MatInv(new_R,nl)==-1||MatInv(new_Px,nx)==-1) return -1;
    int singe_sys_obs[NSYS]={0};
    for(int sys=0,use_sys_=0,n=0;sys<NSYS;sys++){
        for(int f=0;f<NFREQ*2;f++) singe_sys_obs[use_sys_]+=num_obs_[sys][f];
        if(singe_sys_obs[use_sys_]==0) continue;

        VPi.assign(singe_sys_obs[use_sys_],0.0);Pi.assign(singe_sys_obs[use_sys_],0.0);
        APi.assign(nx*singe_sys_obs[use_sys_],0.0);
        Ai.assign(A.begin()+n*nx,A.begin()+(n+singe_sys_obs[use_sys_])*nx);
        Qi[use_sys_].assign(nx*nx,0.0);W.push_back(0.0);

        for(int i=0;i<singe_sys_obs[use_sys_];i++){
            for(int j=0;j<singe_sys_obs[use_sys_];j++){
                Pi[j+i*singe_sys_obs[use_sys_]]=new_R[j+n*(i+n)*nl];
            }
        }
        MatMulVec("NN",1,singe_sys_obs[use_sys_],singe_sys_obs[use_sys_],1.0,v_,Pi,0.0,VPi,n,0,0);
        MatMulVec("NT",1,1,singe_sys_obs[use_sys_],1.0,VPi,v_,0.0,W,0,n,use_sys_);

        MatMulVec("NN",nx,singe_sys_obs[use_sys_],singe_sys_obs[use_sys_],1.0,Ai,Pi,0.0,APi);
        MatMulVec("NT",nx,nx,singe_sys_obs[use_sys_],1.0,APi,Ai,0.0,Qi[use_sys_]);
        for(int i=0;i<nx;i++) for(int j=0;j<nx;j++) QQ[j+i*nx]+=Qi[use_sys_][j+i*nx];
        n+=singe_sys_obs[use_sys_];use_sys_++;
    }

    if(use_sys_<2) return -1;
    for(int i=0;i<nx;i++) for(int j=0;j<nx;j++) QQ[j+i*nx]+=new_Px[j+i*nx];
    if(MatInv(QQ,nx)==-1) return -1;

    sigma.assign(use_sys_,0.0);S.assign(use_sys_*use_sys_,0.0);

    vector<double>QiQ(nx*nx,0.0),QjQ(nx*nx,0.0),QQQQ(nx*nx,0.0);
    for(int i=0;i<use_sys_;i++){
        MatMulVec("NN",nx,nx,nx,1.0,Qi[i],QQ,0.0,QiQ);
        for(int j=0;j<use_sys_;j++){
            MatMulVec("NN",nx,nx,nx,1.0,Qi[j],QQ,0.0,QjQ);
            MatMulVec("NN",nx,nx,nx,1.0,QiQ,QjQ,0.0,QQQQ);
            S[j+i*use_sys_]+=MatTrace(QQQQ.begin(),nx);
            if(j==i) S[j+i*use_sys_]+=singe_sys_obs[i]-2.0*MatTrace(QiQ,nx);
        }
    }

    if(SolverLineEq("T",S,W,sigma,use_sys_,1)) return -1;

    for(int i=use_sys_;i>0;i--) sigma[i-1]/=sigma[0];
    if(hel_iter_==0) sgm2.assign(sigma.begin(),sigma.end());
    else for(int i=0;i<use_sys_;i++) sgm2[i]*=sigma[i];

    for(int sys=0,n=0;sys<use_sys_;sys++){
        for(int i=0;i<singe_sys_obs[sys];i++) for(int j=0;j<singe_sys_obs[sys];j++)
            new_R[j+n+(i+n)*nl]*=sigma[sys];
        n+=singe_sys_obs[sys];
    }

    new_R.clear();new_Px.clear();QQ.clear();APi.clear();Ai.clear();Pi.clear();VPi.clear();
    W.clear();S.clear();sigma.clear();
    for(int i=0;i<NSYS;i++) Qi[i].clear();

    return 0;
}

int HelmertAdj::Adjustment(const vector<double> &L, const vector<double> &A, int nl, int nx, const vector<double> &R,
        vector<double> &X, vector<double> &Px) {

    int stat=0;
    vector<double>ori_Px,ori_X;
    hel_R_.assign(R.begin(),R.end());

    for(hel_iter_=0;hel_iter_<3;hel_iter_++){
        ori_X.assign(X.begin(),X.end());ori_Px.assign(Px.begin(),Px.end());
        if(KFAdj::Adjustment(L,A,nl,nx,R,X,Px)){
            stat=KFAdj::Adjustment(L,A,nl,nx,R,X,Px);
            ori_Px.clear();ori_X.clear();
            return stat;
        }
        if(HelmertEst(A,Px,nl,nx)){
            stat=KFAdj::Adjustment(L,A,nl,nx,R,X,Px);
            ori_Px.clear();ori_X.clear();
            return stat;
        }
    }
    stat=KFAdj::Adjustment(L,A,nl,nx,R,X,Px);
    ori_Px.clear();ori_X.clear();
    return stat;
}