//
// Created by cc on 4/9/20.
//

#include "TrpModel.h"
#include "PPPLibGlo.h"

TrpModel::TrpModel() {
    z_wet_=z_hyd_=0.0;
    map_w_=map_h_=0.0;
    P_=T_=e_=0.0;
    sat_data_= nullptr;
}

TrpModel::~TrpModel() {

}

void TrpModel::StandAtmoPara(const Vec3 blh,const double hgt,const double humi) {
    const double temp0=18.0;                /* temparature at sea level */

    P_=1013.25*pow(1.0-2.2557E-5*hgt,5.225);
    T_=temp0-6.5E-3*hgt+273.16;
    double humidity=0.5*exp(-6.396E-4*hgt);
    if (humi==0.0) humidity=0.0;
    e_=humidity*exp(-37.2465+0.213166*T_-0.000256908*T_*T_);
}

void TrpModel::TrpMapEl(double z) {
    map_h_=map_w_=1.0/cos(z);
}

void TrpModel::TrpMapGMF(Time t,const Vec3 pos,double el) {
    double dfac[20],P[10][10],aP[55],bP[55],tt;
    int i,j,ir,k,n,m,nmax,mmax;
    double phh;
    double ah,bh,ch,aw,bw,cw;
    double ahm,aha,awm,awa;
    double c10h,c11h,c0h;
    double a_ht,b_ht,c_ht;
    double sine,beta,gamma,topcon;
    double hs_km,ht_corr,ht_corr_coef;
    double sum,dt;
    double Pi2=2*PI;
    map_h_=map_w_=0.0;

    static double ah_mean[55] = {
            +1.2517e+02,	+8.503e-01,	+6.936e-02,	-6.760e+00, +1.771e-01,
            +1.130e-02,		+5.963e-01,	+1.808e-02, +2.801e-03, -1.414e-03,
            -1.212e+00,		+9.300e-02,	+3.683e-03, +1.095e-03, +4.671e-05,
            +3.959e-01,		-3.867e-02, +5.413e-03, -5.289e-04, +3.229e-04,
            +2.067e-05,		+3.000e-01, +2.031e-02, +5.900e-03, +4.573e-04,
            -7.619e-05,		+2.327e-06, +3.845e-06, +1.182e-01, +1.158e-02,
            +5.445e-03,		+6.219e-05, +4.204e-06, -2.093e-06, +1.540e-07,
            -4.280e-08,		-4.751e-01, -3.490e-02, +1.758e-03, +4.019e-04,
            -2.799e-06,		-1.287e-06, +5.468e-07, +7.580e-08, -6.300e-09,
            -1.160e-01,		+8.301e-03, +8.771e-04, +9.955e-05, -1.718e-06,
            -2.012e-06,		+1.170e-08, +1.790e-08, -1.300e-09, +1.000e-10
    };

    static double bh_mean[55] = {
            +0.000e+00,	+0.000e+00, +3.249e-02, +0.000e+00, +3.324e-02,
            +1.850e-02, +0.000e+00, -1.115e-01, +2.519e-02, +4.923e-03,
            +0.000e+00, +2.737e-02, +1.595e-02, -7.332e-04, +1.933e-04,
            +0.000e+00, -4.796e-02, +6.381e-03, -1.599e-04, -3.685e-04,
            +1.815e-05, +0.000e+00, +7.033e-02, +2.426e-03, -1.111e-03,
            -1.357e-04, -7.828e-06, +2.547e-06, +0.000e+00, +5.779e-03,
            +3.133e-03, -5.312e-04, -2.028e-05, +2.323e-07, -9.100e-08,
            -1.650e-08, +0.000e+00, +3.688e-02, -8.638e-04, -8.514e-05,
            -2.828e-05, +5.403e-07, +4.390e-07, +1.350e-08, +1.800e-09,
            +0.000e+00, -2.736e-02, -2.977e-04, +8.113e-05, +2.329e-07,
            +8.451e-07, +4.490e-08, -8.100e-09, -1.500e-09, +2.000e-10
    };

    static double ah_amp[55] = {
            -2.738e-01, -2.837e+00, +1.298e-02, -3.588e-01, +2.413e-02,
            +3.427e-02, -7.624e-01, +7.272e-02, +2.160e-02, -3.385e-03,
            +4.424e-01, +3.722e-02, +2.195e-02, -1.503e-03, +2.426e-04,
            +3.013e-01, +5.762e-02, +1.019e-02, -4.476e-04, +6.790e-05,
            +3.227e-05, +3.123e-01, -3.535e-02, +4.840e-03, +3.025e-06,
            -4.363e-05, +2.854e-07, -1.286e-06, -6.725e-01, -3.730e-02,
            +8.964e-04, +1.399e-04, -3.990e-06, +7.431e-06, -2.796e-07,
            -1.601e-07, +4.068e-02, -1.352e-02, +7.282e-04, +9.594e-05,
            +2.070e-06, -9.620e-08, -2.742e-07, -6.370e-08, -6.300e-09,
            +8.625e-02, -5.971e-03, +4.705e-04, +2.335e-05, +4.226e-06,
            +2.475e-07, -8.850e-08, -3.600e-08, -2.900e-09, +0.000e+00
    };

    static double bh_amp[55] = {
            +0.000e+00, +0.000e+00, -1.136e-01, +0.000e+00, -1.868e-01,
            -1.399e-02, +0.000e+00, -1.043e-01, +1.175e-02, -2.240e-03,
            +0.000e+00, -3.222e-02, +1.333e-02, -2.647e-03, -2.316e-05,
            +0.000e+00, +5.339e-02, +1.107e-02, -3.116e-03, -1.079e-04,
            -1.299e-05, +0.000e+00, +4.861e-03, +8.891e-03, -6.448e-04,
            -1.279e-05, +6.358e-06, -1.417e-07, +0.000e+00, +3.041e-02,
            +1.150e-03, -8.743e-04, -2.781e-05, +6.367e-07, -1.140e-08,
            -4.200e-08, +0.000e+00, -2.982e-02, -3.000e-03, +1.394e-05,
            -3.290e-05, -1.705e-07, +7.440e-08, +2.720e-08, -6.600e-09,
            +0.000e+00, +1.236e-02, -9.981e-04, -3.792e-05, -1.355e-05,
            +1.162e-06, -1.789e-07, +1.470e-08, -2.400e-09, -4.000e-10
    };

    static double aw_mean[55] = {
            +5.640e+01, +1.555e+00, -1.011e+00, -3.975e+00, +3.171e-02,
            +1.065e-01, +6.175e-01, +1.376e-01, +4.229e-02, +3.028e-03,
            +1.688e+00, -1.692e-01, +5.478e-02, +2.473e-02, +6.059e-04,
            +2.278e+00, +6.614e-03, -3.505e-04, -6.697e-03, +8.402e-04,
            +7.033e-04, -3.236e+00, +2.184e-01, -4.611e-02, -1.613e-02,
            -1.604e-03, +5.420e-05, +7.922e-05, -2.711e-01, -4.406e-01,
            -3.376e-02, -2.801e-03, -4.090e-04, -2.056e-05, +6.894e-06,
            +2.317e-06, +1.941e+00, -2.562e-01, +1.598e-02, +5.449e-03,
            +3.544e-04, +1.148e-05, +7.503e-06, -5.667e-07, -3.660e-08,
            +8.683e-01, -5.931e-02, -1.864e-03, -1.277e-04, +2.029e-04,
            +1.269e-05, +1.629e-06, +9.660e-08, -1.015e-07, -5.000e-10
    };

    static double bw_mean[55] = {
            +0.000e+00, +0.000e+00, +2.592e-01, +0.000e+00, +2.974e-02,
            -5.471e-01, +0.000e+00, -5.926e-01, -1.030e-01, -1.567e-02,
            +0.000e+00, +1.710e-01, +9.025e-02, +2.689e-02, +2.243e-03,
            +0.000e+00, +3.439e-01, +2.402e-02, +5.410e-03, +1.601e-03,
            +9.669e-05, +0.000e+00, +9.502e-02, -3.063e-02, -1.055e-03,
            -1.067e-04, -1.130e-04, +2.124e-05, +0.000e+00, -3.129e-01,
            +8.463e-03, +2.253e-04, +7.413e-05, -9.376e-05, -1.606e-06,
            +2.060e-06, +0.000e+00, +2.739e-01, +1.167e-03, -2.246e-05,
            -1.287e-04, -2.438e-05, -7.561e-07, +1.158e-06, +4.950e-08,
            +0.000e+00, -1.344e-01, +5.342e-03, +3.775e-04, -6.756e-05,
            -1.686e-06, -1.184e-06, +2.768e-07, +2.730e-08, +5.700e-09
    };

    static double aw_amp[55] = {
            +1.023e-01, -2.695e+00, +3.417e-01, -1.405e-01, +3.175e-01,
            +2.116e-01, +3.536e+00, -1.505e-01, -1.660e-02, +2.967e-02,
            +3.819e-01, -1.695e-01, -7.444e-02, +7.409e-03, -6.262e-03,
            -1.836e+00, -1.759e-02, -6.256e-02, -2.371e-03, +7.947e-04,
            +1.501e-04, -8.603e-01, -1.360e-01, -3.629e-02, -3.706e-03,
            -2.976e-04, +1.857e-05, +3.021e-05, +2.248e+00, -1.178e-01,
            +1.255e-02, +1.134e-03, -2.161e-04, -5.817e-06, +8.836e-07,
            -1.769e-07, +7.313e-01, -1.188e-01, +1.145e-02, +1.011e-03,
            +1.083e-04, +2.570e-06, -2.140e-06, -5.710e-08, +2.000e-08,
            -1.632e+00, -6.948e-03, -3.893e-03, +8.592e-04, +7.577e-05,
            +4.539e-06, -3.852e-07, -2.213e-07, -1.370e-08, +5.800e-09
    };

    static double bw_amp[55] = {
            +0.000e+00, +0.000e+00, -8.865e-02, +0.000e+00, -4.309e-01,
            +6.340e-02, +0.000e+00, +1.162e-01, +6.176e-02, -4.234e-03,
            +0.000e+00, +2.530e-01, +4.017e-02, -6.204e-03, +4.977e-03,
            +0.000e+00, -1.737e-01, -5.638e-03, +1.488e-04, +4.857e-04,
            -1.809e-04, +0.000e+00, -1.514e-01, -1.685e-02, +5.333e-03,
            -7.611e-05, +2.394e-05, +8.195e-06, +0.000e+00, +9.326e-02,
            -1.275e-02, -3.071e-04, +5.374e-05, -3.391e-05, -7.436e-06,
            +6.747e-07, +0.000e+00, -8.637e-02, -3.807e-03, -6.833e-04,
            -3.861e-05, -2.268e-05, +1.454e-06, +3.860e-07, -1.068e-07,
            +0.000e+00, -2.658e-02, -1.947e-03, +7.131e-04, -3.506e-05,
            +1.885e-07, +5.792e-07, +3.990e-08, +2.000e-08, -5.700e-09
    };

    t.Time2Mjd();
    double dmjd=t.mjd_.day+(t.mjd_.sod.sn+t.mjd_.sod.tos)/86400.0;
    double doy=dmjd-44239.0-27;
    double dlat=pos.i_,dlon=pos.j_,dhgt=pos.k_;
    tt=sin(dlat);
    nmax=mmax=9;

    //determine nmax!(faktorielle) moved by 1
    dfac[0]=1;
    for (i=1;i<=2*nmax+1;i++)
        dfac[i]=dfac[i-1]*i;
    //determine Legendre functions (Heiskanen and Moritz, Physical Geodesy, 1967, eq. 1-62)
    for (i=0;i<=nmax;i++) {
        for (j=0;j<=MIN(i,mmax);j++) {
            ir=(int)((i-j)/2);
            sum=0.0;
            for (k=0;k<=ir;k++) {
                sum=sum+pow(-1.0,k)*dfac[2*i-2*k]/dfac[k]/dfac[i-k]/dfac[i-j-2*k]*pow(tt,i-j-2*k);
            }
            //Legender functions moved by 1
            P[i][j]=1.0/pow(2.0,i)*sqrt(pow(1-tt*tt,j))*sum;
        }
    }
    //spherical harmonics
    i=0;
    for (n=0;n<=9;n++) {
        for (m=0;m<=n;m++) {
            i=i+1;
            dt=m*dlon;
            aP[i-1]=P[n][m]*cos(dt);
            bP[i-1]=P[n][m]*sin(dt);
        }
    }
    //hydrostatic
    bh=0.0029;
    c0h=0.062;
    if (dlat<0.0) {	//southern hemisphere
        phh=PI;
        c11h=0.007;
        c10h=0.002;
    }
    else {	//northern hemisphere
        phh=0;
        c11h=0.005;
        c10h=0.001;
    }
    ch=c0h+((cos(doy/365.25*Pi2+phh)+1.0)*c11h/2.0+c10h)*(1.0-cos(dlat));

    ahm=0.0;
    aha=0.0;
    for (i=1;i<=55;i++) {
        ahm=ahm+(ah_mean[i-1]*aP[i-1]+bh_mean[i-1]*bP[i-1])*1.0e-5;
        aha=aha+(ah_amp[i-1] *aP[i-1]+bh_amp[i-1] *bP[i-1])*1.0e-5;
    }
    ah=ahm+aha*cos(doy/365.25*Pi2);
    sine=sin(el);
    beta=bh/(sine+ch);
    gamma=ah/(sine+beta);
    topcon=(1.0+ah/(1.0+bh/(1.0+ch)));
    map_h_=topcon/(sine+gamma);
    //height correction for hydrostatic mapping function from Niell (1996)
    a_ht=2.53e-5;	//2.53 from http://maia.usno.navy.mil/conv2010/chapter9/GMF.F
    b_ht=5.49e-3;
    c_ht=1.14e-3;
    hs_km=dhgt/1000.0;
    beta=b_ht/(sine+c_ht);
    gamma=a_ht/(sine+beta);
    topcon=(1.0+a_ht/(1.0+b_ht/(1.0+c_ht)));
    ht_corr_coef=1.0/sine-topcon/(sine+gamma);
    ht_corr=ht_corr_coef*hs_km;
    map_h_+=ht_corr;

    //wet
    bw=0.00146;
    cw=0.04391;
    awm=0.0;
    awa=0.0;
    for (i=1;i<=55;i++) {
        awm=awm+(aw_mean[i-1]*aP[i-1]+bw_mean[i-1]*bP[i-1])*1e-5;
        awa=awa+(aw_amp[i-1] *aP[i-1]+bw_amp[i-1] *bP[i-1])*1e-5;
    }
    aw=awm+awa*cos(doy/365.25*Pi2);
    beta=bw/(sine+cw);
    gamma=aw/(sine+beta);
    topcon=(1.0+aw/(1.0+bw/(1.0+cw)));
    map_w_=topcon/(sine+gamma);
}

static double Interpc(const double coef[],double lat) {
    int i=(int)(lat/15.0);

    if(i<1) return coef[0];else if(i>4) return coef[4];
    return coef[i-1]*(1.0-lat/15.0+i)+coef[i]*(lat/15.0-i);
}

static double Mapf(double el,double a,double b,double c) {
    double sinel=sin(el);
    return (1.0+a/(1.0+b/(1.0+c)))/(sinel+(a/(sinel+b/(sinel+c))));
}

void TrpModel::TrpMapNeil(Time t, Vec3 pos, double el) {
#if 0
    double lat[5] ={15*D2R,30*D2R,45*D2R,60*D2R,75*D2R};
    double avgad[5] = {1.2769934e-3,1.2683230e-3,1.2465397e-3,1.2196049e-3,1.2045996e-3},
           avgbd[5] = {2.9153695e-3,2.9152299e-3,2.9288445e-3,2.9022565e-3,2.9024912e-3},
           avgcd[5] = {62.620505e-3,62.837393e-3,63.721774e-3,63.824265e-3,64.258455e-3};
    double ampad[5] = {0,1.2709626e-5,2.6523662e-5,3.4000452e-5,4.1202191e-5},
           ampbd[5] ={0,2.1414979e-5,3.0160779e-5,7.2562722e-5,11.723375e-5},
           ampcd[5] = {0,9.0128400e-5,4.3497037e-5,84.795348e-5,170.37206e-5};
    double avgaw[5] = {5.8021879e-4,5.6794847e-4,5.8118019e-4,5.9727542e-4,6.1641693e-4},
           avgbw[5] = {1.4275268e-3,1.5138625e-3,1.4572572e-3,1.5007428e-3,1.7599082e-3},
           avgcw[5] = {4.3472961e-2,4.6729510e-2,4.3908931e-2,4.4626982e-2,5.4736039e-2};

    double aht = 2.53e-5,bht = 5.49e-3,cht = 1.14e-3;
    double ad = 0,bd = 0,cd = 0;//Dry component abc
    double aw = 0,bw = 0,cw = 0;//Wet component abc
    double T0 = 28;
    double doy=t.Time2Doy();
    //Find location
    int flag = 0;
    for(int i = 0;i < 5;i++)
        if (pos.i_>lat[i])
            flag++;
        else
            break;
    //Calculate the dry and wet component abc
    if (pos.i_*R2D>15&&pos.i_*R2D<75)
    {
        ad=avgad[flag-1]+(avgad[flag]-avgad[flag-1])*(pos.i_-lat[flag-1])/(lat[flag]-lat[flag-1])
                +ampad[flag-1]+(ampad[flag]-ampad[flag-1])*(pos.i_-lat[flag-1])*cos(2*PI*(doy-T0)/365.25)/(lat[flag]-lat[flag-1]);
        bd=avgbd[flag-1]+(avgbd[flag]-avgbd[flag-1])*(pos.i_-lat[flag-1])/(lat[flag]-lat[flag-1])
                +ampbd[flag-1]+(ampbd[flag]-ampbd[flag-1])*(pos.i_-lat[flag-1])*cos(2*PI*(doy-T0)/365.25)/(lat[flag]-lat[flag-1]);
        cd=avgcd[flag-1]+(avgcd[flag]-avgcd[flag-1])*(pos.i_-lat[flag-1])/(lat[flag]-lat[flag-1])
                +ampcd[flag-1]+(ampcd[flag]-ampcd[flag-1])*(pos.i_-lat[flag-1])*cos(2*PI*(doy-T0)/365.25)/(lat[flag]-lat[flag-1]);
        aw=avgaw[flag-1]+(avgaw[flag]-avgaw[flag-1])*(pos.i_-lat[flag-1])/(lat[flag]-lat[flag-1]);
        bw=avgbw[flag-1]+(avgbw[flag]-avgbw[flag-1])*(pos.i_-lat[flag-1])/(lat[flag]-lat[flag-1]);
        cw=avgcw[flag-1]+(avgcw[flag]-avgcw[flag-1])*(pos.i_-lat[flag-1])/(lat[flag]-lat[flag-1]);
    }
    else if (pos.i_*R2D<=15)
    {
        ad = avgad[0]+avgad[0]*cos(2*PI*(doy-T0)/365.25);
        bd = avgbd[0]+avgbd[0]*cos(2*PI*(doy-T0)/365.25);
        cd = avgcd[0]+avgcd[0]*cos(2*PI*(doy-T0)/365.25);
        aw = avgaw[0]+avgaw[0]*(pos.i_-lat[flag-1])/(lat[flag]-lat[flag-1]);
        bw = avgbw[0]+avgbw[0]*(pos.i_-lat[flag-1])/(lat[flag]-lat[flag-1]);
        cw = avgcw[0]+avgcw[0]*(pos.i_-lat[flag-1])/(lat[flag]-lat[flag-1]);
    }
    else if (pos.i_*R2D>=75)
    {
        ad=avgad[4]+avgad[4]*cos(2*PI*(doy-T0)/365.25);
        bd=avgbd[4]+avgbd[4]*cos(2*PI*(doy-T0)/365.25);
        cd=avgcd[4]+avgcd[4]*cos(2*PI*(doy-T0)/365.25);
        aw=avgaw[4]+avgaw[4]*(pos.i_-lat[flag-1])/(lat[flag]-lat[flag-1]);
        bw=avgbw[4]+avgbw[4]*(pos.i_-lat[flag-1])/(lat[flag]-lat[flag-1]);
        cw=avgcw[4]+avgcw[4]*(pos.i_-lat[flag-1])/(lat[flag]-lat[flag-1]);
    }
    map_h_=(1+ad/(1+bd/(1+cd)))/(sin(el)+ad/(sin(el)+bd/(sin(el)+cd)))
            +(1/sin(el)-(1+aht/(1+bht/(1+cht)))/(sin(el)+(aht/(sin(el)+bht/(sin(el)+cht)))))*pos.k_/1000;
    map_w_=(1+aw/(1+bw/(1+cw)))/(sin(el)+aw/(sin(el)+bw/(sin(el)+cw)));
#else
    const double coef[][5]={
            { 1.2769934E-3, 1.2683230E-3, 1.2465397E-3, 1.2196049E-3, 1.2045996E-3},
            { 2.9153695E-3, 2.9152299E-3, 2.9288445E-3, 2.9022565E-3, 2.9024912E-3},
            { 62.610505E-3, 62.837393E-3, 63.721774E-3, 63.824265E-3, 64.258455E-3},

            { 0.0000000E-0, 1.2709626E-5, 2.6523662E-5, 3.4000452E-5, 4.1202191E-5},
            { 0.0000000E-0, 2.1414979E-5, 3.0160779E-5, 7.2562722E-5, 11.723375E-5},
            { 0.0000000E-0, 9.0128400E-5, 4.3497037E-5, 84.795348E-5, 170.37206E-5},

            { 5.8021897E-4, 5.6794847E-4, 5.8118019E-4, 5.9727542E-4, 6.1641693E-4},
            { 1.4275268E-3, 1.5138625E-3, 1.4572752E-3, 1.5007428E-3, 1.7599082E-3},
            { 4.3472961E-2, 4.6729510E-2, 4.3908931E-2, 4.4626982E-2, 5.4736038E-2}
    };
    const double aht[]={ 2.53E-5, 5.49E-3, 1.14E-3}; /* height correction */

    double y,cosy,ah[3],aw[3],dm,lat=pos.i_*R2D,hgt=pos.k_;
    int i;

    if(el<=0.0) {map_w_=0.0;return;}
    y=(t.Time2Doy()-28.0)/365.25+(lat<0.0?0.5:0.0);
    cosy=cos(2.0*PI*y);
    lat=fabs(lat);

    for(i=0;i<3;i++){
        ah[i]=Interpc(coef[i],lat)-Interpc(coef[i+3],lat)*cosy;
        aw[i]=Interpc(coef[i+6],lat);
    }
    dm=(1.0/sin(el)-Mapf(el,aht[0],aht[1],aht[2]))*hgt/1E3;
    map_w_=Mapf(el,aw[0],aw[1],aw[2]);
    map_h_=Mapf(el,ah[0],ah[1],ah[2])+dm;

#endif
}

#define ERR_SAAS    0.3         /* saastamoinen model error std (m) */
double TrpModel::SaasModel(Vec3 blh, double humi) {
    int sat=sat_data_->sat_.sat_no_;
#if 0
    double B;
    double hgt=blh.k_<0.0 ? 0.0 : blh.k_;
    if(hgt>15000.0) hgt=15000.0;

    B=1.1561-1.5915E-4*hgt+9.5788E-9*hgt*hgt-2.9393E-13*hgt*hgt*hgt;
    StandAtmoPara(blh,hgt,humi);
    double z=PI/2.0-sat_data_->azel_[1];
    z_hyd_=0.002277*(P_-B*tan(z)*tan(z));
    z_wet_=0.002277*(1255.0/T_+0.05)*e_;
    switch(kPrcOpt.GNSS_opt_.trp_map_opt_){
        case TRPMAP_EL:   TrpMapEl(z);break;
        case TRPMAP_GMF:  TrpMapGMF(sat_data_->sig_recep_,blh,sat_data_->azel_[1]);break;
        case TRPMAP_VMF1:   break;
        case TRPMAP_NEIL: TrpMapNeil(sat_data_->sig_recep_,blh,sat_data_->azel_[1]);break;
    }

    sat_data_->strp_h_=z_hyd_*map_h_;
    sat_data_->strp_w_=z_wet_*map_w_;
    sat_data_->trp_var_=SQR(ERR_SAAS);
#else
    const double temp0=15.0;
    double hgt,pres,temp,e,z;

    if(blh.k_<-100.0||1E6<blh.k_||sat_data_->azel_[1]<=0) return 0.0;
    if(blh.k_>=1.0/2.2557E-5) return 0.0;
    hgt=blh.k_<0.0?0.0:blh.k_;
    if(hgt>15000.0) hgt=15000.0;

    pres=1013.25*pow(1.0-2.2557E-5*hgt,5.2568);
    temp=temp0-6.5E-3*hgt+273.16;
    e=6.108*humi*exp((17.15*temp-4684.0)/(temp-38.45));
    z=PI/2.0-sat_data_->azel_[1];
    z_hyd_=0.0022768*pres/(1.0-0.00266*cos(2.0*blh.i_)-0.00028*hgt/1E3);
    z_wet_=0.0022770*(1255.0/temp+0.05)*e;
    isat_[sat-1].map_trp_h=map_h_=1.0/cos(z);
    isat_[sat-1].map_trp_w=map_w_=1.0/cos(z);
    isat_[sat-1].map_grad_e=isat_[sat-1].map_grad_n=0.0;
    isat_[sat-1].strp_h=z_hyd_*map_h_;
    isat_[sat-1].strp_w=z_wet_*map_w_;
    isat_[sat-1].strp=z_hyd_*map_h_+z_wet_*map_w_;
    isat_[sat-1].var_trp=SQR(ERR_SAAS);
#endif
}

double TrpModel::TrpCorr(const Vec3 pos, const double humi) {
    return 0.0;
}

SaasTrp::SaasTrp() {

}

SaasTrp::~SaasTrp() {

}

double SaasTrp::TrpCorr(const Vec3 pos, const double humi) {
    return SaasModel(pos,humi);
}

EstTrp::EstTrp() {
    
}

EstTrp::~EstTrp() {
    
}

double EstTrp::TrpCorr(Vec3 pos, double humi) {
    int sat=sat_data_->sat_.sat_no_;
    int it=para_.IdxTrpPar();
    // 经过saas模型后，已经获取z_hyd_和z_wet_,map_h_,map_w_
    SaasModel(pos,humi);

    //使用投影函数获取更高精度的map_h_和map_w_;
    TrpMapNeil(sat_data_->sig_recep_,pos,sat_data_->azel_[1]);

    if(kPrcOpt.GNSS_opt_.trp_opt_>TRP_EST_WET&&sat_data_->azel_[1]>0.0){
        double cotz=1.0/tan(sat_data_->azel_[1]);
        // map_w_是投影函数获取的初值
        double grad_n=map_w_*cotz*cos(sat_data_->azel_[0]);
        double grad_e=map_w_*cotz*sin(sat_data_->azel_[0]);
        map_w_+=grad_e*par_X_[it+1]+grad_e*par_X_[it+2];
        // par_X_[it]表示上个历元的对流层湿延迟估计值
        isat_[sat-1].map_grad_n=grad_n*par_X_[it];
        isat_[sat-1].map_grad_e=grad_e*par_X_[it];
    }
    isat_[sat-1].map_trp_w=map_w_;
    isat_[sat-1].map_trp_h=map_h_;
    isat_[sat-1].strp_h=z_hyd_*map_h_;
    isat_[sat-1].strp_w=par_X_[it]*map_w_;
    isat_[sat-1].strp=isat_[sat-1].strp_h+isat_[sat-1].strp_w;
    isat_[sat-1].var_trp=0.0001;
}



