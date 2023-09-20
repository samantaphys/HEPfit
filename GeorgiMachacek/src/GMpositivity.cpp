/* 
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GMpositivity.h"
#include "GMcache.h"
#include "GeorgiMachacek.h"

GMpositivity1::GMpositivity1(const StandardModel& SM_i)
: ThObservable(SM_i), myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double GMpositivity1::computeThValue()
{
    flag_GM = myGM.getFlagGM();
    flag_EGM = myGM.getFlagEGM();
    flag_GGM = myGM.getFlagGGM();
    
    std::cout<< "flag-check4:" << flag_GM << flag_EGM << flag_GGM <<std::endl;
    
    if(flag_GM) {

    double lambda2 = myGM.getMyGMCache()->lambda2;
    double lambda3 = myGM.getMyGMCache()->lambda3;
//    std::cout<<"p1="<<lambda2+lambda3<<std::endl;
    return lambda2+lambda3;
    
    } else if(flag_EGM){
        throw std::runtime_error("WARNING: wrong observables for EGM model");
    } else
        throw std::runtime_error("WARNING: wrong observables for GGM model");
}



GMpositivity2::GMpositivity2(const StandardModel& SM_i)
: ThObservable(SM_i), myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double GMpositivity2::computeThValue()
{
    flag_GM = myGM.getFlagGM();
    flag_EGM = myGM.getFlagEGM();
    flag_GGM = myGM.getFlagGGM();
    
    if(flag_GM) {
         
    double lambda2 = myGM.getMyGMCache()->lambda2;
    double lambda3 = myGM.getMyGMCache()->lambda3;

//    std::cout<<"p2="<<2.0*lambda2+lambda3<<std::endl;
    return 2.0*lambda2+lambda3;
    
    } else if(flag_EGM){
        throw std::runtime_error("WARNING: wrong observables for EGM model");
    } else
        throw std::runtime_error("WARNING: wrong observables for GGM model");
}



GMpositivity3::GMpositivity3(const StandardModel& SM_i)
: ThObservable(SM_i), myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double GMpositivity3::computeThValue()
{
    flag_GM = myGM.getFlagGM();
    flag_EGM = myGM.getFlagEGM();
    flag_GGM = myGM.getFlagGGM();
    
    if(flag_GM) {
         
    double lambda1 = myGM.getMyGMCache()->lambda1;
    double lambda2 = myGM.getMyGMCache()->lambda2;
    double lambda3 = myGM.getMyGMCache()->lambda3;
    double lambda4 = myGM.getMyGMCache()->lambda4;
    double pos3 = -1.0;
    if(lambda1>0 && lambda2+lambda3>0)
    {
        pos3 = -fabs(lambda4)+2.0*sqrt(lambda1*(lambda2+lambda3));
    }
    else
    {
        pos3 = -fabs(lambda1*(lambda2+lambda3));
    }
//    std::cout<<"p3="<<pos3<<std::endl;
    return pos3;
    
    } else if(flag_EGM){
        throw std::runtime_error("WARNING: wrong observables for EGM model");
    } else
        throw std::runtime_error("WARNING: wrong observables for GGM model");
    
}



GMpositivity4::GMpositivity4(const StandardModel& SM_i)
: ThObservable(SM_i), myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double GMpositivity4::computeThValue()
{
    flag_GM = myGM.getFlagGM();
    flag_EGM = myGM.getFlagEGM();
    flag_GGM = myGM.getFlagGGM();
    
    if(flag_GM) {
  
    double lambda1 = myGM.getMyGMCache()->lambda1;
    double lambda2 = myGM.getMyGMCache()->lambda2;
    double lambda3 = myGM.getMyGMCache()->lambda3;
    double lambda4 = myGM.getMyGMCache()->lambda4;
    double lambda5 = myGM.getMyGMCache()->lambda5;
    double pos4 = -1.0;
    if(lambda1>0 && 2.0*lambda2+lambda3>0)
    {
        pos4 = lambda4-0.25*fabs(lambda5)+sqrt(2.0*lambda1*(2.0*lambda2+lambda3));
    }
    else
    {
        pos4 = -fabs(lambda1*(2.0*lambda2+lambda3));
    }
//    std::cout<<"p4="<<pos4<<std::endl;
    return pos4;
    
     } else if(flag_EGM){
        throw std::runtime_error("WARNING: wrong observables for EGM model");
    } else
        throw std::runtime_error("WARNING: wrong observables for GGM model");
}

/**
EGMpositivity1::EGMpositivity1(const StandardModel& SM_i)
: ThObservable(SM_i), myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double EGMpositivity1::computeThValue()
{
    flag_GM = myGM.getFlagGM();
    flag_EGM = myGM.getFlagEGM();
    flag_GGM = myGM.getFlagGGM();
    
    
    if(flag_GM) {
        throw std::runtime_error("WARNING: wrong observables for EGM model");
 
     } else if(flag_EGM){
         
    double lambda1 = myGM.getMyGMCache()->lambda1;
    double lambda2 = myGM.getMyGMCache()->lambda2;
    double lambda3 = myGM.getMyGMCache()->lambda3;
    double lambda4 = myGM.getMyGMCache()->lambda4;
    double lambda5 = myGM.getMyGMCache()->lambda5;
    double pos4 = -1.0;
    if(lambda1>0 && 2.0*lambda2+lambda3>0)
    {
        pos4 = lambda4-0.25*fabs(lambda5)+sqrt(2.0*lambda1*(2.0*lambda2+lambda3));
    }
    else
    {
        pos4 = -fabs(lambda1*(2.0*lambda2+lambda3));
    }
//    std::cout<<"p4="<<pos4<<std::endl;
    return pos4;
        
    } else
        throw std::runtime_error("WARNING: wrong observables for GGM model");
}

EGMpositivity2::EGMpositivity2(const StandardModel& SM_i)
: ThObservable(SM_i), myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double EGMpositivity2::computeThValue()
{
    flagGM = myGM.getflagGM();
    flagEGM = myGM.getflagEGM();
    flagGGM = myGM.getflagGGM();
    flagLO = myGM.getflagLO();
    flagNLO = myGM.getflagNLO();
    
    if(flagGM) {
        throw std::runtime_error("WARNING: wrong observables for EGM model");
 
     } else if(flagEGM){
         
    double lambda1 = myGM.getMyGMCache()->lambda1;
    double lambda2 = myGM.getMyGMCache()->lambda2;
    double lambda3 = myGM.getMyGMCache()->lambda3;
    double lambda4 = myGM.getMyGMCache()->lambda4;
    double lambda5 = myGM.getMyGMCache()->lambda5;
    double pos4 = -1.0;
    if(lambda1>0 && 2.0*lambda2+lambda3>0)
    {
        pos4 = lambda4-0.25*fabs(lambda5)+sqrt(2.0*lambda1*(2.0*lambda2+lambda3));
    }
    else
    {
        pos4 = -fabs(lambda1*(2.0*lambda2+lambda3));
    }
//    std::cout<<"p4="<<pos4<<std::endl;
    return pos4;
        
    } else
        throw std::runtime_error("WARNING: wrong observables for GGM model");
}

EGMpositivity3::EGMpositivity3(const StandardModel& SM_i)
: ThObservable(SM_i), myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double EGMpositivity3::computeThValue()
{
    flagGM = myGM.getflagGM();
    flagEGM = myGM.getflagEGM();
    flagGGM = myGM.getflagGGM();
    flagLO = myGM.getflagLO();
    flagNLO = myGM.getflagNLO();
    
    if(flagGM) {
        throw std::runtime_error("WARNING: wrong observables for EGM model");
 
     } else if(flagEGM){
         
    double lambda1 = myGM.getMyGMCache()->lambda1;
    double lambda2 = myGM.getMyGMCache()->lambda2;
    double lambda3 = myGM.getMyGMCache()->lambda3;
    double lambda4 = myGM.getMyGMCache()->lambda4;
    double lambda5 = myGM.getMyGMCache()->lambda5;
    double pos4 = -1.0;
    if(lambda1>0 && 2.0*lambda2+lambda3>0)
    {
        pos4 = lambda4-0.25*fabs(lambda5)+sqrt(2.0*lambda1*(2.0*lambda2+lambda3));
    }
    else
    {
        pos4 = -fabs(lambda1*(2.0*lambda2+lambda3));
    }
//    std::cout<<"p4="<<pos4<<std::endl;
    return pos4;
        
    } else
        throw std::runtime_error("WARNING: wrong observables for GGM model");
}

EGMpositivity4::EGMpositivity4(const StandardModel& SM_i)
: ThObservable(SM_i), myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double EGMpositivity4::computeThValue()
{
    flagGM = myGM.getflagGM();
    flagEGM = myGM.getflagEGM();
    flagGGM = myGM.getflagGGM();
    flagLO = myGM.getflagLO();
    flagNLO = myGM.getflagNLO();
    
    if(flagGM) {
        throw std::runtime_error("WARNING: wrong observables for EGM model");
 
     } else if(flagEGM){
         
    double lambda1 = myGM.getMyGMCache()->lambda1;
    double lambda2 = myGM.getMyGMCache()->lambda2;
    double lambda3 = myGM.getMyGMCache()->lambda3;
    double lambda4 = myGM.getMyGMCache()->lambda4;
    double lambda5 = myGM.getMyGMCache()->lambda5;
    double pos4 = -1.0;
    if(lambda1>0 && 2.0*lambda2+lambda3>0)
    {
        pos4 = lambda4-0.25*fabs(lambda5)+sqrt(2.0*lambda1*(2.0*lambda2+lambda3));
    }
    else
    {
        pos4 = -fabs(lambda1*(2.0*lambda2+lambda3));
    }
//    std::cout<<"p4="<<pos4<<std::endl;
    return pos4;
        
    } else
        throw std::runtime_error("WARNING: wrong observables for GGM model");
}

EGMpositivity5::EGMpositivity5(const StandardModel& SM_i)
: ThObservable(SM_i), myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double EGMpositivity5::computeThValue()
{
    flagGM = myGM.getflagGM();
    flagEGM = myGM.getflagEGM();
    flagGGM = myGM.getflagGGM();
    flagLO = myGM.getflagLO();
    flagNLO = myGM.getflagNLO();
    
    if(flagGM) {
        throw std::runtime_error("WARNING: wrong observables for EGM model");
 
     } else if(flagEGM){
         
    double lambda1 = myGM.getMyGMCache()->lambda1;
    double lambda2 = myGM.getMyGMCache()->lambda2;
    double lambda3 = myGM.getMyGMCache()->lambda3;
    double lambda4 = myGM.getMyGMCache()->lambda4;
    double lambda5 = myGM.getMyGMCache()->lambda5;
    double pos4 = -1.0;
    if(lambda1>0 && 2.0*lambda2+lambda3>0)
    {
        pos4 = lambda4-0.25*fabs(lambda5)+sqrt(2.0*lambda1*(2.0*lambda2+lambda3));
    }
    else
    {
        pos4 = -fabs(lambda1*(2.0*lambda2+lambda3));
    }
//    std::cout<<"p4="<<pos4<<std::endl;
    return pos4;
        
    } else
        throw std::runtime_error("WARNING: wrong observables for GGM model");
}
 */