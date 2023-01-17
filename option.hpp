//
//  option.hpp
//  Option
//
//  Created by SenChen on 2023/1/16.
//

#ifndef option_hpp
#define option_hpp
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <cmath>
#include <random>

using namespace std;



double normalCDF(double x) // Phi(-∞, x) aka N(x)
{
    return std::erfc(-x/std::sqrt(2))/2;
}
class Option{
protected:
    double Strike;
    double Sigma;
    double T;
    double RiskFreeRate;
    virtual double getCurrentExerciseValue(double spot, double t) = 0;
    virtual double calcBinomialNodeValue(double spot, double t, double time_val);

public:
    Option(double strike, double sigma, double t, double riskFreeRate) : Strike(strike), Sigma(sigma), T(t), RiskFreeRate(riskFreeRate){}
                                                 
    virtual double getValue(double spot) = 0;
    virtual double getIntrinsicValue(double spot) = 0;
    virtual double calcBlackScholesValue(double spot)=0;
    double getDelta(double spot);
    double calcBinomialTreeValue(double spot, int treeDepth);

};

class EuropeanOption : public Option{
public:
    EuropeanOption(double strike, double sigma, double t, double riskFreeRate) : Option(strike, sigma, t, riskFreeRate){}
    virtual double getValue(double spot);
    double getCurrentExerciseValue(double spot, double t);
};

class EuropeanCall : public EuropeanOption{
public:
    EuropeanCall(double strike, double sigma, double t, double riskFreeRate) : EuropeanOption(strike, sigma, t, riskFreeRate){}
    virtual double getIntrinsicValue(double spot) ;
    virtual double calcBlackScholesValue(double spot);
};
class EuropeanCallKnockout : public EuropeanCall{
protected:
    double Barrier;
public:
    EuropeanCallKnockout(double strike,double barrier,double sigma,double expiry,double riskFreeRate):EuropeanCall(strike,sigma,expiry,riskFreeRate),Barrier(barrier){}
    double getIntrinsicValue(double spot);
    double getValue(double spot);
    double calcBinomialNodeValue(double spot, double t, double time_val);
};

class EuropeanPut : public EuropeanOption{
public:
    EuropeanPut(double strike, double sigma, double t, double riskFreeRate) : EuropeanOption(strike, sigma, t, riskFreeRate){}
    virtual double getIntrinsicValue(double spot);
    virtual double calcBlackScholesValue(double spot);
};

class AmericanOption : public Option{
public:
    AmericanOption(double strike, double sigma, double t, double riskFreeRate) : Option(strike, sigma, t, riskFreeRate){}
    virtual double getValue(double spot);
    double getCurrentExerciseValue(double spot, double t);
};

class AmericanCall : public AmericanOption{
public:
    AmericanCall(double strike, double sigma, double t, double riskFreeRate) : AmericanOption(strike, sigma, t, riskFreeRate){}
    virtual double getIntrinsicValue(double spot) ;
    virtual double calcBlackScholesValue(double spot);
};

class AmericanPut : public AmericanOption{
public:
    AmericanPut(double strike, double sigma, double t, double riskFreeRate) : AmericanOption(strike, sigma, t, riskFreeRate){}
    virtual double getIntrinsicValue(double spot);
    virtual double calcBlackScholesValue(double spot);
};
class Position {
public:
    double weight;
    Option* option;
    Position() : weight(0), option(NULL) {}
    Position(double w, Option* o) : weight(w), option(o) {}
    double getweight()const;

    Option*getoption()const;
   
};

#endif /* option_hpp */
