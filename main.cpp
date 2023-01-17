//
//  main.cpp
//  Option
//
//  Created by SenChen on 2022/3/16.
//
#include <cmath>
#include "option.hpp"
#include <numeric>
#include <iostream>
#include <random>

std::mt19937 mt_rand(12345);
std::normal_distribution<double> dis_normal(0.0, 1.0);

using namespace std;\

bool myfunction (const pair<int, double>&i,const pair<int, double>&j){
    return (i.second>j.second); }


bool myfunction2 (int i,int j) { return (i>j); }

double Option::calcBinomialNodeValue(double spot, double t, double time_val){
    if (time_val < getCurrentExerciseValue(spot, t)){
        return getCurrentExerciseValue(spot, t);
    }
    else{
        return time_val;
    }
}

double Option::getDelta(double spot){
    return getValue(spot*1.01) - getValue(spot);
}

double Option::calcBinomialTreeValue(double spot, int treeDepth){
    double dt = T/treeDepth;
    double u = exp(Sigma * sqrt(dt));
    double d = 1/u;
    double p = (exp(RiskFreeRate*dt) - d)/(u-d);
    vector<double> vals;
    for (int i = 0; i<=treeDepth; i++){
        double s = spot * pow(u, (2*i-treeDepth));
        vals. push_back (getCurrentExerciseValue(s, T));
    }
    double t =T;
    for (int j = treeDepth-1; j >=0; j--){
        t = t - dt;
        for (int i = 0; i <=j; i++){
            double s = spot * pow(u, (2*i-j));
            double v_time = exp(-RiskFreeRate*dt)*(p*vals[i+1] + (1-p)*vals[i]);
            double v = calcBinomialNodeValue(s, t, v_time);
            vals[i] = v;
        }
        }
            return vals[0];
}

//EuroOption
double EuropeanOption::getValue(double spot){
    return calcBlackScholesValue(spot);
}
double EuropeanOption::getCurrentExerciseValue(double spot, double t){
    if (t==T) {
        return getIntrinsicValue(spot);
    }
    else{
        return 0;
    }
}

///EuroCall
double EuropeanCall::getIntrinsicValue(double spot) {
    return max(spot - Strike,0.0);
}
double EuropeanCall::calcBlackScholesValue(double spot){
    double d1 = ((1 / (sqrt(T) * Sigma))*(log( spot / Strike )+(RiskFreeRate + Sigma * Sigma /2) * T));
    double d2 = d1 - Sigma * sqrt(T);
    return normalCDF(d1) * spot - normalCDF(d2) * Strike * exp(-RiskFreeRate * T);
}
double EuropeanCallKnockout::getIntrinsicValue(double spot){
    if (spot < Barrier) {
        return max(spot - Strike, 0.0);
    }
    else{
        return 0;
    }
}
double EuropeanCallKnockout::getValue(double spot){
    return calcBinomialTreeValue (spot, 500);
}
double EuropeanCallKnockout::calcBinomialNodeValue(double spot, double t, double time_val){
    if (spot >= Barrier){
        return 0;
    }
    else{
        return Option::calcBinomialNodeValue(spot, t, time_val);
    }
}

//EuroPut
double EuropeanPut::getIntrinsicValue(double spot) {
    return max(Strike - spot,0.0);
}
double EuropeanPut::calcBlackScholesValue(double spot){
    double d1 = (log(spot/Strike) + (RiskFreeRate + 0.5 * Sigma * Sigma) * T)/(Sigma * sqrt(T));
    double d2 = d1 - (Sigma * sqrt(T));
    return ( (1-normalCDF(d2))*Strike*exp(-RiskFreeRate*T) - (1 - normalCDF(d1)) *spot);
}

//AmericanOption
double AmericanOption::getValue(double spot){
    return calcBinomialTreeValue (spot, 500);
}
double AmericanOption::getCurrentExerciseValue(double spot, double t){
    return getIntrinsicValue(spot);
}
//AmericanCall
double AmericanCall::getIntrinsicValue(double spot) {
    return max(spot - Strike,0.0);
}
double AmericanCall::calcBlackScholesValue(double spot){
    double d1 = ((1 / (sqrt(T) * Sigma))*(log( spot / Strike )+(RiskFreeRate + Sigma * Sigma /2) * T));
    double d2 = d1 - Sigma * sqrt(T);
    return normalCDF(d1) * spot - normalCDF(d2) * Strike * exp(-RiskFreeRate * T);
}
//AmericanPut
double AmericanPut::getIntrinsicValue(double spot) {
    return max(Strike - spot, 0.0);
}
double AmericanPut::calcBlackScholesValue(double spot){
    double d1 = (log(spot/Strike) + (RiskFreeRate + 0.5 * Sigma * Sigma) * T)/(Sigma * sqrt(T));
    double d2 = d1 - (Sigma * sqrt(T));
    return ( (1-normalCDF(d2))*Strike*exp(-RiskFreeRate*T) - (1 - normalCDF(d1)) *spot);
}

//Position
double Position::getweight()const
{
    return weight;
}
Option*Position::getoption()const
{
    return option;
}


void print_values(double spot, Option& opt, string desc) {
    cout << desc << opt.getValue(spot) << " binomial_tree(depth=100): " << opt.calcBinomialTreeValue(spot, 100) << " Black Scholes: " << opt.calcBlackScholesValue(spot) << endl;

}

int main() {



    double spot = 100;

    double r = 0.05;
    double t = 1;
    double sigma = 0.25;

    EuropeanCall ec1(spot, sigma, t, r);
    AmericanCall ac1(spot, sigma, t, r);
    EuropeanCallKnockout ecko1(spot, 125, sigma, t, r);
    EuropeanPut ep1(spot, sigma, t, r);
    AmericanPut ap1(spot, sigma, t, r);
 


    print_values(spot, ec1, "EURO CALL value: ");
    print_values(spot, ac1, "AMERICAN CALL value: ");
    print_values(spot, ecko1, "KNOCKOUT CALL value: ");
    print_values(spot, ep1, "EURO PUT value: ");
    print_values(spot, ap1, "AMERICAN PUT value: ");



    return 0;

}
