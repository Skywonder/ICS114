//
//  main.cpp
//  ICS114-P2-Part1
//
//  Created by Kuan-Ping Chang  on 4/21/17.
//  Copyright © 2017 Kuan-Ping Chang . All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <random>

using namespace std;


double PI = 3.14159265359;
//part 1-1
vector<pair<double, double> > sm_d;
//part 1-2
vector<pair<double, double> > sm_d_w1;
vector<pair<double, double> > sm_d_w2;
vector<pair<double, double> > sm_d_w3;
//part 2-1
vector<pair<double, double> > sm_d_w4;
//part 2-2
vector<pair<double, double> > sm_d_w5;
//part 3
vector<pair<double, double> > sm_d_w6;

//1-1
void generator_1()
{
    float total = 0.0;
    float p = 0.25;
    const int N = 100000;
    vector<double> random_x {};
    vector<double> i_hat_list {};
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(-2,2);
    
    //get N number of x
    for (int i = 0; i < N; i++){
        double r_x = dist(gen);
        random_x.push_back(r_x);
    }
    //find sum of i_hat
    for (int i =0; i < N; i++)
    {
        double x_2 = pow(random_x[i], 2);
        double e = exp(-x_2/2);
        double i_hat = (e/p);
        i_hat_list.push_back(i_hat);
        total = total + i_hat;
     
    }
    double sample_mean = total/N;
    //std::cout << "Sample Mean: " << sample_mean << std::endl;
    
    float s_total = 0.0;
    for(int i = 0; i < N; i++)
    {
        s_total = s_total + pow((i_hat_list[i] - sample_mean),2);
    }
    double stad = sqrt(s_total/(N-1));
    double deviation = (2 * stad) / sqrt(N);
    //std::cout << "Deviation: " << deviation << std::endl;
    
    sm_d.push_back(make_pair(sample_mean, deviation));
    
}

//1-2
void generator_2(){
    
    const int N = 100000;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0,1);
    vector<float> freq = {0.1, 1.0, 10.0};
    int count = 0;
    for (int j = 0; j < freq.size(); j++)
    {
        vector<double> random_x {}; //reset random_x
        vector<double> i_hat_list {};//reset i_hat
        float total = 0.0; //reset total
        //get N number of x
        for (int i = 0; i < N; i++){
            double r_x = dist(gen);
            random_x.push_back(1 - (log(r_x)/freq[j]));
        }
        
        //find sum of i_hat
        for (int i =0; i < N; i++)
        {
            double x_2 = pow(random_x[i], 2);
            double e = exp(-x_2/2);
            double p = freq[j] * exp(-freq[j]*(random_x[i] -1));
            double i_hat = (e/p);
            i_hat_list.push_back(i_hat);
            total = total + i_hat;
            
        }
        double sample_mean = total/N;
        //std::cout << "Sample Mean: " << sample_mean << std::endl;
    
        float s_total = 0.0;
        for(int i = 0; i < N; i++)
        {
            s_total = s_total + pow((i_hat_list[i] - sample_mean),2);
        }
        double stad = sqrt(s_total/(N-1));
        double deviation = (2 * stad) / sqrt(N);
        //std::cout << "Deviation: " << deviation << std::endl;
    
        if(count == 0)
            sm_d_w1.push_back(make_pair(sample_mean, deviation));
        else if(count == 1)
            sm_d_w2.push_back(make_pair(sample_mean, deviation));
        else if(count == 2)
            sm_d_w3.push_back(make_pair(sample_mean, deviation));
        count++;
    }
}

//2-1
void Integrals_UnitSphere(){
    const int N = 100000;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> rnd(0,1);
    double total = 0.0;
    vector<double> i_hat_list;
    vector<double> x_1;
    //get all the A
    for (int i = 0; i < N; i++)
    {
        x_1.push_back((PI/2)*rnd(gen));
    }
    //find sum of i_hat
    for (int i =0; i < N; i++)
    {
        double p = 1/pow(PI,2);
        double e = sin(x_1[i]);
        double i_hat = (e/p);
        i_hat_list.push_back(i_hat);
        total = total + i_hat;
        
    }
    double sample_mean = total/N;
    //std::cout << "Sample Mean: " << sample_mean << std::endl;
    
    float s_total = 0.0;
    for(int i = 0; i < N; i++)
    {
        s_total = s_total + pow((i_hat_list[i] - sample_mean),2);
    }
    double stad = sqrt(s_total/(N-1));
    double deviation = (2 * stad) / sqrt(N);
    //std::cout << "Deviation: " << deviation << std::endl;
    
    sm_d_w4.push_back(make_pair(sample_mean, deviation));
    
}
//2-2
void ArbitrarySphereical()
{
    const int N = 100000;
    std::random_device rd;
    std::random_device rd2;
    std::mt19937 gen(rd());
    std::mt19937 gen2(rd2());
    std::uniform_real_distribution<double> rnd(0,1);
    double total = 0.0;
    vector<double> random_theta;
    vector<double> random_phi;
    vector<double> i_hat_list;
    //get all the A
    for (int i = 0; i < N; i++)
    {
        random_theta.push_back(acos(rnd(gen)));
        random_phi.push_back(2*PI*rnd(gen2));
    }
    //find sum of i_hat
    for (int i =0; i < N; i++)
    {
        double p = 1/(2*PI);
        double e = pow(cos(random_theta[i])*sin(random_theta[i])*cos(random_phi[i]),2) * sin(PI/2);
        double i_hat = (e/p);
        i_hat_list.push_back(i_hat);
        total = total + i_hat;
        
    }
    double sample_mean = total/N;
    //std::cout << "Sample Mean: " << sample_mean << std::endl;
    
    float s_total = 0.0;
    for(int i = 0; i < N; i++)
    {
        s_total = s_total + pow((i_hat_list[i] - sample_mean),2);
    }
    double stad = sqrt(s_total/(N-1));
    double deviation = (2 * stad) / sqrt(N);
    //std::cout << "Deviation: " << deviation << std::endl;
    
    sm_d_w5.push_back(make_pair(sample_mean, deviation));

}

//Estimating Radiant Flux
void radiantflux(){
    const int N = 100000;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::mt19937 gen2(rd());
    std::mt19937 gen3(rd());
    std::mt19937 gen4(rd());
    std::uniform_real_distribution<double> rng(0,1);
    std::uniform_real_distribution<double> rng2(-0.5, 0.5);

    vector<pair<double, double> > X;
    vector<double> discriminant;
    vector<double> i_hat_list;
    double total = 0.0;
    vector<double> l_theta;
    for (int i =0; i< N; i++)
    {
        //w - x = direction = ray
        double theta = acos(rng(gen)); //theta
        double phi = 2*PI*rng(gen2); //phi
        l_theta.push_back(theta);
        double wx = theta;
        double wy = phi;
        double dx = rng2(gen3);
        double dy = rng2(gen4);
        double sphere[3] = {1.0,1.0,5.0};
        double radius = 1.0;
        double translated[3];
        double l = sphere[2]/cos(wx);
        double l2 = l * sin(wx);
        double l3 = l2 * sin(wy);
        double l4 = l2 * cos(wy);
        translated[0] = 1.0 - l3;
        translated[1] = 1.0 - l4;
        translated[2] = 0.0;
        double check = pow(pow(dx-translated[0], 2.0) + pow(dy - translated[1], 2.0), 0.5);
        if (check <= pow(radius, 2.0))
            discriminant.push_back(0);
        else
            discriminant.push_back(1);
    }
    
    //corss product of 2 vector
    
    //find sum of i_hat
    for (int i =0; i < N; i++)
    {
        double p = 1.0/(2.0*PI);
        double e = 100 * discriminant[i] * cos(l_theta[i])*sin(PI/2);
        double i_hat = (e/p);
        i_hat_list.push_back(i_hat);
        total = total + i_hat;
        
    }
    double sample_mean = total/N;
    //std::cout << "Sample Mean: " << sample_mean << std::endl;
    
    float s_total = 0.0;
    for(int i = 0; i < N; i++)
    {
        s_total = s_total + pow((i_hat_list[i] - sample_mean),2);
    }
    double stad = sqrt(s_total/(N-1));
    double deviation = (2 * stad) / sqrt(N);
    //std::cout << "Deviation: " << deviation << std::endl;
    
    sm_d_w6.push_back(make_pair(sample_mean, deviation));

}

double cross_product(double a1, double a2, double a3, double b1, double b2, double b3)
{
    return sqrt(pow(a2*b3 - a3*b2,2)+pow(a3*b1 - a1*b3,2)+ pow(a1*b2 - a2*b1,2));

}


int main(int argc, const char * argv[]) {

    for (int i = 0; i < 10; i++)
    {
        generator_1();
        generator_2();
        Integrals_UnitSphere();
        ArbitrarySphereical();
        radiantflux();
    }
    for (int i = 0; i < sm_d.size(); i++)
    {
        cout << "Sample Mean: " << sm_d[i].first << " Deviation: " << sm_d[i].second << endl;
    }

    //cout << sm_d_w1.size() << endl;
    
    cout << "For freq 0.1" << endl;
    for(int i = 0; i < sm_d_w1.size(); i++)
    {
        cout << "Sample Mean: " << sm_d_w1[i].first << " Deviation: " << sm_d_w1[i].second << endl;
    }
    cout << endl;
    cout << "For freq 1.0" << endl;
    for(int i = 0; i < sm_d_w2.size(); i++)
    {
        cout << "Sample Mean: " << sm_d_w2[i].first << " Deviation: " << sm_d_w2[i].second << endl;
    }
    cout << endl;
    cout << "For freq 10.0" << endl;
    for(int i = 0; i < sm_d_w3.size(); i++)
    {
        cout << "Sample Mean: " << sm_d_w3[i].first << " Deviation: " << sm_d_w3[i].second << endl;
    }
    cout << endl;
    cout << "For Unit Sphere" << endl;
    for(int i = 0; i < sm_d_w4.size(); i++)
    {
        cout << "Sample Mean: " << sm_d_w4[i].first << " Deviation: " << sm_d_w4[i].second << endl;
    }
    cout << endl;
    cout << "For Arbitary Sphere" << endl;
    for(int i = 0; i < sm_d_w5.size(); i++)
    {
        cout << "Sample Mean: " << sm_d_w5[i].first << " Deviation: " << sm_d_w5[i].second << endl;
    }
    cout << endl;
    cout << "Radian Flux Estimate" << endl;
    for(int i = 0; i < sm_d_w6.size(); i++)
    {
        cout << "Sample Mean: " << sm_d_w6[i].first << " Deviation: " << sm_d_w6[i].second << endl;
    }
    cout << endl;
}
