/* 
==================================================
 Bot
 
 An original simulation of a non-traditional 
 biology-inspired neural network evolving in 
 a naturally selective environment to demonstrate
 the emergence of directed survival behavior.
 
Copyright (C) 11 November 2013 Souvik Das
ALL RIGHTS RESERVED
=================================================
*/

#pragma once

#include "Brain.h"
#include "Food.h"

#include "TLine.h"
#include<vector>
using namespace std;

class Bot: public Entity
{
  public:
  
    Brain *brain_;
    std::string name_;
    TEllipse *visRange1_, *visRange2_, *visRange3_;
    TLine *visPeriphery1_, *visPeriphery2_;
    double visualAngle_;
    int bodyColor_;
    double speed_;
    double worldSize_;
    float avgTimeToFood_;
    float avgTimeToFoodTemp_;
    int nEatenFoods_;
    int birthTime_;
    int penaltyTerm_;
    int eatenInGeneration_;
    int timeLastEaten_;
    int debug_;
    int kids_;
    
  public:
  
    Bot(std::string type, double x, double y, double theta, double visualAngle, int brainSize, int bodyColor, double speed, std::string name, double worldSize, int birthTime, int debug = 1);
    Bot(Bot *parentBot, double mu_newNeuron, double mu_newConnection, double mu_modConnection, double mu_visualAngle, int birthTime);
    typedef std::vector<vector<double> >Vec2d;
    Bot(Bot *parentBot, Vec2d *dist_matrix, int birthTime);
    ~Bot();
    
    void draw();
    void moveForward();
    void moveBackward();
    void seeEntity(Entity *entity);
    void seeFoods(std::vector<Food*> *foods);
    void seeBots(std::vector<Bot*> *bots);
    
    void stepInTime();
    void incrementTimeToFood(int time);
    void updateAvgTimeToEat();
    void updatePenaltyTerm(int time);
    void mutate(double mu_modConnection, double mu_visualAngle,int generations);
    
    void printBrain();
};  
