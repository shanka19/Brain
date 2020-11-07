/* 
==================================================
 BrainInWorld
 
 An original simulation of a non-traditional 
 biology-inspired neural network evolving in 
 a naturally selective environment to demonstrate
 the emergence of directed survival behavior.
 
Copyright (C) 11 November 2013 Souvik Das
ALL RIGHTS RESERVED
=================================================
*/

#include <vector>
#include <TApplication.h>
#include <TROOT.h>
#include <TStyle.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <string.h>
#include <TFile.h>

#include <TCanvas.h>
#include <TEllipse.h>
#include <TLine.h>
#include <TText.h>
#include <TRandom3.h>
#include <TTimer.h>

#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>

#include "interface/Bot.h"
#include "interface/Neuron.h"
#include "interface/Brain.h"
#include "interface/Fire.h"
#include "interface/Food.h"
#include "interface/ToolBox.h"
#include "interface/CommandLineArguments.h"
#include "interface/Crossover.h"
#include "interface/Sort.h"

using namespace std;
int skipGenerations=0;
int endGeneration=100;
int timeStep=200;
double worldSize=100;
double regenFood=1.0;
int seed=100;

unsigned int nFoods=5;
unsigned int nBots=10;
// unsigned int nPredators=5;

// Mutation parameters
double mu_newNeuron=0; // 0.001;
double mu_newConnection=0.05;
double mu_modConnection=0.05;
double mu_visualAngle=0.05;

// Debug Levels
// bits: xxxx
// bit 0 = TCanvas visualization
// bit 1 = Verbalization
// bit 2 = Fill histograms
// bit 3 = Draw the histograms

int debug = 0x2;

// Reproduction Types
// 0 = Only Mutation
// 1 = Only Crossover
// 2 = Mutation and Crossover

int reproductionType;

// Reproduction Types
// 0 = Only Mutation
// 1 = Only Crossover
// 2 = Mutation and Crossover

int crossoverType;
// 1 = Crossover Type 1 
// 2 = Crossover Type 2
// 3 = Crossover Type 3
// 4 = Crossover Type 4

int main(int argc, char *argv[])
{
  // Get command line arguments
  std::map<std::string, float> cmdMap=commandLineArguments(argc, argv);
  if (cmdMap.find("-debug")!=cmdMap.end())           debug=(unsigned int)cmdMap["-debug"];
  if (cmdMap.find("-skipGenerations")!=cmdMap.end()) skipGenerations=(unsigned int)cmdMap["-skipGenerations"];
  if (cmdMap.find("-endGeneration")!=cmdMap.end())   endGeneration=(unsigned int)cmdMap["-endGeneration"];
  if (cmdMap.find("-timeStep")!=cmdMap.end())        timeStep=(unsigned int)cmdMap["-timeStep"];
  if (cmdMap.find("-worldSize")!=cmdMap.end())       worldSize=(unsigned int)cmdMap["-worldSize"];
  if (cmdMap.find("-nBots")!=cmdMap.end())           nBots=(unsigned int)cmdMap["-nBots"];
  if (cmdMap.find("-nFoods")!=cmdMap.end())          nFoods=(unsigned int)cmdMap["-nFoods"];
  if (cmdMap.find("-seed")!=cmdMap.end())            seed=(unsigned int)cmdMap["-seed"];
  if (cmdMap.find("-mu_modConnection")!=cmdMap.end()) mu_modConnection=cmdMap["-mu_modConnection"];
  if (cmdMap.find("-mu_visualAngle")!=cmdMap.end())   mu_visualAngle=cmdMap["-mu_visualAngle"];
  if (cmdMap.find("-reproductionType")!=cmdMap.end()) reproductionType=cmdMap["-reproductionType"];
  if (cmdMap.find("-crossoverType")!=cmdMap.end())    crossoverType=cmdMap["-crossoverType"];

  r3->SetSeed(seed);
  
  std::cout<<"debug = "<<debug<<std::endl;
  std::cout<<"visualization = "<<decodeDebug(debug, 0)<<std::endl;
  std::cout<<"mu_modConnection = "<< mu_modConnection << std::endl;
  std::cout<<"mu_visualAngle = "<< mu_visualAngle << std::endl;
  std::cout<<"reproductionType = "<< reproductionType << std::endl;
  std::cout<<"crossoverType = "<< crossoverType << std::endl;
  std::cout<<"seed = "<< seed << std::endl;
  TApplication *myapp=new TApplication("myapp",0,0);
  gStyle->SetCanvasPreferGL(true);
  gStyle->SetPalette(1);
  
  typedef std::vector<Bot*> Bots;
  Bots bots;
  for (unsigned int i=0; i<nBots; ++i)
  {
    Bot *bot=new Bot("Bot", r3->Rndm()*worldSize, r3->Rndm()*worldSize, r3->Rndm()*2.*pi, pi/4., 30, kBlue, 1.0, "Bot_"+itoa(i), worldSize, 0, debug);
    //cout<< "Initial Bot address : " << bot << endl;
    bots.push_back(bot);
  }
  std::cout<<"Instantiated bots."<<std::endl;
  
  typedef std::vector<Food*> Foods;
  Foods foods;
  for (unsigned int i=0; i<nFoods; ++i)
  {
    Food *food=new Food(r3->Rndm()*worldSize, r3->Rndm()*worldSize, r3->Rndm()*2.*pi, worldSize);
    foods.push_back(food);
  }
  std::cout<<"Instantiated food."<<std::endl;
  
  std::vector <double> time_vector;
  std::vector <double> avgBrainSize_vector;
  std::vector <double> generation_vector;
  std::vector <double> dtime_vector;
  
  TCanvas *c_World;
  TText *text=new TText(0.01, 0.01, "Generation 0");
  text=new TText(0.01, 0.01, "Generation 0");
  text->SetNDC();
  text->SetTextFont(42);
  if (decodeDebug(debug, 0)==1)
  {
    c_World=new TCanvas("c_World", "Natural Neural Network in Genetic Algorithm", 500, 500);
    // Safety Circle
    // TEllipse *e_safe=new TEllipse(worldSize/2., worldSize/2., 70, 70);
    // e_safe->Draw();
    c_World->Range(0,0,worldSize,worldSize);
  }
  
  TCanvas *c_Potential_Histograms;
  TCanvas *c_SynapticStrength_Histograms;
  TCanvas *c_Distance_Histograms;
  if (decodeDebug(debug, 3)==1)
  {
    c_Potential_Histograms=new TCanvas("c_Potential_Histograms", "Brain Data - Neural Potentials", 700, 700);
    // c_SynapticStrength_Histograms=new TCanvas("c_SynapticStrength_Histograms", "Brain Data - Synaptic Strengths", 700, 700);
    c_Distance_Histograms=new TCanvas("c_Distance_Histograms", "Brain Data - Neural Distances", 700, 700);
    c_Potential_Histograms->Divide(ceil(bots.size()/3.), 3);
    c_SynapticStrength_Histograms->Divide(ceil(bots.size()/3.), 3);
    c_Distance_Histograms->Divide(ceil(bots.size()/3.), 3);
  }
    
  int time=0;
  int generations=0;
  //int generations_predator=0;
  int oldGeneration=generations;
  int dtime=0;
  
  typedef vector<int> Vec1d;
  typedef vector<vector<double> > Vec2d;
   // Create a user defined datatype 
  typedef vector<Vec2d> Vec3d; // 3d Vector to store the two parent matrices 
  
  // Create and initialize Parent and Child Matrices to zeroes 
  Vec2d parent1;
  Vec2d parent2;
  Vec2d child1;
  Vec2d child2;
  
  Vec3d storage; // 3d vector to store the 2 parent matrices (vector of vectors)
  
  Vec2d vv_tempDistanceMatrix; // Vector of vector to store the distance matrix 
  int matrix_size = bots.at(0)->brain_->neurons_.size();
  std::vector<double> v_zeroVector(matrix_size,0);
  
  // For Loop to initialize matrices to all zeroes 
  for(unsigned int i = 0;i < matrix_size ; i++)
  {
     parent1.push_back(v_zeroVector);
     parent2.push_back(v_zeroVector);
     child1.push_back(v_zeroVector);
     child2.push_back(v_zeroVector);
     vv_tempDistanceMatrix.push_back(v_zeroVector);
  }
  
  Crossover *crossover_ = new Crossover();
  
  Sort *obj = new Sort();
  std::vector<int> parent_index;
  std::vector<int> avgTimeToFoods(bots.size(),0); // Vector with same size as bots to keep a record of average time taken by each bot to eat the foods
  std::vector<int> avgTimeToFoodsTemp(bots.size(),0); // Vector with same size as bots to keep a record of average time taken by each bot to eat the foods
  std::vector<int> delIndex; // vector to store indexes to be deleted
  int counter = 0; // Variable to keep a counter of how many bots have eaten food 
  bool start_mate = false; // Temp Variable to come out of the loop 
  //int dtime_predator=0;
  // Time loop
  while (foods.size()>0 && generations<endGeneration)
  {
    
    ++time;
    ++dtime;
    for (unsigned int i=0; i<bots.size(); ++i)
    {
      bots.at(i)->seeFoods(&foods);
      bots.at(i)->stepInTime();
    }
    
    for (unsigned int i=0; i<foods.size(); ++i)
    {
      foods.at(i)->moveForward();
    } 
    
    // check for bots eating food
    int nEaten=0;
    int nBots=bots.size();
    for (unsigned int i=0; i<nBots; ++i)
    {
      int eatenFood=-1;
      
      for (unsigned int j=0; j<foods.size(); ++j)
      { 
        double d2=pow(bots.at(i)->x_-foods.at(j)->x_, 2)+pow(bots.at(i)->y_-foods.at(j)->y_, 2);
        if (d2<13)
        {
          eatenFood=j;
	  ++counter;
	  bots.at(i)->incrementTimeToFood(time);
	  bots.at(i)->eatenInGeneration_ = 1;
	  //bots.at(i)->timeLastEaten_ = time;
	  parent_index.push_back(i);
	  for(unsigned int k=0; k<bots.at(i)->brain_->neurons_.size();k++)
	  {
	    for(unsigned int l=0; l<bots.at(i)->brain_->neurons_.at(k)->neuralRelations_.size();l++)
	    {
	      int index = bots.at(i)->brain_->neurons_.at(k)->neuralRelations_.at(l)->index;
              double distance_temp = bots.at(i)->brain_->neurons_.at(k)->neuralRelations_.at(l)->distance;
	      vv_tempDistanceMatrix[k][index] = distance_temp; // Distance Matrix
            }
	   }
	   storage.push_back(vv_tempDistanceMatrix); // Store parent in 3D vector
	  
	   vv_tempDistanceMatrix.clear(); // Clear and then reinitialize it with zeros for the next parent
	   for(unsigned int i = 0;i < matrix_size ; i++) vv_tempDistanceMatrix.push_back(v_zeroVector);
           if (counter==2)
           {
             start_mate = true;
             counter = 0;
           }
          if (decodeDebug(debug, 1)==1) std::cout<<"Bot "<<bots.at(i)->name_<<" ate food "<<j<<std::endl;
        }
      }
      
      if (eatenFood!=-1) 
      {
        ++nEaten;
        
        if (decodeDebug(debug, 1)==1) std::cout<<"foods.size() = "<<foods.size()<<" and eatenFood = "<<eatenFood<<std::endl;
        delete *(foods.begin()+eatenFood);
        foods.erase(foods.begin()+eatenFood);
        
        if (r3->Rndm()<regenFood)
        {
          Food *food=new Food(r3->Rndm()*worldSize, r3->Rndm()*worldSize, r3->Rndm()*(2.*pi), worldSize);
          foods.push_back(food);
        }       

        if (decodeDebug(debug, 1)==1) std::cout<<"removed from vector, foods.size() = "<<foods.size()<<std::endl;
      }
    }
     
    // Mating of two bots that ate food begins here
    if (start_mate)
    {
      if(reproductionType == 0)
      {
        //cout << "Type-0" << endl;
	int p1_index = parent_index[0];                            
        int p2_index = parent_index[1];                           
        parent_index.clear();                                   
                                                                
        for(unsigned int k = 0;k<matrix_size;k++)                  
        {                                                         
          for(unsigned int l = 0;l<matrix_size;l++)               
          {                                                       
            //cout << "Before " << endl;                           
            parent1[k][l] = storage[0][k][l];                     
            parent2[k][l] = storage[1][k][l];                      
            //cout << "After" << endl;                            
          }
        } 
	
	Bot *bot1 = new Bot(bots.at(p1_index),&parent1,time);
        Bot *bot2 = new Bot(bots.at(p2_index),&parent2,time);
	
        bot1->mutate(mu_modConnection,mu_visualAngle,generations);
        bot2->mutate(mu_modConnection,mu_visualAngle,generations);
	
	bots.push_back(bot1);
        bots.push_back(bot2);
      }
      
      else if(reproductionType == 1)
      {
        //cout << "Type-1" << endl;                      
        int p1_index = parent_index[0];                            
        int p2_index = parent_index[1];                         
        parent_index.clear();                                   
                                                                
        for(unsigned int k = 0;k<matrix_size;k++)                  
        {                                                         
          for(unsigned int l = 0;l<matrix_size;l++)               
          {                                                       
            //cout << "Before " << endl;                           
            parent1[k][l] = storage[0][k][l];                     
            parent2[k][l] = storage[1][k][l];                      
            //cout << "After" << endl;                            
          }
        }   
        if(crossoverType == 1)
	{
	  crossover_->crossover_type1(&parent1,&parent2,&child1,&child2);
	}
	else if(crossoverType == 2)
	{ 
	  crossover_->crossover_type2(&parent1,&parent2,&child1,&child2);
	}
	else if(crossoverType == 3)
	{
	  crossover_->crossover_type3(&parent1,&parent2,&child1,&child2);
	}
	else
	{
	  crossover_->crossover_type4(&parent1,&parent2,&child1,&child2);
	}
        Bot *bot1 = new Bot(bots.at(p1_index),&child1,time);
        Bot *bot2 = new Bot(bots.at(p2_index),&child2,time);
	
	bots.push_back(bot1);
        bots.push_back(bot2);
      }
      else
      {
        //cout << "Type-2" << endl; 
        //cout<< "Creating Children" << endl;                      
        int p1_index = parent_index[0];                            
        int p2_index = parent_index[1];                          
        parent_index.clear();                                   
                                                                
        for(unsigned int k = 0;k<matrix_size;k++)                  
        {                                                         
          for(unsigned int l = 0;l<matrix_size;l++)               
          {                                                                                
            parent1[k][l] = storage[0][k][l];                     
            parent2[k][l] = storage[1][k][l];                                                  
          }
        }   
        if(crossoverType == 1)
	{
	  crossover_->crossover_type1(&parent1,&parent2,&child1,&child2);
	}
	else if(crossoverType == 2)
	{ 
	  crossover_->crossover_type2(&parent1,&parent2,&child1,&child2);
	}
	else if(crossoverType == 3)
	{
	  crossover_->crossover_type3(&parent1,&parent2,&child1,&child2);
	}
	else
	{
	  crossover_->crossover_type4(&parent1,&parent2,&child1,&child2);
	}
        Bot *bot1 = new Bot(bots.at(p1_index),&child1,time);
        Bot *bot2 = new Bot(bots.at(p2_index),&child2,time);
	
	bot1->mutate(mu_modConnection,mu_visualAngle,generations);
        bot2->mutate(mu_modConnection,mu_visualAngle,generations);
	
	bots.push_back(bot1);
        bots.push_back(bot2);
	
      }
     
      /*cout << " " << endl;
      for(unsigned int i=0;i<bots.size();i++)
      {
        cout << "avgTimeToFoods : " << i << "Bot Name " << bots.at(i)->name_ << " " << bots.at(i)->avgTimeToFood_ << endl;
      }
      cout << " " << endl;*/

      Vec1d indexes(bots.size()-2,0);
      int countDelete = 0;
      int del_index1 = 0;
      int del_index2 = 0;
      bool delete_index = false;
      
      // Updating the Avg Time to eat foods for bots that did not eat in a generation. We will add a penalty term to them
      
      int numberOfBots = bots.size()-2;

      for(unsigned int i=0;i<numberOfBots;i++)
      {
        bots.at(i)->updatePenaltyTerm(time,bots.at(i)->birthTime_,bots.at(i)->timeLastEaten_,bots.at(i)->eatenInGeneration_);
	bots.at(i)->updateAvgTimeToEat(bots.at(i)->penaltyTerm_,bots.at(i)->eatenInGeneration_,bots.at(i)->nEatenFoods_);
	avgTimeToFoods[i] = bots.at(i)->avgTimeToFood_;
	avgTimeToFoodsTemp[i] = bots.at(i)->avgTimeToFoodTemp_;
	bots.at(i)->eatenInGeneration_ = 0;
	//cout << "avgTimeToFoodsTemp : " << i << "Bot Name " << bots.at(i)->name_ << " " << bots.at(i)->avgTimeToFoodTemp_ << endl;
      }
      
      //cout << " " << endl;
      // We now sort our avgTimeToFoods vector and get their indexes
      obj->sort_with_index(&avgTimeToFoodsTemp,&indexes,bots.size()-2);
      
      /* We first check if index to be deleted is smaller than the other. For instance we need to delete bot-0 and bot-1 present at
      indexes 0 and 1 respectively. If we first delete bot-0, the vector changes and indexes no. 1 from bots.begin() would be bot-2 and not bot-1*/
      
      del_index1 = indexes[indexes.size()-1];
      del_index2 = indexes[indexes.size()-2];

      if(del_index1<del_index2)
      {
        delete *(bots.begin()+ del_index2);
        bots.erase(bots.begin() + del_index2);
        delete *(bots.begin() + del_index1);
        bots.erase(bots.begin() + del_index1);
      }
      else
      {
        delete *(bots.begin()+ del_index1);
        bots.erase(bots.begin() + del_index1);
        delete *(bots.begin() + del_index2);
        bots.erase(bots.begin() + del_index2);
      }

      /*cout << "After Delete " << endl;
      cout << " " << endl;
      
      for(unsigned int i=0;i<bots.size();i++)
      {
        cout << "avgTimeToFoods : " << i << "Bot Name " << bots.at(i)->name_ << " " << bots.at(i)->avgTimeToFood_ << endl;
      }*/
      
      //cout << " " << endl;
      delIndex.clear();
      start_mate = false;
      ++generations;
      generation_vector.push_back(generations);
      dtime_vector.push_back(dtime);
      dtime=0;
    }
    
    // This code has been replaced by mating code above
    /*for (unsigned int i=0; i<nEaten; ++i)
    {
      delete *(bots.begin());
      bots.erase(bots.begin());
    }*/
    
    // Draw visualization
    if (decodeDebug(debug, 0)==1 && ((time%timeStep==0 && generations>=skipGenerations) || time==1))
    {
      //cout << "x_,y_ : " << bots.at(0)->x_ << " " << bots.at(0)->y_  << endl;
      c_World->cd();
      for (unsigned int i=0; i<bots.size(); ++i) bots.at(i)->draw();
      for (unsigned int i=0; i<foods.size(); ++i) foods.at(i)->draw();
      // for (unsigned int i=0; i<predators.size(); ++i) predators.at(i)->draw();
      text->SetText(0.01, 0.01, ("Generation "+itoa(generations)).c_str());
      text->Draw();
      c_World->Update();
      //cout<< "Reaching End of Draw Visualization" << endl;
      // c_World->SaveAs(("Movie/c_World_"+itoa(time)+".png").c_str());
      // c_World->Print("Movie/Movie_basic.gif+");
    }
    
    if (decodeDebug(debug, 3)==1 && time%timeStep==0 && generations>skipGenerations) // Flash histograms
    {
      for (unsigned int i=0; i<bots.size(); ++i)
      {
        c_Potential_Histograms->cd(i+1);
        bots.at(i)->brain_->drawPotentials();
        c_Distance_Histograms->cd(i+1);
        bots.at(i)->brain_->drawDistances();
        // c_SynapticStrength_Histograms->cd(i+1);
        // bots.at(i)->brain_->drawSynapticStrengths();
      }
      c_Potential_Histograms->Modified();
      c_Potential_Histograms->Update();
      // c_SynapticStrength_Histograms->Modified();
      // c_SynapticStrength_Histograms->Update();
      c_Distance_Histograms->Modified();
      c_Distance_Histograms->Update();
      
      // c_Potential_Histograms->SaveAs("c_Potential_Histograms.png");
      // c_SynapticStrength_Histograms->SaveAs("c_SynapticStrength_Histograms.png");
      // c_Distance_Histograms->SaveAs("c_Distance_Histograms.png");
    }
    
    if (generations%100==0 && generations!=oldGeneration) 
    {
      std::cout<<"Generation "<<generations<<std::endl;
      oldGeneration=generations;
      
      TGraph *g_avgBrainSize_time=new TGraph(avgBrainSize_vector.size(), &time_vector[0], &avgBrainSize_vector[0]);
      g_avgBrainSize_time->SetName("g_avgBrainSize_time");
      g_avgBrainSize_time->SetTitle("; time steps; Average size of brains");

      TGraph *g_avgBrainSize_generation=new TGraph(avgBrainSize_vector.size(), &generation_vector[0], &avgBrainSize_vector[0]);
      g_avgBrainSize_generation->SetName("g_avgBrainSize_generation");
      g_avgBrainSize_generation->SetTitle("; generations; Average size of brains");
      
      TGraph *g_dtime_generation=new TGraph(dtime_vector.size(), &generation_vector[0], &dtime_vector[0]);
      g_dtime_generation->SetName("g_dtime_generation");
      g_dtime_generation->SetTitle("; generations; Time to next meal");
      
      TGraph *g_dtime_time=new TGraph(dtime_vector.size(), &time_vector[0], &dtime_vector[0]);
      g_dtime_time->SetName("g_dtime_time");
      g_dtime_time->SetTitle("; time steps; Time to next meal");
      
      int nSizeMatrix=bots.at(0)->brain_->neurons_.size();
      TH1F *h_spontaneity=new TH1F(("h_spontaneity_"+itoa(generations)).c_str(), "; i^{th} Neuron", nSizeMatrix, 0, nSizeMatrix);
      TH2F *h_distances_matrix=new TH2F(("h_distances_matrix_"+itoa(generations)).c_str(), "; i^{th} Neuron; j^{th} Neuron", nSizeMatrix, 0, nSizeMatrix, nSizeMatrix, 0, nSizeMatrix);
      TH1F *h_distances=new TH1F(("h_distances_Generation_"+itoa(generations)).c_str(), "; distance", 50, 0, 1.0);
      for (unsigned int i=0; i<nBots; ++i)
      {
        // file->mkdir(("Bot_brain_"+itoa(i)).c_str());
        Brain *brain=bots.at(i)->brain_;
        for (unsigned int j=0; j<brain->neurons_.size(); ++j)
        {
          Neuron *neuron=brain->neurons_.at(j);
          h_spontaneity->Fill(j, neuron->spontaneousRate_/double(nBots));
          double sumDistance=0;
          for (unsigned int k=0; k<neuron->neuralRelations_.size(); ++k)
          {
            sumDistance+=neuron->neuralRelations_.at(k)->distance;
            h_distances_matrix->Fill(j, neuron->neuralRelations_.at(k)->index, (neuron->neuralRelations_.at(k)->distance)/double(nSizeMatrix*nBots));
          }
          h_distances->Fill(sumDistance/double(brain->neurons_.size()), 1./double(nBots));
        }
      }
      
      TFile *file;
      if(generations==100)
      {
        file=new TFile("AnalyzeThis.root", "recreate");
        file->mkdir("Brain");
      }
      else file=new TFile("AnalyzeThis.root", "update");
      g_dtime_generation->Write(g_dtime_generation->GetName(), 5 );
      g_dtime_time->Write(g_dtime_time->GetName(), 5 );
      file->cd("Brain");
      h_spontaneity->Write();
      h_distances->Write();
      h_distances_matrix->Write();
      file->Close();
      
      delete g_dtime_generation;
      delete g_dtime_time;
      delete h_spontaneity;
      delete h_distances_matrix;
      delete h_distances;
    }
  }
  
  /*if (decodeDebug(debug, 0)==1) 
  {
    c_World->Print("Movie/Movie_basic.gif++");
    delete c_World;
  }*/
  
  std::cout<<"Exited program after "<<endGeneration<<" generations as requested."<<std::endl;
  
  delete text;
  delete myapp;
  return 0;
}
    
