#include<iostream>
#include<vector>

#include "../interface/Crossover.h"
#include "../interface/ToolBox.h"
using namespace std;

typedef vector<vector<double> >Vec2d;
Crossover::Crossover()
{

}
/*void Crossover::crossover_type1(Vec2d *A,Vec2d *B,Vec2d *child1,Vec2d *child2)
{
    int n = A->size();
    Vec2d child1_temp;
    Vec2d child2_temp;
    for(int i=0;i<n;i++)
    {
        vector<double> vec1(n, 0);
        child1_temp.push_back(vec1);
        child2_temp.push_back(vec1);
        vec1.clear();
    }
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            if(j<n/2)
            {
                child1_temp[i][j] = (*A)[i][j];
                child2_temp[j][i] = (*A)[j][i];
            }
            else
            {
                child1_temp[i][j] = (*B)[i][j];
                child2_temp[j][i] = (*B)[j][i];
            }
        }
    }
    *child1 = child1_temp;
    *child2 = child2_temp;
}*/

void Crossover::crossover_type1(Vec2d *A,Vec2d *B,Vec2d *child1,Vec2d *child2)
{
    int n = A->size();
    Vec2d child1_temp;
    Vec2d child2_temp;
    for(int i=0;i<n;i++)
    {
        vector<double> vec1(n, 0);
        child1_temp.push_back(vec1);
        child2_temp.push_back(vec1);
        vec1.clear();
    }
    double prob = r3->Uniform();
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            if(prob<=0.5)
	    {
              if(j<n/2)
              {
                  child1_temp[j][i] = (*A)[j][i];
                  child2_temp[j][i] = (*B)[j][i];
              }
              else
              {
                  child1_temp[j][i] = (*B)[j][i];
                  child2_temp[j][i] = (*A)[j][i];
              }
	    }
	    else
	    {
              if(j<n/2)
              {
                child1_temp[i][j] = (*A)[i][j];
                child2_temp[i][j] = (*B)[i][j];
              }
              else
              {
                child1_temp[i][j] = (*B)[i][j];
                child2_temp[i][j] = (*A)[i][j];
              }

            }
        }
    }
    *child1 = child1_temp;
    *child2 = child2_temp;
}

void Crossover::crossover_type2(Vec2d *A,Vec2d *B,Vec2d *child1,Vec2d *child2)
{
    int n = A->size();
    Vec2d child1_temp;
    Vec2d child2_temp;
    for(int i=0;i<n;i++)
    {
        vector<double> vec1(n, 0);
        child1_temp.push_back(vec1);
        child2_temp.push_back(vec1);
        vec1.clear();
    }
    int crossover_point = r3->Integer(n);
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            if(j<crossover_point)
            {
                child1_temp[i][j] = (*A)[i][j];
                child2_temp[j][i] = (*B)[j][i];
            }
            else
            {
                child1_temp[i][j] = (*B)[i][j];
                child2_temp[j][i] = (*A)[j][i];
            }
        }
    }
    *child1 = child1_temp;
    *child2 = child2_temp;
}

void Crossover::crossover_type3(Vec2d *A,Vec2d *B,Vec2d *child1,Vec2d *child2)
{
    int n = A->size();
    Vec2d child1_temp;
    Vec2d child2_temp;
    for(int i=0;i<n;i++)
    {
        vector<double> vec1(n, 0);
        child1_temp.push_back(vec1);
        child2_temp.push_back(vec1);
        vec1.clear();
    }
    int crossoverPoint = r3->Integer(n);
    double alphaCross = r3->Uniform();
    for(int i=0;i<n;i++)
    {
      for(int j=0;j<n;j++)
      {
	child1_temp[i][j] = (1-alphaCross)*(*A)[i][j] + alphaCross*(*B)[i][j];
	child2_temp[i][j] = (1-alphaCross)*(*B)[i][j] + alphaCross*(*A)[i][j];             
      }
    }
    *child1 = child1_temp;
    *child2 = child2_temp;
}

void Crossover::crossover_type4(Vec2d *A,Vec2d *B,Vec2d *child1,Vec2d *child2)
{
    int n = A->size();
    Vec2d child1_temp;
    Vec2d child2_temp;
    for(int i=0;i<n;i++)
    {
        vector<double> vec1(n, 0);
        child1_temp.push_back(vec1);
        child2_temp.push_back(vec1);
        vec1.clear();
    }
    for(unsigned int i=0;i<n;i++)
    {
      for(unsigned int j=0;j<n;j++)
      {
        child1_temp[i][j] = (*A)[i][j];
	child2_temp[i][j] = (*B)[i][j];
      }
    }
    int startIndex = 11;
    int endIndex = 30;
    int endIndexRow = 16;
    
    for(unsigned int i=startIndex;i<endIndexRow;i++)
    {
      for(unsigned int j=0;j<endIndex;j++)
      {
        child1_temp[i][j] = (*B)[i][j];
	child2_temp[i][j] = (*A)[i][j];
      }
    }
    
    *child1 = child1_temp;
    *child2 = child2_temp;
}

