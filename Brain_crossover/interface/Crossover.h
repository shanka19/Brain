#include<vector>
#include<iostream>

using namespace std;

class Crossover
{
   public:
   
   typedef vector<vector<double> > Vec2d;
   Crossover();
   void crossover_type1(Vec2d *A,Vec2d *B,Vec2d *child1,Vec2d *child2);
   void crossover_type2(Vec2d *A,Vec2d *B,Vec2d *child1,Vec2d *child2);
   void crossover_type3(Vec2d *A,Vec2d *B,Vec2d *child1,Vec2d *child2);
   void crossover_type4(Vec2d *A,Vec2d *B,Vec2d *child1,Vec2d *child2);
};
