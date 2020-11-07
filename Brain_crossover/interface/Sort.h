#include<vector>
#include<iostream>
#include<numeric>

using namespace std; 

class Sort
{
  public:
  
  typedef vector<int> Vec1d;
  Sort();
  void sort_with_index(Vec1d *vec,Vec1d *indexes,int n);
};
