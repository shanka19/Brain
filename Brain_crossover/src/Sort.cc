#include<iostream>
#include<vector>
#include<algorithm>
#include <numeric>

#include "../interface/Sort.h"
#include "../interface/ToolBox.h"
using namespace std;

typedef vector<int> Vec1d;
Sort::Sort()
{

}
void Sort::sort_with_index(Vec1d *vec,Vec1d *indexes,int n)
{
  Vec1d vec_temp(n,0);
  for(unsigned int i=0;i<n;i++)
  {
    vec_temp[i] = (*vec)[i];
  }
  std::vector<int> y(n);
  std::iota(y.begin(), y.end(), 0);
  auto comparator = [&vec_temp](int a, int b){ return vec_temp[a] < vec_temp[b]; };
  std::sort(y.begin(), y.end(), comparator);
  *indexes = y;
}
