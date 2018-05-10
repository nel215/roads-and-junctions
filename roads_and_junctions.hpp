#ifndef ROADS_AND_JUNCTIONS_HPP_
#define ROADS_AND_JUNCTIONS_HPP_
#include <vector>
#include <iostream>
using namespace std;

struct P {
  int y, x;
  P(int _y, int _x): y(_y), x(_x) {}
};


class RoadsAndJunctions {
  int S;
  int NC;
  vector<P> cities;
  double junctionCost;
  double failureProbability;

 public:
  vector<int> buildJunctions(int _S, vector<int> _cities, double _junctionCost, double _failureProbability) {
    S = _S;
    NC = _cities.size() / 2;
    cities.assign(NC, P(0, 0));
    for (int i=0; i < NC; i++) {
      int x = _cities[2*i+0];
      int y = _cities[2*i+1];
      cities[i] = P(y, x);
    }
    junctionCost = _junctionCost;
    failureProbability = _failureProbability;
    return vector<int>();
  }
  vector<int> buildRoads(vector<int> junctionStatus) {
    return vector<int>();
  }
};
#endif  // ROADS_AND_JUNCTIONS_HPP_
