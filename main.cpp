#include "roads_and_junctions.hpp"

template<class T> void getVector(vector<T>& v) {
    for (int i = 0; i < v.size(); ++i)
        cin >> v[i];
}

int main() {
  RoadsAndJunctions rj;
  int S, C;
  cin >> S >> C;
  vector<int> cities(C);
  getVector(cities);
  double junctionCost, failureProbability;
  cin >> junctionCost >> failureProbability;

  vector<int> ret = rj.buildJunctions(S, cities, junctionCost, failureProbability);
  cout << ret.size() << endl;
  for (int i = 0; i < ret.size(); ++i) {
    cout << ret[i] << endl;
  }
  cout.flush();

  int J;
  cin >> J;
  vector<int> junctionStatus(J);
  getVector(junctionStatus);

  ret = rj.buildRoads(junctionStatus);
  cout << ret.size() << endl;
  for (int i = 0; i < ret.size(); ++i) {
    cout << ret[i] << endl;
  }
  cout.flush();
}

