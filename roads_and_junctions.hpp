#ifndef ROADS_AND_JUNCTIONS_HPP_
#define ROADS_AND_JUNCTIONS_HPP_
#include <cmath>
#include <iostream>
#include <vector>
#include <queue>
#include <utility>
using namespace std;

const int dy[4] = {0, 1, 0, -1};
const int dx[4] = {1, 0, -1, 0};

struct P {
  int y, x;
  P(int _y, int _x): y(_y), x(_x) {}
  bool operator<(const P &p)const{
    return y < p.y;
  }
};


class RoadsAndJunctions {
  int S;
  int NC;
  vector<P> cities;
  double junctionCost;
  double failureProbability;
  vector<vector<int> > nearestCity;
  vector<P> candidates;

 public:
  bool valid(int y, int x)const {
    return 0 <= y && y < S && 0 <= x && x < S;
  }
  void buildRegions() {
    nearestCity.assign(S, vector<int>(S, -1));
    priority_queue<pair<double, pair<P, int>>> que;
    for (int i=0; i < NC; i++) {
      que.push(make_pair(0, make_pair(cities[i], i)));
    }
    while (!que.empty()) {
      auto top = que.top().second;
      auto pos = top.first;
      auto parent = top.second;
      que.pop();
      if (nearestCity[pos.y][pos.x] != -1) continue;
      nearestCity[pos.y][pos.x] = parent;
      for (int i=0; i < 4; i++) {
        int ny = pos.y + dy[i];
        int nx = pos.x + dx[i];
        if (!valid(ny, nx)) continue;
        if (nearestCity[ny][nx] != -1) continue;
        int cy = ny-cities[parent].y;
        int cx = nx-cities[parent].x;
        double d = cy*cy+cx*cx;
        que.push(make_pair(-d, make_pair(P(ny, nx), parent)));
      }
    }
  }
  void buildCandidates() {
    candidates.clear();
    vector<int> used(NC, 0);
    const int ey[4] = {0, 0, 1, 1};
    const int ex[4] = {0, 1, 0, 1};
    for (int y=0; y < S-1; y++) {
      for (int x=0; x < S-1; x++) {
        int d = 0;
        for (int i=0; i < 4; i++) {
          d += used[nearestCity[y+ey[i]][x+ex[i]]] == 0;
          used[nearestCity[y+ey[i]][x+ex[i]]]++;
        }
        for (int i=0; i < 4; i++) used[nearestCity[y+ey[i]][x+ex[i]]]--;
        if (d >= 3) candidates.push_back(P(y, x));
      }
    }
  }
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
    cerr << "S:" << S << "\tNC:" << NC << endl;
    buildRegions();
    buildCandidates();
    vector<int> res;
    for (int i=0; i < candidates.size(); i++) {
      res.push_back(candidates[i].x);
      res.push_back(candidates[i].y);
      if (res.size() >= NC*2*2) break;
    }
    cerr << "JS:" << res.size()/2 << endl;
    return res;
  }
  vector<int> buildRoads(vector<int> junctionStatus) {
    vector<P> points = cities;
    vector<int> status(NC, 1);
    for (int i=0; i < junctionStatus.size(); i++) {
      points.push_back(candidates[i]);
      status.push_back(junctionStatus[i]);
    }
    priority_queue<pair<int, pair<int, int>>> que;
    que.push(make_pair(0, make_pair(0, -1)));
    vector<int> res;
    vector<int> used(points.size(), 0);
    while (!que.empty()) {
      auto top = que.top();
      auto idx = top.second.first;
      auto prev = top.second.second;
      que.pop();
      if (used[idx]) continue;
      used[idx] = 1;
      if (prev != -1) {
        res.push_back(prev);
        res.push_back(idx);
      }
      for (int i=0; i < points.size(); i++) {
        if (used[i] || status[i] == 0) continue;
        int cy = points[idx].y - points[i].y;
        int cx = points[idx].x - points[i].x;
        int d = cy*cy+cx*cx;
        que.push(make_pair(-d, make_pair(i, idx)));
      }
    }
    return res;
  }
};
#endif  // ROADS_AND_JUNCTIONS_HPP_
