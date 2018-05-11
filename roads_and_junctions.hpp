#ifndef ROADS_AND_JUNCTIONS_HPP_
#define ROADS_AND_JUNCTIONS_HPP_
#include <cmath>
#include <iostream>
#include <vector>
#include <queue>
#include <utility>
#include <algorithm>
using namespace std;

#ifdef LOCAL
const double ticks_per_sec = 3200000000;
const double timeLimit = 6.0;
#else
const double ticks_per_sec = 2800000000;
const double timeLimit = 9.0;
#endif  // LOCAL
inline double getTime() {
    uint32_t lo, hi;
    asm volatile ("rdtsc" : "=a" (lo), "=d" (hi));
    return (((uint64_t)hi << 32) | lo) / ticks_per_sec;
}

// rng
class XorShift {
  uint32_t x;
  uint32_t y;
  uint32_t z;
  uint32_t w;
  uint64_t max_uint32 = static_cast<uint32_t>(-1);
 public:
  explicit XorShift(int seed) {
    std::srand(seed);
    x = std::rand();
    y = std::rand();
    z = std::rand();
    w = std::rand();
  }
  uint32_t rand() {
    uint32_t t = x ^ (x << 11);
    x = y; y = z; z = w;
    return w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
  }
  double uniform() {
    double a = rand();
    return a/(max_uint32+1);
  }
};
XorShift rng(215);

const int dy[4] = {0, 1, 0, -1};
const int dx[4] = {1, 0, -1, 0};
const int ey[4] = {0, 0, 1, 1};
const int ex[4] = {0, 1, 0, 1};

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
  double startTime;
  vector<int> bestNumJunctions;
  vector<P> junctions;

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
  double solve(const vector<int> &numJunctions) {
    double res = 0;
    vector<P> junctions;
    for (int i=0; i < numJunctions.size(); i++) {
      if (numJunctions[i] == 0) continue;
      // TODO: consider probablity
      res += numJunctions[i] * junctionCost;
      for (int j=0; j < min(4, numJunctions[i]); j++) {
        int y = candidates[i].y + ey[j];
        int x = candidates[i].x + ex[j];
        junctions.push_back(P(y, x));
      }
    }
    auto points = cities;
    points.insert(points.end(), junctions.begin(), junctions.end());
    priority_queue<pair<int, int>> que;
    que.push(make_pair(0, 0));
    vector<int> used(points.size(), 0);
    while (!que.empty()) {
      auto top = que.top();
      auto idx = top.second;
      que.pop();
      if (used[idx]) continue;
      used[idx] = 1;
      res += sqrt(-top.first);
      for (int i=0; i < points.size(); i++) {
        if (used[i]) continue;
        int cy = points[idx].y - points[i].y;
        int cx = points[idx].x - points[i].x;
        int d = cy*cy+cx*cx;
        que.push(make_pair(-d, i));
      }
    }
    return res;
  }
  void anneal() {
    vector<int> numJunctions(candidates.size(), 0);
    bestNumJunctions = numJunctions;
    int sumNumber = 0;
    double bestScore = 1e10;
    double currentScore = 1e10;
    while (getTime() - startTime < timeLimit) {
      int idx = rng.rand() % candidates.size();
      int inst = 1 - 2*(rng.rand() & 1);
      if (sumNumber == 0) inst = 1;
      if (sumNumber >= 2*NC) inst = -1;
      if (numJunctions[idx] == 0 && inst == -1) continue;
      numJunctions[idx] += inst;
      sumNumber += inst;
      double score = solve(numJunctions);
      // cerr << score << endl;
      if (score < currentScore) {
        currentScore = score;
      } else {
        numJunctions[idx] -= inst;
        sumNumber -= inst;
      }
      if (currentScore < bestScore) {
        for (int i=0; i < numJunctions.size(); i++) {
          cerr << numJunctions[i] << " ";
        }
        cerr << endl;
        bestScore = currentScore;
        bestNumJunctions = numJunctions;
        cerr << bestScore << endl;
      }
    }
  }
  vector<int> buildJunctions(int _S, vector<int> _cities, double _junctionCost, double _failureProbability) {
    startTime = getTime();
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
    anneal();
    junctions.clear();
    vector<int> res;
    for (int i=0; i < bestNumJunctions.size(); i++) {
      if (bestNumJunctions[i] == 0) continue;
      for (int j=0; j < min(4, bestNumJunctions[i]); j++) {
        int y = candidates[i].y+ey[j];
        int x = candidates[i].x+ex[j];
        res.push_back(x);
        res.push_back(y);
        junctions.push_back(P(y, x));
      }
    }
    cerr << "JS:" << junctions.size() << endl;
    return res;
  }
  vector<int> buildRoads(vector<int> junctionStatus) {
    vector<P> points = cities;
    vector<int> status(NC, 1);
    for (int i=0; i < junctionStatus.size(); i++) {
      points.push_back(junctions[i]);
      status.push_back(junctionStatus[i]);
    }
    priority_queue<pair<int, pair<int, int>>> que;
    que.push(make_pair(0, make_pair(0, -1)));
    vector<int> res;
    vector<int> used(points.size(), 0);
    double c = 0;
    while (!que.empty()) {
      auto top = que.top();
      auto idx = top.second.first;
      auto prev = top.second.second;
      que.pop();
      if (used[idx]) continue;
      used[idx] = 1;
      if (prev != -1) {
        c += sqrt(-top.first);
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
    cerr << c << endl;
    return res;
  }
};
#endif  // ROADS_AND_JUNCTIONS_HPP_
