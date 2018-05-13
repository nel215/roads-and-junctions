#ifndef ROADS_AND_JUNCTIONS_HPP_
#define ROADS_AND_JUNCTIONS_HPP_
#include <cmath>
#include <set>
#include <stack>
#include <iostream>
#include <vector>
#include <queue>
#include <utility>
#include <complex>
#include <algorithm>
using namespace std;
const double pi = acos(-1);

#ifdef LOCAL
const double ticks_per_sec = 3200000000;
const double timeLimit = 6.0;
#else
const double ticks_per_sec = 2800000000;
const double timeLimit = 9.1;
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

struct Point {
  int y, x;
  Point(): y(0), x(0) {}
  Point(int _y, int _x): y(_y), x(_x) {}
  bool operator<(const Point &p)const{
    return y*2000+x < p.y*2000+p.x;
  }
};

// geometry
typedef complex<double> P;
namespace std {
  bool operator < (const P& a, const P& b) {
    return real(a) != real(b) ? real(a) < real(b) : imag(a) < imag(b);
  }
}

double cross(const P& a, const P& b) {
  return imag(conj(a)*b);
}

double dot(const P& a, const P& b) {
  return real(conj(a)*b);
}

P rot(const P &p, double theta) {
	double s = sin(theta);
	double c = cos(theta);
	return P(c, -s)*p;
}

int ccw(P a, P b, P c) {
  b -= a; c -= a;
  if (cross(b, c) > 0)   return +1;       // counter clockwise
  if (cross(b, c) < 0)   return -1;       // clockwise
  if (dot(b, c) < 0)     return +2;       // c--a--b on line
  if (norm(b) < norm(c)) return -2;       // a--b--c on line
  return 0;
}

P crossLL(const P &l0, const P &l1, const P &m0, const P &m1) {
    double num = cross(m1-m0, m0-l0);
    double denom = cross(m1-m0, l1-l0);
    return l0+(l1-l0)*num/denom;
}

struct Cluster {
  Point p;
  int cnt;
  Cluster(Point _p, int _cnt) {
    p = _p;
    cnt = _cnt;
  }
};

class DelaunayTriangulation {
  const vector<Point> &origin;
  vector<P> ps;
  stack<pair<int, int>> S;
  vector<vector<int>> em;
  vector<set<int>> E;
  void setTriangle(int i, int j, int r) {
    // cerr << "set " << i << " " << j << " " << r << endl;
    E[i].insert(j);
    E[j].insert(r);
    E[r].insert(i);
    em[i][j] = r;
    em[j][r] = i;
    em[r][i] = j;
    S.push(make_pair(i, j));
  }
  void removeEdge(int i, int j) {
    // cerr << "remove " << i << " " << j << endl;
    E[i].erase(j);
    E[j].erase(i);
    em[i][j] = em[j][i] = -1;
  }
  void decomposeOn(int i, int j, int k, int r) {
    // cerr << "on" << i << " " << j << " " << k << " " << r << endl;
    int m = em[j][i];
    removeEdge(j, i);
    setTriangle(i, m, r);
    setTriangle(m, j, r);
    setTriangle(j, k, r);
    setTriangle(k, i, r);
  }
  void decomposeIn(int i, int j, int k, int r) {
    // cerr << "in" << i << " " << j << " " << k << " " << r << endl;
    setTriangle(i, j, r);
    setTriangle(j, k, r);
    setTriangle(k, i, r);
  }

  void flipEdge(int i, int j, int r) {
    // cerr << "flip " << i << " " << j << " " << r << endl;
    int k = em[j][i];
    removeEdge(i, j);
    setTriangle(i, k, r);
    setTriangle(k, j, r);
  }
  bool isLegal(int i, int j) {
    // cerr << "legal " << i << " " << j << endl;
    return (em[i][j] < 0 || em[j][i] < 0 || !incircle(ps[i], ps[j], ps[em[i][j]], ps[em[j][i]]));
  }
  bool incircle(P a, P b, P c, P p) {
    a -= p; b -= p; c -= p;
    return norm(a) * cross(b, c)
         + norm(b) * cross(c, a)
         + norm(c) * cross(a, b) >= 0; // < : inside, = cocircular, > outside
  }

 public:
  explicit DelaunayTriangulation(const vector<Point> &_origin): origin(_origin) {
  }
  vector<set<int>> solve() {
    const int n = origin.size();
    const int m = n + 3;
    ps.resize(n);
    for (int i=0; i < n; i++) ps[i] = P(origin[i].x, origin[i].y);
    ps.push_back(P(-1e7, -1e7));
    ps.push_back(P(+1e7, -1e7));
    ps.push_back(P(0, 1e7));
    em.assign(m, vector<int>(m, -1));
    E.assign(m, set<int>());
    S = stack<pair<int, int>>();
    setTriangle(n, n+1, n+2);
    for (int r=0; r < n; r++) {
      int i = n, j = n+1, k;
      while (1) {
        k = em[i][j];
        if (ccw(ps[i], ps[em[i][j]], ps[r]) == +1) {
          j = k;
        } else if (ccw(ps[j], ps[em[i][j]], ps[r]) == -1) {
          i = k;
        } else { break; }
      }
      if (ccw(ps[i], ps[j], ps[r]) != +1) {
        decomposeOn(i, j, k, r);
      } else if (ccw(ps[j], ps[k], ps[r]) != +1) {
        decomposeOn(j, k, i, r);
      } else if (ccw(ps[k], ps[i], ps[r]) != +1) {
        decomposeOn(k, i, j, r);
      } else {
        decomposeIn(i, j, k, r);
      }
      while (!S.empty()) {
        int u = S.top().first, v = S.top().second;
        S.pop();
        if (!isLegal(u, v)) flipEdge(u, v, r);
      }
    }
    for (int r=0; r < n; r++) {
      E[r].erase(n);
      E[r].erase(n+1);
      E[r].erase(n+2);
    }
    return E;
  }
};

class RoadsAndJunctions {
  int S;
  int NC;
  vector<Point> cities;
  double junctionCost;
  double failureProbability;
  vector<vector<int> > nearestCity;
  vector<Cluster> candidates;
  double startTime;
  vector<int> bestNumJunctions;
  vector<Cluster> bestCandidates;
  vector<Point> junctions;
  vector<vector<uint8_t>> occupied;
  struct Node {
    int idx;
    int prev;
    int cost;
    Node(int _idx, int _prev, int _cost) {
      idx = _idx;
      prev = _prev;
      cost = _cost;
    }
    bool operator<(const Node &n)const {
      return cost > n.cost;
    }
  };
  struct Edge {
    int f, t;
    int cost;
    Edge(int _f, int _t, int _cost) {
      f = _f;
      t = _t;
      cost = _cost;
    }
  };

 public:
  bool valid(int y, int x)const {
    return 0 <= y && y < S && 0 <= x && x < S;
  }
  inline int mypow(int a) {
    return a*a;
  }
  int calcSteinerCost(int i, int j, int k, int y, int x, const vector<Point> &ps) {
    return mypow(ps[i].y-y) + mypow(ps[j].y-y) + mypow(ps[k].y-y) +
           mypow(ps[i].x-x) + mypow(ps[j].x-x) + mypow(ps[k].x-x);
  }
  void prepareBaseScore() {
#ifdef LOCAL
    double res = 0;
    priority_queue<pair<int, int>> que;
    vector<int> used(NC, 0);
    que.push(make_pair(0, 0));
    while (!que.empty()) {
      double ccost = -que.top().first;
      auto idx = que.top().second;
      que.pop();
      if (used[idx]) continue;
      res += sqrt(ccost);
      int cy = cities[idx].y;
      int cx = cities[idx].x;
      used[idx] = 1;
      for (int i=0; i < NC; i++) {
        if (used[i]) continue;
        int yy = cy-cities[i].y;
        int xx = cx-cities[i].x;
        int ncost = yy*yy+xx*xx;
        que.push(make_pair(-ncost, i));
      }
    }
    cerr << "base_score:" << res << endl;
#endif
  }

  void buildInitialMST() {
    vector<Point> ps = cities;
    auto candOccupied = occupied;
    while (1) {
      const int n = ps.size();
      priority_queue<Node> que;
      vector<int> used(n, 0);
      vector<vector<pair<int, int>>> link(n, vector<pair<int, int>>());
      vector<Edge> edge;
      que.push(Node(0, -1, 0));
      while (!que.empty()) {
        auto cn = que.top();
        que.pop();
        int cy = ps[cn.idx].y;
        int cx = ps[cn.idx].x;
        if (used[cn.idx]) continue;
        used[cn.idx] = 1;
        if (cn.prev != -1) {
          edge.push_back(Edge(cn.prev, cn.idx, cn.cost));
          link[cn.prev].push_back(make_pair(cn.idx, cn.cost));
          link[cn.idx].push_back(make_pair(cn.prev, cn.cost));
        }
        for (int i=0; i < n; i++) {
          if (used[i]) continue;
          int yy = cy-ps[i].y;
          int xx = cx-ps[i].x;
          int ncost = yy*yy+xx*xx;
          que.push(Node(i, cn.idx, ncost));
        }
      }
      const double thres = 2.0943951023931953;
      auto bestPoint = make_pair(make_pair(0.0, 0), Point(-1, -1));
      for (int idx=0; idx < edge.size(); idx++) {
        auto e = edge[idx];
        for (auto tmp : link[e.f]) {
          int v = tmp.first;
          if (v == e.t) continue;
          double oldCost = sqrt(e.cost + tmp.second);
          // v -> f -> t
          P a(ps[v].x, ps[v].y);
          P b(ps[e.f].x, ps[e.f].y);
          P c(ps[e.t].x, ps[e.t].y);
          double d = 0;
          d = max(d, abs(arg((c-b)/(a-b))));
          d = max(d, abs(arg((b-a)/(c-a))));
          d = max(d, abs(arg((a-c)/(b-c))));
          if (d >= thres) continue;
          P p = a + b + c;
          int y = p.imag() / 3;
          int x = p.real() / 3;
          int cost = calcSteinerCost(v, e.f, e.t, y, x, ps);
          while (1) {
            int bestDir = 0;
            int bestCost = 1<<30;
            for (int i=0; i < 4; i++) {
              if (!valid(y+dy[i], x+dx[i])) continue;
              if (candOccupied[y+dy[i]][x+dx[i]]) continue;
              int newCost = calcSteinerCost(v, e.f, e.t, y+dy[i], x+dx[i], ps);
              if (newCost < bestCost) {
                bestDir = i;
                bestCost = newCost;
              }
            }
            if (bestCost < cost) {
              cost = bestCost;
              y += dy[bestDir];
              x += dx[bestDir];
            } else {
              break;
            }
          }
          if (candOccupied[y][x]) continue;
          auto maxDiff = make_pair(0.0, 0);
          for (int cnt=1; cnt < multipleCost.size(); cnt++) {
            double sucP = multipleCost[cnt].second;
            if (sucP < 0.70) continue;
            double d = sucP * (oldCost - (sqrt(cost) + multipleCost[cnt].first + 2*(cnt-1)));
            if (maxDiff.first < d) maxDiff = make_pair(d, cnt);
          }
          if (maxDiff.first > bestPoint.first.first) {
            bestPoint = make_pair(maxDiff, Point(y, x));
          }
        }
      }
      if (bestPoint.first.first <= 0) break;
      cerr << "improve:" << bestPoint.first.first << " " << bestPoint.first.second << " " << bestPoint.second.x << " " << bestPoint.second.y << endl;
      auto p = bestPoint.second;
      if (candOccupied[p.y][p.x] == 0) {
        candOccupied[p.y][p.x] = 1;
        candidates.push_back(Cluster(p, bestPoint.first.second));
        ps.push_back(p);
      }
    }
  }
  double solve(const vector<int> &numJunctions) {
    double res = 0;
    vector<Point> junctions;
    for (int i=0; i < numJunctions.size(); i++) {
      if (numJunctions[i] == 0) continue;
      for (int j=0; j < min(1, numJunctions[i]); j++) {
        // if (rng.uniform() < failureProbability) continue;
        int c = candidates[i].cnt;
        res += multipleCost[c].first;
        int y = candidates[i].p.y + ey[j];
        int x = candidates[i].p.x + ex[j];
        junctions.push_back(Point(y, x));
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
  bool movable(const double score, const double currentScore, const double temperatue) {
    double p = exp(min(2.0, currentScore - score)) - temperatue;
    return p > rng.uniform();
  }
  void anneal() {
    vector<int> actions = {
      0,  // up/down
      1   // move
    };
    cerr << "msg:start annealing" << endl;
    if (candidates.empty()) return;
    vector<int> numJunctions(candidates.size(), 0);
    bestNumJunctions = numJunctions;
    bestCandidates = bestCandidates;
    int sumNumber = 0;
    double bestScore = 1e10;
    double currentScore = 1e10;
    const int trial = 1;
    while (1) {
      double elapsed = getTime() - startTime;
      if (elapsed >= timeLimit) break;
      double temperatue = 1. - elapsed / timeLimit;  // 1->0
      int idx = rng.rand() % candidates.size();
      int act = actions[rng.rand() % 2];
      if (act == 0) {
        int inst = 1 - 2*(rng.rand() & 1);
        if (sumNumber == 0) inst = 1;
        if (sumNumber >= 2*NC) inst = -1;
        if (numJunctions[idx] == 0 && inst == -1) continue;
        if (numJunctions[idx] == 1 && inst == +1) continue;
        numJunctions[idx] += inst;
        sumNumber += inst;
        double score = 0;
        for (int i=0; i < trial; i++) {
          score += solve(numJunctions);;
        }
        score /= trial;
        if (movable(score, currentScore, temperatue)) {
          currentScore = score;
        } else {
          numJunctions[idx] -= inst;
          sumNumber -= inst;
        }
      }
      if (act == 1) {
        int xy = rng.rand() & 1;
        int inst = 1 - 2*(rng.rand() & 1);
        if (xy) {
          candidates[idx].p.x += inst;
        } else {
          candidates[idx].p.y += inst;
        }
        double score = 1e10;
        if (numJunctions[idx] == 0) {
          score = currentScore;
        } else {
          if (valid(candidates[idx].p.y, candidates[idx].p.x)) {
            score = 0;
            for (int i=0; i < trial; i++) {
              score += solve(numJunctions);
            }
            score /= trial;
          }
        }
        if (movable(score, currentScore, temperatue)) {
          currentScore = score;
        } else {
          if (xy) {
            candidates[idx].p.x -= inst;
          } else {
            candidates[idx].p.y -= inst;
          }
        }
      }
      if (currentScore < bestScore) {
        // for (int i=0; i < numJunctions.size(); i++) {
        //   cerr << numJunctions[i] << " ";
        // }
        // cerr << endl;
        bestScore = currentScore;
        bestNumJunctions = numJunctions;
        bestCandidates = candidates;
        cerr << bestScore << endl;
      }
    }
  }
  vector<pair<double, double>> multipleCost;
  void prepareProbability() {
    const int n = 20;
    vector<vector<double> > dp(n, vector<double>(n, 0));
    dp[0][0] = 1.0;
    for (int i=0; i+1 < n; i++) {
      for (int j=0; j+1 < n; j++) {
        dp[i+1][j] += failureProbability * dp[i][j];
        dp[i+1][j+1] += (1.-failureProbability) * dp[i][j];
      }
    }
    multipleCost.resize(n);
    for (int i=0; i < n; i++) {
      multipleCost[i] = make_pair(i*junctionCost, 1.-dp[i][0]);
    }
  }
  vector<int> buildJunctions(int _S, vector<int> _cities, double _junctionCost, double _failureProbability) {
    startTime = getTime();
    S = _S+1;
    NC = _cities.size() / 2;
    cities.assign(NC, Point(0, 0));
    occupied.assign(S, vector<uint8_t>(S, 0));
    for (int i=0; i < NC; i++) {
      int x = _cities[2*i+0];
      int y = _cities[2*i+1];
      cities[i] = Point(y, x);
      occupied[y][x] = 1;
    }
    junctionCost = _junctionCost;
    failureProbability = _failureProbability;
    cerr << "S:" << S << "\tNC:" << NC << "\tFP:" << failureProbability << "\tJC:" << junctionCost << endl;
    prepareBaseScore();
    prepareProbability();
    buildInitialMST();
    anneal();
    junctions.clear();
    vector<int> res;
    cerr << "msg:build multiple" << endl;
    for (int i=0; i < bestNumJunctions.size(); i++) {
      // bestNumJunctions[i] = 1;
      if (bestNumJunctions[i] == 0) continue;
      int c = bestCandidates[i].cnt;
      if (junctions.size()+c >= 2*NC) continue;
      auto sp = bestCandidates[i].p;
      priority_queue<pair<int, Point>> que;
      que.push(make_pair(0, sp));
      while (!que.empty()) {
        int cd = -que.top().first;
        auto p = que.top().second;
        que.pop();
        if (occupied[p.y][p.x] == 0) {
          occupied[p.y][p.x] = 1;
          res.push_back(p.x);
          res.push_back(p.y);
          junctions.push_back(p);
          c--;
          if (c == 0) break;
        }
        for (int j=0; j < 4; j++) {
          int y = p.y + dy[j];
          int x = p.x + dx[j];
          if (!valid(y, x)) continue;
          int d = mypow(y-sp.y) + mypow(x-sp.x);
          if (d <= cd) continue;
          que.push(make_pair(-d, Point(y, x)));
        }
      }
    }
    cerr << "JS:" << junctions.size() << endl;
    return res;
  }
  vector<int> buildRoads(vector<int> junctionStatus) {
    double score = 0;
    vector<Point> points = cities;
    vector<int> status(NC, 1);
    for (int i=0; i < junctionStatus.size(); i++) {
      score += junctionCost;
      points.push_back(junctions[i]);
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
        score += sqrt(-top.first);
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
    cerr << "score:" << score << endl;
    return res;
  }
};
#endif  // ROADS_AND_JUNCTIONS_HPP_
