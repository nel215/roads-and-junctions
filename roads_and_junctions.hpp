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
const double timeLimit = 9.6;
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

const int dy[8] = {0, 1, 0, -1, -1, -1, 1, 1};
const int dx[8] = {1, 0, -1, 0, -1, 1, -1, 1};
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
  double calcSteinerCostDouble(int i, int j, int k, int y, int x, const vector<Point> &ps) {
    return sqrt(mypow(ps[i].y-y) + mypow(ps[i].x-x)) +
           sqrt(mypow(ps[j].y-y) + mypow(ps[j].x-x)) +
           sqrt(mypow(ps[k].y-y) + mypow(ps[k].x-x));
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

  double calcSteinerCost2Double(int a, int b, int c, int d, int y1, int x1, int y2, int x2, const vector<Point> &ps) {
    double ayy = ps[a].y - y1, axx = ps[a].x - x1;
    double byy = ps[b].y - y1, bxx = ps[b].x - x1;
    double cyy = ps[c].y - y2, cxx = ps[c].x - x2;
    double dyy = ps[d].y - y2, dxx = ps[d].x - x2;
    double yy = y1 - y2, xx = x1 - x2;
    return sqrt(mypow(ayy) + mypow(axx)) +
           sqrt(mypow(byy) + mypow(bxx)) +
           sqrt(mypow(cyy) + mypow(cxx)) +
           sqrt(mypow(dyy) + mypow(dxx)) +
           sqrt(mypow(yy) + mypow(xx));
  }
  double solveForBS(vector<Point> ps, vector<set<int>> &link, vector<Edge> &edge) {
    const int n = ps.size();
    priority_queue<Node> que;
    que.push(Node(0, -1, 0));
    vector<uint8_t> used(n, 0);
    link.assign(n, set<int>());
    edge.clear();
    double res = 0;
    while (!que.empty()) {
      auto cn = que.top();
      que.pop();
      auto p = ps[cn.idx];
      // if (cn.cost == 0) return 1e10;
      if (used[cn.idx]) continue;
      used[cn.idx] = 1;
      res += sqrt(cn.cost);
      if (cn.prev != -1) {
        edge.push_back(Edge(cn.prev, cn.idx, cn.cost));
        link[cn.prev].insert(cn.idx);
        link[cn.idx].insert(cn.prev);
      }
      for (int i=0; i < n; i++) {
        if (used[i]) continue;
        int yy = p.y-ps[i].y;
        int xx = p.x-ps[i].x;
        int ncost = yy*yy+xx*xx;
        que.push(Node(i, cn.idx, ncost));
      }
    }
    while (1) {
      bool deleted = false;
      for (int i=NC; i < n; i++) {
        if (link[i].size() == 1) {
          auto it = link[i].begin();
          int a = *it;
          double ai = sqrt(dist2(ps[a], ps[i]));
          link[a].erase(i);
          link[i].clear();
          res -= ai;
          deleted = true;
        }
      }
      for (int i=NC; i < n; i++) {
        if (link[i].size() == 2) {
          auto it = link[i].begin();
          int a = *it;
          it++;
          int b = *it;
          double ab = sqrt(dist2(ps[a], ps[b]));
          double ai = sqrt(dist2(ps[a], ps[i]));
          double bi = sqrt(dist2(ps[b], ps[i]));
          link[a].erase(i);
          link[b].erase(i);
          link[a].insert(b);
          link[b].insert(a);
          link[i].clear();
          res += ab - (ai+bi);
          deleted = true;
        }
      }
      if (deleted) continue;
      for (int i=NC; i < n; i++) {
        if (link[i].size() == 0) continue;
        if (link[i].size() == 3) {
          auto it = link[i].begin();
          int a = *it;
          it++;
          int b = *it;
          it++;
          int c = *it;
          double ab = sqrt(dist2(ps[a], ps[b]));
          double bc = sqrt(dist2(ps[b], ps[c]));
          double ca = sqrt(dist2(ps[c], ps[a]));
          double ai = sqrt(dist2(ps[a], ps[i]));
          double bi = sqrt(dist2(ps[b], ps[i]));
          double ci = sqrt(dist2(ps[c], ps[i]));
          double maxEdge = max(ab, max(bc, ca));
          double baseCost = ab + bc + ca - maxEdge;
          double bestCost = baseCost;
          int bestCnt = 0;
          for (int cnt=1; cnt < 10; cnt++) {
            double sucP = multipleCost[cnt].second;
            double cost = sucP*(ai+bi+ci) + (1.-sucP)*baseCost + cnt*junctionCost;
            if (bestCost > cost) {
              bestCost = cost;
              bestCnt = cnt;
            }
          }
          if (bestCnt > 0) continue;
          link[a].erase(i);
          link[b].erase(i);
          link[c].erase(i);
          link[i].clear();
          if (maxEdge == ab) {
            link[b].insert(c);
            link[c].insert(b);
            link[c].insert(a);
            link[a].insert(c);
          } else if (maxEdge == bc) {
            link[a].insert(b);
            link[b].insert(a);
            link[c].insert(a);
            link[a].insert(c);
          } else {
            link[a].insert(b);
            link[b].insert(a);
            link[b].insert(c);
            link[c].insert(b);
          }
          res += baseCost - (ai+bi+ci);
          deleted = true;
          break;
        }
      }
      if (!deleted) break;
    }
    for (int i=NC; i < n; i++) {
      if (link[i].size() == 0) continue;
      if (link[i].size() == 3) {
        auto it = link[i].begin();
        int a = *it;
        it++;
        int b = *it;
        it++;
        int c = *it;
        double ab = sqrt(dist2(ps[a], ps[b]));
        double bc = sqrt(dist2(ps[b], ps[c]));
        double ca = sqrt(dist2(ps[c], ps[a]));
        double ai = sqrt(dist2(ps[a], ps[i]));
        double bi = sqrt(dist2(ps[b], ps[i]));
        double ci = sqrt(dist2(ps[c], ps[i]));
        double baseCost = ab + bc + ca - max(ab, max(bc, ca));
        double bestCost = baseCost;
        int bestCnt = 0;
        for (int cnt=1; cnt < 10; cnt++) {
          double sucP = multipleCost[cnt].second;
          double cost = sucP*(ai+bi+ci) + (1.-sucP)*baseCost + cnt*junctionCost;
          if (bestCost > cost) {
            bestCost = cost;
            bestCnt = cnt;
          }
        }
        res += bestCost - (ai+bi+ci);
      } else {
        res += junctionCost;
      }
    }
    return res;
  }
  struct BeamNode {
    vector<Point> ps;
    double cost;
    BeamNode() {
      cost = 1e10;
    }
    BeamNode(vector<Point> _ps, double _cost) {
      ps = _ps;
      cost = _cost;
    }
    bool operator<(const BeamNode &n)const {
      return cost < n.cost;
    }
  };
  inline int dist2(const Point &a, const Point &b) {
    return mypow(a.y-b.y) + mypow(a.x-b.x);
  }
  void buildInitialMST() {
    const double thres = 2.0943951023931953;
    vector<vector<uint32_t> > hashSeed(S, vector<uint32_t>(S));
    for (int i=0; i < S; i++) {
      for (int j=0; j < S; j++) {
        hashSeed[i][j] = rng.rand();
      }
    }
    int beamSize = 1 * 100 / NC;
    if (NC < 50) beamSize *= 3;
    if (NC < 30) beamSize *= 3;
    vector<BeamNode> curQue;
    vector<set<int>> link;
    vector<Edge> edge;
    vector<set<int>> tmpLink;
    vector<Edge> tmpEdge;
    curQue.push_back(BeamNode(cities, solveForBS(cities, tmpLink, tmpEdge)));
    auto best = curQue[0];
    for (int q=0; q < 30; q++) {
      vector<BeamNode> nextQue;
      set<uint32_t> hashSet;
      for (int bi=0; bi < curQue.size(); bi++) {
        // construct
        auto bn = curQue[bi];
        // cerr << q << " " << bn.cost << endl;
        vector<Point> ps = bn.ps;
        auto candOccupied = occupied;
        uint32_t hash = 0;
        for (int i=NC; i < ps.size(); i++) {
          candOccupied[ps[i].y][ps[i].x] = 1;
          hash ^= hashSeed[ps[i].y][ps[i].x];
        }
        // solve current graph
        double currentScore = solveForBS(ps, link, edge);
        DelaunayTriangulation tri(ps);
        auto _graph = tri.solve();
        vector<vector<int>> graph(NC, vector<int>());
        for (int i=0; i < NC; i++) graph[i] = vector<int>(_graph[i].begin(), _graph[i].end());
        for (int i=0; i < NC; i++) {
          for (int j=0; j < graph[i].size(); j++) {
            int v = graph[i][j];
            if (i > v) continue;
            for (int k=j+1; k < graph[i].size(); k++) {
              int u = graph[i][k];
              // a -> b -> c
              P a(ps[i].x, ps[i].y);
              P b(ps[v].x, ps[v].y);
              P c(ps[u].x, ps[u].y);
              double oldCost = abs(a-b)+abs(b-c)+abs(c-a)-max(abs(a-b), max(abs(b-c), abs(c-a)));
              double d = 0;
              d = max(d, abs(arg((c-b)/(a-b))));
              d = max(d, abs(arg((b-a)/(c-a))));
              d = max(d, abs(arg((a-c)/(b-c))));
              if (d >= thres) continue;
              P p = a + b + c;
              int y = p.imag() / 3;
              int x = p.real() / 3;
              double cost = calcSteinerCostDouble(i, u, v, y, x, ps);
              while (1) {
                int bestDir = 0;
                double bestCost = 1<<30;
                for (int i=0; i < 8; i++) {
                  if (!valid(y+dy[i], x+dx[i])) continue;
                  if (candOccupied[y+dy[i]][x+dx[i]]) continue;
                  double newCost = calcSteinerCostDouble(i, u, v, y+dy[i], x+dx[i], ps);
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
              if (oldCost < cost) continue;
              if (candOccupied[y][x]) continue;
              vector<Point> jcs = ps;
              jcs.push_back(Point(y, x));
              double score = solveForBS(jcs, tmpLink, tmpEdge);
              uint32_t nextHash = hash ^ hashSeed[y][x];
              if (hashSet.count(nextHash)) continue;
              hashSet.insert(nextHash);
              nextQue.push_back(BeamNode(jcs, score));
            }
          }
        }

        // find point
        for (int idx=0; idx < edge.size(); idx++) {
          auto e = edge[idx];
          for (auto v : link[e.f]) {
            // int v = tmp.first;
            if (v == e.t) continue;
            double oldCost = sqrt(e.cost) + sqrt(dist2(ps[v], ps[e.f]));
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
            double cost = calcSteinerCostDouble(v, e.f, e.t, y, x, ps);
            while (1) {
              int bestDir = 0;
              double bestCost = 1<<30;
              for (int i=0; i < 8; i++) {
                if (!valid(y+dy[i], x+dx[i])) continue;
                if (candOccupied[y+dy[i]][x+dx[i]]) continue;
                double newCost = calcSteinerCostDouble(v, e.f, e.t, y+dy[i], x+dx[i], ps);
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
            if (oldCost < cost) continue;
            if (candOccupied[y][x]) continue;
            vector<Point> jcs = ps;
            jcs.push_back(Point(y, x));

            double score = solveForBS(jcs, tmpLink, tmpEdge);
            uint32_t nextHash = hash ^ hashSeed[y][x];
            if (hashSet.count(nextHash)) continue;
            hashSet.insert(nextHash);
            nextQue.push_back(BeamNode(jcs, score));
          }
        }
      }
      if (nextQue.empty()) break;
      sort(nextQue.begin(), nextQue.end());
      if (nextQue.size() > beamSize) nextQue.resize(beamSize);
      curQue = move(nextQue);
      best = min(best, curQue[0]);
    }
    cerr << best.cost << endl;
    for (int i=NC; i < best.ps.size(); i++) {
      // cerr << best.ps[i].x << " " << best.ps[i].y << endl;
      candidates.push_back(Cluster(best.ps[i], 1));
    }
  }
  double solve(const vector<int> &numJunctions) {
    vector<Point> junctions;
    for (int i=0; i < numJunctions.size(); i++) {
      auto a = candidates[i].p;
      if (!valid(a.y, a.x)) return 1e10;
      if (occupied[a.y][a.x]) return 1e10;
      if (numJunctions[i] > 0) junctions.push_back(a);
      for (int j=i+1; j < numJunctions.size(); j++) {
        auto b = candidates[j].p;
        if (a.y == b.y && a.x == b.x) return 1e10;
      }
    }
    auto ps = cities;
    ps.insert(ps.end(), junctions.begin(), junctions.end());
    vector<set<int>> link;
    vector<Edge> edge;
    double res = solveForBS(ps, link, edge);
    return res;
  }
  bool movable(const double score, const double currentScore, const double temperatue) {
    double p = exp(min(2.0, (currentScore - score))) + temperatue;
    return p > rng.uniform();
  }
  void anneal() {
    vector<int> actions = {
      // 1,  // up/down
      // 1,  // move
      // 1,  // move
      // 1,  // move
      // 1,  // move
      // 1,  // move
      // 1,  // move
      // 1,  // move
      1   // move
    };
    cerr << "msg:start annealing" << endl;
    if (candidates.empty()) return;
    vector<int> numJunctions(candidates.size(), 1);
    bestNumJunctions = numJunctions;
    bestCandidates = candidates;
    // return;
    double bestScore = solve(numJunctions);
    double currentScore = bestScore;
    int annealCount = 0;
    while (1) {
      double elapsed = getTime() - startTime;
      if (elapsed >= timeLimit) break;
      annealCount++;
      // double temperatue = 1. - elapsed / timeLimit;  // 1->0
      double temperatue = 0.6;
      int idx = rng.rand() % candidates.size();
      int act = actions[rng.rand() % actions.size()];
      if (act == 0) {
        numJunctions[idx] = 1 - numJunctions[idx];
        double score = solve(numJunctions);
        if (movable(score, currentScore, temperatue)) {
          currentScore = score;
        } else {
          numJunctions[idx] = 1 - numJunctions[idx];
        }
      }
      if (act == 1) {
        int xy = rng.rand() & 1;
        int inst = 1 - rng.rand() % 3;
        if (xy) {
          candidates[idx].p.x += inst;
        } else {
          candidates[idx].p.y += inst;
        }
        double score = solve(numJunctions);
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
        bestScore = currentScore;
        bestNumJunctions = numJunctions;
        bestCandidates = candidates;
        cerr << bestScore << endl;
      }
    }
    cerr << "anneal:" << annealCount << endl;
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
      // multipleCost[i] = make_pair(i*junctionCost, min(1.0, (1.-dp[i][0])*1.20));
      multipleCost[i] = make_pair(i*junctionCost, min(1.0, (1.-dp[i][0])*1.00));
      // multipleCost[i] = make_pair(i*junctionCost, min(1.0, (1.-dp[i][0])*1.00));
    }
  }
  struct Cand {
    double diff;
    int count;
    int a, b, c;
    Point p;
    Cand() {
      diff = -1e9;
    }
    Cand(double _diff, int _count, int _a, int _b, int _c, Point _p) {
      diff = _diff;
      count = _count;
      a = _a;
      b = _b;
      c = _c;
      p = _p;
    }
  };
  double evaluate(const set<int> &index, const vector<Point> &ps, const vector<int> jcCount) {
    double sum = 0;
    priority_queue<Node> que;
    vector<uint8_t> used(3*NC, 0);
    que.push(Node(0, -1, 0));
    while (!que.empty()) {
      auto cn = que.top();
      que.pop();
      if (used[cn.idx]) continue;
      if (cn.idx >= NC) sum += junctionCost * jcCount[cn.idx];
      sum += sqrt(cn.cost);
      used[cn.idx] = 1;
      for (int i : index) {
        if (used[i]) continue;
        int d = mypow(ps[i].y-ps[cn.idx].y) + mypow(ps[i].x-ps[cn.idx].x);
        que.push(Node(i, cn.idx, d));
      }
    }
    return sum;
  }
  void setBestClusters(const set<int> &index, const vector<Point> &ps, const vector<int> &jcCount, vector<Cluster> &bestCluster) {
    bestCluster.clear();
    for (int i=NC; i < 3*NC; i++) {
      if (index.count(i) == 0) continue;
      bestCluster.push_back(Cluster(ps[i], jcCount[i]));
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
    vector<Point> ps = cities;
    for (int i=0; i < bestCandidates.size() ; i++) {
      if (bestNumJunctions[i] == 0) continue;
      ps.push_back(bestCandidates[i].p);
    }
    const int n = ps.size();
    bestNumJunctions.clear();
    vector<set<int>> link;
    vector<Edge> edge;
    double tmpScore = solveForBS(ps, link, edge);
    for (int i=NC; i < n; i++) {
      if (link[i].size() == 0) {
        bestNumJunctions.push_back(0);
      } else if (link[i].size() == 3) {
        auto it = link[i].begin();
        int a = *it;
        it++;
        int b = *it;
        it++;
        int c = *it;
        double ab = sqrt(dist2(ps[a], ps[b]));
        double bc = sqrt(dist2(ps[b], ps[c]));
        double ca = sqrt(dist2(ps[c], ps[a]));
        double ai = sqrt(dist2(ps[a], ps[i]));
        double bi = sqrt(dist2(ps[b], ps[i]));
        double ci = sqrt(dist2(ps[c], ps[i]));
        double baseCost = ab + bc + ca - max(ab, max(bc, ca));
        double bestCost = baseCost;
        int bestCnt = 0;
        for (int cnt=1; cnt < 10; cnt++) {
          double sucP = multipleCost[cnt].second;
          double cost = sucP*(ai+bi+ci) + (1.-sucP)*baseCost + cnt*junctionCost;
          if (bestCost > cost) {
            bestCost = cost;
            bestCnt = cnt;
          }
        }
        bestNumJunctions.push_back(bestCnt);
      } else {
        bestNumJunctions.push_back(1);
      }
    }
    cerr << tmpScore << endl;

    for (int i=0; i < bestNumJunctions.size(); i++) {
      if (bestNumJunctions[i] == 0) continue;
      int c = bestNumJunctions[i];
      if (junctions.size()+c >= 2*NC) continue;
      auto sp = ps[i+NC];
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
