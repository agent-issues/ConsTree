#include <Rcpp.h>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <vector>

// CCD0 maximum a posteriori (CCD0-MAP) tree -- Berling, Klawitter, Bouckaert,
// Xie, Gavryushkin & Drummond (2025), PLOS Comput. Biol. 21:e1012789.
//
// Definition.  Let C be the set of clades observed in k rooted trees, with
// clade credibility cc(X) = f(X) / k (leaves and the root have cc = 1).  The
// CCD0 forest network contains EVERY clade split {A, B} of a clade X in C for
// which A and B are also in C (A u B = X, A n B = 0), observed together or not.
// A tree T assembled from such splits has probability
//   Pr(T) = prod_{X in T} cc(X) / Z,  Z = sum over all such trees of the product,
// and the CCD0-MAP tree maximizes prod cc(X), i.e. sum log cc(X).
//
// Reference: the BEAST 2 'CCD' package (CompEvol/CCD v1.1.0, used by
// TreeAnnotator 2.7.7 -topology CCD0): CCD0.expand() adds the splits,
// setPartitionProbabilities() sets CCP(A|B) = V(A) V(B) / W(X) with
// V(X) = cc(X) W(X), W(X) = sum_splits V(A) V(B), and
// Clade.computeMaxSubtreeLogCCP() runs the MAP dynamic programme.  Maximizing
// the sum of log CCPs is equivalent to maximizing sum log cc because the W
// terms telescope.  Ties: the more balanced split (larger smaller child) wins;
// then the split whose smaller child is lexicographically first, i.e. holds the
// lowest-indexed taxon of the symmetric difference (Clade.java,
// BitSetUtil.getLexicographicFirst).  BEAST numbers leaves in alphabetical
// (Java String.compareTo, i.e. ordinal) order of their labels when it reads a
// NEXUS translate block (TreeParser, sortNodesAlphabetically); the caller
// passes that rank as `tieRank`, so ties are broken identically.
//
// Split enumeration (differs from, but is equivalent to, the reference).
// Taxa are mapped to bit positions in the preorder tip order of the first
// tree, so most clades occupy short position ranges.  Exactly one child of any
// split of X contains min(X), the lowest set position of X; so the candidates
// for that child are the observed clades A with min(A) = min(X),
// |A| < |X| and max(A) <= max(X).  For each candidate that is a subset of X
// the partner X \ A is looked up in a hash table keyed by a Zobrist hash
// (XOR of per-taxon random keys), which for A subset of X is hash(X) ^ hash(A).
// Every split is therefore visited exactly once.  Splits are never stored:
// the bottom-up pass (MAP + normalizer) and the top-down pass (CCD0 clade
// probabilities) each re-enumerate them, so memory is O(m * n / 64) words for
// m clades on n taxa.

namespace {

using std::vector;
typedef uint64_t u64;

const double NEG_INF = -std::numeric_limits<double>::infinity();

inline u64 splitmix64(u64& state) {
  u64 z = (state += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}

inline int popcount64(u64 x) {
  int n = 0;
  while (x) {
    x &= x - 1;
    ++n;
  }
  return n;
}

inline int lowestBit(u64 x) {
  int i = 0;
  while (!(x & 1ULL)) {
    x >>= 1;
    ++i;
  }
  return i;
}

inline int highestBit(u64 x) {
  int i = 63;
  while (!(x & (1ULL << i))) --i;
  return i;
}

struct CladeTable {
  int nWord;
  vector<u64> bits;        // m * nWord, bit positions (not tip indices)
  vector<u64> hash;        // Zobrist hash
  vector<int> count;       // number of trees containing the clade
  vector<int> size, minPos, maxPos;
  vector<int> lastTree;    // guards against counting a clade twice per tree
  vector<int> next;        // hash-chain successor (-1 = none)
  std::unordered_map<u64, int> head;

  explicit CladeTable(int nw) : nWord(nw) {}

  int nClade() const { return static_cast<int>(count.size()); }

  const u64* get(int i) const {
    return bits.data() + static_cast<size_t>(i) * nWord;
  }

  // Index of the clade with these bits, or -1.
  int find(const u64* b, u64 h, int sz, int mn, int mx) const {
    auto it = head.find(h);
    if (it == head.end()) return -1;
    const int w0 = mn >> 6, w1 = mx >> 6;
    for (int i = it->second; i >= 0; i = next[i]) {
      if (size[i] != sz || minPos[i] != mn || maxPos[i] != mx) continue;
      const u64* c = get(i);
      bool same = true;
      for (int w = w0; w <= w1; ++w) {
        if (c[w] != b[w]) {
          same = false;
          break;
        }
      }
      if (same) return i;
    }
    return -1;
  }

  int add(const u64* b, u64 h, int sz, int mn, int mx) {
    const int idx = nClade();
    bits.insert(bits.end(), b, b + nWord);
    hash.push_back(h);
    count.push_back(0);
    size.push_back(sz);
    minPos.push_back(mn);
    maxPos.push_back(mx);
    lastTree.push_back(-1);
    auto it = head.find(h);
    if (it == head.end()) {
      next.push_back(-1);
      head.emplace(h, idx);
    } else {
      next.push_back(it->second);
      it->second = idx;
    }
    return idx;
  }
};

struct Ccd0 {
  int nTip, nTree, nWord;
  vector<int> posOfTip;    // 0-based tip index -> bit position
  vector<int> tipOfPos;    // bit position -> 0-based tip index
  vector<int> tieRank;     // 0-based tip index -> tie-breaking rank
  vector<u64> key;         // Zobrist key per bit position
  CladeTable tab;
  int root;
  vector<vector<int> > byMin;  // clades with a given min position, by size
  vector<int> bySize;          // all clades, increasing size

  Ccd0(int n, int k) : nTip(n), nTree(k), nWord((n + 63) / 64),
    posOfTip(n, -1), tipOfPos(n, -1), key(n), tab((n + 63) / 64), root(-1) {
    u64 state = 0x5DEECE66DULL ^ static_cast<u64>(n);
    for (int i = 0; i < n; ++i) key[i] = splitmix64(state);
  }

  // Lexicographically first of two distinct, equal-sized clades: the one
  // holding the lowest-ranked taxon of their symmetric difference.
  int lexFirst(int a, int b) const {
    const u64* x = tab.get(a);
    const u64* y = tab.get(b);
    int bestRank = nTip;
    int owner = a;
    for (int w = 0; w < nWord; ++w) {
      u64 d = x[w] ^ y[w];
      while (d) {
        const int bit = lowestBit(d);
        d &= d - 1;
        const int rank = tieRank[tipOfPos[(w << 6) + bit]];
        if (rank < bestRank) {
          bestRank = rank;
          owner = (x[w] >> bit) & 1ULL ? a : b;
        }
      }
    }
    return owner;
  }

  int smallerChild(int a, int b) const {
    if (tab.size[a] < tab.size[b]) return a;
    if (tab.size[a] > tab.size[b]) return b;
    return lexFirst(a, b);
  }

  // TRUE if split (a2, b2) beats (a1, b1) under the reference tie rule.
  bool tieWins(int a1, int b1, int a2, int b2) const {
    const int s1 = smallerChild(a1, b1);
    const int s2 = smallerChild(a2, b2);
    if (tab.size[s2] != tab.size[s1]) return tab.size[s2] > tab.size[s1];
    return lexFirst(s1, s2) == s2;
  }

  // Call f(a, b) for every split {a, b} of clade x whose children are both
  // observed clades; `a` is the child containing min(x).
  template <typename F>
  void forEachSplit(int x, F f) const {
    const int sx = tab.size[x];
    const int mx = tab.maxPos[x];
    const u64* bx = tab.get(x);
    const u64 hx = tab.hash[x];
    const vector<int>& cand = byMin[tab.minPos[x]];
    const int w0x = tab.minPos[x] >> 6, w1x = mx >> 6;
    for (const int a : cand) {
      const int sa = tab.size[a];
      if (sa >= sx) break;
      if (tab.maxPos[a] > mx) continue;
      const u64* ba = tab.get(a);
      const int w1a = tab.maxPos[a] >> 6;
      bool subset = true;
      for (int w = w0x; w <= w1a; ++w) {
        if (ba[w] & ~bx[w]) {
          subset = false;
          break;
        }
      }
      if (!subset) continue;
      const u64 hb = hx ^ tab.hash[a];
      auto it = tab.head.find(hb);
      if (it == tab.head.end()) continue;
      const int sb = sx - sa;
      for (int b = it->second; b >= 0; b = tab.next[b]) {
        if (tab.size[b] != sb || tab.minPos[b] <= tab.minPos[x] ||
            tab.maxPos[b] > mx) {
          continue;
        }
        const u64* bb = tab.get(b);
        bool match = true;
        for (int w = w0x; w <= w1x; ++w) {
          if (bb[w] != (bx[w] & ~ba[w])) {
            match = false;
            break;
          }
        }
        if (match) {
          f(a, b);
          break;
        }
      }
    }
  }
};

// Add one tree (ape edge matrix, 1-based, tips 1..nTip) to the clade table.
void addTree(Ccd0& g, const Rcpp::IntegerMatrix& edge, int treeIdx,
             vector<u64>& nodeBits, vector<u64>& nodeHash,
             vector<int>& nodeSize, vector<int>& nodeMin, vector<int>& nodeMax,
             vector<int>& nChild, vector<int>& firstChild,
             vector<int>& sibling, vector<int>& parent, vector<int>& stack,
             vector<int>& order) {
  const int nTip = g.nTip, nWord = g.nWord;
  const int nEdge = edge.nrow();
  int maxNode = nTip;
  for (int e = 0; e < nEdge; ++e) {
    maxNode = std::max(maxNode, std::max(edge(e, 0), edge(e, 1)));
  }
  const int nNode = maxNode + 1;
  if (static_cast<int>(parent.size()) < nNode) {
    nodeBits.resize(static_cast<size_t>(nNode) * nWord);
    nodeHash.resize(nNode);
    nodeSize.resize(nNode);
    nodeMin.resize(nNode);
    nodeMax.resize(nNode);
    nChild.resize(nNode);
    firstChild.resize(nNode);
    sibling.resize(nNode);
    parent.resize(nNode);
  }
  std::fill(parent.begin(), parent.begin() + nNode, 0);
  std::fill(nChild.begin(), nChild.begin() + nNode, 0);
  std::fill(firstChild.begin(), firstChild.begin() + nNode, -1);
  for (int e = 0; e < nEdge; ++e) {
    const int p = edge(e, 0), ch = edge(e, 1);
    if (p < 1 || ch < 1) Rcpp::stop("Invalid edge matrix in tree %i.", treeIdx + 1);
    if (parent[ch]) Rcpp::stop("Node with two parents in tree %i.", treeIdx + 1);
    parent[ch] = p;
    sibling[ch] = firstChild[p];
    firstChild[p] = ch;
    ++nChild[p];
  }
  int rootNode = -1;
  for (int v = nTip + 1; v < nNode; ++v) {
    if (nChild[v] && !parent[v]) {
      if (rootNode != -1) Rcpp::stop("Tree %i has several roots.", treeIdx + 1);
      rootNode = v;
    }
  }
  if (rootNode == -1) Rcpp::stop("Tree %i has no root.", treeIdx + 1);
  // Iterative preorder; process in reverse for a postorder.
  order.clear();
  stack.clear();
  stack.push_back(rootNode);
  while (!stack.empty()) {
    const int v = stack.back();
    stack.pop_back();
    order.push_back(v);
    for (int ch = firstChild[v]; ch != -1; ch = sibling[ch]) stack.push_back(ch);
  }
  int nTipSeen = 0;
  for (auto it = order.rbegin(); it != order.rend(); ++it) {
    const int v = *it;
    u64* bv = nodeBits.data() + static_cast<size_t>(v) * nWord;
    if (v <= nTip) {
      std::fill(bv, bv + nWord, 0ULL);
      const int pos = g.posOfTip[v - 1];
      bv[pos >> 6] |= 1ULL << (pos & 63);
      nodeHash[v] = g.key[pos];
      nodeSize[v] = 1;
      nodeMin[v] = pos;
      nodeMax[v] = pos;
      ++nTipSeen;
      continue;
    }
    std::fill(bv, bv + nWord, 0ULL);
    u64 h = 0;
    int sz = 0, mn = nTip, mx = -1;
    for (int ch = firstChild[v]; ch != -1; ch = sibling[ch]) {
      const u64* bc = nodeBits.data() + static_cast<size_t>(ch) * nWord;
      for (int w = nodeMin[ch] >> 6; w <= (nodeMax[ch] >> 6); ++w) bv[w] |= bc[w];
      h ^= nodeHash[ch];
      sz += nodeSize[ch];
      mn = std::min(mn, nodeMin[ch]);
      mx = std::max(mx, nodeMax[ch]);
    }
    nodeHash[v] = h;
    nodeSize[v] = sz;
    nodeMin[v] = mn;
    nodeMax[v] = mx;
    if (v == rootNode || sz < 2) continue;  // root & leaves counted globally
    int idx = g.tab.find(bv, h, sz, mn, mx);
    if (idx < 0) idx = g.tab.add(bv, h, sz, mn, mx);
    if (g.tab.lastTree[idx] != treeIdx) {
      g.tab.lastTree[idx] = treeIdx;
      ++g.tab.count[idx];
    }
  }
  if (nTipSeen != nTip || nodeSize[rootNode] != nTip) {
    Rcpp::stop("Tree %i does not contain each of the %i leaves exactly once.",
               treeIdx + 1, nTip);
  }
}

}  // namespace

// CCD0-MAP tree.  `edgeList`: one ape edge matrix per tree, tips numbered
// 1..nTip in a common label order, each on its own root.  `tieRank`: 0-based
// rank of each tip for lexicographic tie-breaking.  Returns the MAP tree
// as a cladewise edge matrix (tips 1..nTip, internal nodes nTip+1.. in
// preorder) with per-internal-node clade counts and CCD0 clade probabilities,
// plus summary statistics of the CCD0.
// [[Rcpp::export]]
Rcpp::List ccd0MapCpp(Rcpp::List edgeList, int nTip, Rcpp::IntegerVector tieRank,
                      bool cladeProbability) {
  const int k = edgeList.size();
  if (k < 1) Rcpp::stop("No trees supplied.");
  if (nTip < 2) Rcpp::stop("Need at least two leaves.");
  if (tieRank.size() != nTip) Rcpp::stop("`tieRank` must have one entry per leaf.");
  Ccd0 g(nTip, k);
  g.tieRank.assign(tieRank.begin(), tieRank.end());
  for (int t = 0; t < nTip; ++t) {
    if (g.tieRank[t] < 0 || g.tieRank[t] >= nTip) Rcpp::stop("Invalid `tieRank`.");
  }
  const int nWord = g.nWord;

  // Bit positions: preorder tip order of the first tree, for locality.
  {
    Rcpp::IntegerMatrix e0 = edgeList[0];
    const int nEdge = e0.nrow();
    int maxNode = nTip;
    for (int e = 0; e < nEdge; ++e) maxNode = std::max(maxNode, std::max(e0(e, 0), e0(e, 1)));
    vector<vector<int> > kids(maxNode + 1);
    vector<char> hasParent(maxNode + 1, 0);
    for (int e = 0; e < nEdge; ++e) {
      kids[e0(e, 0)].push_back(e0(e, 1));
      hasParent[e0(e, 1)] = 1;
    }
    int rootNode = -1;
    for (int v = nTip + 1; v <= maxNode; ++v) {
      if (!kids[v].empty() && !hasParent[v]) rootNode = v;
    }
    int next = 0;
    if (rootNode > 0) {
      vector<int> st(1, rootNode);
      while (!st.empty()) {
        const int v = st.back();
        st.pop_back();
        if (v <= nTip) {
          if (v >= 1 && g.posOfTip[v - 1] < 0) g.posOfTip[v - 1] = next++;
          continue;
        }
        for (auto it = kids[v].rbegin(); it != kids[v].rend(); ++it) st.push_back(*it);
      }
    }
    for (int t = 0; t < nTip; ++t) {
      if (g.posOfTip[t] < 0) g.posOfTip[t] = next++;  // malformed: caught below
    }
    for (int t = 0; t < nTip; ++t) g.tipOfPos[g.posOfTip[t]] = t;
  }

  // Leaves (indices 0..nTip-1, by position) and the root (index nTip).
  {
    vector<u64> b(nWord);
    u64 hAll = 0;
    for (int p = 0; p < nTip; ++p) {
      std::fill(b.begin(), b.end(), 0ULL);
      b[p >> 6] |= 1ULL << (p & 63);
      g.tab.add(b.data(), g.key[p], 1, p, p);
      hAll ^= g.key[p];
    }
    std::fill(b.begin(), b.end(), 0ULL);
    for (int p = 0; p < nTip; ++p) b[p >> 6] |= 1ULL << (p & 63);
    g.root = g.tab.add(b.data(), hAll, nTip, 0, nTip - 1);
  }

  {
    vector<u64> nodeBits, nodeHash;
    vector<int> nodeSize, nodeMin, nodeMax, nChild, firstChild, sibling, parent,
      stack, order;
    for (int t = 0; t < k; ++t) {
      Rcpp::IntegerMatrix edge = edgeList[t];
      addTree(g, edge, t, nodeBits, nodeHash, nodeSize, nodeMin, nodeMax,
              nChild, firstChild, sibling, parent, stack, order);
      if ((t & 255) == 255) Rcpp::checkUserInterrupt();
    }
  }
  const int m = g.tab.nClade();
  for (int p = 0; p < nTip; ++p) g.tab.count[p] = k;
  g.tab.count[g.root] = k;
  g.tab.lastTree.clear();
  g.tab.lastTree.shrink_to_fit();

  // Order clades by size; index candidate children by min position.
  {
    vector<int> nOfSize(nTip + 2, 0);
    for (int i = 0; i < m; ++i) ++nOfSize[g.tab.size[i]];
    vector<int> start(nTip + 2, 0);
    for (int s = 1; s <= nTip + 1; ++s) start[s] = start[s - 1] + nOfSize[s - 1];
    g.bySize.assign(m, 0);
    for (int i = 0; i < m; ++i) g.bySize[start[g.tab.size[i]]++] = i;
    g.byMin.assign(nTip, vector<int>());
    for (const int i : g.bySize) g.byMin[g.tab.minPos[i]].push_back(i);
  }

  // Bottom-up: MAP value (sum of log cc), log normalizer W, log #trees.
  const double logK = std::log(static_cast<double>(k));
  vector<double> best(m, NEG_INF), logV(m, NEG_INF), logW(m, NEG_INF),
    logN(m, NEG_INF);
  vector<int> bestA(m, -1), bestB(m, -1);
  double nSplit = 0;
  for (int i = 0; i < nTip; ++i) {
    best[i] = 0;
    logV[i] = 0;
    logW[i] = 0;
    logN[i] = 0;
  }
  int done = 0;
  for (const int x : g.bySize) {
    if (g.tab.size[x] < 2) continue;
    double bx = NEG_INF, wMax = NEG_INF, wSum = 0, nMax = NEG_INF, nSum = 0;
    int ba = -1, bb = -1;
    g.forEachSplit(x, [&](int a, int b) {
      nSplit += 1;
      if (best[a] == NEG_INF || best[b] == NEG_INF) return;
      const double val = best[a] + best[b];
      if (ba < 0 || val > bx + 1e-12 * (1.0 + std::fabs(bx))) {
        bx = val;
        ba = a;
        bb = b;
      } else if (val >= bx - 1e-12 * (1.0 + std::fabs(bx))) {
        if (g.tieWins(ba, bb, a, b)) {
          bx = std::max(bx, val);
          ba = a;
          bb = b;
        }
      }
      const double lw = logV[a] + logV[b];
      if (lw > wMax) {
        wSum = wSum * std::exp(wMax - lw) + 1.0;
        wMax = lw;
      } else {
        wSum += std::exp(lw - wMax);
      }
      const double ln = logN[a] + logN[b];
      if (ln > nMax) {
        nSum = nSum * std::exp(nMax - ln) + 1.0;
        nMax = ln;
      } else {
        nSum += std::exp(ln - nMax);
      }
    });
    if (ba >= 0) {
      const double lcc = std::log(static_cast<double>(g.tab.count[x])) - logK;
      best[x] = bx + lcc;
      bestA[x] = ba;
      bestB[x] = bb;
      logW[x] = wMax + std::log(wSum);
      logV[x] = logW[x] + lcc;
      logN[x] = nMax + std::log(nSum);
    }
    if ((++done & 1023) == 0) Rcpp::checkUserInterrupt();
  }
  if (bestA[g.root] < 0) {
    Rcpp::stop("No fully resolved tree can be assembled from the observed clades.");
  }

  // Top-down: probability of each clade under the CCD0.
  vector<double> prob;
  if (cladeProbability) {
    prob.assign(m, 0.0);
    prob[g.root] = 1.0;
    for (auto it = g.bySize.rbegin(); it != g.bySize.rend(); ++it) {
      const int x = *it;
      if (g.tab.size[x] < 2 || prob[x] == 0 || logW[x] == NEG_INF) continue;
      const double px = prob[x], lw = logW[x];
      g.forEachSplit(x, [&](int a, int b) {
        if (logV[a] == NEG_INF || logV[b] == NEG_INF) return;
        const double q = px * std::exp(logV[a] + logV[b] - lw);
        prob[a] += q;
        prob[b] += q;
      });
      if ((++done & 1023) == 0) Rcpp::checkUserInterrupt();
    }
  }

  // Emit the MAP tree, cladewise.
  Rcpp::IntegerMatrix edge(2 * nTip - 2, 2);
  Rcpp::IntegerVector nodeCount(nTip - 1);
  Rcpp::NumericVector nodeProb(nTip - 1);
  {
    // Depth-first: each edge is followed by the edges of its subtree, and
    // internal nodes are numbered in the order they are reached.  The child
    // containing min(x) -- the first taxon in the first tree's preorder -- is
    // visited first.
    int nextNode = nTip + 1, nextEdge = 0;
    const int rootNode = nextNode++;
    nodeCount[0] = g.tab.count[g.root];
    nodeProb[0] = cladeProbability ? prob[g.root] : NA_REAL;
    vector<std::pair<int, int> > st;  // (clade, parent node number)
    st.emplace_back(bestB[g.root], rootNode);
    st.emplace_back(bestA[g.root], rootNode);
    while (!st.empty()) {
      const int c = st.back().first;
      const int par = st.back().second;
      st.pop_back();
      int node;
      if (g.tab.size[c] == 1) {
        node = g.tipOfPos[g.tab.minPos[c]] + 1;
      } else {
        node = nextNode++;
        nodeCount[node - nTip - 1] = g.tab.count[c];
        nodeProb[node - nTip - 1] = cladeProbability ? prob[c] : NA_REAL;
        st.emplace_back(bestB[c], node);
        st.emplace_back(bestA[c], node);
      }
      edge(nextEdge, 0) = par;
      edge(nextEdge, 1) = node;
      ++nextEdge;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("edge") = edge,
    Rcpp::Named("count") = nodeCount,
    Rcpp::Named("probability") = nodeProb,
    Rcpp::Named("logCredibility") = best[g.root],
    Rcpp::Named("logNormalizer") = logW[g.root],
    Rcpp::Named("logNTrees") = logN[g.root],
    Rcpp::Named("nClades") = m,
    Rcpp::Named("nSplits") = nSplit
  );
}
