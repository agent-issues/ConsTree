#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// HIPSTR / MrHIPSTR summary trees -- Baele et al. (2025), "HIPSTR: highest
// independent posterior subtree reconstruction in TreeAnnotator X",
// Bioinformatics 41(10):btaf488.  Implemented from the paper and checked
// line-by-line against the reference implementations (TreeAnnotator X,
// BEAST X v10.5.0, dr.app.tools.treeannotator.HIPSTRTreeBuilder; and
// phylotreelib 2.3.1 SummaryTreeBuilder.hipstr_tree(), as used by sumt 4.1.5).
//
// Pass 1 (one sweep over the input) records every observed clade, the number
// of trees containing it, and every observed pair of child clades.  Pass 2 is a
// dynamic programme over clades in increasing size:
//   score(C) = log CF(C) + max over observed child pairs (A, B) of
//              score(A) + score(B),        score(tip) = 0,
// and the tree is read off the argmax pairs from the root down.  MrHIPSTR
// maximises (number of clades with CF > 1/2, score) lexicographically; this is
// the limit of the references' "+1e10 per majority clade" reward, without the
// loss of precision that the reward incurs in double arithmetic.
//
// Clades are identified by a 128-bit Zobrist hash: each tip carries two
// independent pseudo-random 64-bit words (fixed seed, so runs are
// reproducible) and a clade's key is their sum modulo 2^64.  Two distinct leaf
// sets collide with probability 2^-128 per pair, so identification costs O(1)
// per node and the whole computation is linear in the total size of the input
// (expected time, via hash tables), rather than O(n / 64) per node for bitsets.
//
// Ties between child pairs (equal objective, within a relative tolerance of
// 1e-12 that absorbs rounding in sums of logarithms) are broken, in order, by
//   (i)  the larger summed count of the two child clades -- the rule
//        phylotreelib applies, and the rule TreeAnnotator X evidently intends
//        (its `breakTies` block never fires; it keeps the first pair in Java
//        HashSet iteration order);
//   (ii) the pair whose child containing the lowest-numbered tip of the clade
//        has the lexicographically smallest sorted vector of tip numbers
//        (a proper prefix counts as smaller).
// Tip numbers follow the label order supplied by the R wrapper, so the result
// does not depend on the order of the input trees.
//
// Reentrant: no file-scope mutable state.

namespace {

using std::uint64_t;
using std::vector;

inline uint64_t splitMix64(uint64_t& state) {
  uint64_t z = (state += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}

struct CladeKey {
  uint64_t h1;
  uint64_t h2;
  bool operator==(const CladeKey& o) const noexcept {
    return h1 == o.h1 && h2 == o.h2;
  }
};

struct CladeKeyHash {
  std::size_t operator()(const CladeKey& k) const noexcept {
    uint64_t x = k.h1 ^ (k.h2 * 0x9E3779B97F4A7C15ULL);
    x ^= x >> 32;
    return static_cast<std::size_t>(x);
  }
};

struct PairHash {
  std::size_t operator()(uint64_t k) const noexcept {
    k ^= k >> 33;
    k *= 0xFF51AFD7ED558CCDULL;
    k ^= k >> 33;
    k *= 0xC4CEB9FE1A85EC53ULL;
    k ^= k >> 33;
    return static_cast<std::size_t>(k);
  }
};

struct CladeTable {
  int nTip = 0;
  int nTree = 0;
  vector<int> count;   // trees containing the clade
  vector<int> size;    // leaves in the clade
  vector<int> minTip;  // lowest tip number (0-based) in the clade
  vector<int> pairParent, pairA, pairB;  // observed child pairs; A holds minTip
};

// Tie-break (ii): compare the sorted tip vectors of clades x and y.
void collectTips(int clade, int nTip, const vector<int>& bestA,
                 const vector<int>& bestB, vector<int>& out) {
  vector<int> stack(1, clade);
  while (!stack.empty()) {
    const int c = stack.back();
    stack.pop_back();
    if (c < nTip) {
      out.push_back(c);
    } else {
      stack.push_back(bestA[c]);
      stack.push_back(bestB[c]);
    }
  }
  std::sort(out.begin(), out.end());
}

bool lexSmaller(int x, int y, int nTip, const vector<int>& bestA,
                const vector<int>& bestB) {
  vector<int> tx, ty;
  collectTips(x, nTip, bestA, bestB, tx);
  collectTips(y, nTip, bestA, bestB, ty);
  return std::lexicographical_compare(tx.begin(), tx.end(),
                                      ty.begin(), ty.end());
}

inline bool nearlyEqual(double a, double b) {
  const double scale = std::max(1.0, std::max(std::fabs(a), std::fabs(b)));
  return std::fabs(a - b) <= 1e-12 * scale;
}

}  // namespace

// Internal: HIPSTR / MrHIPSTR tree from a list of rooted binary trees.
//
// `edges[[i]]` is the ape edge matrix of tree i (any row order); `tipIndex[[i]]`
// maps the tree's tip j to its label number (1..nTip) in the shared label
// order.  The R wrapper has already checked labels; the topology checks
// (rooted, strictly bifurcating, connected) happen here.
//
// Returns the result tree's edge matrix (cladewise order: root nTip + 1,
// internal nodes numbered in the order visited, the child holding the
// lowest-numbered tip first), each internal node's clade count and
// credibility, the log clade credibility, and diagnostics used in validation.
// [[Rcpp::export]]
Rcpp::List hipstrConsensus(const Rcpp::List edges, const Rcpp::List tipIndex,
                           const int nTip, const bool majority) {
  const int nTree = edges.size();
  if (nTree < 1) Rcpp::stop("No trees supplied.");
  if (tipIndex.size() != nTree) Rcpp::stop("`tipIndex` must match `edges`.");
  if (nTip < 2) Rcpp::stop("At least two tips are required.");
  const int nNode = 2 * nTip - 1;  // nodes in a rooted binary tree
  const int nEdge = nNode - 1;

  // Tip hashes, fixed seed.
  vector<CladeKey> tipKey(nTip);
  {
    uint64_t state = 0x48495053545221ULL;  // "HIPSTR!"
    for (int i = 0; i < nTip; ++i) {
      tipKey[i].h1 = splitMix64(state);
      tipKey[i].h2 = splitMix64(state);
    }
  }

  CladeTable tab;
  tab.nTip = nTip;
  tab.nTree = nTree;
  tab.count.assign(nTip, nTree);
  tab.size.assign(nTip, 1);
  tab.minTip.resize(nTip);
  for (int i = 0; i < nTip; ++i) tab.minTip[i] = i;

  std::unordered_map<CladeKey, int, CladeKeyHash> cladeIndex;
  std::unordered_set<uint64_t, PairHash> pairSeen;
  cladeIndex.reserve(static_cast<std::size_t>(nTip) * 4);
  pairSeen.reserve(static_cast<std::size_t>(nTip) * 4);

  // Per-tree scratch, indexed by ape node number (1..nNode).
  vector<int> parentOf(nNode + 1), nKids(nNode + 1), kids(2 * (nNode + 1));
  vector<int> nodeClade(nNode + 1), order;
  vector<CladeKey> nodeKey(nNode + 1);
  order.reserve(nNode);

  for (int t = 0; t < nTree; ++t) {
    if ((t & 0xFF) == 0) Rcpp::checkUserInterrupt();
    const Rcpp::IntegerMatrix edge = Rcpp::as<Rcpp::IntegerMatrix>(edges[t]);
    const Rcpp::IntegerVector tips = Rcpp::as<Rcpp::IntegerVector>(tipIndex[t]);
    if (edge.ncol() != 2 || edge.nrow() != nEdge || tips.size() != nTip) {
      Rcpp::stop("Tree %d is not rooted and fully resolved: rooted, strictly "
                 "bifurcating trees on %d leaves have %d edges, not %d.",
                 t + 1, nTip, nEdge, edge.nrow());
    }
    std::fill(parentOf.begin(), parentOf.end(), 0);
    std::fill(nKids.begin(), nKids.end(), 0);
    for (int r = 0; r < nEdge; ++r) {
      const int p = edge(r, 0);
      const int c = edge(r, 1);
      if (p <= nTip || p > nNode || c < 1 || c > nNode) {
        Rcpp::stop("Tree %d has an invalid edge matrix.", t + 1);
      }
      if (parentOf[c] != 0) {
        Rcpp::stop("Tree %d: node %d has more than one parent.", t + 1, c);
      }
      if (nKids[p] == 2) {
        Rcpp::stop("Tree %d is not fully resolved (node %d has more than two "
                   "children); HIPSTR requires rooted, strictly bifurcating "
                   "trees.", t + 1, p);
      }
      parentOf[c] = p;
      kids[2 * p + nKids[p]++] = c;
    }
    // With nEdge edges and at most two children per internal node, every
    // internal node has exactly two children, and exactly one node lacks a
    // parent; if that node is a tip, the internal nodes contain a cycle.
    int root = 0;
    for (int v = nTip + 1; v <= nNode; ++v) {
      if (parentOf[v] == 0) root = v;
    }
    if (!root) Rcpp::stop("Tree %d has no root.", t + 1);

    // Preorder from the root; reversed, it visits children before parents.
    order.clear();
    order.push_back(root);
    for (std::size_t i = 0; i < order.size(); ++i) {
      const int v = order[i];
      if (v > nTip) {
        order.push_back(kids[2 * v]);
        order.push_back(kids[2 * v + 1]);
      }
    }
    if (static_cast<int>(order.size()) != nNode) {
      Rcpp::stop("Tree %d is not connected.", t + 1);
    }

    for (int i = nNode - 1; i >= 0; --i) {
      const int v = order[i];
      if (v <= nTip) {
        const int label = tips[v - 1] - 1;
        if (label < 0 || label >= nTip) {
          Rcpp::stop("Tree %d: invalid tip index.", t + 1);
        }
        nodeClade[v] = label;
        nodeKey[v] = tipKey[label];
        continue;
      }
      const int k1 = kids[2 * v], k2 = kids[2 * v + 1];
      CladeKey key{nodeKey[k1].h1 + nodeKey[k2].h1,
                   nodeKey[k1].h2 + nodeKey[k2].h2};
      nodeKey[v] = key;
      int a = nodeClade[k1], b = nodeClade[k2];
      auto found = cladeIndex.find(key);
      int id;
      if (found == cladeIndex.end()) {
        id = static_cast<int>(tab.count.size());
        cladeIndex.emplace(key, id);
        tab.count.push_back(0);
        tab.size.push_back(tab.size[a] + tab.size[b]);
        tab.minTip.push_back(std::min(tab.minTip[a], tab.minTip[b]));
      } else {
        id = found->second;
      }
      ++tab.count[id];
      nodeClade[v] = id;
      if (tab.minTip[b] < tab.minTip[a]) std::swap(a, b);
      const uint64_t pairKey = (static_cast<uint64_t>(a) << 32) |
        static_cast<uint64_t>(static_cast<std::uint32_t>(b));
      if (pairSeen.insert(pairKey).second) {
        tab.pairParent.push_back(id);
        tab.pairA.push_back(a);
        tab.pairB.push_back(b);
      }
    }
  }

  const int nClade = static_cast<int>(tab.count.size());
  const int nPair = static_cast<int>(tab.pairParent.size());

  // Clades by increasing size (counting sort), and pairs grouped by parent.
  vector<int> bySize;
  {
    vector<int> start(nTip + 2, 0);
    for (int c = 0; c < nClade; ++c) ++start[tab.size[c] + 1];
    for (int s = 1; s <= nTip + 1; ++s) start[s] += start[s - 1];
    bySize.resize(nClade);
    for (int c = 0; c < nClade; ++c) bySize[start[tab.size[c]]++] = c;
  }
  vector<int> pairStart(nClade + 1, 0), pairList(nPair);
  for (int p = 0; p < nPair; ++p) ++pairStart[tab.pairParent[p] + 1];
  for (int c = 0; c < nClade; ++c) pairStart[c + 1] += pairStart[c];
  {
    vector<int> fill(pairStart.begin(), pairStart.end() - 1);
    for (int p = 0; p < nPair; ++p) pairList[fill[tab.pairParent[p]]++] = p;
  }

  const double logN = std::log(static_cast<double>(nTree));
  vector<double> score(nClade, 0.0);
  vector<int> nMaj(nClade, 0), bestA(nClade, -1), bestB(nClade, -1);
  vector<unsigned char> tiedScore(nClade, 0), tiedCount(nClade, 0);

  for (const int c : bySize) {
    if (tab.size[c] < 2) continue;
    int best = -1;
    int bestMaj = 0;
    double bestScore = 0;
    long long bestCount = 0;
    for (int i = pairStart[c]; i < pairStart[c + 1]; ++i) {
      const int p = pairList[i];
      const int a = tab.pairA[p], b = tab.pairB[p];
      const int m = majority ? nMaj[a] + nMaj[b] : 0;
      const double sc = score[a] + score[b];
      const long long cnt = static_cast<long long>(tab.count[a]) + tab.count[b];
      bool better;
      if (best < 0) {
        better = true;
      } else if (m != bestMaj) {
        better = m > bestMaj;
      } else if (!nearlyEqual(sc, bestScore)) {
        better = sc > bestScore;
      } else if (cnt != bestCount) {
        better = cnt > bestCount;
      } else {
        better = lexSmaller(a, tab.pairA[best], nTip, bestA, bestB);
      }
      if (better) {
        best = p;
        bestMaj = m;
        bestScore = sc;
        bestCount = cnt;
      }
    }
    // Diagnostics: was the choice tied on the objective / also on counts?
    int nTiedScore = 0, nTiedCount = 0;
    for (int i = pairStart[c]; i < pairStart[c + 1]; ++i) {
      const int p = pairList[i];
      const int a = tab.pairA[p], b = tab.pairB[p];
      const int m = majority ? nMaj[a] + nMaj[b] : 0;
      if (m == bestMaj && nearlyEqual(score[a] + score[b], bestScore)) {
        ++nTiedScore;
        if (static_cast<long long>(tab.count[a]) + tab.count[b] == bestCount) {
          ++nTiedCount;
        }
      }
    }
    tiedScore[c] = nTiedScore > 1;
    tiedCount[c] = nTiedCount > 1;
    bestA[c] = tab.pairA[best];
    bestB[c] = tab.pairB[best];
    score[c] = std::log(static_cast<double>(tab.count[c])) - logN + bestScore;
    nMaj[c] = bestMaj + ((majority && 2LL * tab.count[c] > nTree) ? 1 : 0);
  }

  // Root clade: the clade containing every tip (the last, largest clade).
  const int rootClade = bySize.back();

  // Read the tree off the argmax pairs, depth first, child holding the lowest
  // tip first; internal nodes numbered in visiting order (cladewise).
  Rcpp::IntegerMatrix outEdge(nEdge, 2);
  Rcpp::IntegerVector outCount(nTip - 1);
  Rcpp::NumericVector outCred(nTip - 1);
  int nextNode = nTip + 1, row = 0;
  int nTieScore = 0, nTieCount = 0, nMajOut = 0;
  double logCred = 0;
  vector<std::pair<int, int>> stack;  // (clade, parent node)
  stack.emplace_back(rootClade, 0);
  while (!stack.empty()) {
    const int c = stack.back().first;
    const int parent = stack.back().second;
    stack.pop_back();
    int node;
    if (c < nTip) {
      node = c + 1;
    } else {
      node = nextNode++;
      outCount[node - nTip - 1] = tab.count[c];
      outCred[node - nTip - 1] = static_cast<double>(tab.count[c]) / nTree;
      logCred += std::log(static_cast<double>(tab.count[c])) - logN;
      nTieScore += tiedScore[c];
      nTieCount += tiedCount[c];
      nMajOut += 2LL * tab.count[c] > nTree;
      stack.emplace_back(bestB[c], node);
      stack.emplace_back(bestA[c], node);
    }
    if (parent) {
      outEdge(row, 0) = parent;
      outEdge(row, 1) = node;
      ++row;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("edge") = outEdge,
    Rcpp::Named("count") = outCount,
    Rcpp::Named("credibility") = outCred,
    Rcpp::Named("logCredibility") = logCred,
    Rcpp::Named("score") = score[rootClade],
    Rcpp::Named("nMajority") = nMajOut,
    Rcpp::Named("nClade") = nClade - nTip,
    Rcpp::Named("nPair") = nPair,
    Rcpp::Named("nTiedScore") = nTieScore,
    Rcpp::Named("nTiedCount") = nTieCount
  );
}
