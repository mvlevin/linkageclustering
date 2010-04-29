// Copyright Michael Levin 2010.
//
// Example program that uses compute_linkage_sequence.h.
//
// Takes input from file "input.txt" in format
//
// n = <Number of elements>
// k_0( = Number of i such that i >= 0, dissimilarity_{0, i) != 0) index_01 dissimilarity_{0, index_01} .. index_0k0 dissimilarity_{0, index_0k0}
// k_1( = Number of i such that i >= 1, dissimilarity_{1, i} != 0) index_11 dissimilarity_{1, index_11} .. index_1k1 dissimilarity_{1, index_1k1}
// \dots
// k_{n - 1}( = Number of i such that i >= n - 1, dissimilarity_{n - 1, i} != 0) index_n-10 dissimilarity_{n-1, index_n-10} .. index_n-1kn-1 dissimilarity_{n-1, index_n-1kn-1}
//
// Outputs sequence like
//    i = 8 F_i = 1
//    i = 5 F_i = 1
//    i = 2 F_i = 1
//    i = 0 F_i = 2
//    i = 1 F_i = 2
//    i = 3 F_i = 2
//    i = 6 F_i = 1
//    i = 4 F_i = 1
//    i = 7 F_i = 0

#include "compute_linkage_sequence.h"

#include <cassert>
#include <cstdio>
#include <iostream>
#include <map>
#include <set>
#include <utility>
#include <vector>


using std::cout;
using std::endl;
using std::pair;
using std::vector;

int main() {
  freopen("input.txt", "r", stdin);
  freopen("output.txt", "w", stdout);
  size_t element_count;
  scanf("%u", &element_count);
  vector< vector< pair<size_t, long long> > > dissimilarities(element_count);
  for (size_t element_index = 0; element_index < element_count; ++element_index) {
    size_t dissimilarities_count;
    scanf("%d", &dissimilarities_count);
    for (size_t dissimilarity_index = 0;
         dissimilarity_index < dissimilarities_count;
         ++dissimilarity_index) {
      size_t to;
      size_t dissimilarity;
      scanf("%u%u", 
            &to,
            &dissimilarity);
      dissimilarities[element_index].push_back(make_pair(to, dissimilarity));
      dissimilarities[to].push_back(make_pair(element_index, dissimilarity));
    }
  }
  vector< pair<size_t, long long> > linkage_sequence = GetLinkageSequence(dissimilarities);
  for (size_t index = 0; index < linkage_sequence.size(); ++index) {
    cout << "i = " << linkage_sequence[index].first 
         << " F_i = " << linkage_sequence[index].second << endl;
  }
  return 0;
}
