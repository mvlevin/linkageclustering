#include "linkage_sequence_computer.h"
#include "memory_based_dissimilarity_matrix.h"

#include <cassert>
#include <cstdio>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <utility>
#include <vector>

using std::cout;
using std::endl;
using std::numeric_limits;
using std::pair;
using std::vector;

const int MAX_ELEMENT_COUNT = 100;
const int MAX_DISSIMILARITY = 1000000000;

map<size_t, long long> GetLinkageValues(    
    const set<int>& current_set,
    const vector< vector< pair<size_t, size_t> > >& dissimilarities) {
  map<size_t, long long> result;
  for (set<int>::const_iterator it = current_set.begin();
       it != current_set.end();
       ++it) {
    int element_index = *it;
    result[element_index] = 0;
    for (int index = 0; index < dissimilarities[element_index].size(); ++index) {
      int other_element_index = dissimilarities[element_index][index].first;
      if (current_set.find(other_element_index) != current_set.end()) {
        result[element_index] += dissimilarities[element_index][index].second;
      }
    }
  }
  return result;
}

pair<size_t, long long> GetMinPair(const map<size_t, long long>& linkage_values) {
  pair<size_t, long long> result(
      numeric_limits<size_t>::max(), 
      numeric_limits<long long>::max());
  for (map<size_t, long long>::const_iterator it = linkage_values.begin();
       it != linkage_values.end();
       ++it) {
    if (it->second < result.second) {
      result = *it;
    }
  }
  return result;
}

vector< pair<size_t, long long> > GetLinkageSequenceSimple(
    const vector< vector< pair<size_t, size_t> > >& dissimilarities) {
  set<int> current_set;
  for (int i = 0; i < dissimilarities.size(); ++i) {
    current_set.insert(i);
  }
  map<size_t, long long> linkage_values = GetLinkageValues(current_set, dissimilarities);
  vector< pair<size_t, long long> > linkage_sequence;
  for (int step = 0; step < dissimilarities.size(); ++step) {
    pair<size_t, long long> min_pair = GetMinPair(linkage_values);
    linkage_sequence.push_back(min_pair);
    current_set.erase(min_pair.first);
    linkage_values = GetLinkageValues(current_set, dissimilarities);
  }
  return linkage_sequence;
}

int main() {
  srand(239);
  while (true) {
    int element_count = rand() % MAX_ELEMENT_COUNT + 1;
    vector< vector< pair<size_t, size_t> > > dissimilarities(element_count);
    for (int element_index = 0; element_index < element_count; ++element_index) {
      for (int column_index = element_index + 1; column_index < element_count; ++column_index) {
        if (rand() % 5 == 0) {
          int dissimilarity = rand() % MAX_DISSIMILARITY + 1;
          dissimilarities[element_index].push_back(make_pair(column_index, dissimilarity));
          dissimilarities[column_index].push_back(make_pair(element_index, dissimilarity));
        }
      }
    }
    vector< pair<size_t, long long> > linkage_sequence = 
        LinkageSequenceComputer<MemoryBasedDissimilarityMatrix<size_t>, long long>::GetLinkageSequence(
            MemoryBasedDissimilarityMatrix<size_t>(dissimilarities));
    vector< pair<size_t, long long> > linkage_sequence_correct = GetLinkageSequenceSimple(dissimilarities);
    if (linkage_sequence == linkage_sequence_correct) {
      cerr << "OK" << endl;
    } else {
      cerr << "Wrong answer!" << endl;
      cerr << dissimilarities.size() << endl;
      for (int i = 0; i < dissimilarities.size(); ++i) {
        cerr << dissimilarities[i].size();
        for (int j =  0; j < dissimilarities[i].size(); ++j) {
          cerr << " (" << dissimilarities[i][j].first << ", " << dissimilarities[i][j].second << ")";
        }
        cerr << endl;
      }
      for (int i = 0; i < linkage_sequence.size(); ++i) {
        cerr << "(" << linkage_sequence[i].first << ", " << linkage_sequence[i].second << ") ";
      }
      cerr << endl;
      for (int i = 0; i < linkage_sequence_correct.size(); ++i) {
        cerr << "(" << linkage_sequence_correct[i].first << ", " << linkage_sequence_correct[i].second << ") ";
      }
      cerr << endl;
      break;
    }
  }
  return 0;
}
