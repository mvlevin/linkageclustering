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

#include "linkage_sequence_computer.h"
#include "file_set_dissimilarity_matrix.h"
#include "text_file_reader.h"

#include <cassert>
#include <cstdio>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <utility>
#include <vector>


using std::cout;
using std::endl;
using std::ostringstream;
using std::pair;
using std::vector;
using LinkageClustering::FileSetDissimilarityMatrix;
using LinkageClustering::LinkageSequenceComputer;

FileSetDissimilarityMatrix<size_t, TextFileReader> InitializeMatrix(const string& basename, int rowcount) {
  vector<string> filenames;
  for (int i = 0; i < rowcount; ++i) {
    ostringstream out;
    out << basename << "-" << i << ".txt";
    filenames.push_back(out.str());
  }
  FileSetDissimilarityMatrix<size_t, TextFileReader> matrix(filenames);
  return matrix;
}

int main(int argc, char** argv) {
  if (argc != 3) {
    cerr << "Usage: calculate_using_files basename rowcount" << endl;
    return -1;
  }
  string basename = argv[1];
  size_t rowcount;
  int result = sscanf(argv[2], "%u", &rowcount);
  if (result != 1) {
    cerr << "rowcount must be unsigned int" << endl;
    return -2;
  }
  freopen("output.txt", "w", stdout);
  vector< pair<size_t, long long> > linkage_sequence = 
    LinkageSequenceComputer<FileSetDissimilarityMatrix<size_t, TextFileReader>, long long>::GetLinkageSequence(InitializeMatrix(basename, rowcount));
  for (size_t index = 0; index < linkage_sequence.size(); ++index) {
    cout << "i = " << linkage_sequence[index].first 
         << " F_i = " << linkage_sequence[index].second << endl;
  }
  return 0;
}
