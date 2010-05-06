// Copyring 2010 Michael Levin
//
// Implementation of Linkage Clustering procedure
// based on "MULTIPARTITE GRAPH CLUSTERING FOR
//           STRUCTURED DATASETS AND AUTOMATING
//           ORTHOLOG EXTRACTION"
// by Akshay Vashisht.
//
// Assumes linkage function L(i, X) is based on 
// the dissimilarity matrix a_ij.
//
// Instead of the best cluster, produces sequence
// (i_1, F_{i_1}), (i_2, F_{i_2}), \dots, (i_n, F_{i_n})
// of the elements removed in the best possible order
// and the linkage values of the sets remaining after
// each of the consecutive removes.
//
// The running time of the algorithm is O(M + n log n),
// where M is the number of non-zero elements in a_ij,
// and n is the the number of elements in the whole set
// V.
#ifndef LINKAGE_SEQUENCE_COMPUTER_H_
#define LINKAGE_SEQUENCE_COMPUTER_H_

#include "linkage_container.h"

#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <utility>
#include <vector>

using std::cerr;
using std::endl;
using std::make_pair;
using std::map;
using std::pair;
using std::set;
using std::vector;

namespace LinkageClustering {

// Computes the linkage sequence for a set provided a dissimilarity
// matrix.
//
// The DissimilarityMatrix class must define:
//   The following types:
//
//     IndexType - The type to represent an index of an element in the set.
//                 Can be size_t, long long, double, etc.
//     DissimilarityValueType - The type to represent an element of the matrix.
//                              Can be size_t, long long, double, etc.
//     RowIterator - The type to pass one-way through a row of the matrix, with
//                   the following interface:
//         IndexType GetToIndex() - returns the index of the current column.
//         DissimilarityValueType GetDissimilarity() - returns the value of the current matrix cell.
//         bool Done() - returns true if we passed through the whole row.
//         void Next() - moves iterator to the next column.
//
//   And the following interface:
//
//     IndexType GetRowsCount() - returns the number of rows in the matrix
//     RowIterator GetRowIterator(IndexType index) - returns the RowIterator for index-th row
//
// The LinkageValueType is the type to store the values of linkages L(i, X) of an
// element to a set. It can be a different type from the type for matrix elements
// because the linkage values may be significantly bigger than dissimilarity values.
// For example, if DissimilarityValueType is unsigned int32, then LinkageValueType may be
// unsigned int64 to fit the maximum possible value of the linkage function.
// It is required that DissimilarityValueType can be converted to LinkageValueType automatically
// and without loss of precision.
template<class DissimilarityMatrix, class LinkageValueType>
class LinkageSequenceComputer {
 public:
  typedef typename DissimilarityMatrix::IndexType IndexType;
  typedef typename DissimilarityMatrix::DissimilarityValueType DissimilarityValueType;
  typedef typename DissimilarityMatrix::RowIterator RowIterator;

  // Returns sequence { (i_1, F_{i_1}), (i_2, F_{i_2}), \dots,
  // (i_n, F_{i_n}) }, where i_1 is the element with the least
  // linkage L(i_1, V), F_{i_1} is the linkage of V \ {i_1},
  // i_2 is the element with the least linkage L(i_2, V \ {i_1}),
  // F_{i_2} is the linkage of V \ {i_1} \ {i_2}, and so on.
  static vector< pair<IndexType, LinkageValueType> > GetLinkageSequence(
      const DissimilarityMatrix& dissimilarity_matrix);
 private:
  // Initializes linkage_container by adding all elements of V,
  // along with their linkages L(i, V) to the whole set.
  static void InitializeLinkageContainerWithWholeSet(
      const DissimilarityMatrix& dissimilarity_matrix,
      LinkageContainer<IndexType, LinkageValueType>* linkage_container) {
    for (size_t index = 0; index < dissimilarity_matrix.GetRowCount(); ++index) {
      LinkageValueType linkage = 0;
      for (RowIterator dissimilarity_iterator = 
               dissimilarity_matrix.GetRowIterator(index);
           !dissimilarity_iterator.Done();
           dissimilarity_iterator.Next()) {
        linkage += dissimilarity_iterator.GetDissimilarity();
      }
      linkage_container->AddElement(index, linkage);
    }
  }
};

template<class DissimilarityMatrix, class LinkageValueType>
vector< pair<typename DissimilarityMatrix::IndexType, LinkageValueType> > 
LinkageSequenceComputer<DissimilarityMatrix, LinkageValueType>::GetLinkageSequence(
    const DissimilarityMatrix& dissimilarities) {
  LinkageContainer<IndexType, LinkageValueType> linkage_container;
  InitializeLinkageContainerWithWholeSet(dissimilarities, &linkage_container);
  vector< pair<IndexType, LinkageValueType> > linkage_sequence;
  for (IndexType step = 0; step < dissimilarities.GetRowCount(); ++step) {
    LinkageValueType min_linkage = linkage_container.GetMinLinkage();
    IndexType removed_element_index = linkage_container.RemoveMinLinkageElement();
    linkage_sequence.push_back(make_pair(removed_element_index, min_linkage));
    for (RowIterator dissimilarity_iterator = 
             dissimilarities.GetRowIterator(removed_element_index);
         !dissimilarity_iterator.Done();
         dissimilarity_iterator.Next()) {
      IndexType to_index = dissimilarity_iterator.GetToIndex();
      LinkageValueType dissimilarity = dissimilarity_iterator.GetDissimilarity();
      if (linkage_container.Contains(to_index)) {
        linkage_container.DecreaseLinkage(to_index, dissimilarity);
      }
    }    
  }
  return linkage_sequence;
}

}  // namespace LinkageClustering

#endif  // LINKAGE_SEQUENCE_COMPUTER_H_
