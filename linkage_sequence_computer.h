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

// An auxilliary data structure to store elements
// of the current set along with their linkages
// to the current set, and retrieve and delete
// the element with the lowest linkage value fast.
template<class IndexType, class LinkageValueType>
class LinkageContainer;

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

// Contains subset X of the initial set V of elements, along
// with linkage L(i, X) of each element i of the current
// subset X to the whole subset. In this implementation,
// it is assumed that L(i, X) = sum for j in X a_ij, where
// a_ij is the dissimilarity between i and j.
// 
// Initialize it with the whole initial set V and correct
// linkage values L(i, V) for each i in V.
//
// After initialization, LinkageContainer is able to 
//   * find the element with the least linkage value
//   * remove the element with the least linkage value
//   * check if current subset contains a given element
//   * decrease the linkage of a given element by some value
// Each of the operations runs in O(log n) time, where |X| = n.
//
// LinkageValueType is the type in which to accumulate linkages.
//
// !!!IMPORTANT!!!
// Use either long long or double for the LinkageValueType if you have a big
// matrix. Don't use int (int32), as if there are
// a lot of elements, even if each a_ij fits in int32, L(i, X) can
// be out of the int32 range.
template<class IndexType, class LinkageValueType>
class LinkageContainer {
public:
  LinkageContainer() {}

  // Adds element index with L(index, V) = linkage to the current subset X.
  //
  // Aborts on trying to add the same element twice.
  void AddElement(IndexType index, LinkageValueType linkage);

  // Returns true if the current subset X is empty.
  bool Empty() const;

  // Returns minimum linkage value L(i, X) for i in the current subset X.
  //
  // Aborts if called when X is empty.
  LinkageValueType GetMinLinkage() const;

  // Removes the element i with the minimum linkage value L(i, X) to the
  // current subset X.
  //
  // Returns the index of i (from 0 to n - 1).
  //
  // Aborts if called when X is empty.
  IndexType RemoveMinLinkageElement();

  // Returns true if the current subset X contains element index.
  bool Contains(IndexType index) const;

  // Decreases the linkage value of element index.
  //
  // Aborts if element index is not in the current subset X.
  void DecreaseLinkage(IndexType index, LinkageValueType decrease);
private:
  // Contains pairs (L(i, X), i) for each element i in the current subset X.
  set< pair<LinkageValueType, IndexType> > elements_;

  // Contains map from i to L(i, X) for each element i in the current subset X.
  map<IndexType, LinkageValueType> linkage_;
};

template<class IndexType, class LinkageValueType>
void LinkageContainer<IndexType, LinkageValueType>::AddElement(
    IndexType index, 
    LinkageValueType linkage) {
  if (linkage_.find(index) != linkage_.end()) {
    cerr << "Trying to add an existing element to the current subset" << endl;
    abort();
  }
  elements_.insert(make_pair(linkage, index));
  linkage_[index] = linkage;
  assert(elements_.size() == linkage_.size());
}

template<class IndexType, class LinkageValueType>
bool LinkageContainer<IndexType, LinkageValueType>::Empty() const {
  return elements_.empty();
}

template<class IndexType, class LinkageValueType>
LinkageValueType LinkageContainer<IndexType, LinkageValueType>::GetMinLinkage() const {
  if (elements_.empty()) {
    cerr << "Can't get min linkage of an empty set" << endl;
    abort();
  }
  return elements_.begin()->first;
}

template<class IndexType, class LinkageValueType>
IndexType LinkageContainer<IndexType, LinkageValueType>::RemoveMinLinkageElement() {
  if (elements_.empty()) {
    cerr << "Can't remove min linkage element from an empty set" << endl;
    abort();
  }
  IndexType index = elements_.begin()->second;
  linkage_.erase(index);
  elements_.erase(elements_.begin());
  assert(linkage_.size() == elements_.size());
  return index;
}

template<class IndexType, class LinkageValueType>
bool LinkageContainer<IndexType, LinkageValueType>::Contains(IndexType index) const {
  if (linkage_.find(index) == linkage_.end()) {
    return false;
  } else {
    return true;
  }
}

template<class IndexType, class LinkageValueType>
void LinkageContainer<IndexType, LinkageValueType>::DecreaseLinkage(IndexType index, LinkageValueType decrease) {
  if (linkage_.find(index) == linkage_.end()) {
    cerr << "There is no element " 
         << index 
         << " in the current set" 
         << endl;
    abort();
  }
  if (decrease > linkage_[index]) {
    cerr << "Trying to make linkage negative" << endl;
    abort();
  }
  pair<LinkageValueType, IndexType> element = make_pair(linkage_[index], index);
  if (elements_.find(element) == elements_.end()) {
    cerr << "There is no element with index = " 
         << index 
         << " in the current set" 
         << endl;
    abort();
  }
  elements_.erase(element);
  linkage_[index] -= decrease;
  elements_.insert(make_pair(linkage_[index], index));
  assert(linkage_.size() == elements_.size());
}

template<class DissimilarityMatrix, class LinkageValueType>
vector< pair<typename DissimilarityMatrix::IndexType, LinkageValueType> > 
LinkageSequenceComputer<DissimilarityMatrix, LinkageValueType>::GetLinkageSequence(
    const DissimilarityMatrix& dissimilarities) {
  LinkageContainer<IndexType, LinkageValueType> linkage_container;
  InitializeLinkageContainerWithWholeSet(dissimilarities, &linkage_container);
  vector< pair<IndexType, LinkageValueType> > linkage_sequence;
  for (size_t step = 0; step < dissimilarities.GetRowCount(); ++step) {
    LinkageValueType min_linkage = linkage_container.GetMinLinkage();
    size_t removed_element_index = linkage_container.RemoveMinLinkageElement();
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

#endif  // LINKAGE_SEQUENCE_COMPUTER_H_
