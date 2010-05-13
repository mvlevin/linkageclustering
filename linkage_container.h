// Copyright 2010 Michael Levin
//
// Data structure to effectively store a subset of a
// set of elements along with their linkages to the
// current set. Provides fast way to add an element
// along with its linkage value, to reduce the linkage
// value of an element, to find the element with the
// minimum linkage value and to remove this element
// from the current subset.
#ifndef LINKAGE_CONTAINER_H_
#define LINKAGE_CONTAINER_H_

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

}  // namespace LinkageClustering

#endif  // LINKAGE_CONTAINER_H_

