// Copyright 2010 Michael Levin
//
// Example class representing a sparse
// square dissimilarity matrix.

#include <vector>

using std::vector;

// Represents a sparse square dissimilarity matrix.
//
// Stored as a vector of rows.
// Each row is stored as a vector of pairs (index, dissimilarity)
// for all such indices that dissimilarity is non-zero.
//
// param DissimilarityValueType_ - type to use for storage of matrix
//       cells value. Can be, for example, size_t, unsigned long long
//       or double.
template<class DissimilarityValueType_>
class MemoryBasedDissimilarityMatrix {
 public:
  // The following typedefs, iterator class and methods are required to use the GetLinkageSequence
  // algorithm with a dissimilarity matrix.
  // See compute_linkage_sequence.h for details.
  typedef typename size_t IndexType;
  typedef typename DissimilarityValueType_ DissimilarityValueType;

  // Represents an iterator over a row of the matrix.
  // 
  // Fulfills the interface needed for a DissimilarityIterator 
  // in the GetLinkageSequence algorithm.
  // See compute_linkage_sequence.h for details.
  class RowIterator {
   public:
    typedef typename vector<pair<IndexType, DissimilarityValueType> >::const_iterator Iterator;

    // Initialize with iterators to begin and end of the row's underlying vector.
    RowIterator(Iterator begin, Iterator end) : iterator_(begin), end_(end) {}

    // Returns the index of the current column of the matrix.
    IndexType GetToIndex() const { return iterator_->first; }

    // Returns the value of the current cell of the matrix.
    DissimilarityValueType GetDissimilarity() const { return iterator_->second; }

    // Returns true if all the non-zero elements are passed.
    bool Done() const { return iterator_ == end_; }

    // Moves to the next non-zero element in the row.
    void Next() { ++iterator_; }
   private:
    Iterator iterator_;
    Iterator end_;
  };

  MemoryBasedDissimilarityMatrix(
    const vector< vector<pair<IndexType, DissimilarityValueType> > >& sparse_matrix)
  : sparse_matrix_(sparse_matrix) {}

  // Returns number of rows in the matrix.
  IndexType GetRowCount() const { return sparse_matrix_.size(); }

  // Returns 
  RowIterator GetRowIterator(IndexType index) const { 
    return RowIterator(
        sparse_matrix_[index].begin(),
        sparse_matrix_[index].end());
  }

 private:
  // The underlying sparse matrix is stored as a vector of sparse rows.
  // Each row consists of an ordered by index list of pairs (index, dissimilarity)
  // for each such index that dissimilarity is non-zero.
  vector< vector<pair<IndexType, DissimilarityValueType> > > sparse_matrix_;
};
