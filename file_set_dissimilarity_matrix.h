// Copyright 2010 Michael Levin
//
// Example class representing a
// square dissimilarity matrix as a set of files.

#ifndef FILE_SET_DISSIMILARITY_MATRIX_H_
#define FILE_SET_DISSIMILARITY_MATRIX_H_

#include <string>
#include <vector>

using std::string;
using std::vector;

// Represents a square dissimilarity matrix in a set of files.
//
// Stored as a vector of file names. Each file contains
// the description of one row of the matrix, and can be read
// correctly by a FileReader.
//
// File must be organized as
// <number of non-zero elements> <index> <dissimilarity> <index> <dissimilarity>
// in some format, understood by the FileReader template argument provided.
//
// param DissimilarityValueType_ - type to use for storage of matrix
//       cells value. Can be, for example, size_t, unsigned long long
//       or double.
// 
// param FileReader - class that reads the files representing matrix's rows.
template<class DissimilarityValueType_, class FileReader>
class FileSetDissimilarityMatrix {
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
  //
  // Doesn't store the whole row in memory, reads the data from file on request.
  class RowIterator {
   public:
    // Initialize with a pointer to a FileReader object.
    RowIterator(const string& filename) {
      file_reader_ = new FileReader(filename);
      records_count_ = file_reader_->Read<IndexType>();
      records_read_ = 0;
      ReadNextPair();
    }

    ~RowIterator() {
      delete file_reader_;
    }

    // Returns the index of the current column of the matrix.
    IndexType GetToIndex() const { return index_; }

    // Returns the value of the current cell of the matrix.
    DissimilarityValueType GetDissimilarity() const { return dissimilarity_; }

    // Returns true if all the non-zero elements are passed.
    bool Done() const { return records_read_ > records_count_; }

    // Moves to the next non-zero element in the row.
    void Next() { 
      ReadNextPair();
    }
   private:
    // Reads next pair (index, dissimilarity) from the file.
    void ReadNextPair() {
      if (records_read_ < records_count_) {
        index_ = file_reader_->Read<IndexType>();
        dissimilarity_ = file_reader_->Read<DissimilarityValueType>();
      }
      ++records_read_;
    }

    FileReader* file_reader_;
    IndexType records_count_;
    IndexType records_read_;
    IndexType index_;
    DissimilarityValueType dissimilarity_;
  };

  FileSetDissimilarityMatrix(
    const vector<string>& filenames)
  : filenames_(filenames) {}

  // Returns number of rows in the matrix.
  IndexType GetRowCount() const { return filenames_.size(); }

  // Returns iterator over one row of the matrix. Uses file-based iterator.
  RowIterator GetRowIterator(IndexType index) const { 
    return RowIterator(filenames_[index]);
  }

 private:
  // Vector of file names containing rows of the dissimilarity matrix.
  vector<string> filenames_;
};

#endif  // FILE_SET_DISSIMILARITY_MATRIX_H_
