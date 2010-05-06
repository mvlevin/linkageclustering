#ifndef TEXT_FILE_READER_H_
#define TEXT_FILE_READER_H_

#include <fstream>
#include <string>

using std::ifstream;
using std::string;

namespace LinkageClustering {

class TextFileReader {
 public:
   TextFileReader(const string& filename) : filename_(filename), file_stream_(filename_.c_str(), ifstream::in) {
    if (!file_stream_.is_open()) {
      cerr << "Couldn't open file " << filename_ << endl;
      abort();
    }
  }

  ~TextFileReader() {
    file_stream_.close();
  }

  template<class ValueType> 
  ValueType Read() {
    if (!file_stream_.good()) {
      cerr << "Can't read from file " << filename_ << endl;
      abort();
    }
    ValueType value;
    file_stream_ >> value;
    return value;
  }
  
 private:
  string filename_;
  ifstream file_stream_;
};

}  // namespace LinkageClustering

#endif  // TEXT_FILE_READER_H_