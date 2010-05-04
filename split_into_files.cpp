#include <cmath>
#include <cstdio>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

using namespace std;

void AddName(const string& name, map<string, size_t>* ids, size_t* id_count) {
  if (ids->find(name) == ids->end()) {
    (*ids)[name] = *id_count;
    ++(*id_count);
  }
}

void ComputeIdsAndDegrees(
    const char* input_filename, 
    map<string, size_t>* ids, 
    size_t* id_count, 
    map<string, size_t>* degrees) {
  freopen(input_filename, "r", stdin);
  *id_count = 0;
  while (!feof(stdin)) {
    char buffer[1000];
    // From
    scanf("%s", buffer);
    AddName(buffer, ids, id_count);
    ++(*degrees)[buffer];
    // To
    scanf("%s", buffer);
    AddName(buffer, ids, id_count);
    ++(*degrees)[buffer];
    // Edge weight
    scanf("%s\n", buffer);
  }
}

size_t GetId(const map<string, size_t>& ids) {
  char buffer[1000];
  scanf("%s", buffer);
  map<string, size_t>::const_iterator it = ids.find(buffer);
  if (it == ids.end()) {
    cerr << "Couldn't find name " << buffer << " in ids" << endl;
    abort();
  } else {
    return it->second;
  }
}

void WriteToFiles(
    const char* input_filename, 
    const string& output_basename, 
    const map<string, size_t>& ids, 
    size_t id_count,
    const map<string, size_t>& degrees) {
  freopen(input_filename, "r", stdin);
  vector<FILE*> output_files(id_count);
  for (int i = 0; i < id_count; ++i) {
    ostringstream out;
    out << output_basename << "-" << i << ".txt";
    output_files[i] = fopen(out.str().c_str(), "w");
    if (!output_files[i]) {
      cerr << "Couldn't open file " << out.str() << " for writing" << endl;
      abort();
    }
  }
  for (map<string, size_t>::const_iterator ids_it = ids.begin(), degrees_it = degrees.begin();
       ids_it != ids.end() && degrees_it != degrees.end();
       ++ids_it, ++degrees_it) {
    fprintf(output_files[ids_it->second], "%u", degrees_it->second);
  }
  while (!feof(stdin)) {
    size_t from = GetId(ids);
    if (from >= id_count) {
      cerr << "from is " << from << ", but must be from 0 to " << id_count << endl;
      abort();
    }
    size_t to = GetId(ids);
    if (to >= id_count) {
      cerr << "to is " << to << ", but must be from 0 to " << id_count << endl;
      abort();
    }
    double initial_weight;
    scanf("%lf\n", &initial_weight);
    size_t weight = floor(initial_weight * 10000 + 1e-8);
    fprintf(output_files[from], " %u %u", to, weight);
    fprintf(output_files[to], " %u %u", from, weight);
  }
  for (int i = 0; i < id_count; ++i) {
    fclose(output_files[i]);
  }
}

int main(int argc, char** argv) {
  if (argc != 3) {
    cerr << "Usage: split_into_files input_filename basename" << endl;
    return -1;
  }
  map<string, size_t> ids;
  size_t id_count = 0;
  map<string, size_t> degrees;
  ComputeIdsAndDegrees(argv[1], &ids, &id_count, &degrees);
  WriteToFiles(argv[1], argv[2], ids, id_count, degrees);
  return 0;
}