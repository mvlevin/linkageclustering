#include <cmath>
#include <cstdio>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

using namespace std;

size_t AddName(const string& name, map<string, size_t>* ids, size_t* id_count, FILE* ids_file) {
  map<string, size_t>::const_iterator it = ids->find(name);
  if (it == ids->end()) {
    size_t new_id = *id_count;
    (*ids)[name] = new_id;
    ++(*id_count);
    fprintf(ids_file, "%s %u\n", name.c_str(), new_id);
    return new_id;
  } else {
    return it->second;
  }
}

void ConvertToAdjLists(
    const char* input_filename, 
    map<string, size_t>* ids, 
    size_t* id_count, 
    vector< vector< pair<size_t, size_t> > >* adj_lists) {
  freopen(input_filename, "r", stdin);
  FILE* ids_file = fopen("ids.txt", "w");
  if (!ids_file) {
    cerr << "Couldn't open file ids.txt to write down ids" << endl;
    abort();
  }
  *id_count = 0;
  while (!feof(stdin)) {
    char buffer[1000];
    // From
    scanf("%s", buffer);
    size_t from = AddName(buffer, ids, id_count, ids_file);
    // To
    scanf("%s", buffer);
    size_t to = AddName(buffer, ids, id_count, ids_file);
    double initial_weight;
    scanf("%lf\n", &initial_weight);
    size_t weight = floor(initial_weight * 10000 + 1e-8);
    if (adj_lists->size() < *id_count) {
      adj_lists->resize(*id_count);
    }
    adj_lists->at(from).push_back(make_pair(to, weight));
    adj_lists->at(to).push_back(make_pair(from, weight));
  }
  fclose(ids_file);
}

void WriteDown(const char* filename, const vector< vector< pair<size_t, size_t> > >& adj_lists) {
  FILE* output = fopen(filename, "w");
  if (!output) {
    cerr << "Couldn't open file " << filename << " to write down the adjacency lists" << endl;
    abort();
  }
  fprintf(output, "%u\n", adj_lists.size());
  for (size_t row = 0; row < adj_lists.size(); ++row) {
    fprintf(output, "%u", adj_lists[row].size());
    for (size_t column = 0; column < adj_lists[row].size(); ++column) {
      fprintf(output, " %u %u", adj_lists[row][column].first, adj_lists[row][column].second);
    }
    fprintf(output, "\n");
  }
  fclose(output);
}

int main(int argc, char** argv) {
  if (argc != 3) {
    cerr << "Usage: convert_to_adj_lists input_filename output_filename" << endl;
    return -1;
  }
  map<string, size_t> ids;
  size_t id_count = 0;
  vector< vector< pair<size_t, size_t> > > adj_lists;
  ConvertToAdjLists(argv[1], &ids, &id_count, &adj_lists);
  WriteDown(argv[2], adj_lists);
  return 0;
}