#include <iostream>
#include <fstream>
#include <map>

using namespace std;

int main(int argc, char **argv) {
  // Read filename from argv
  if (argc != 2) {
    cerr << "Usage: " << argv[0] << " idxfile";
  }
  map<unsigned long long, unsigned long long> mp;
  unsigned long long i;
  ifstream f(argv[1]);
  while (f) {
    f >> i;
    mp[i]++;
  }
  i = 0;
  for (auto it : mp) {
    if (it.second > 10) {
      if (it.first - i > 1000)
	cout << endl;
      cout << it.first << '\t' << it.second << endl;
      i = it.first;
    }
  }
  return 0;
}
