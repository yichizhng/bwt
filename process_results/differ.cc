#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char **argv) {
  if (argc != 3) {
    cerr << "Usage: " << argv[1] << " file1 file2" << endl;
  }
  ifstream file1(argv[1]), file2(argv[2]);
  int i = 0;
  int pos1, pos2;
  while (file1 && file2) {
    ++i;
    file1 >> pos1;
    file2 >> pos2;
    if (pos1 != pos2)
      cout << i << ' ' << pos1 << ' ' << pos2 << endl;
  }
}
