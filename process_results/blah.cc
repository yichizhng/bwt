#include <iostream>

using namespace std;

int blah(int &a) {
  a = 5;
  return 9;
}

int bleh (int &a) {
  a = 7;
  return 8;
}

int main() {
  int x = 10;
  (x = blah(x)) = bleh(x);
  cerr << x << endl;
  return 0;
}
