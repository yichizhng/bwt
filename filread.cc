// Strips away the stuff I don't feel like reading
// Replaces N with G so that matches
// are not done against them (unless you somehow
// found a nucleotide sequence containing a lot of G's
// in a row; unlikely)

// Reminder to self: HG19 is found at
// http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/

#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <limits>
#include <cctype>


using namespace std;

int main(int argc, char **argv) {
  if (argc != 3) {
    cerr << "Usage: " << argv[0] << " seq.fa out.sq" << endl;
    exit(1);
  }
  ifstream inf(argv[1]);
  if (!inf) {
    cerr << "Could not open input file" << endl;
    exit(1);
  }
  ofstream outf(argv[2], ofstream::out);
  if (!outf) {
    cerr << "Could not write to output file" << endl;
    exit(1);
  }

  inf >> skipws;
  char a;

  while(!inf.eof()) {
    inf >> a;
    switch(a) {
    case 'A':
    case 'a':
    case 'C':
    case 'c':
    case 'T':
    case 't':
      outf << (char)toupper(a);
      break;
    case '>':
      // Skip this entire line
      inf.ignore(numeric_limits<streamsize>::max(), '\n');
      // Replace it with a newline
      outf << '\n';
      break;
    default: /* Any other character is replaced with 'G' */
      outf << 'G';
      break;
    }
  }
}
