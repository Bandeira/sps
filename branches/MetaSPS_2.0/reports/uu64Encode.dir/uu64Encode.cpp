#include <iostream>
#include <fstream>

#include "base64.h"


using namespace std;


int main(int argc, char **argv)
{
  // test for empty parameters
  if(argc != 2) {
    cout << "uuEncode64 - a uu64 file encoding tool." << endl;
    cout << "Usage: uuEncode64 <input file>" << endl;
    return 0;
  }

  // get the filename
  string fn = argv[1];

  // load the file (image), and store it in internal image
  int length;
  char * buffer;

  ifstream is;
  is.open (fn.c_str(), ios::binary );

  // check for file open
  if(!is.is_open()) {
    cout << "Unable to open file " << fn << endl;
    return false;
  }

  // get length of file:
  is.seekg (0, ios::end);
  length = is.tellg();
  is.seekg (0, ios::beg);

  // allocate memory:
  buffer = new char [length+1];

  // read data as a block:
  is.read(buffer, length);
  is.close();
  // terminator
  buffer[length] = 0;

  // if uuencoded, encode
  string m_image = base64_encode(reinterpret_cast<const unsigned char*>(buffer), length);

  // dump encoded file to stdout
  cout <<  m_image;

  // return ok
  return 0;
}
