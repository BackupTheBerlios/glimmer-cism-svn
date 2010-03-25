#include <cstdio>
#include <cstdlib>
#include <iostream>

using namespace std;

extern "C" {

  void helloworld_() {
    cout << "Hello World!" << endl;
    //return 0;
  }

} //extern "C"
