/* This program takes in the output of ZTOPtchanPDF.x in the delta file
   and computes the PDF error. */

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cmath>

int main(int argc, char* argv[])
{
  using namespace std;

  int sets = 59;
  double tolerance = 10.;

  vector<double> nums;
  
  for(int i=0; i < sets; ++i) {
    double vals[10];
    for (int j=0; j<10; ++j) {
      cin >> vals[j];
    }
    nums.push_back(vals[8]);
  }

  /*
  for(int i=0; i<=58; ++i) {
    cout << nums[i] << " ";
  }
  cout << endl;
  */

  double center = nums[0];
  vector<double> errplus, errminus;
  
  for(int i=1; i <= (sets-1)/2; ++i) {
    // Look at pairs for asymmetric errors
    double errA = nums[2*i-1] - center;
    double errB = nums[2*i] - center;
    if (errA >= 0. && errB >= 0.) {
      errplus.push_back(max(errA, errB));
      errminus.push_back(0.);
    }
    else if (errA <= 0. && errB <= 0.) {
      errplus.push_back(0.);
      errminus.push_back(min(errA, errB));
    }
    else {
      errplus.push_back(max(errA, errB));
      errminus.push_back(min(errA, errB));
    }
  }

  // Get tolerance error
  double possum = 0.;
  double negsum = 0.;
  for (int i=0; i < (sets-1)/2; ++i) {
    possum += errplus[i]*errplus[i];
    negsum += errminus[i]*errminus[i];
  }
  //  double poserr = sqrt(tolerance * possum/(sets-1.)*2.);
  //  double negerr = sqrt(tolerance * negsum/(sets-1.)*2.);
  double poserr = sqrt(possum);
  double negerr = sqrt(negsum);

  cout << "Result:" << endl;
  cout << center << " + " << poserr << " - " << negerr << endl;
  
  return 0;
}
