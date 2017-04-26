#include <math.h>
#include <assert.h>
#include <string.h>
#include "mex.h"
#include "MxArray.hpp"
#include <vector>

#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

/*
 *	Use opencv function to resample image quickly
 */

// matlab entry point
// dst = resize(src, scale)
// image should be color with double values
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { 
  if (nrhs < 2){
  	mexErrMsgTxt("Wrong number of inputs");
  }
  if (nlhs != 1){
  	mexErrMsgTxt("Wrong number of outputs");
  }

  vector<MxArray> rhs(prhs,prhs+nrhs);

  //convert input data to opencv matrix
  Mat img = rhs[0].toMat();
  Mat imgr;
  Size s = rhs[1].toSize();
  Size newSize = Size(s.height,s.width);
  Size oldSize = img.size();
  //interpolation method
  int interpolation = INTER_LINEAR;

  //if interpolation method provided set it
  if(nrhs == 3){
    string interp = rhs[2].toString();
    if(interp.compare("antialias") == 0){
      interpolation = INTER_AREA;
    }else if(interp.compare("linear") == 0){
      interpolation = INTER_LINEAR;
    }else if(interp.compare("auto") == 0){ //if we are zooming, use linear else use area interpolation
      //old array has width and height swapped, newArray does not
      if(newSize.width > oldSize.height){
        interpolation = INTER_LINEAR;
      }else{
        interpolation = INTER_AREA;
      }
    }else{
      mexErrMsgTxt("Invalid interpolation provided, valid is linear (default), antialias, auto");
    }
  }

  //use opencv resize function
  resize(img,imgr,newSize,0,0,interpolation);
  //convert back to matlab representation
  plhs[0] = MxArray(imgr);
}



