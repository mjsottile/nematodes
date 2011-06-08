#include "mex.h"

void mexFunction (int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[]) 
{
  double *xs, *ys, *zs, *segs;
  mwSize numPts, numSegs, nzmax;
  double percentSparse = 0.2;
  int k;
  int curcol;

  if (nrhs != 4) {
    mexErrMsgTxt("Four inputs required.  xs, ys, zs, tets.");
  }
  if (nlhs != 1) {
    mexErrMsgTxt("one output required.");
  }

  xs = mxGetPr(prhs[0]);
  ys = mxGetPr(prhs[1]);
  zs = mxGetPr(prhs[2]);
  segs = mxGetPr(prhs[3]);

  numPts = mxGetM(prhs[0]);
  numSegs = mxGetM(prhs[3]);

  nzmax = (mwSize)ceil((double)numPts * (double)numPts * percentSparse);

  plhs[0] = mxCreateSparse(numPts, numPts, nzmax, mxREAL);

  /* get components of sparse array */
  sr = mxGetPr(plhs[0]);
  irs = mxGetIr(plhs[0]);
  jcs = mxGetJc(plhs[0]);

  k = 0;
  curcol = 1;
  end1_idx = 0;
  end2_idx = numSegs;
  
  for (int i = 0; i < numPts; i++) {
    while (segs[end1_idx] == curcol && end1_idx < numSegs) {
      int p1 = segs[end1_idx];
      int p2 = segs[end2_idx];

      p1 = p1-1; p2 = p2-1; /* adjust to 0 based indices */
      if (abs(zs[p1]-zs[p2]) > 0.0000001) {
	if (zs[p1] > zs[p2]) {
	  
	} else {

	}
      }

      end1_idx++;				
      end2_idx++;
    }

    if (k >= nzmax) {
      mwSize oldnzmax = nzmax;
      percent_sparse += 0.1;
      nzmax= (mwSize)ceil((double)numPts *(double)numPts * percentSparse);
      if (oldnzmax == nzmax) nzmax++;
      mxSetNzmax(plhs[0],nzmax);
      mxSetPr(plhs[0], mxRealloc(sr, nzmax*sizeof(double)));
      mxSetIr(plhs[0], mxRealloc(irs, nzmax*sizeof(mwIndex)));
      sr = mxGetPr(plhs[0]);
      irs = mxGetIr(plhs[0]);
    }
    
  }
}
