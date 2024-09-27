/* Convert fortran function to cuda
   What variables do we need?
   Textures for X, T/Q, and PDF   
*/

#include "defns.cuh"

void bin_search (double * data, double val, int * lo, int * hi) {
  /* perform binary search to find an index bin pair [lo,hi] 
     in data closest to value of val */
   while(*hi - *lo > 1) {
     int mid = (*hi + *lo) / 2;
     if(val >= data[mid]) {
       *lo = mid;
     } else {
       *hi = mid;
    }
  }
}

void bin_search_check_x(int lo, int upper, int size, int * ret) {
  /* ensure that nothing went wrong when binary searching on x */
  if (lo <= -1) {
    // X cannot be less than '0'
    printf("X:Severe error in Binary Search! Lower than min\n");
    exit(20);
  } else if(lo == 0)
    // actual value is in first bin, closest to min
    *ret = 0;
  else if(lo <= size-2)
    // clip if high val
    *ret = lo - 1;
  else if(lo == size-1)
    *ret = lo - 2;
  else {
    // x is not bigger than 1
    printf("X:Severe error in Binary Search! Higher than max\n");
    exit(21);
  }
}

void bin_search_check_q(int lower, int upper, int size, int * ret) {
  /* ensure that nothing went wrong when binary searching on Q */
  if (lower <= 0) {
    *ret = 0;
  } else if(lower <= size-2) {
    *ret = lower - 1;
  } else {
    *ret = size - 3;
  }
}

extern "C" void partonx12_wrapper_
(int * iSetch, int *iParton, int * nx, int * nt, int * npts,
 int * NfMx, int * MxVal,
 double * XX, double* QQ, double* qB,
 cudaTextureObject_t* xTex, cudaTextureObject_t* tTex, cudaTextureObject_t* pTex,
 double * ret) {

  // DECLARATIONS
  double X = *XX, Q = *QQ, qBase = *qB; 
  double svec1, svec2, svec3, svec4, s12, s13, s23, s24, s34, sy2, sy3,
    const1, const2, const3, const4, const5, const6, s1213, s2434, sdet, tmp,
    tvec1, tvec2, tvec3, tvec4, t12, t13, t23, t24, t34, ty2, ty3, tmp1, tmp2,
    tdet, *Fx, ss, tt;

  double Fij[4], Fvec[4];

  double * xv, * tv, * UPD, * xvpow, xpow = 0.3;

  int ip;
  // END DECLARATIONS
  
  // allocate size for x array, t array, and table array
  xv  = (double *) malloc(sizeof(double) * (*nx  ));
  tv  = (double *) malloc(sizeof(double) * (*nt  ));
  UPD = (double *) malloc(sizeof(double) * (*npts));

  // Place the values in the arrays from the textures
  texturevals_kernel_ (xTex, *nx  , xv); 
  texturevals_kernel_ (tTex, *nt  , tv);   
  texturevals_kernel_ (pTex, *npts, UPD);  

  // allocate xvpow array
  xvpow = (double*)malloc(sizeof(double) * (*nx));
  
  xvpow[0] = 0.0;
  // index at 0 
  for(int i = 0; i < *nx; i++) {
    xvpow[i] = pow(xv[i], xpow);
  }
    
  tt = log(log(Q/qBase));
  
  // binary search on x
  int lowerX = -1, upperX = *nx + 1, idxX;
  // search for most appropriate bin [lowerX, upperX] for x vals
  bin_search(xv, X, &lowerX, &upperX);
  // perform checks and find actual index value
  bin_search_check_x(lowerX, upperX, *nx, &idxX);

  ss = pow(X,xpow);
  
  if(lowerX >= 2 && lowerX <= *nx-2) {
    
    svec1 = xvpow[idxX];
    svec2 = xvpow[idxX+1];
    svec3 = xvpow[idxX+2];
    svec4 = xvpow[idxX+3];
    
    s12 = svec1 - svec2;
    s13 = svec1 - svec3;
    s23 = svec2 - svec3;
    s24 = svec2 - svec4;
    s34 = svec3 - svec4;
    
    sy2 = ss - svec2;
    sy3 = ss - svec3;
    
    const1 = s13/s23;
    const2 = s12/s23;
    const3 = s34/s23;
    const4 = s24/s23;
    s1213 = s12 + s13;
    s2434 = s24 + s34;
    sdet = s12*s34 - s1213*s2434;
    tmp = sy2*sy3/sdet;
    const5 = (s34*sy2-s2434*sy3)*tmp/s12;
    const6 = (s1213*sy2-s12*sy3)*tmp/s34;

  }
  
  // binary search on t (q)
  int lowerQ = -1, upperQ = *nt + 1, idxQ;
  bin_search(tv, tt, &lowerQ, &upperQ);
  bin_search_check_q(lowerQ, upperQ, *nt, &idxQ);

  if(lowerQ >= 1 && lowerQ <= *nt - 2) {
    tvec1 = tv[idxQ];
    tvec2 = tv[idxQ+1];
    tvec3 = tv[idxQ+2];
    tvec4 = tv[idxQ+3];
    
    t12 = tvec1 - tvec2;
    t13 = tvec1 - tvec3;
    t23 = tvec2 - tvec3;
    t24 = tvec2 - tvec4;
    t34 = tvec3 - tvec4;
    
    ty2 = tt - tvec2;
    ty3 = tt - tvec3;
    
    tmp1 = t12 + t13;
    tmp2 = t24 + t34;
    
    tdet = t12*t34 - tmp1*tmp2;
  }
  
  // done with first time setup

  if(*iParton > *MxVal) 
    ip = -*iParton;
  else
    ip = *iParton;

  int jtmp = ((ip + *NfMx)*(*nt+1)+(idxQ-1))*(*nx+1)+idxX+1;
  
  for(int i = 1; i <= nqvec; i++) {
    int tempIdx = jtmp + (i * (*nx+1));

    if(idxX == 0) {

      Fij[0] = 0;
      Fij[1] = UPD[tempIdx+1] * (pow(xv[1],2));
      Fij[2] = UPD[tempIdx+2] * (pow(xv[2],2));
      Fij[3] = UPD[tempIdx+3] * (pow(xv[3],2));

      polint_wrapper_(xvpow, Fij, &ss, Fx);

      if (X > 0.0)
	Fvec[i] = *Fx / pow(X,2.0);
      
    } else if(lowerX == *nx - 1) {

      polint_wrapper_(&xvpow[*nx-3], &UPD[tempIdx-1], &ss, Fx);
      Fvec[i] = *Fx;
      
    } else {
      
      double sf2, sf3, g1, g4;
      sf2 = UPD[tempIdx+0];
      sf3 = UPD[tempIdx+1];

      g1 =  sf2*const1 - sf3*const2;
      g4 = -sf2*const3 + sf3*const4;

      Fvec[i-1] = (const5 * (UPD[tempIdx-1]-g1)
	       + const6 * (UPD[tempIdx+2]-g4)
	       + sf2 * sy3 - sf3 * sy2) / s23;

    }
  }

  if(lowerQ <= 0) 
    polint_wrapper_(tv, Fvec, &tt, ret);
  else if (lowerQ >= *nt-1)
    polint_wrapper_(&tv[*nt-3], Fvec, &tt, ret);
  else {

    double tf2, tf3, g1,g4,h00;

    tf2 = Fvec[1];
    tf3 = Fvec[2];

    g1 = ( tf2*t13 - tf3*t12) / t23;
    g4 = (-tf2*t34 + tf3*t24) / t23;

    h00 = ((t34*ty2-tmp2*ty3)*(Fvec[0]-g1)/t12 + (tmp1*ty2-t12*ty3)*(Fvec[3]-g4)/t34);

    *ret = (h00*ty2*ty3/tdet + tf2*ty3 - tf3*ty2) / t23;
    
  }
  
  free(xv); free(tv); free(UPD);
  return;
  
}


