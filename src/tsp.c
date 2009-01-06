#include <R.h>
#include <Rinternals.h>
#include <Rmath.h> 
#include <R_ext/Utils.h> 
#define MAXTSP 2500


SEXP ctspair(SEXP Rdat, SEXP Rgrp)
{
  int i,j,k;
  double score; 
  int nProtected = 0; 

  int probsum[2]; probsum[0] = 0; probsum[1] = 0;
  
  double *grp, *dat;
  grp = REAL(Rgrp);
  dat = REAL(Rdat); 
  
  int n, m;
  n = LENGTH(Rgrp);
  m = LENGTH(Rdat)/n;
  

  double n0, n1;
  n0 = 0; 
  n1 = 0;
  for(i = 0; i < n; i++){
    if(grp[i] == 0){n0 = n0 + 1;}
    if(grp[i] == 1){n1 = n1 + 1;}
  }
  
  SEXP maxval, ind, ntsp, list; 
  int *xind, *xntsp;
  double *xmaxval;

  int mcount; mcount = 1;  
  
  PROTECT(maxval = allocVector(REALSXP,1));
  nProtected++;
  PROTECT(ntsp = allocVector(INTSXP,1));
  nProtected++;
  PROTECT(ind = allocVector(INTSXP,m*m));
  nProtected++;
   

  xmaxval = REAL(maxval);
  xntsp = INTEGER(ntsp);
  xind = INTEGER(ind);

  xmaxval[0] = 0; 
  xntsp[0] = 0; 

  for(i = 0; i < m; i++){
    for(j = 0; j < i; j++){
      for(k = 0; k < n; k++){
	if(grp[k] == 1){
	  probsum[1] += (dat[i + k * m] < dat[j + k * m]);
	}
	if(grp[k] == 0){
	  probsum[0] += (dat[i + k * m] < dat[j + k * m]);
	}
      }
      score = fabs(probsum[1]/n1 - probsum[0]/n0);
      if(score == xmaxval[0]){
	xntsp[0] = xntsp[0] + 1;
	if(xntsp[0] >= (mcount*2500 - 1)){
 	  mcount++;
 	  ind = lengthgets(ind, (mcount*2*MAXTSP));
	  UNPROTECT(1);
 	  PROTECT(ind);
	  xind = INTEGER(ind);
	  printf("Mcount is %d",mcount);
 	}
	xind[(2*(xntsp[0]-1))] = i + 1;
	xind[(2*(xntsp[0]-1) + 1)] = j + 1;
	
      }
      if(score > xmaxval[0]){
	UNPROTECT(1);
	PROTECT(ind = allocVector(INTSXP,(2*MAXTSP)));
	xind = INTEGER(ind);
	xntsp[0] = 1;
	xind[0] = i + 1;
	xind[1] = j + 1;
	xmaxval[0] = score;
	mcount = 1; 
      }
      probsum[0] = 0;
      probsum[1] = 0;
      void R_CheckUserInterrupt(void); 
    }
  }

  ind = lengthgets(ind, (2*xntsp[0]));
  UNPROTECT(1);
  PROTECT(ind);
  
  PROTECT(list=allocVector(VECSXP,3));
  ++nProtected;
  SET_VECTOR_ELT(list, 0, ind);
  SET_VECTOR_ELT(list, 1, maxval);
  SET_VECTOR_ELT(list, 2, ntsp);
   
  UNPROTECT(nProtected);
  return(list);  

}


