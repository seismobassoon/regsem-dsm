/* b-splines on arbitrary spacing 
   CHM 12/97
*/
#include <math.h>
#include<stdio.h>

#define U (unsigned)
#define TOL 1.5

extern float *farray1();
extern void free_farray1();
extern void stop(); 
void fill_hh_A3d();
float spl_A3d();

/*wrapper for call from fortran 90 code*/
#ifdef _IBM_
void fspl(ord,nknots,knot,xi,rho)
#else
void fspl_(ord,nknots,knot,xi,rho)
#endif  
int *ord,*nknots;
double *knot,*xi,*rho;
{
  float *knot_s,xi_s,val;
  int ord_s,i;

  ord_s=*ord-1;/*change from fortran to c array convention*/
  xi_s=(float)(*xi);
  knot_s=farray1(0,*nknots-1);
  for(i=0;i<*nknots;i++)
    knot_s[i]=(float)knot[i];
  val=spl_A3d(ord_s,*nknots,knot_s,xi_s);
  *rho=(double)val;
  free_farray1(knot_s,0);
  return;
}

float spl_A3d(ord,nknots,knot,xi)
int ord,nknots;
float *knot,xi;
/* ord: number of rho(x)
   nknots : # of knkots = Nx+1 (Nx=index of highest spline)
   xi: point of interest
   splines defined as :
   f_i(x) = a_i(x-x_i)^3 + b_i(x-x_i)^2 + c_i(x-x_i) + d_i
*/
{
int ii,Nx;
float *hh,rho_x;
float coefa,coefb,coefc,coefd;
  
  Nx=nknots-1;
  /* Compute vector hh of spacings */
  hh=farray1(0,Nx-1);
  fill_hh_A3d(hh,knot,Nx);
  
  /* Consistency checks */
  if((xi-(float)TOL)>knot[Nx]){
    printf("xi=%g / knot[%d]=%g",xi,Nx,knot[Nx]);
    stop("spl: xi>knot[Nx]");
  }
  else if((xi+(float)TOL)<knot[0]){
    printf("xi=%g / knot[0]=%g",xi,knot[0]);
    stop("spl: xi<knot[0]");
  }
  else if(ord>Nx)
    stop("spl: order > Nx");
  
  if(ord==0){	/* LHS */
  float denom;
    denom=3.*hh[ord]*hh[ord]+3.*hh[ord]*hh[ord+1]+hh[ord+1]*hh[ord+1];
    if(xi>=knot[ord]&&xi<=knot[ord+1]){			/* x0<=x<=x1 */
      coefa=4./(hh[ord]*(hh[ord]+hh[ord+1])*denom);
      coefb=0.0;
      coefc=-12/denom;
      coefd=4*(2*hh[ord]+hh[ord+1])/denom;
      
      rho_x=coefa*(xi-knot[ord])*(xi-knot[ord])*(xi-knot[ord]);
      rho_x+=coefb*(xi-knot[ord])*(xi-knot[ord]);
      rho_x+=coefc*(xi-knot[ord]);
      rho_x+=coefd;
    }
    else if(xi>knot[ord+1]&xi<=knot[ord+2]){		/* x1<=x<=x2 */
      coefa=-4./(hh[ord+1]*(hh[ord]+hh[ord+1])*denom);
      coefb=12/((hh[ord]+hh[ord+1])*denom);
      coefc=-12.*hh[ord+1]/((hh[ord]+hh[ord+1])*denom);
      coefd=4.*hh[ord+1]*hh[ord+1]/((hh[ord]+hh[ord+1])*denom);
      
      rho_x=coefa*(xi-knot[ord+1])*(xi-knot[ord+1])*(xi-knot[ord+1]);
      rho_x+=coefb*(xi-knot[ord+1])*(xi-knot[ord+1]);
      rho_x+=coefc*(xi-knot[ord+1]);
      rho_x+=coefd;
    }
    else						/* x>x2 */
      rho_x=0.0; 
  }
  
  else if(ord==1){	/* LHS+1 */
  float denom,denomsum,dd;
    denom=(3.*hh[ord-1]*hh[ord-1]+4.*hh[ord-1]*hh[ord]+hh[ord]*hh[ord]+
           2.*hh[ord-1]*hh[ord+1]+hh[ord]*hh[ord+1]);
    denomsum=hh[ord-1]+hh[ord]+hh[ord+1];
    dd=denomsum*denom;
    if(xi>=knot[ord-1]&&xi<=knot[ord]){			/* x0<=x<=x1 */
      coefa=-4.*(3.*hh[ord-1]+2.*hh[ord]+hh[ord+1])/
      	    (hh[ord-1]*(hh[ord-1]+hh[ord])*dd);
      coefb=0.;
      coefc=12./denom;
      coefd=0.;
      
      rho_x =coefa*(xi-knot[ord-1])*(xi-knot[ord-1])*(xi-knot[ord-1]);
      rho_x+=coefb*(xi-knot[ord-1])*(xi-knot[ord-1]);
      rho_x+=coefc*(xi-knot[ord-1]);
      rho_x+=coefd;
    }
    
    else if(xi>=knot[ord]&&xi<=knot[ord+1]){			/* x1<=x<=x2 */
      coefa=4.*(2.*hh[ord-1]*hh[ord-1]+6.*hh[ord-1]*hh[ord]+3.*hh[ord]*hh[ord]+3.*hh[ord-1]*hh[ord+1]+
            3.*hh[ord]*hh[ord+1]+hh[ord+1]*hh[ord+1])/
            (hh[ord]*(hh[ord-1]+hh[ord])*(hh[ord]+hh[ord+1])*dd);
      coefb=-12.*(3.*hh[ord-1]+2.*hh[ord]+hh[ord+1])/
      	    ((hh[ord-1]+hh[ord])*dd);
      coefc=12.*(-2.*hh[ord-1]*hh[ord-1]+hh[ord]*hh[ord]+hh[ord]*hh[ord+1])/
      	    ((hh[ord-1]+hh[ord])*dd);
      coefd=4.*hh[ord-1]*(4.*hh[ord-1]*hh[ord]+3.*hh[ord]*hh[ord]+2.*hh[ord-1]*hh[ord+1]+3.*hh[ord]*hh[ord+1])/
      	    ((hh[ord-1]+hh[ord])*dd);
      	    
      rho_x=coefa*(xi-knot[ord])*(xi-knot[ord])*(xi-knot[ord]);
      rho_x+=coefb*(xi-knot[ord])*(xi-knot[ord]);
      rho_x+=coefc*(xi-knot[ord]);
      rho_x+=coefd;
    }
    
    else if(xi>=knot[ord+1]&&xi<=knot[ord+2]){			/* x2<=x<=x3 */
      dd*=(hh[ord]+hh[ord+1]);
      coefa=-4.*(2.*hh[ord-1]+hh[ord])/(hh[ord+1]*dd);
      coefb=12.*(2.*hh[ord-1]+hh[ord])/dd;
      coefc=-12.*(2.*hh[ord-1]+hh[ord])*hh[ord+1]/dd;
      coefd=4.*(2.*hh[ord-1]+hh[ord])*hh[ord+1]*hh[ord+1]/dd;
      
      rho_x=coefa*(xi-knot[ord+1])*(xi-knot[ord+1])*(xi-knot[ord+1]);
      rho_x+=coefb*(xi-knot[ord+1])*(xi-knot[ord+1]);
      rho_x+=coefc*(xi-knot[ord+1]);
      rho_x+=coefd;
    }
    else						/* x>x3 */
      rho_x=0.0; 
  }
  
  else if(ord==Nx-1){		/* RHS-1 */
  float denom,denomsum,dd;
    denom=hh[ord-2]*hh[ord-1]+hh[ord-1]*hh[ord-1]+2.*hh[ord-2]*hh[ord]+4.*hh[ord-1]*hh[ord]+3.*hh[ord]*hh[ord];
    denomsum=hh[ord-2]+hh[ord-1]+hh[ord];
    dd=denomsum*denom;
    if(xi>=knot[ord-2]&&xi<=knot[ord-1]){	/* x0<=x<=x1 */
      coefa=4.*(hh[ord-1]+2.*hh[ord])/(hh[ord-2]*(hh[ord-2]+hh[ord-1])*dd);
      coefb=coefc=coefd=0.0;
      
      rho_x =coefa*(xi-knot[ord-2])*(xi-knot[ord-2])*(xi-knot[ord-2]);
      rho_x+=coefb*(xi-knot[ord-2])*(xi-knot[ord-2]);
      rho_x+=coefc*(xi-knot[ord-2]);
      rho_x+=coefd;
    }
    
    else if(xi>=knot[ord-1]&&xi<=knot[ord]){	/* x1<=x<=x2 */
      coefa=-4.*(hh[ord-2]*hh[ord-2]+3.*hh[ord-2]*hh[ord-1]+3.*hh[ord-1]*hh[ord-1]+3.*hh[ord-2]*hh[ord]+
             6.*hh[ord-1]*hh[ord]+2.*hh[ord]*hh[ord])/
             (hh[ord-1]*(hh[ord-2]+hh[ord-1])*(hh[ord-1]+hh[ord])*dd);
      coefb=12.*(hh[ord-1]+2.*hh[ord])/((hh[ord-2]+hh[ord-1])*dd);
      coefc=12.*hh[ord-2]*(hh[ord-1]+2.*hh[ord])/((hh[ord-2]+hh[ord-1])*dd);
      coefd=4.*hh[ord-2]*hh[ord-2]*(hh[ord-1]+2.*hh[ord])/((hh[ord-2]+hh[ord-1])*dd);
      
      rho_x =coefa*(xi-knot[ord-1])*(xi-knot[ord-1])*(xi-knot[ord-1]);
      rho_x+=coefb*(xi-knot[ord-1])*(xi-knot[ord-1]);
      rho_x+=coefc*(xi-knot[ord-1]);
      rho_x+=coefd;
    }
    
    else if(xi>=knot[ord]&&xi<=knot[ord+1]){	/* x2<=x<=x3 */
      dd*=(hh[ord-1]+hh[ord]);
      coefa=4.*(hh[ord-2]+2.*hh[ord-1]+3.*hh[ord])/(hh[ord]*dd);
      coefb=-12.*(hh[ord-2]+2.*hh[ord-1]+3.*hh[ord])/dd;
      coefc=12.*(-hh[ord-2]*hh[ord-1]-hh[ord-1]*hh[ord-1]+2.*hh[ord]*hh[ord])/dd;
      coefd=4.*hh[ord]*(3.*hh[ord-2]*hh[ord-1]+3.*hh[ord-1]*hh[ord-1]+2.*hh[ord-2]*hh[ord]+4.*hh[ord-1]*hh[ord])/dd;
      
      rho_x=coefa*(xi-knot[ord])*(xi-knot[ord])*(xi-knot[ord]);
      rho_x+=coefb*(xi-knot[ord])*(xi-knot[ord]);
      rho_x+=coefc*(xi-knot[ord]);
      rho_x+=coefd;
    }
    else						/* x>x4 */
      rho_x=0.0;
  }
  
  else if(ord==Nx){		/* RHS */
  float denom;
    denom=(hh[ord-2]+hh[ord-1])*(hh[ord-2]*hh[ord-2]+3.*hh[ord-2]*hh[ord-1]+3.*hh[ord-1]*hh[ord-1]);
    if(xi>=knot[ord-2]&&xi<=knot[ord-1]){	/* x0<=x<=x1 */
      coefa=4./(hh[ord-2]*denom);
      coefb=coefc=coefd=0.0;
      
      rho_x =coefa*(xi-knot[ord-2])*(xi-knot[ord-2])*(xi-knot[ord-2]);
      rho_x+=coefb*(xi-knot[ord-2])*(xi-knot[ord-2]);
      rho_x+=coefc*(xi-knot[ord-2]);
      rho_x+=coefd;
    }
    
    else if(xi>=knot[ord-1]&&xi<=knot[ord]){	/* x1<=x<=x2 */
      coefa=-4./(hh[ord-1]*denom);
      coefb=12/denom;
      coefc=12*hh[ord-2]/denom;
      coefd=4.*hh[ord-2]*hh[ord-2]/denom;
      
      rho_x =coefa*(xi-knot[ord-1])*(xi-knot[ord-1])*(xi-knot[ord-1]);
      rho_x+=coefb*(xi-knot[ord-1])*(xi-knot[ord-1]);
      rho_x+=coefc*(xi-knot[ord-1]);
      rho_x+=coefd;
    }
    
    else						/* x>x2 */
      rho_x=0.0;
  }
  
  else{			/* Away from borders */
  float denom1,denom2,denom;
    denom1=hh[ord-2]+hh[ord-1]+hh[ord]+hh[ord+1];
    if(xi>=knot[ord-2]&&xi<=knot[ord-1]){	/* x0<=x<=x1 */
      coefa=4./(hh[ord-2]*(hh[ord-2]+hh[ord-1])*(hh[ord-2]+hh[ord-1]+hh[ord])*denom1);
      coefb=coefc=coefd=0.;
      
      rho_x =coefa*(xi-knot[ord-2])*(xi-knot[ord-2])*(xi-knot[ord-2]);
      rho_x+=coefb*(xi-knot[ord-2])*(xi-knot[ord-2]);
      rho_x+=coefc*(xi-knot[ord-2]);
      rho_x+=coefd;
    }
    else if(xi>=knot[ord-1]&&xi<=knot[ord]){	/* x1<=x<=x2 */
      denom2=(hh[ord-2]+hh[ord-1])*(hh[ord-2]+hh[ord-1]+hh[ord]);
      denom=denom1*denom2;
      
      coefa=-4.*(hh[ord-2]*hh[ord-2]+3.*hh[ord-2]*hh[ord-1]+3.*hh[ord-1]*hh[ord-1]+2.*hh[ord-2]*hh[ord]+
      4.*hh[ord-1]*hh[ord]+hh[ord]*hh[ord]+hh[ord-2]*hh[ord+1]+2.*hh[ord-1]*hh[ord+1]+hh[ord]*hh[ord+1])/
      (hh[ord-1]*(hh[ord-1]+hh[ord])*(hh[ord-1]+hh[ord]+hh[ord+1])*denom);
      coefb=12./denom;
      coefc=12.*hh[ord-2]/denom;
      coefd=4.*hh[ord-2]*hh[ord-2]/denom;
      
      rho_x =coefa*(xi-knot[ord-1])*(xi-knot[ord-1])*(xi-knot[ord-1]);
      rho_x+=coefb*(xi-knot[ord-1])*(xi-knot[ord-1]);
      rho_x+=coefc*(xi-knot[ord-1]);
      rho_x+=coefd;
    }
    
    else if(xi>=knot[ord]&&xi<=knot[ord+1]){	/* x2<=x<=x3 */
      denom2=(hh[ord-1]+hh[ord])*(hh[ord-2]+hh[ord-1]+hh[ord])*(hh[ord-1]+hh[ord]+hh[ord+1]);
      denom=denom1*denom2;
      
      coefa=4.*(hh[ord-2]*hh[ord-1]+hh[ord-1]*hh[ord-1]+2.*hh[ord-2]*hh[ord]+4.*hh[ord-1]*hh[ord]+3.*hh[ord]*hh[ord]+
      	    hh[ord-2]*hh[ord+1]+2.*hh[ord-1]*hh[ord+1]+3.*hh[ord]*hh[ord+1]+hh[ord+1]*hh[ord+1])/
      	    (hh[ord]*(hh[ord]+hh[ord+1])*denom);
      coefb=-12.*(hh[ord-2]+2.*hh[ord-1]+2.*hh[ord]+hh[ord+1])/denom;
      coefc=12.*(-hh[ord-2]*hh[ord-1]-hh[ord-1]*hh[ord-1]+hh[ord]*hh[ord]+hh[ord]*hh[ord+1])/denom;
      coefd=4.*(2.*hh[ord-2]*hh[ord-1]*hh[ord]+2.*hh[ord-1]*hh[ord-1]*hh[ord]+hh[ord-2]*hh[ord]*hh[ord]+
      	    2.*hh[ord-1]*hh[ord]*hh[ord]+hh[ord-2]*hh[ord-1]*hh[ord+1]+hh[ord-1]*hh[ord-1]*hh[ord+1]+
      	    hh[ord-2]*hh[ord]*hh[ord+1]+2.*hh[ord-1]*hh[ord]*hh[ord+1])/denom;
     
      rho_x=coefa*(xi-knot[ord])*(xi-knot[ord])*(xi-knot[ord]);
      rho_x+=coefb*(xi-knot[ord])*(xi-knot[ord]);
      rho_x+=coefc*(xi-knot[ord]);
      rho_x+=coefd;
    }
             
    else if(xi>=knot[ord+1]&&xi<=knot[ord+2]){	/* x3<=x<=x4 */
      denom2=(hh[ord]+hh[ord+1])*(hh[ord-1]+hh[ord]+hh[ord+1]);
      denom=denom1*denom2;
      
      coefa=-4./(hh[ord+1]*denom);
      coefb=12/denom;
      coefc=-12*hh[ord+1]/denom;
      coefd=4.*hh[ord+1]*hh[ord+1]/denom;
      
      rho_x=coefa*(xi-knot[ord+1])*(xi-knot[ord+1])*(xi-knot[ord+1]);
      rho_x+=coefb*(xi-knot[ord+1])*(xi-knot[ord+1]);
      rho_x+=coefc*(xi-knot[ord+1]);
      rho_x+=coefd;
    }
    
    else						/* x>x4 */
      rho_x=0.0;
  }
  free_farray1(hh,0);
  return(rho_x);
}

void fill_hh_A3d(hh,knot,Nx)
float *hh,*knot;
int Nx;
{
int ii;

  for(ii=0;ii<Nx;ii++)
    hh[ii]=knot[ii+1]-knot[ii];
}

float *farray1(n11,n12)
int n11,n12;
{
float *m;
        m=(float *)malloc( U (n12-n11+1)*sizeof(float) );
        if(!m) {
          stop("allocation error in farray1");
        }
        return(m-n11);
}

void free_farray1(a,n11)
float *a;
int n11;
{
        free(&a[n11]);
}

void stop(char *message){
  printf("\n\a <error> %s\n",message);
  exit(-1);
} 
