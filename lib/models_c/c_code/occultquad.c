/*
C  This routine computes the lightcurve for occultation
C  of a quadratically limb-darkened source without microlensing.
C  Please cite Mandel & Agol (2002) if you make use of this routine
C  in your research.  Please report errors or bugs to agol@tapir.caltech.edu

	This code is a translation of the Mandel/Agol Fortran code into C.
	LK 9/13/12
*/
#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>
#include <omp.h>

#define IND(a,i) *((double *)(a->data+i*a->strides[0]))

double rc(double x, double y);
double rj(double x, double y, double z, double p);
double ellec(double k);
double ellk(double k);
double rf(double x, double y, double z);
double ellpic_bulirsch(double n, double k);
double max(double a, double b);
double min(double a, double b);
static PyObject *occultquad(PyObject *self, PyObject *args);

static PyObject *occultquad(PyObject *self, PyObject *args)
{
/*	Input: *************************************
 
	 z0   impact parameter in units of rs
	 u1   linear    limb-darkening coefficient (gamma_1 in paper)
	 u2   quadratic limb-darkening coefficient (gamma_2 in paper)
	 p    occulting star size in units of rs
	 nz   number of points in the light curve (# elements in z0)
	
	 Output: ***********************************
	
	 default:  muo1, fraction of flux at each z0 for a limb-darkened source
	 optional:  mu0, fraction of flux at each z0 for a uniform source
	
	 Limb darkening has the form:
	 I(r)=[1-u1*(1-sqrt(1-(r/rs)^2))-u2*(1-sqrt(1-(r/rs)^2))^2]/(1-u1/3-u2/6)/pi
*/
	int i, nz;
	double u1, u2, p, *mu, *lambdad, *etad, \
		*lambdae, pi, x1, x2, x3, z, omega, kap0 = 0.0, kap1 = 0.0, \
		q, Kk, Ek, Pk, n;
	PyArrayObject *z0, *muo1, *mu0;
	npy_intp dims[1];

  	if(!PyArg_ParseTuple(args,"Odddi", &z0, &u1, &u2, &p, &nz))
   	 {
      		return NULL;
   	 }

	lambdad = (double *)malloc(nz*sizeof(double));
	lambdae = (double *)malloc(nz*sizeof(double));
	etad = (double *)malloc(nz*sizeof(double));
	mu = (double *)malloc(nz*sizeof(double));

	if(fabs(p - 0.5) < 1.0e-3) p = 0.5;
	dims[0] = z0->dimensions[0];
	mu0 = (PyArrayObject *) PyArray_SimpleNew(1, dims, PyArray_DOUBLE);
	muo1 = (PyArrayObject *) PyArray_SimpleNew(1, dims, PyArray_DOUBLE);

	omega=1.0-u1/3.0-u2/6.0;
	pi=acos(-1.0);
	
	#pragma omp parallel for private(z,x1,x2,x3,n,q,Kk,Ek,Pk,kap0,kap1)
	for(i = 0; i < nz; i++)
	{	
		z  = IND(z0, i);
		x1 = pow((p - z), 2.0);
		x2 = pow((p + z), 2.0);
		x3 = p*p - z*z;

		//source is unocculted:
		if(z > 1.0 + p)					
		{
	//		printf("zone 1\n");
			lambdad[i] = 0.0;
			etad[i] = 0.0;
			lambdae[i] = 0.0;
			IND(muo1, i) = 1.0-((1.0-u1-2.0*u2)*lambdae[i]+(u1+2.0*u2)*lambdad[i]+u2*etad[i])/omega;
			IND(mu0, i) = 1.0-lambdae[i];
			continue;
		}
		//source is completely occulted:
		if(p > 1.0 && z < p - 1.0)			
		{
	//		printf("zone 2\n");
			lambdad[i] = 0.0;
			etad[i] = 0.5;		//error in fortran code corrected here, following Eastman's python code
			lambdae[i] = 1.0;
			IND(muo1, i) = 1.0-((1.0-u1-2.0*u2)*lambdae[i]+(u1+2.0*u2)*(lambdad[i] + 2.0/3.0)+u2*etad[i])/omega;
			IND(mu0, i) = 1.0-lambdae[i];
			continue;
		}
		//source is partly occulted and occulting object crosses the limb:
		if(z > fabs(1.0 - p) && z < 1.0 + p)	
		{				
			//printf("zone 3\n");
			kap1 = acos(min((1.0-p*p+z*z)/2.0/z,1.0));
			kap0 = acos(min((p*p+z*z-1.0)/2.0/p/z,1.0));
			lambdae[i] = p*p*kap0+kap1;
			lambdae[i] = (lambdae[i] - 0.50*sqrt(max(4.0*z*z-pow((1.0+z*z-p*p), 2.0),0.0)))/pi;
		}
		//occulting object transits the source but doesn't completely cover it:
		if(z < 1.0 - p)				
		{					
			//printf("zone 4\n");
			lambdae[i] = p*p;
		}
		//edge of the occulting star lies at the origin
		if(fabs(z-p) < 1.0e-4*(z+p))			
		{
			//printf("zone 5\n");
			if(z > 0.5)
			{
				//printf("zone 5.1\n");
				//lam = 0.5*pi;
				q = 0.5/p;
				Kk = ellk(q);
				Ek = ellec(q);
				lambdad[i] = 1.0/3.0+16.0*p/9.0/pi*(2.0*p*p-1.0)*Ek- \
				 	 	(32.0*pow(p, 4.0)-20.0*p*p+3.0)/9.0/pi/p*Kk;
				etad[i] = 1.0/2.0/pi*(kap1+p*p*(p*p+2.0*z*z)*kap0- \
				              	(1.0+5.0*p*p+z*z)/4.0*sqrt((1.0-x1)*(x2-1.0)));
				if(p == 0.5)
				{
					lambdad[i] = 1.0/3.0-4.0/pi/9.0;
					etad[i] = 3.0/32.0;
				}
				IND(muo1, i) = 1.0-((1.0-u1-2.0*u2)*lambdae[i]+(u1+2.0*u2)*lambdad[i]+u2*etad[i])/omega;
				IND(mu0, i) = 1.0-lambdae[i];
				continue;
			}
			else
			{
				//printf("zone 5.2\n");
				//lam = 0.50*pi;
				q = 2.0*p;
				Kk = ellk(q);
				Ek = ellec(q);
				lambdad[i] = 1.0/3.0+2.0/9.0/pi*(4.0*(2.0*p*p-1.0)*Ek + (1.0-4.0*p*p)*Kk);
				etad[i] = p*p/2.0*(p*p+2.0*z*z);
				IND(muo1, i) = 1.0-((1.0-u1-2.0*u2)*lambdae[i]+(u1+2.0*u2)*lambdad[i]+u2*etad[i])/omega;
				IND(mu0, i) = 1.0-lambdae[i];
				continue;
			}
		}
		//occulting star partly occults the source and crosses the limb:
		if((z > 0.5 + fabs(p -0.5) && z < 1.0 + p) || (p > 0.5 && z > fabs(1.0-p)*1.0001 \
			&& z < p))
		{
			//printf("zone 6\n");
			//lam = 0.50*pi;
			q = sqrt((1.0-pow((p-z), 2.0))/4.0/z/p);
			Kk = ellk(q);
			Ek = ellec(q);
			n = 1.0/x1-1.0;
			//Pk = Kk-n/3.0*rj(0.0,1.0-q*q,1.0,1.0+n);
			Pk = ellpic_bulirsch(n,q);
			lambdad[i] = 1.0/9.0/pi/sqrt(p*z)*(((1.0-x2)*(2.0*x2+ \
			        x1-3.0)-3.0*x3*(x2-2.0))*Kk+4.0*p*z*(z*z+ \
			        7.0*p*p-4.0)*Ek-3.0*x3/x1*Pk);
			if(z < p) lambdad[i] += 2.0/3.0;
			etad[i] = 1.0/2.0/pi*(kap1+p*p*(p*p+2.0*z*z)*kap0- \
				(1.0+5.0*p*p+z*z)/4.0*sqrt((1.0-x1)*(x2-1.0)));
			IND(muo1, i) = 1.0-((1.0-u1-2.0*u2)*lambdae[i]+(u1+2.0*u2)*lambdad[i]+u2*etad[i])/omega;
			IND(mu0, i) = 1.0-lambdae[i];
			continue;
		}
		//occulting star transits the source:
		if(p < 1.0  && z < (1.0 - p)*1.0001)
		{
			//printf("zone 7\n");
			//lam = 0.50*pi;
			q = sqrt((x2-x1)/(1.0-x1));
			Kk = ellk(q);
			Ek = ellec(q);
			n = x2/x1-1.0;
			//Pk = Kk-n/3.0*rj(0.0,1.0-q*q,1.0,1.0+n);
			Pk = ellpic_bulirsch(n,q);
			lambdad[i] = 2.0/9.0/pi/sqrt(1.0-x1)*((1.0-5.0*z*z+p*p+ \
			         x3*x3)*Kk+(1.0-x1)*(z*z+7.0*p*p-4.0)*Ek-3.0*x3/x1*Pk);
	                    //ellpic_bulirsch(n,q)) FINDME: Check calc
            //printf("%14.10f\n",Pk);
			if(z < p) lambdad[i] += 2.0/3.0;
			if(fabs(p+z-1.0) < 1.0e-4)
			{
				lambdad[i] = 2.0/3.0/pi*acos(1.0-2.0*p)-4.0/9.0/pi* \
				            sqrt(p*(1.0-p))*(3.0+2.0*p-8.0*p*p);
			}
			etad[i] = p*p/2.0*(p*p+2.0*z*z);
		}
		IND(muo1, i) = 1.0-((1.0-u1-2.0*u2)*lambdae[i]+(u1+2.0*u2)*lambdad[i]+u2*etad[i])/omega;
		IND(mu0, i) = 1.0-lambdae[i];
	}
	
	free(lambdad);
	free(lambdae);
	free(etad);
	free(mu);
	Py_XDECREF(mu0);

	return PyArray_Return(muo1);
}
double min(double a, double b)
{
	if(a < b) return a;
	else return b;
}

double max(double a, double b)
{
	if(a > b) return a;
	else return b;
}
double rc(double x, double y)
{
	double rc, ERRTOL,TINY,SQRTNY,BIG,TNBG,COMP1,COMP2,THIRD,C1,C2, C3,C4;
	ERRTOL=0.04; TINY=1.69e-38; SQRTNY=1.3e-19; BIG=3.0e37;
	TNBG=TINY*BIG; COMP1=2.236/SQRTNY; COMP2=TNBG*TNBG/25.0;
	THIRD=1.0/3.0; C1=0.3; C2=1.0/7.0; C3=0.375; C4=9.0/22.0;

	double alamb,ave,s,w,xt,yt;
	if(x < 0.0 || y == 0.0 || (x+fabs(y)) < TINY || (x+fabs(y)) > BIG || (y < -COMP1 && x > 0 && x < COMP2)){
		printf("Invalid argument(s) in rc\n");
		return 0;
	}
	if(y > 0.0)
	{
		xt=x;
		yt=y;
		w=1.0;
	}
	else
	{
		xt=x-y;
		yt=-y;
		w=sqrt(x)/sqrt(xt);
	}
	s = ERRTOL*10.0;
	while(fabs(s) > ERRTOL)
	{
		alamb = 2.0*sqrt(xt)*sqrt(yt)+yt;
		xt = 0.25*(xt+alamb);
		yt  =0.25*(yt+alamb);
		ave = THIRD*(xt+yt+yt);
		s = (yt-ave)/ave;
	}
	rc = w*(1.0+s*s*(C1+s*(C2+s*(C3+s*C4))))/sqrt(ave);
	return rc;
}

double rj(double x, double y, double z, double p)
{
	double rj, ERRTOL,TINY,BIG,C1,C2,C3,C4,C5,C6,C7,C8, tempmax;
	
	ERRTOL=0.05; TINY=2.5e-13; BIG=9.0e11; C1=3.0/14.0;
	C2=1.0/3.0; C3=3.0/22.0; C4=3.0/26.0; C5=.750*C3;
     	C6=1.50*C4; C7=.50*C2; C8=C3+C3;
	
	double  a = 0.0,alamb,alpha,ave,b = 0.0,beta,delp,delx,dely,delz,ea,eb,ec,ed,ee, \
     		fac,pt,rcx = 0.0,rho,sqrtx,sqrty,sqrtz,sum,tau,xt,yt,zt;
      
	if(x < 0.0 || y < 0.0 || z < 0.0 || (x+y) < TINY || (x+z) < TINY || (y+z) < TINY || fabs(p) < TINY \
		|| x > BIG || y > BIG || z > BIG || fabs(p) > BIG)
	{
		printf("Invalid argument(s) in rj\n");
		return 0;
	}
	sum=0.0;
	fac=1.0;
	if(p > 0.0)
	{
		xt=x;
		yt=y;
		zt=z;
		pt=p;
	}
	else
	{
		xt = min(x, y);
		xt = min(xt,z);
		zt = max(x, y);
		zt = max(zt, z);
		yt = x+y+z-xt-zt;
		a = 1.0/(yt-p);
		b = a*(zt-yt)*(yt-xt);
		pt = yt+b;
		rho = xt*zt/yt;
		tau = p*pt/yt;
		rcx = rc(rho,tau);
	}
	tempmax = ERRTOL*10.0;
	while(tempmax > ERRTOL)
	{
		sqrtx = sqrt(xt);
		sqrty = sqrt(yt);
		sqrtz = sqrt(zt);
		alamb = sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
		alpha = pow((pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz), 2.0);
		beta = pt*(pt+alamb)*(pt + alamb);
		sum = sum+fac*rc(alpha,beta);
		fac = 0.25*fac;
		xt = 0.25*(xt+alamb);
		yt = 0.25*(yt+alamb);
		zt = 0.250*(zt+alamb);
		pt = 0.25*(pt+alamb);
		ave = 0.2*(xt+yt+zt+pt+pt);
		delx = (ave-xt)/ave;
		dely = (ave-yt)/ave;
		delz = (ave-zt)/ave;
		delp = (ave-pt)/ave;
		tempmax = max(fabs(delx), fabs(dely));
		tempmax = max(tempmax, fabs(delz));
		tempmax = max(tempmax, fabs(delp));
	}
	ea = delx*(dely+delz)+dely*delz;
	eb = delx*dely*delz;
	ec = delp*delp;
	ed = ea-3.0*ec;
	ee = eb+2.0*delp*(ea-ec);
	rj = 3.0*sum+fac*(1.0+ed*(-C1+C5*ed-C6*ee)+eb*(C7+delp*(-C8+delp*C4)) +\
		delp*ea*(C2-delp*C3)-C2*delp*ec)/(ave*sqrt(ave));
	if(p < 0.0) rj=a*(b*rj+3.0*(rcx-rf(xt,yt,zt)));
	return rj;  
}
	
double ellec(double k)
{
	double m1,a1,a2,a3,a4,b1,b2,b3,b4,ee1,ee2,ellec;
	// Computes polynomial approximation for the complete elliptic
	// integral of the second kind (Hasting's approximation):
	m1=1.0-k*k;
	a1=0.44325141463;
	a2=0.06260601220;
	a3=0.04757383546;
	a4=0.01736506451;
	b1=0.24998368310;
	b2=0.09200180037;
	b3=0.04069697526;
	b4=0.00526449639;
	ee1=1.0+m1*(a1+m1*(a2+m1*(a3+m1*a4)));
	ee2=m1*(b1+m1*(b2+m1*(b3+m1*b4)))*log(1.0/m1);
	ellec=ee1+ee2;
	return ellec;
}

double ellk(double k)
{
	double a0,a1,a2,a3,a4,b0,b1,b2,b3,b4,ellk, ek1,ek2,m1;
	// Computes polynomial approximation for the complete elliptic
	// integral of the first kind (Hasting's approximation):
	m1=1.0-k*k;
	a0=1.38629436112;
	a1=0.09666344259;
	a2=0.03590092383;
	a3=0.03742563713;
	a4=0.01451196212;
	b0=0.5;
	b1=0.12498593597;
	b2=0.06880248576;
	b3=0.03328355346;
	b4=0.00441787012;
	ek1=a0+m1*(a1+m1*(a2+m1*(a3+m1*a4)));
	ek2=(b0+m1*(b1+m1*(b2+m1*(b3+m1*b4))))*log(m1);
	ellk=ek1-ek2;
	return ellk;
}

double ellpic_bulirsch(double n, double k)
{
  int isdone;
  double output,kc,p,c,d,e,f,g,m0,pi;
  
  kc = sqrt(1.-pow(k,2)); 
  e  = kc;
  p  = sqrt(n+1.);
  d  = 1./p;
  c  = 1.;
  m0 = 1.;
  pi = acos(-1.0);
  
  isdone = 0;
  while (isdone == 0)
    {
      f  = c;
      c  = d/p+c;
      g  = e/p; 
      d  = 2.*(f*g+d);
      p  = g + p; 
      g  = m0; 
      m0 = kc + m0;
      
      if (fabs(1.-kc/g) > 1.e-8)
        {
          kc = 2*sqrt(e);
          e  = kc*m0;
        }
      else
        {
          output = 0.5*pi*(c*m0+d)/(m0*(m0+p));
          isdone = 1;
        }
    }
  
  return output;
}

double rf(double x, double y, double z)
{
	double rf, ERRTOL,TINY,BIG,THIRD,C1,C2,C3,C4, tempmax;
	
	ERRTOL=0.08; TINY=1.5e-38; BIG=3.0e37; THIRD=1.0/3.0;
	C1=1.0/24.0; C2=0.1; C3=3.0/44.0; C4=1.0/14.0;
	
	double alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt;

	if(min(x,y) < 0.0 || z < 0.0 || min(x+y, x+z) < TINY || y+z < TINY || max(x,y) > BIG || z > BIG)
	{
		printf("Invalid argument(s) in rf\n");
		return 0;
	}
	xt=x;
	yt=y;
	zt=z;
	tempmax = ERRTOL*10.0;
	while(tempmax > ERRTOL)
	{
		sqrtx=sqrt(xt);
		sqrty=sqrt(yt);
		sqrtz=sqrt(zt);
		alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
		xt=0.25*(xt+alamb);
		yt=0.25*(yt+alamb);
		zt=0.25*(zt+alamb);
		ave=THIRD*(xt+yt+zt);
		delx=(ave-xt)/ave;
		dely=(ave-yt)/ave;
		delz=(ave-zt)/ave;
		tempmax = max(fabs(delx), fabs(dely));
		tempmax = max(tempmax, fabs(delz));
	}
	e2=delx*dely-delz*delz;
	e3=delx*dely*delz;
	rf=(1.0+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave);
	return rf;
}

static char module_docstring[]="\
  FINDME: This function NEEDS OF DOC_STRING\n\
\n\
  Parameters\n\
  ----------\n\
        \n\
  Returns\n\
  -------\n\
\n\
  Revisions\n\
  ---------\n\
    Created: 2012-09-12  Laura Kreidberg \n\
    Updated: 2018-11-22  Jonathan Fraine Updated for Python, with support for Python2.7\n\
";

static PyMethodDef module_methods[] = {
  {"occultquad", occultquad,METH_VARARGS,module_docstring},{NULL}};

PyMODINIT_FUNC
#if PY_MAJOR_VERSION >= 3
    PyInit_occultquad(void)
#else
    initoccultquad(void)
#endif
{
    #if PY_MAJOR_VERSION >= 3
        PyObject *module;
        static struct PyModuleDef moduledef = {
            PyModuleDef_HEAD_INIT,
            "occultquad",             /* m_name */
            module_docstring,    /* m_doc */
            -1,                  /* m_size */
            module_methods,      /* m_methods */
            NULL,                /* m_reload */
            NULL,                /* m_traverse */
            NULL,                /* m_clear */
            NULL,                /* m_free */
        };
    #endif

    #if PY_MAJOR_VERSION >= 3
        module = PyModule_Create(&moduledef);
        if (!module)
            return NULL;
        /* Load `numpy` functionality. */
        import_array();
        return module;
    #else
        PyObject *m = Py_InitModule3("occultquad", module_methods, module_docstring);
        if (m == NULL)
            return;
        /* Load `numpy` functionality. */
        import_array();
    #endif
}
