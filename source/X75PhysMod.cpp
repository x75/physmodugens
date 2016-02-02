/*
 * X75PhysMod UGens
 * Oswald Berthold 20071107
 * GPL
 *
 * see Experiments in computer controlled acoustic modelling (a
 * step backward?) by Ron Berry
 * paper suggested by Manfred Hild / MatKlang 07/08
 *
 * XXX: make waveshaping coeffs settable via array and modulatable
 * XXX: loose string
 */

#include "SC_PlugIn.h"

#define PI 3.141592653589793
#define PI2 2*M_PI
#define MAX_DELAY_LENGTH SAMPLERATE // 1 second max. period duration

static InterfaceTable *ft;

// declare structs to hold unit generator state
// simple pluck unit (KS)
struct X75_Pluck1 : public Unit
{
  uint32 dlen;
  uint32 mPosition;
  float freq;
  float dur;
  float Q;
  float rho, S; // Lowpass
  float lastz;
  float lastf1, lastf2;
  float C;      // Allpass
  float *mData, *mImp; // delay buffer
};

// waveguide 1, abstracted delay module
struct X75_WG1 : public Unit
{
  uint32 dlen;
  uint32 mPosition;
  float *mData; // delay buffer
  // more state
  float p, p2;
  uint32 p1;
  float rho, fc, sigma;
  float lastz; //, lastzs;
  float lpx; // low-pass x
};

// 2D-Waveguide, 2-dimensional (plate) model
struct X75_WG2D : public Unit
{
  uint32 dlen;
  uint32 pos1, pos2;
  float *data1, *data2; // delay buffer
  float p1, p2;
  uint32 p1_1, p2_1;
  float p1_2, p2_2;
  float rho1, rho2;
  float fc1, sigma1;
  float fc2, sigma2;
  float lastz1; // , lastzs1;
  float lastz2; // , lastzs2; 
  float lpx1; // low-pass coeff
  float lpx2; // low-pass coeff
};

// 2Dplus-Waveguide
struct X75_WG2DPlus : public Unit
{
  uint32 dlen;
  uint32 pos1, pos2, pos3;
  float *data1, *data2, *data3; // delay buffer
  float p1, p2, p3; // loop length float
  uint32 p1_1, p2_1, p3_1; // loop length integer part
  float p1_2, p2_2, p3_2; // loop length float residual
  float rho1, rho2, rho3; // feedback, feedback 3 local
  float rho1l, rho2l, rho3l;
  float fc1, sigma1; // filter, inversion
  float fc2, sigma2; // filter, inversion
  float fc3, sigma3; // filter, inversion
  float lastz1, lastz2, lastz3; // , lastzs1;
  float lpx1, lpx2, lpx3; // low-pass coeff
};

// waveguide 3: trumpet model
// see R.W. Berry1988, icmc proceedings 1988, p.340
struct X75_WG3_Trumpet : public Unit
{
  uint32 dlen; // maximum delay length
  uint32 pos;  // current position
  float *data; // delay buffer
  // more state
  float p, p2; // float p, p2 is float residue floor(p)
  uint32 p1;   // integer part of p
  float rho, fc1, fc2, sigma; // feedback, LP-filter cutoff, HP-filter cutoff, inversion
  float lastz, lasty; // filter memory
  float lpx, hpx; // low-pass x, clalculated from cutoff
  float squamp, squoff; // squaring amplitude, squaring offset
};

// waveguide 4: saxophone/clarinet/flute model
// see R.W. Berry, icmc proceedings 1988, p.342
struct X75_WG4_Saxophone : public Unit
{
  uint32 dlen; // maximum delay length
  uint32 pos;  // current position
  float *data; // delay buffer
  // more state
  float p, p2; // float p, p2 is float residue floor(p)
  uint32 p1;   // integer part of p
  float rho1, rho2, sigma; // feedback, LP-filter cutoff, HP-filter
			   // cutoff, inversion
  float fc1, q1, a0, b1, b2; /* filter freq, Q and coeffs for resonant
				LP filter */
  float z1, y1, y2; // filter memory
  // float lpx, hpx; // low-pass x, clalculated from cutoff
  float squamp, squoff; // squaring amplitude, squaring offset
};

////////////////////////////////////////////////////////////
// waveshaping functions

// waveshaper
double X75_chebyshev(int n, double x)
{
  // XXX: nr, p. 237
  double T_n;       
  if (n == 0) T_n = 1.f;
  if (n == 1) T_n = x;
  if (n > 1) T_n = 2*x*X75_chebyshev((n-1),x) - X75_chebyshev((n-2),x);
  return T_n;
}

// WaveShaping
double X75_waveshape(double x, double *coeffs, int numcoeffs)
{
  double out = 0.0;
  int i;
  //out = chebyshev(1,x) + 0.8*chebyshev(2,x) + 0.6*chebyshev(3,x) + 0.4*chebyshev(4,x) + 0.2*chebyshev(5,x);
  // some sigmoid, static coefficients used below from X75PhysModTest.m
  if(x > 1.0) x = 1.0;
  if(x < -1.0) x = -1.0;
  for(i=0; i<numcoeffs; i++)
    out += coeffs[i] * X75_chebyshev(i,x);
  return out;
}

// declare unit generator functions
extern "C"
{
  void load(InterfaceTable *inTable);
  // simple pluck
  void X75_Pluck1_next(X75_Pluck1 *unit, int inNumSamples);
  float X75_Pluck1_Gf(X75_Pluck1* unit); // required gain for specified decay by X dB in tau
                                                                                 // seconds
  float X75_Pluck1_Gnorm(X75_Pluck1* unit); // actual gain at f
  void X75_Pluck1_rho_S(X75_Pluck1* unit); // shorten or lengthen?
  void X75_Pluck1_Ctor(X75_Pluck1* unit);
  void X75_Pluck1_Dtor(X75_Pluck1* unit);
  // generalized waveguide
  void X75_WG1_next(X75_WG1 *unit, int inNumSamples);
  void X75_WG1_Ctor(X75_WG1* unit);
  void X75_WG1_Dtor(X75_WG1* unit);
  float X75_WG1_Softlim(X75_WG1* unit, float z);
  void X75_WG1_setP(X75_WG1* unit, float p);
  void X75_WG1_setLpx(X75_WG1* unit, float fc);
  // 2D WaveGuide
  void X75_WG2D_next(X75_WG2D *unit, int inNumSamples);
  void X75_WG2D_Ctor(X75_WG2D* unit);
  void X75_WG2D_Dtor(X75_WG2D* unit);
  void X75_WG2D_setP(X75_WG2D* unit);
  void X75_WG2D_setLpx(X75_WG2D* unit);
  // 2Dplus WaveGuide
  void X75_WG2DPlus_next(X75_WG2DPlus *unit, int inNumSamples);
  void X75_WG2DPlus_Ctor(X75_WG2DPlus* unit);
  void X75_WG2DPlus_Dtor(X75_WG2DPlus* unit);
  void X75_WG2DPlus_setP(X75_WG2DPlus* unit);
  void X75_WG2DPlus_setLpx(X75_WG2DPlus* unit);
  // WG3_Trumpet model
  void X75_WG3_Trumpet_next(X75_WG3_Trumpet *unit, int inNumSamples);
  void X75_WG3_Trumpet_Ctor(X75_WG3_Trumpet* unit);
  void X75_WG3_Trumpet_Dtor(X75_WG3_Trumpet* unit);
  void X75_WG3_Trumpet_setP(X75_WG3_Trumpet* unit, float p);
  void X75_WG3_Trumpet_setLpx(X75_WG3_Trumpet* unit, float fc);
  void X75_WG3_Trumpet_setHpx(X75_WG3_Trumpet* unit, float fc);
  // WG4_Saxophone model
  void X75_WG4_Saxophone_next(X75_WG4_Saxophone *unit, int inNumSamples);
  void X75_WG4_Saxophone_Ctor(X75_WG4_Saxophone* unit);
  void X75_WG4_Saxophone_Dtor(X75_WG4_Saxophone* unit);
  void X75_WG4_Saxophone_setP(X75_WG4_Saxophone* unit, float p);
  void X75_WG4_Saxophone_setRLPFcoeffs(X75_WG4_Saxophone* unit, float fc, float q);
  //void X75_WG4_Saxophone_setHpx(X75_WG4_Saxophone* unit, float fc);
};

//////////////////////////////////////////////////////////////////
void X75_Pluck1_Ctor(X75_Pluck1* unit)
{
  float ex_mu = 0.f; // excitation mean
  float ex_max = 0;
  float loop, delay, D;
  int p;
  // seed PRNG
  //srand(); XXX: seed mit gettimeofday, oder random()
  // 1. set the calculation function.
  SETCALC(X75_Pluck1_next);
  // 1.a set up buffer
  //float fbufnum  = ZIN0(1);
  //uint32 bufnum = (int)fbufnum; 

  //World *world = unit->mWorld; 
  //if (bufnum >= world->mNumSndBufs) bufnum = 0; 
  //SndBuf *buf = world->mSndBufs + bufnum;
  
  // 2. initialize the unit generator state variables.
  // get the delay length
  unit->freq = ZIN0(0);
  unit->dlen = (uint32)(SAMPLERATE/unit->freq);// (uint32)ZIN0(0); // buf->frames; 
  unit->dur = ZIN0(1);
  unit->Q = 40.f; // decay by Q dB in dur secs
  unit->rho = ZIN0(2);
  unit->S = ZIN0(3);
  unit->C = ZIN0(4);

  X75_Pluck1_rho_S(unit);
  //  calculate C
  loop = SAMPLERATE/unit->freq;
  p = (int)loop;
  delay = p + unit->S;
  if(delay > loop)
    delay = --p + unit->S;
  D = loop - delay;
  unit->C = (1. - D)/(1. + D);

  unit->dlen = p;

  // memory
  unit->lastz = 0.f;
  unit->lastf1 = unit->lastf2 = 0.f;

  // allocate the buffer
  unit->mData = (float*)RTAlloc(unit->mWorld, unit->dlen * sizeof(float));
  // delay buffer initialisieren / anregung
  // anregung mittelwertfrei machen und normalisieren
  LOOP(unit->dlen,
       //unit->mData[xxi] = rand()/(RAND_MAX/2) - 1.f;
       unit->mData[xxi] = (rand()/(RAND_MAX+1.0)) * 2.0 - 1.0;
       ex_max = ex_max >= fabs(unit->mData[xxi]) ? ex_max : fabs(unit->mData[xxi]);
       ex_mu += unit->mData[xxi]);
  ex_mu = ex_mu / (unit->dlen + 0.f);
  ex_max = ex_max - ex_mu;
  LOOP(unit->dlen,
       //printf("u->mD: %f, ex_mu: %f, ex_max: %f\n", unit->mData[xxi], ex_mu, ex_max);
       unit->mData[xxi] = (unit->mData[xxi] - ex_mu);
       );
    
  unit->mPosition = 0;
  
  // debug out
  /* 
  printf("X75_Pluck1 Init: delbuf len: %d, rho: %f, S: %f, C: %f, excit.mean: %f, excit. max: %f\n",
	 unit->dlen, unit->rho,
	 unit->S, unit->C,
	 ex_mu, ex_max);
  */

  // calculate one sample of output.
  X75_Pluck1_next(unit, 1);
}

//////////////////////////////////////////////////////////////////
// Dtor is called to perform any clean up for the unit generator. 
void X75_Pluck1_Dtor(X75_Pluck1* unit)
{
  // free the buffer
  RTFree(unit->mWorld, unit->mData);
}

//////////////////////////////////////////////////////////////////
// calculate the required filter gain for f and tau
float X75_Pluck1_Gf(X75_Pluck1* unit)
{
  return(pow(10, -(unit->Q)/(20*unit->freq*unit->dur)));
}

//////////////////////////////////////////////////////////////////
// calculate the actual filter gain for f and tau
float X75_Pluck1_Gnorm(X75_Pluck1* unit)
{
  return(cos(PI*unit->freq/SAMPLERATE));
}

//////////////////////////////////////////////////////////////////
// calculate the actual filter gain for f and tau
void X75_Pluck1_rho_S(X75_Pluck1* unit)
{
  float gf, gnorm;
  float cosf1, a, b, c, D, a2, S1, S2;

  gf = X75_Pluck1_Gf(unit);
  gnorm = X75_Pluck1_Gnorm(unit);

  if(gnorm >= gf) {
    unit->rho = gf/gnorm;
    unit->S = 0.5;
  }
  else {
    unit->rho = 1.0;
    cosf1 = cos(2.*PI*unit->freq/SAMPLERATE);
    a = 2. - 2.*cosf1;

    b = 2.*cosf1 - 2.;
    c = 1. - gf*gf;
    D = sqrt(b*b - 4.*a*c);
    a2 = 2.*a;
    S1 = (-b + D)/a2;
    S2 = (-b - D)/a2;
    unit->S = S1; // XXX decide which one
  }
}

//////////////////////////////////////////////////////////////////
// calculation function when the buffer has been filled
void X75_Pluck1_next(X75_Pluck1 *unit, int inNumSamples)
{
  // get the pointer to the output buffer
  float *out = OUT(0);
  // get the pointer to the input buffer
  //float *in = IN(0);
  // get values from struct and store them in local variables.
  // The optimizer will cause them to be loaded it into a register.
  float *data = unit->mData;
  uint32 length = unit->dlen;
  uint32 position = unit->mPosition;
  float rho = unit->rho;
  float S = unit->S;
  float C = unit->C;
  // short-term memory
  float lastz = unit->lastz;
  float lastf1 = unit->lastf1;
  float lastf2 = unit->lastf2;
  float f1, f2;

  // perform a loop for the number of samples in the control period.
  // If this unit is audio rate then inNumSamples will be 64 or whatever
  // the block size is. If this unit is control rate then inNumSamples will
  // be 1.
  // float z_sum = 0.f; // sum of samples * coeffs
  for (int i=0; i < inNumSamples; ++i)
    {
      // get old value in delay line
      float z = data[position];
      int idx = 0;
      // store new value in delay line
      idx = position ? position : length-1;
      // Lowpass Filter mit \rho und S
      // data[position] = rho * ((1-S) * z + (data[idx-1]) * S); // in[i];
      f1 = rho * ((1-S) * z + (S * lastz)); // in[i];
      lastz = z;
      // Allpass Filter fuer das loop delay tuning
      f2 = C*f1 + lastf1 - C*lastf2;
      
      data[position] = f2;
      lastf1 = f1;
      lastf2 = f2;

      // debug out
      // printf("next: z: %f, z_sum: %f, cursamp: %f\n", z, z_sum, data[position]);

      // see if the position went to the end of the buffer 
      if (++position >= length) {
	position = 0; // go back to beginning
      }

      out[i] = z; // ; // z_sum;
    }
  // store the position back to the struct
  unit->mPosition = position;
  unit->lastz = lastz;
  unit->lastf1 = lastf1;
  unit->lastf2 = lastf2;
}

//////////////////////////////////////////////////////////////////
// WaveGuide 1
// - anregung durch input, nicht intern
void X75_WG1_Ctor(X75_WG1* unit)
{
  // 1. set the calculation function.
  SETCALC(X75_WG1_next);
  // 1.a set up buffer

  // 2. initialize the unit generator state variables.
  // get the delay length
  unit->dlen = (uint32)MAX_DELAY_LENGTH;
  // allocate the buffer
  unit->mData = (float*)RTAlloc(unit->mWorld, unit->dlen * sizeof(float));
  LOOP(unit->dlen, unit->mData[xxi] = 0.f;);

  unit->p = ZIN0(1);
  X75_WG1_setP(unit, unit->p);
  unit->rho = ZIN0(2);
  //unit->fc = ZIN0(3);
  X75_WG1_setLpx(unit, ZIN0(3));
  unit->sigma = ZIN0(4);


  unit->mPosition = 0;
  unit->lastz = 0.f;
  
  // debug out
  /* 
  printf("X75_WG1 Init: delbuf maxlen: %d, eff. len: %f, p1: %d, p2: %f, rho: %f, fc: %f, lpx: %f, sigma: %f\n",
	 unit->dlen, unit->p,
	 unit->p1, unit->p2, unit->rho, unit->fc, unit->lpx, unit->sigma);
  */
  // calculate one sample of output.
  X75_WG1_next(unit, 1);
}

//////////////////////////////////////////////////////////////////
// Dtor is called to perform any clean up for the unit generator. 
void X75_WG1_Dtor(X75_WG1* unit)
{
  // free the buffer
  RTFree(unit->mWorld, unit->mData);
}

//////////////////////////////////////////////////////////////////
// calculate softlim function
// XXX: should be turned into full waveshaping module
float X75_WG1_Softlim(X75_WG1 *unit, float z)
{
  float slope, intercept;
  float b[6];
  float lim = 1.2;
  float s[5]; // slopes
  float absz = fabs(z);
  float ret;
  b[0] = 0.7;
  float step = (lim - b[0])/5;
  for(int i = 1; i < 6; i++) {
    b[i] = b[i-1] + 0.1 * s[i-1];
  }
  s[0] = 0.9; s[1] = 0.8; s[2] = 0.6; s[3] = 0.4; s[4] = 0.3;

  if(absz <= b[0])
    ret = absz;
  else {
    if(b[0] < absz && absz <= b[1])
      ret = b[0] + (absz-b[0]) * s[0];
    else if(b[1] < absz && absz <= b[2])
      ret = b[1] + (absz-b[1]) * s[1];
    else if(b[2] < absz && absz <= b[3])
      ret = b[2] + (absz-b[2]) * s[2];
    else if(b[3] < absz && absz <= b[4])
      ret = b[3] + (absz-b[3]) * s[3];
    else if(b[4] < absz && absz <= b[5])
      ret = b[4] + (absz-b[4]) * s[4];
    else
      ret = 1.;
  }
  if(z < 0)
    return(-ret);
  return(ret);
}

//////////////////////////////////////////////////////////////////
// calculate p1,p2 from p
void X75_WG1_setP(X75_WG1* unit, float p)
{
  unit->p1 = (uint32)p;
  unit->p2 = p - unit->p1;
}

//////////////////////////////////////////////////////////////////
// calculate lpx from fc
void X75_WG1_setLpx(X75_WG1* unit, float fc)
{
  unit->fc = fc;
  unit->lpx = exp(-2 * PI * (unit->fc/SAMPLERATE));
}

//////////////////////////////////////////////////////////////////
// calculation function when the buffer has been filled
void X75_WG1_next(X75_WG1 *unit, int inNumSamples)
{
  // get the pointer to the output buffer
  float *out = OUT(0);
  // get the pointer to the input buffer
  float *in = IN(0);
  // get values from struct and store them in local variables.
  // The optimizer will cause them to be loaded it into a register.
  // get new values
  float p = ZIN0(1);
  if(p != unit->p) {
    X75_WG1_setP(unit, p);
    unit->p = p;
  }

  float rho = ZIN0(2);
  if(rho != unit->rho)
    unit->rho = rho;

  float fc = ZIN0(3);
  if(fc != unit->fc)
    X75_WG1_setLpx(unit, fc);

  float sigma = ZIN0(4); //unit->sigma;
  if(sigma != unit->sigma)
    unit->sigma = sigma;

  float *data = unit->mData;
  uint32 length = unit->p1; // unit->dlen;
  uint32 pos2 = unit->mPosition;
  uint32 pos1 = pos2 + 1;
  float p2 = unit->p2;
  float lastz = unit->lastz;
  float z; // z -> y_n
  float zs; // zs -> y_n schlange
  float lpx = unit->lpx;

  // test changes

  // perform a loop for the number of samples in the control period.
  for (int i=0; i < inNumSamples; ++i)
    {
      // x_n-p1-1
      if(pos1 >= length)
      //if(pos2 < 0)
	pos1 = 0;

      data[pos1] += in[i]; // get excitation into the loop

      // FIR interpolation
      zs = (1-p2) * data[pos2] + p2 * data[pos1];
      // Low Pass Filter
      z = (1-lpx)*zs + lpx * lastz; // lowpass aus dspguide, ch19, S.323, lastz statt lastzs
      lastz = z;
      // store new value in delay line
      data[pos2] = X75_WG1_Softlim(unit, (rho * lastz)); // mixer
      //data[pos1] = in[i] + (rho * lastz); // mixer

      // debug out
      //printf("next: z: %f, z_sum: %f, cursamp: %f\n", z, z_sum, data[position]);

      // see if the position went to the end of the buffer
      if (++pos2 >= length) {
	pos2 = 0; // go back to beginning
      }
      pos1 = pos2+1;
      
      // debug out
      //printf("next: z_sum: %f\n", z_sum);
      if(sigma > 0.0)
	z = -z;
      out[i] = z;
    }
  // store the position back to the struct
  unit->mPosition = pos2;
  unit->lastz = lastz;
}

//////////////////////////////////////////////////////////////////
// 2D WaveGuide
// excitation from external input
void X75_WG2D_Ctor(X75_WG2D* unit)
{
  SETCALC(X75_WG2D_next);

  // set max delay length
  unit->dlen = (uint32)MAX_DELAY_LENGTH;
  // allocate the buffers and initialize them to 0
  unit->data1 = (float*)RTAlloc(unit->mWorld, unit->dlen * sizeof(float));
  LOOP(unit->dlen, unit->data1[xxi] = 0.f;);
  unit->data2 = (float*)RTAlloc(unit->mWorld, unit->dlen * sizeof(float));
  LOOP(unit->dlen, unit->data2[xxi] = 0.f;);

  // base cycle length
  unit->p1 = ZIN0(1);
  unit->p2 = ZIN0(5);
  X75_WG2D_setP(unit);

  // feeback
  unit->rho1 = ZIN0(2);
  unit->rho2 = ZIN0(6);

  // LP cutoff in Hz
  unit->fc1 = ZIN0(3);
  unit->fc2 = ZIN0(7);
  X75_WG2D_setLpx(unit);

  // invert delay output
  unit->sigma1 = ZIN0(4);
  unit->sigma2 = ZIN0(8);

  // delay line management
  unit->pos1 = unit->pos2 = 0;
  // signals
  unit->lastz1 = 0.f;
  unit->lastz2 = 0.f;
  
  // debug out
  /* 
  printf("X75_WG2D Init: delbuf maxlen: %d, eff-len1: %f, efflen2: %f, p1_1: %d, p1_2: %f, p2_1: %d, p2_2: %f, fc1: %f, fc2: %f, lpx1: %f, lpx2: %f,s, sigma1: %f, sigma2: %f\n",
	 unit->dlen, unit->p1, unit->p2,
	 unit->p1_1, unit->p1_2, unit->p2_1, unit->p2_2, unit->fc1, unit->fc2, unit->lpx1, unit->lpx2, unit->sigma1, unit->sigma2);
  */
  // init
  X75_WG2D_next(unit, 1);
}

// Dtroy everything
void X75_WG2D_Dtor(X75_WG2D* unit)
{
  // free the buffers
  RTFree(unit->mWorld, unit->data1);
  RTFree(unit->mWorld, unit->data2);
}

// calculate p1,p2 from p
void X75_WG2D_setP(X75_WG2D* unit)
{
  uint32 tmp;
  //  printf("p1: %f, p2: %f\n", unit->p1, unit->p2);
  tmp = (uint32)unit->p1;
  unit->p1_1 = tmp;
  unit->p1_2 = unit->p1 - tmp;
  tmp = (uint32)unit->p2;
  unit->p2_1 = tmp;
  unit->p2_2 = unit->p2 - tmp;
  printf("p1: %f, p2: %f\n", unit->p1, unit->p2);
}

// calculate lpx from fc
void X75_WG2D_setLpx(X75_WG2D* unit)
{
  //unit->fc1 = fc;
  unit->lpx1 = exp(-2 * PI * (unit->fc1/SAMPLERATE));
  unit->lpx2 = exp(-2 * PI * (unit->fc2/SAMPLERATE));
}

// calculation function when the buffer has been filled
void X75_WG2D_next(X75_WG2D *unit, int inNumSamples)
{
  // get the pointer to the output buffer
  float *out = OUT(0);
  // get the pointer to the input buffer
  float *in = IN(0);
  // get values from struct and store them in local variables.
  // The optimizer will cause them to be loaded it into a register.
  // get new values
  float p1 = ZIN0(1);
  float p2 = ZIN0(5);

  if((p1 != unit->p1) || (p2 != unit->p2)) { // check if p arguments changed
    unit->p1 = p1;
    unit->p2 = p2;
    X75_WG2D_setP(unit);
  }

  float rho1 = ZIN0(2);
  float rho2 = ZIN0(6);
  if(rho1 != unit->rho1 || rho2 != unit->rho2) {  // update sigma
    unit->rho1 = rho1;
    unit->rho2 = rho2;
  }

  float fc1 = ZIN0(3);
  float fc2 = ZIN0(7);
  if(fc1 != unit->fc1 || fc2 != unit->fc2) { // update cutoff freq(s)
    unit->fc1 = fc1;
    unit->fc2 = fc2;
    X75_WG2D_setLpx(unit);
  }

  float sigma1 = ZIN0(4);
  float sigma2 = ZIN0(8);
  if(sigma1 != unit->sigma1 || sigma2 != unit->sigma2) {  // update sigma
    unit->sigma1 = sigma1;
    unit->sigma2 = sigma2;
  }

  float *data1 = unit->data1;
  float *data2 = unit->data2;

  uint32 length1 = unit->p1_1; // unit->dlen;
  uint32 length2 = unit->p2_1; // unit->dlen;

  float p1_2 = unit->p1_2;
  float p2_2 = unit->p2_2;

  uint32 pos1_1 = unit->pos1;
  uint32 pos1_2 = pos1_1 + 1;
  uint32 pos2_1 = unit->pos2;
  uint32 pos2_2 = pos2_1 + 1;

  //float rho1 = unit->rho1;
  //float rho2 = unit->rho2;

  float lastz1 = unit->lastz1;
  //float lastzs1 = unit->lastzs1;
  float lastz2 = unit->lastz2;
  //float lastzs2 = unit->lastzs2;

  float z1, zs1; // z -> y_n
  float z2, zs2; // zs -> y_n schlange
  float mix = 0.f;
  float lpx1 = unit->lpx1;
  float lpx2 = unit->lpx2;

  int ws_numcoeffs = 4;
  double ws_coeffs[ws_numcoeffs];
  ws_coeffs[0] = 0;
  ws_coeffs[1] = -1.011; // -1.0109;
  ws_coeffs[2] = 0;
  ws_coeffs[3] = 0.122; //1.2221e-01;

  // test changes

  // perform a loop for the number of samples in the control period.
  for (int i=0; i < inNumSamples; ++i)
    {
      // x_n-p1-1
      if(pos1_2 >= length1)
	pos1_2 = 0;
      if(pos2_2 >= length2)
	pos2_2 = 0;

      data1[pos1_1] += in[i];
      data2[pos2_1] += in[i];
      // FIR interpolation
      zs1 = (1-p1_2) * data1[pos1_1] + p1_2 * data1[pos1_2];
      zs2 = (1-p2_2) * data2[pos2_1] + p2_2 * data2[pos2_2];
      // LowPass Filter
      z1 = (1-lpx1)*zs1 + lpx1 * lastz1; // lowpass aus dspguide, ch19, S.323, lastz statt lastzs
      z2 = (1-lpx2)*zs2 + lpx2 * lastz2; // lowpass aus dspguide, ch19, S.323, lastz statt lastzs

      // invert or not invert
      if(sigma1 > 0.0)
	z1 = -z1;
      if(sigma2 > 0.0)
	z2 = -z2;

      // store
      lastz1 = z1;
      lastz2 = z2;

      //mix = in[i] + rho1 * lastzs1 + rho2 * lastzs2;
      mix = X75_waveshape(rho1 * z1 + rho2 * z2, ws_coeffs, ws_numcoeffs);
      // store new value in delay line
      data1[pos1_1] = mix; //
      data2[pos2_1] = mix; // 
      //data[pos1] = in[i] + (rho * lastz); // mixer
      // softlim: hm ..

      // debug out
      //printf("next: z: %f, z_sum: %f, cursamp: %f\n", z, z_sum, data[position]);

      // see if the position went to the end of the buffer 
      if (++pos1_1 >= length1) {
	pos1_1 = 0; // go back to beginning
      }
      pos1_2 = pos1_1+1;
      // same for 2nd dimension
      if (++pos2_1 >= length2) {
	pos2_1 = 0; // go back to beginning
      }
      pos2_2 = pos2_1+1;
      
      // debug out
      //printf("next: z_sum: %f\n", z_sum);
      out[i] = mix; // z1 + z2; // write it out
      //lastzs1 = z1;
      //lastzs2 = z2;
    }
  // store the position back to the struct
  unit->pos1 = pos1_1;
  unit->pos2 = pos2_1;
  unit->lastz1 = lastz1;
  //  unit->lastzs1 = lastzs1;
  unit->lastz2 = lastz2;
  // unit->lastzs2 = lastzs2;
}

//////////////////////////////////////////////////////////////////
// 2DPlus WaveGuide
// excited by external signal
void X75_WG2DPlus_Ctor(X75_WG2DPlus* unit)
{
  // set the calculation function.
  SETCALC(X75_WG2DPlus_next);

  // init stuff
  // get the delay length
  unit->dlen = (uint32)MAX_DELAY_LENGTH;
  // allocate the buffers and initialize them to 0.0
  unit->data1 = (float*)RTAlloc(unit->mWorld, unit->dlen * sizeof(float));
  LOOP(unit->dlen, unit->data1[xxi] = 0.f;);
  unit->data2 = (float*)RTAlloc(unit->mWorld, unit->dlen * sizeof(float));
  LOOP(unit->dlen, unit->data2[xxi] = 0.f;);
  unit->data3 = (float*)RTAlloc(unit->mWorld, unit->dlen * sizeof(float));
  LOOP(unit->dlen, unit->data3[xxi] = 0.f;);

  unit->p1 = ZIN0(1);
  unit->p2 = ZIN0(5);
  unit->p3 = ZIN0(9);
  X75_WG2DPlus_setP(unit);

  unit->rho1 = ZIN0(2);
  unit->rho2 = ZIN0(6);
  unit->rho3 = ZIN0(10);

  unit->fc1 = ZIN0(3);
  unit->fc2 = ZIN0(7);
  unit->fc3 = ZIN0(11);
  X75_WG2DPlus_setLpx(unit);

  unit->sigma1 = ZIN0(4);
  unit->sigma2 = ZIN0(8);
  unit->sigma3 = ZIN0(12);

  unit->rho1l = ZIN0(13);
  unit->rho2l = ZIN0(14);
  unit->rho3l = ZIN0(15);

  unit->pos1 = unit->pos2 = unit->pos3 = 0;
  unit->lastz1 = 0.f; // unit->lastzs1 = 0.f;
  unit->lastz2 = 0.f; // unit->lastzs2 = 0.f;
  unit->lastz3 = 0.f; // unit->lastzs2 = 0.f;
  
  // debug out
//   printf("X75_WG2DPlus Init: delbuf maxlen: %d, eff-len1: %f, efflen2: %f, p1_1: %d, p1_2: %f, p2_1: %d, p2_2: %f, fc1: %f, fc2: %f, lpx1: %f, lpx2: %f,s, sigma1: %f, sigma2: %f, etc ..\n",
// 	 unit->dlen, unit->p1, unit->p2,
// 	 unit->p1_1, unit->p1_2, unit->p2_1, unit->p2_2, unit->fc1, unit->fc2, unit->lpx1, unit->lpx2, unit->sigma1, unit->sigma2);
  /* for(int i = 0; i < unit->dlen; i++)
    printf("%f, %f, ", unit->mImp[i], buf->data[i]);
    printf("\n"); */
  // calculate one sample of output.
  X75_WG2DPlus_next(unit, 1);
}

//////////////////////////////////////////////////////////////////
// Dtor is called to perform any clean up for the unit generator. 
void X75_WG2DPlus_Dtor(X75_WG2DPlus* unit)
{
  // free the buffer
  RTFree(unit->mWorld, unit->data1);
  RTFree(unit->mWorld, unit->data2);
  RTFree(unit->mWorld, unit->data3);
}

//////////////////////////////////////////////////////////////////
// calculate p1,p2 from p
void X75_WG2DPlus_setP(X75_WG2DPlus* unit)
{
  uint32 tmp;
  //  printf("p1: %f, p2: %f\n", unit->p1, unit->p2);
  tmp = (uint32)unit->p1;
  unit->p1_1 = tmp;
  unit->p1_2 = unit->p1 - tmp;
  tmp = (uint32)unit->p2;
  unit->p2_1 = tmp;
  unit->p2_2 = unit->p2 - tmp;
  tmp = (uint32)unit->p3;
  unit->p3_1 = tmp;
  unit->p3_2 = unit->p3 - tmp;
  //  printf("p1: %f, p2: %f, p3: %f\n", unit->p1, unit->p2, unit->p3);
}

//////////////////////////////////////////////////////////////////
// calculate lpx from fc
void X75_WG2DPlus_setLpx(X75_WG2DPlus* unit)
{
  //unit->fc1 = fc;
  unit->lpx1 = exp(-2 * PI * (unit->fc1/SAMPLERATE));
  unit->lpx2 = exp(-2 * PI * (unit->fc2/SAMPLERATE));
  unit->lpx3 = exp(-2 * PI * (unit->fc3/SAMPLERATE));
}

//////////////////////////////////////////////////////////////////
// calculation function when the buffer has been filled
void X75_WG2DPlus_next(X75_WG2DPlus *unit, int inNumSamples)
{
  // get the pointer to the output buffer
  float *out = OUT(0);
  // get the pointer to the input buffer
  float *in = IN(0);
  // get values from struct and store them in local variables.
  // The optimizer will cause them to be loaded it into a register.
  // get new values
  float p1 = ZIN0(1);
  float p2 = ZIN0(5);
  float p3 = ZIN0(9);

  if((p1 != unit->p1) || (p2 != unit->p2) || (p3 != unit->p3)) { // check if p arguments changed
    unit->p1 = p1;
    unit->p2 = p2;
    unit->p3 = p3;
    X75_WG2DPlus_setP(unit);
  }

  float rho1 = ZIN0(2);
  float rho2 = ZIN0(6);
  float rho3 = ZIN0(10);
  if(rho1 != unit->rho1 || rho2 != unit->rho2 || rho3 != unit->rho3) {  // update sigma
    unit->rho1 = rho1;
    unit->rho2 = rho2;
    unit->rho3 = rho3;
  }

  float fc1 = ZIN0(3);
  float fc2 = ZIN0(7);
  float fc3 = ZIN0(11);
  if(fc1 != unit->fc1 || fc2 != unit->fc2 || fc3 != unit->fc3) { // update cutoff freq(s)
    unit->fc1 = fc1;
    unit->fc2 = fc2;
    unit->fc3 = fc3;
    X75_WG2DPlus_setLpx(unit);
  }

  float sigma1 = ZIN0(4);
  float sigma2 = ZIN0(8);
  float sigma3 = ZIN0(12);
  if(sigma1 != unit->sigma1 || sigma2 != unit->sigma2 || sigma3 != unit->sigma3) {  // update sigma
    unit->sigma1 = sigma1;
    unit->sigma2 = sigma2;
    unit->sigma3 = sigma3;
  }

  float rho1l = ZIN0(13);
  float rho2l = ZIN0(14);
  float rho3l = ZIN0(15);
  if(rho1l != unit->rho1l || rho2l != unit->rho2l || rho3l != unit->rho3l) {  // update sigma
    unit->rho1l = rho1l;
    unit->rho2l = rho2l;
    unit->rho3l = rho3l;
  }

  float *data1 = unit->data1;
  float *data2 = unit->data2;
  float *data3 = unit->data3;

  uint32 length1 = unit->p1_1; // unit->dlen;
  uint32 length2 = unit->p2_1; // unit->dlen;
  uint32 length3 = unit->p3_1; // unit->dlen;

  float p1_2 = unit->p1_2;
  float p2_2 = unit->p2_2;
  float p3_2 = unit->p3_2;

  uint32 pos1_1 = unit->pos1;
  uint32 pos1_2 = pos1_1 + 1;
  uint32 pos2_1 = unit->pos2;
  uint32 pos2_2 = pos2_1 + 1;
  uint32 pos3_1 = unit->pos3;
  uint32 pos3_2 = pos3_1 + 1;

  //float rho1 = unit->rho1;
  //float rho2 = unit->rho2;

  float lastz1 = unit->lastz1;
  //float lastzs1 = unit->lastzs1;
  float lastz2 = unit->lastz2;
  float lastz3 = unit->lastz3;
  //float lastzs2 = unit->lastzs2;

  float z1, zs1; // z -> y_n
  float z2, zs2; // zs -> y_n schlange
  float z3, zs3; // zs -> y_n schlange

  float lpx1 = unit->lpx1;
  float lpx2 = unit->lpx2;
  float lpx3 = unit->lpx3;

  float mix = 0.f;
  float mix1 = 0.f;
  float mix2 = 0.f;

  // perform a loop for the number of samples in the control period.
  for (int i=0; i < inNumSamples; ++i)
    {
      // x_n-p1-1
      if(pos1_2 >= length1)
	pos1_2 = 0;
      if(pos2_2 >= length2)
	pos2_2 = 0;
      if(pos3_2 >= length3)
	pos3_2 = 0;

      mix1 = in[i]; //  + rho1 * lastz1 + rho2 * lastz2 + rho3 * lastz3;
      mix2 = rho1l * lastz1 + rho2l * lastz2; //  + rho3l * lastz3;

      data1[pos1_1] += mix1;
      data2[pos2_1] += mix1;
      data3[pos3_1] += mix2; // no
      // FIR interpolation
      zs1 = (1-p1_2) * data1[pos1_1] + p1_2 * data1[pos1_2];
      zs2 = (1-p2_2) * data2[pos2_1] + p2_2 * data2[pos2_2];
      zs3 = (1-p3_2) * data3[pos3_1] + p3_2 * data3[pos3_2];
      // LowPass Filter
      z1 = (1-lpx1)*zs1 + lpx1 * lastz1; // lowpass from dspguide, ch19, p.323, lastz statt lastzs
      z2 = (1-lpx2)*zs2 + lpx2 * lastz2; // 
      z3 = (1-lpx3)*zs3 + lpx3 * lastz3; // 

      // invert or not invert
      if(sigma1 > 0.0)
	z1 = -z1;
      if(sigma2 > 0.0)
	z2 = -z2;
      if(sigma3 > 0.0)
	z3 = -z3;

      // remember
      lastz1 = z1;
      lastz2 = z2;
      lastz3 = z3;

      //mix = in[i] + rho1 * lastzs1 + rho2 * lastzs2;
      mix1 = rho1 * z1 + rho2 * z2 + rho3 * z3; // XXX: X75_waveshape, set coeffs first
      mix2 = rho3l * z3;
      // store new value in delay line
      data1[pos1_1] = mix1; //
      data2[pos2_1] = mix1; // 
      data3[pos3_1] = mix2; // 

      //data[pos1] = in[i] + (rho * lastz); // mixer
      // softlim: hm ..

      // debug out
      //printf("next: z: %f, z_sum: %f, cursamp: %f\n", z, z_sum, data[position]);

      // see if the position went to the end of the buffer 
      if (++pos1_1 >= length1) {
	pos1_1 = 0; // go back to beginning
      }
      pos1_2 = pos1_1+1;
      // same for 2nd dimension
      if (++pos2_1 >= length2) {
	pos2_1 = 0; // go back to beginning
      }
      pos2_2 = pos2_1+1;
      // same for 2nd dimension
      if (++pos3_1 >= length3) {
	pos3_1 = 0; // go back to beginning
      }
      pos3_2 = pos3_1+1;
      
      // debug out
      //printf("next: z_sum: %f\n", z_sum);
      out[i] = mix2; // z1 + z2; // write it out
      //lastzs1 = z1;
      //lastzs2 = z2;
    }
  // store the position back to the struct
  unit->pos1 = pos1_1;
  unit->pos2 = pos2_1;
  unit->pos3 = pos3_1;
  unit->lastz1 = lastz1;
  //  unit->lastzs1 = lastzs1;
  unit->lastz2 = lastz2;
  unit->lastz3 = lastz3;
  // unit->lastzs2 = lastzs2;
}

//////////////////////////////////////////////////////////////////
// WaveGuide 3, Trumpet model
// - excitation through external enveloped, modulated oscillator input
void X75_WG3_Trumpet_Ctor(X75_WG3_Trumpet* unit)
{
  SETCALC(X75_WG3_Trumpet_next);

  // get the delay length
  unit->dlen = (uint32)MAX_DELAY_LENGTH;
  // allocate the buffer and init to zero
  unit->data = (float*)RTAlloc(unit->mWorld, unit->dlen * sizeof(float));
  LOOP(unit->dlen, unit->data[xxi] = 0.f;);

  unit->p = ZIN0(1);
  X75_WG3_Trumpet_setP(unit, unit->p);
  unit->rho = ZIN0(2);
  //unit->fc = ZIN0(3);
  X75_WG3_Trumpet_setLpx(unit, ZIN0(3));
  X75_WG3_Trumpet_setHpx(unit, ZIN0(4));
  unit->squamp = ZIN0(5);
  unit->squoff = ZIN0(6);
  unit->sigma = ZIN0(7);


  unit->pos = 0;
  unit->lastz = unit->lasty = 0.f;
  
  // debug out
  /* 
  printf("X75_WG3_Trumpet Init: delbuf maxlen: %d, eff. len: %f, p1: %d, p2: %f, rho: %f, fc: %f, lpx: %f, sigma: %f\n",
	 unit->dlen, unit->p,
	 unit->p1, unit->p2, unit->rho, unit->fc, unit->lpx, unit->sigma);
  */
  // calculate one sample of output.
  X75_WG3_Trumpet_next(unit, 1);
}

//////////////////////////////////////////////////////////////////
// Dtor is called to perform any clean up for the unit generator. 
void X75_WG3_Trumpet_Dtor(X75_WG3_Trumpet* unit)
{
  // free the buffer
  RTFree(unit->mWorld, unit->data);
}

//////////////////////////////////////////////////////////////////
// calculate p1,p2 from p
void X75_WG3_Trumpet_setP(X75_WG3_Trumpet* unit, float p)
{
  unit->p1 = (uint32)p;
  unit->p2 = p - unit->p1;
}

// calculate lpx from fc1
void X75_WG3_Trumpet_setLpx(X75_WG3_Trumpet* unit, float fc)
{
  unit->fc1 = fc;
  unit->lpx = exp(-2 * PI * (unit->fc1/SAMPLERATE));
}

// calculate hpx from fc2
void X75_WG3_Trumpet_setHpx(X75_WG3_Trumpet* unit, float fc)
{
  unit->fc2 = fc;
  unit->hpx = exp(-2 * PI * (unit->fc2/SAMPLERATE));
}

//////////////////////////////////////////////////////////////////
// calculation function when the buffer has been filled
void X75_WG3_Trumpet_next(X75_WG3_Trumpet *unit, int inNumSamples)
{
  // get output pointer
  float *out = OUT(0);
  // get input pointer
  float *in = IN(0);

  // localize and update state values
  float p = ZIN0(1);
  if(p != unit->p) {
    X75_WG3_Trumpet_setP(unit, p);
    unit->p = p;
  }

  float rho = ZIN0(2);
  if(rho != unit->rho)
    unit->rho = rho;

  float fc1 = ZIN0(3);
  if(fc1 != unit->fc1)
    X75_WG3_Trumpet_setLpx(unit, fc1);

  float fc2 = ZIN0(4);
  if(fc2 != unit->fc2)
    X75_WG3_Trumpet_setHpx(unit, fc2);

  float squamp = ZIN0(5); //unit->squamp;
  if(squamp != unit->squamp)
    unit->squamp = squamp;

  float squoff = ZIN0(6); //unit->squoff;
  if(squoff != unit->squoff)
    unit->squoff = squoff;

  float sigma = ZIN0(7); //unit->sigma;
  if(sigma != unit->sigma)
    unit->sigma = sigma;

  float *data = unit->data;
  uint32 length = unit->p1; // unit->dlen;
  uint32 pos1 = unit->pos;
  uint32 pos2 = pos1 + 1;
  float p2 = unit->p2;
  float lastz = unit->lastz;
  float lasty = unit->lasty;
  float z; // z -> y_n
  float zs; // zs -> y_n schlange
  float lpx = unit->lpx;
  float hpx = unit->hpx;
  float a0, a1, mix;

  int ws_numcoeffs = 4;
  double ws_coeffs[ws_numcoeffs];
  ws_coeffs[0] = 0;
  ws_coeffs[1] = -1.011; // -1.0109;
  ws_coeffs[2] = 0;
  ws_coeffs[3] = 0.122; //1.2221e-01;

  a0 = a1 = 0.5*(1+hpx); // XXX calc only when changed

  // perform dsp block loop
  for (int i=0; i < inNumSamples; ++i)
    {
      // x_n-p1-1
      if(pos2 >= length)
      //if(pos2 < 0)
	pos2 = 0;

      mix = data[pos1] + in[i];
      data[pos1] = X75_waveshape(mix * (mix * squamp + squoff), ws_coeffs, ws_numcoeffs         ); // get excitation into the loop

      // FIR interpolation
      zs = (1-p2) * data[pos1] + p2 * data[pos2];
      // Low Pass Filter
      z = (1-lpx) * zs + lpx * lastz; // lowpass aus dspguide, ch19, S.323, lastz statt lastzs
      // store new value in delay line
      //data[pos1] = X75_waveshape(rho * lastz); // mixer
      data[pos1] = rho * z; // mixer
      //data[pos1] = in[i] + (rho * lastz); // mixer

      // debug out
      //printf("next: z: %f, z_sum: %f, cursamp: %f\n", z, z_sum, data[position]);

      // see if the position went to the end of the buffer
      if (++pos1 >= length) {
	pos1 = 0; // go back to beginning
      }
      pos2 = pos1 + 1;
      
      // debug out
      //printf("next: z_sum: %f\n", z_sum);
      if(sigma > 0.0)
	z = -z;
      out[i] = lasty = a0 * z - a1 * lastz + hpx * lasty;
      lastz = z;
    }
  // store the position back to the struct
  unit->pos = pos1;
  unit->lastz = lastz;
  unit->lasty = lasty;
}

//////////////////////////////////////////////////////////////////
// WaveGuide #4, Saxophone model
// - excitation through external enveloped, modulated oscillator input
void X75_WG4_Saxophone_Ctor(X75_WG4_Saxophone* unit)
{
  SETCALC(X75_WG4_Saxophone_next);

  // get the delay length
  unit->dlen = (uint32)MAX_DELAY_LENGTH;
  unit->pos = 0; // initial pointer position
  // initialize the buffer
  unit->data = (float*)RTAlloc(unit->mWorld, unit->dlen * sizeof(float));
  LOOP(unit->dlen, unit->data[xxi] = 0.f;);

  unit->p = ZIN0(2); // get p
  X75_WG4_Saxophone_setP(unit, unit->p); // set p1,p2
  unit->rho1 = ZIN0(3);
  unit->rho2 = ZIN0(4);
  X75_WG4_Saxophone_setRLPFcoeffs(unit, ZIN0(5), ZIN0(6));
  // p and fc1 should the same according to schema
  //X75_WG4_Saxophone_setHpx(unit, ZIN0(4));
  unit->squamp = ZIN0(8);
  unit->squoff = ZIN0(9);
  unit->sigma = ZIN0(10);


  unit->z1 = unit->y1 = unit->y2 = 0.f;
  
  // debug out
  /* 
  printf("X75_WG4_Saxophone Init: delbuf maxlen: %d, eff. len: %f, p1: %d, p2: %f, rho: %f, fc: %f, lpx: %f, sigma: %f\n",
	 unit->dlen, unit->p,
	 unit->p1, unit->p2, unit->rho, unit->fc, unit->lpx, unit->sigma);
  */
  // calculate one sample of output.
  X75_WG4_Saxophone_next(unit, 1);
}

// Dtor is called to perform any clean up for the unit generator. 
void X75_WG4_Saxophone_Dtor(X75_WG4_Saxophone* unit)
{
  // free the buffer
  RTFree(unit->mWorld, unit->data);
}

// calculate p1,p2 from p
void X75_WG4_Saxophone_setP(X75_WG4_Saxophone* unit, float p)
{
  unit->p1 = (uint32)p;
  unit->p2 = p - unit->p1;
}

// calculate RLPF coefficients from fc1, q
void X75_WG4_Saxophone_setRLPFcoeffs(X75_WG4_Saxophone* unit, float fc, float q)
{
  float qres = sc_max(0.001, q);
  float pfreq = fc * unit->mRate->mRadiansPerSample;
  
  float D = tan(pfreq * qres * 0.5);
  float C = ((1.f-D)/(1.f+D));
  float cosf = cos(pfreq);
  
  unit->b1 = (1.f + C) * cosf;
  unit->b2 = -C;
  unit->a0 = (1.f + C - unit->b1) * .25;
#ifdef DEBUG
  printf("%g %g %g %g   %g %g %g  %g %g %g\n", fc, q, pfreq, qres, D, C, cosf, unit->b1, unit->b2, unit->a0);
#endif

  unit->fc1 = fc;
  unit->q1 = q;
  //unit->lpx = exp(-2 * PI * (unit->fc1/SAMPLERATE));
}

//////////////////////////////////////////////////////////////////
// calculation function when the buffer has been filled
void X75_WG4_Saxophone_next(X75_WG4_Saxophone *unit, int inNumSamples)
{
  float *out = OUT(0);
  float *in0 = IN(0);
  float *in1 = IN(1);

  // localize and update state values
  float p = ZIN0(2);
  if(p != unit->p) {
    X75_WG4_Saxophone_setP(unit, p);
    unit->p = p;
  }

  float rho1 = ZIN0(3);
  if(rho1 != unit->rho1)
    unit->rho1 = rho1;

  float rho2 = ZIN0(4);
  if(rho2 != unit->rho2)
    unit->rho2 = rho2;

  float fc1 = ZIN0(5);
  float q1 = ZIN0(6);
  if(fc1 != unit->fc1 || q1 != unit->q1)
    X75_WG4_Saxophone_setRLPFcoeffs(unit, fc1, q1);

  /* float fc2 = ZIN0(7);
  if(fc2 != unit->fc2)
  X75_WG4_Saxophone_setHpx(unit, fc2); */

  float squamp = ZIN0(8); //unit->squamp;
  if(squamp != unit->squamp)
    unit->squamp = squamp;

  float squoff = ZIN0(9); //unit->squoff;
  if(squoff != unit->squoff)
    unit->squoff = squoff;

  float sigma = ZIN0(10); //unit->sigma;
  if(sigma != unit->sigma)
    unit->sigma = sigma;

  // more local values
  float *data = unit->data;
  uint32 length = unit->p1;
  uint32 pos1 = unit->pos;
  uint32 pos2 = pos1 + 1;
  float p2 = unit->p2;
  float z; // unit output
  float zs, y0; // sub-signals
  float z1 = unit->z1; // unit mem
  float y1 = unit->y1; // filter mem
  float y2 = unit->y2; // filter mem
  float a0 = unit->a0; // filter coeffs
  float b1 = unit->b1;
  float b2 = unit->b2;
  float mix;

  int ws_numcoeffs = 4;
  double ws_coeffs[ws_numcoeffs];
  // cheby coeffs for a symmetrically limiting transfer function
  ws_coeffs[0] = 0;
  ws_coeffs[1] = -1.011; // -1.0109;
  ws_coeffs[2] = 0;
  ws_coeffs[3] = 0.122; //1.2221e-01;

  // los gehts
  for (int i=0; i < inNumSamples; ++i)
    {
      // x_n-p1-1
      if(pos2 >= length)
      //if(pos2 < 0)
	pos2 = 0;

      zs = in1[i] * z1; // excitation mix
 
      y0 = a0 * zs + b1 * y1 + b2 * y2; // resonant low pass filter
      //y0 = y0 + 2.f * y1 + y2; // 

      // mixer config with completed filter calculation
      mix = in0[i] + rho1 * (y0 + 2.f * y1 + y2) + rho2 * z1;

      y2 = y1; y1 = y0; // 2nd order filter memory

      // FIR interpolation, unit output
      z = (1-p2) * data[pos1] + p2 * data[pos2];
      if(sigma > 0.0)
	z = -z;

      // write input and feedback into delay line
      data[pos1] = X75_waveshape(mix * (mix * squamp + squoff), ws_coeffs, ws_numcoeffs); // get excitation into the loop

      // debug out
      //printf("next: z: %g, z_sum: %f, cursamp: %f\n", z, z_sum, data[position]);

      // see if the position went to the end of the buffer
      if (++pos1 >= length) {
	pos1 = 0; // go back to beginning
      }
      pos2 = pos1 + 1;
      
      // debug out
      //printf("next: z_sum: %f\n", z_sum);
      out[i] = z; // y1 = a0 * z - a1 * z1 + hpx * y1;
      z1 = z;
    }
  // store the position back to the struct
  unit->pos = pos1;
  unit->z1 = z1;
  unit->y1 = y1;
  unit->y2 = y2;
}

////////////////////////////////////////////////////////////////////
// the load function is called by the host when the plug-in is loaded
// void load(InterfaceTable *inTable)
PluginLoad(X75PhysMod)
{//
  ft = inTable;
  DefineDtorUnit(X75_Pluck1);
  DefineDtorUnit(X75_WG1);
  DefineDtorUnit(X75_WG2D);
  DefineDtorUnit(X75_WG2DPlus);
  DefineDtorUnit(X75_WG3_Trumpet);
  DefineDtorUnit(X75_WG4_Saxophone);
}
