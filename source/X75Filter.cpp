/*
 * X75FilterUGens
 * Oswald Berthold 20071107
 * GPL
 *
 * inspired by MatKlang 07/08
 */

#include "SC_PlugIn.h"

#define PI 3.141592653589793
#define PI2 2*M_PI

static InterfaceTable *ft;

// considerations:
// an element for use in a filterbank
// bank is constructed in sclang
// 
struct X75BP2 : public Unit
{
  int inp, oup; // input/output pointers
  float x_n[3]; // input at t_0, t_1, t_2
  float y_n[3]; // output at t_0, t_1, t_2
  float c_n[3]; // numerator coefficients
  float d_n[3]; // denominator coeffs
  float Q, freq, R, omega, amp; // parameter
};

// Amplitude Follower
struct X75AmpF : public Unit
{
  // y_n := \max{|x_n|, \beta y_{n-1}}
  float ym1; // y_{n-1}
  float alpha; // beta :)
};

extern "C"
{
  void load(InterfaceTable *inTable);
  void X75BP2_next(X75BP2 *unit, int inNumSamples);
  void X75BP2_Ctor(X75BP2* unit);
  void X75BP2_setCoeffs(X75BP2* unit);
  void X75AmpF_next(X75AmpF *unit, int inNumSamples);
  void X75AmpF_Ctor(X75AmpF *unit);
}

//////////////////////////////////////////////////////////////////////
//
void X75BP2_Ctor(X75BP2* unit)
{	
  //  printf("X75BP2_Reset\n");
  SETCALC(X75BP2_next);
  unit->inp = 0;
  unit->oup = 0;
  // init storage
  LOOP(3, unit->x_n[xxi] = 0.f);
  LOOP(3, unit->y_n[xxi] = 0.f);

  unit->freq = ZIN0(1);
  unit->Q = ZIN0(2);
  unit->amp = 1.0;

  X75BP2_setCoeffs(unit);
  X75BP2_next(unit, 1);
}

void X75BP2_setCoeffs(X75BP2* unit) {
  unit->omega = unit->freq * unit->mRate->mRadiansPerSample; // frequenz als kreisfrequenz
  unit->R = exp(-0.5*unit->omega/unit->Q); // bandbreite approx. als -2 ln R, moore S.134
  unit->amp = 1-unit->R;
  //  printf("freq: %f, kfreq: %f, Q: %f, R: %f, amp: %f\n", unit->freq, unit->omega, unit->Q, unit->R, unit->amp);
  unit->c_n[0] = 1.f;
  unit->c_n[1] = 0.f; // 1-unit->R;
  unit->c_n[2] = -(unit->R); // 3-unit->R; // 
  unit->d_n[0] = 1;
  unit->d_n[1] = 2 * unit->R * cos(unit->omega);
  unit->d_n[2] = -(pow(unit->R, 2));
//   LOOP(3,
//        printf("c_%d: %f d_%d: %f\n", xxi, unit->c_n[xxi], xxi, unit->d_n[xxi]);
//        );
}

void X75BP2_next(X75BP2* unit, int inNumSamples)
{
  //printf("X75BP2_next_a\n");

  float *out = ZOUT(0);
  float *in = ZIN(0);

  float freq = ZIN0(1);
  float Q = ZIN0(2);

  //  update coeffs
  if(freq != unit->freq || Q != unit->Q) {
    unit->freq = freq;
    unit->Q = Q;
    X75BP2_setCoeffs(unit); }

  // restore current pointers
  int inp = unit->inp;
  int oup = unit->oup;

  LOOP(inNumSamples,
       // printf("delp: %d, orig: %f, -4-samps: delp4: %d, val: %f\n", delp, unit->m_x[delp], delp4, unit->m_x[delp4]);
       unit->x_n[inp] = ZXP(in);
       // band 1
       unit->y_n[oup] = (unit->c_n[0] * unit->x_n[inp]) + (unit->c_n[1] * unit->x_n[(inp+2)%3]) + (unit->c_n[2] * unit->x_n[(inp+1)%3])
       + (unit->d_n[1] * unit->y_n[(oup+2)%3]) + (unit->d_n[2] * unit->y_n[(oup+1)%3]);
       // unit->y_n[oup] = unit->y_n[oup] * unit->amp;       
       ZXP(out) = unit->y_n[oup] * unit->amp;
       inp = (inp+1)%3;
       oup = (inp+1)%3;
       );
  // preserve current pointer
  unit->inp = inp;
  unit->oup = oup;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

void X75AmpF_Ctor(X75AmpF* unit)
{	
  //  printf("X75AmpF_Reset\n");
  SETCALC(X75AmpF_next);
  //unit->m_x = ZIN0(0);
  unit->alpha = ZIN0(1); // 0.15;
  unit->ym1 = 0.f;
  X75AmpF_next(unit, 1);
}


void X75AmpF_next(X75AmpF* unit, int inNumSamples)
{
  //printf("X75AmpF_next_a\n");

  float *out = ZOUT(0);
  float *in = ZIN(0);

  float alpha = ZIN0(1);
  float x;
  float ym1 = unit->ym1;

  if(unit->alpha != alpha) unit->alpha = alpha;

  //printf("inNumSamples: %d\ninNumSamples / 4.floor: %d\n", inNumSamples, inNumSamples >> 2);
  //printf("inNumSamples + 3: %d\n", inNumSamples & 3);

  LOOP(inNumSamples,
       x = ZXP(in);
       //printf("delp: %d, orig: %f, -4-samps: delp4: %d, val: %f\n", delp, unit->m_x[delp], delp4, unit->m_x[delp4]);
       ym1 = (x >= ym1 ? fabs(x) : alpha * ym1);
       ZXP(out) = ym1;
       );
  unit->ym1 = ym1;
}

//////////////////////////////////////////////////////////////////////
//void load(InterfaceTable *inTable)
PluginLoad(X75PhysMod)
{
  ft = inTable;

  DefineSimpleUnit(X75BP2);
  DefineSimpleUnit(X75AmpF);
}
