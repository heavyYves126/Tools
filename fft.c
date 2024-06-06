#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fft.h"

#define PI 3.14159265358979323846

unsigned long * do_decimation(
  unsigned long *N,
  unsigned short *l
);
void initialize_w(
  double *WNReal,
  double *WNImag,
  unsigned long *N
);
void initialize_xfft(
  double *xn,
  double *xfftReal,
  double *xfftImag,
  unsigned long *N,
  unsigned long *M
);
void initialize_xp(
  double *xfftReal,
  double *xfftImag,
  double *xpReal,
  double *xpImag,
  unsigned long *indexes,
  unsigned long *N
);
void calculate_wn(
  double *WNReal,
  double *WNImag,
  unsigned long *N,
  unsigned short *nPf,
  unsigned short fc
);
void do_butterfly_operation(
  double *xpReal,
  double *xpImag,
  double *xfftReal,
  double *xfftImag,
  double *WNReal,
  double *WNImag,
  unsigned short *BFc,
  unsigned short *nFc
);
double calculate_angle(
  double *a,
  double *b
);

/**
 * Calculate the Fast Fourier Transform of a given signal.
 *
 * @param * xn Discrete signal.
 * @param * N  Size of the signal.
 * @return    The Fourier Transform.
 */
struct DFT calculate_fft(double *xn, unsigned long *M)
{
  struct DFT fft;
  double *xpReal, *xpImag, *xfftReal, *xfftImag, *WNReal, *WNImag;
  unsigned long *indexes;
  unsigned long i, N;
  unsigned short l = 0, sc = 1, nPf = 1, BFc = 0;

  l = (unsigned short) ceil((log((double) *M) / log(2.0)));
  N = (unsigned long) pow(2, (double) l);
  xpReal = calloc(N, sizeof(double));
  xpImag = calloc(N, sizeof(double));
  xfftReal = calloc(N, sizeof(double));
  xfftImag = calloc(N, sizeof(double));
  WNReal = calloc(N, sizeof(double));
  WNImag = calloc(N, sizeof(double));

  if (
    NULL == xpReal
    || NULL == xpImag
    || NULL == xfftReal
    || NULL == xfftImag
    || NULL == WNReal
    || NULL == WNImag
  ) {
    printf("Unable to allocate a variable for FFT\n");
		exit(0);
  }

  indexes = do_decimation(&N, &l); //Calculo de indeces en orden binario inverso
  initialize_xfft(xn, xfftReal, xfftImag, &N, M); //Inicializacion y Zero-padding
  initialize_xp(xfftReal, xfftImag, xpReal, xpImag, indexes, &N); //Decimacion en el tiempo y preparacion de X(w)

  for (sc = 1; sc < l + 1; sc++) {
    nPf = 1 << (sc - 1); // Calcula # de factores de fase diferentes por etapa
    BFc = 0; //Contador de operaciones mariposa/Indice para resultados (BFc)
    initialize_w(WNReal, WNImag, &N);
    calculate_wn(WNReal, WNImag, &N, &nPf, (l - sc));

    for (i = 0; i< (unsigned long) (N / 2); i++) {
      do_butterfly_operation(xpReal, xpImag, xfftReal, xfftImag, WNReal, WNImag, &BFc, &nPf);
      // Si BFc es multiplo de nPF(2^sc-1), avanza nPf elementos en arreglo de
      // operandos para realizar mariposa, sino avanza 1 elemento.
      BFc += 0 == (BFc + 1) % nPf ? nPf + 1 : 1;
    }

    for(i = 0; i < N; i++){ // Inicializa valores para la siguiente etapa
      xpReal[i] = xfftReal[i];
    	xpImag[i] = xfftImag[i];
    }
  }

  for(i = 0; i < N; i++) {
    //Calcula modulo de los valores de la FFT
    xpReal[i] = sqrt(pow(xfftReal[i], 2) + pow(xfftImag[i], 2));
    //Calcula angulo de los valores de la FFT
    xpImag[i] = calculate_angle(&xfftReal[i], &xfftImag[i]);
  }

  free(indexes);
  free(xfftReal);
  free(xfftImag);
  free(WNReal);
  free(WNImag);

  fft.magnitude = xpReal;
  fft.phase = xpImag;
  fft.size = N;

  return fft;
}

/**
 * Initialize Xfft values.
 *
 * @param *xn       The signal Xn.
 * @param *xfftReal Real part of Xfft.
 * @param *xfftImag Imaginary part of Xfft.
 * @param *N        Size of Xfft.
 * @param *N        Size of Xn.
 */
void initialize_xfft(
  double *xn,
  double *xfftReal,
  double *xfftImag,
  unsigned long *N,
  unsigned long *M
)
{
    unsigned long i;

    for (i = 0; i < *N; i++) {
      xfftReal[i] = i >= *M ? 0 : xn[i];
      xfftImag[i] = 0;
    }
}

/**
 * Initialize Xp by rearranging the xn input.
 *
 * @param *xfftReal Real part of Xfft.
 * @param *xfftImag Imaginary part of Xfft.
 * @param *xpReal   Real part of Xp.
 * @param *xpImag   Imaginary part of Xp.
 * @param *indexes  Decimation indexes.
 * @param *N        Size of the signal.
 */
void initialize_xp(
  double *xfftReal,
  double *xfftImag,
  double *xpReal,
  double *xpImag,
  unsigned long *indexes,
  unsigned long *N
)
{
    unsigned long i = 0, j = 0;

    for (i = 0; i < *N; i++) {
      j = indexes[i];
      xpReal[i] = xfftReal[j];
      xpImag[i] = xfftImag[j];
    }
}

/**
 * Do the decimation of the signal.
 *
 * @param *N  The size of the signal.
 * @param *l  The l value.
 * @return   Array of decimation indexes.
 */
unsigned long * do_decimation(unsigned long *N, unsigned short *l)
{
  unsigned long *indexes;
  unsigned long i = 0;
  unsigned short j = 0, temp = 0, reverse_num = 0;

  indexes = calloc(*N, sizeof(unsigned long));

  if ( NULL == indexes) {
    printf("Unable to allocate indexes\n");
    exit(0);
  }


  for (i = 0; i < *N; i++) {
    for (j = 0; j < *l; j++) {
      temp = (i & (1 << j));
      if(0 != temp) {
        reverse_num |= (1 << ((*l - 1) - j));
      }
    }

    indexes[i] = reverse_num;
    temp = 0;
    reverse_num = 0;
  }

  return indexes;
}

/**
 * Initialize Wn on each phase.
 *
 * @param *WNReal Real part of Wn.
 * @param *WNImag Imaginary part of Wn.
 * @param *N      Size of the signal.
 */
void initialize_w(double *WNReal, double *WNImag, unsigned long *N)
{
  unsigned long i = 0;

  // Inicializa factores de fase
  for(i = 0; i < *N; i++) {
    WNReal[i] = 1;
    WNImag[i] = 0;
  }
}

/**
 * Calculate the values of Wn on each phase.
 *
 * @param *WNReal Real part of Wn.
 * @param *WNImag Imaginary part of Wn.
 * @param *N      Size of the signal.
 * @param *nPf    Phase factor of the current phase.
 * @param fc
 */
void calculate_wn(
  double *WNReal,
  double *WNImag,
  unsigned long *N,
  unsigned short *nPf,
  unsigned short fc
)
{
  unsigned long i = 0, j = 0;
  unsigned short r = 0;
  double WNrReal = 0, WNrImag = 0;

  for (i = 0; i < *nPf; i++) {
    WNrReal = cos(2 * PI * (r) / *N); // Calcula factores de fase por etapa
    WNrImag = sin(-2 * PI * (r) / *N);
    j = *nPf + i; // Evita factores de fase WN^0 al inicio de cada etapa

    while(j < *N) { // Acomoda factores de fase en arreglo
      WNReal[j] = WNrReal;
      WNImag[j] = WNrImag;
      j += *nPf << 1;
    }

    r += 1 << fc; // Modifica exponente para calcular siguiente factor de fase
  }
}

/**
 * Reallocate vectors to a power of 2.
 *
 * @param *xpReal   Real part of Xp..
 * @param *xpImag   Imaginary part of Xp.
 * @param *xfftReal Real part of Xfft.
 * @param *xfftImag Imaginary part of Xfft.
 * @param *WNReal   Real part of Wn.
 * @param *WNImag   Imaginary part of Xn.
 * @param *nFc      Current index.
 * @param *BFc      Next index.
 */
void do_butterfly_operation(
  double *xpReal,
  double *xpImag,
  double *xfftReal,
  double *xfftImag,
  double *WNReal,
  double *WNImag,
  unsigned short *BFc,
  unsigned short *nPf
)
{
  // Obtiene el elemento superior resultado de la meriposa
  xfftReal[*BFc] = xpReal[*BFc]
                  + (xpReal[*BFc + *nPf] * WNReal[*BFc + *nPf])
                  - (xpImag[*BFc + *nPf] * WNImag[*BFc + *nPf]);
  xfftImag[*BFc] = xpImag[*BFc]
                  + (xpReal[*BFc + *nPf] * WNImag[*BFc + *nPf])
                  + (xpImag[*BFc + *nPf] * WNReal[*BFc + *nPf]);
  // Obtiene elemento inferior resultado de la mariposa
  xfftReal[*BFc + *nPf] = xpReal[*BFc]
                  - (xpReal[*BFc + *nPf] * WNReal[*BFc + *nPf])
                  + (xpImag[*BFc + *nPf] * WNImag[*BFc + *nPf]);
  xfftImag[*BFc + *nPf] = xpImag[*BFc]
                  - (xpReal[*BFc + *nPf] * WNImag[*BFc + *nPf])
                  - (xpImag[*BFc + *nPf] * WNReal[*BFc + *nPf]);

}

/**
 * Calculate the angle of a complex number.
 *
 * @param  a Real part.
 * @param  b Imaginary part.
 * @return   Angle
 */
double calculate_angle(double *a, double *b)
{
  double fi = 0.0;

  if (0.0 == fabs(*a)) {
    return *b >= 0.0 ? 90.0 : -90.0;
  } else if (0.0 == fabs(*b)) {
    return *a >= 0.0 ? 0.0 : 180.0;
  }

  fi = fabs(atan(*b / *a) * 180 / PI);

  if (0.0 > *a) {
    return 180 + (*b >= 0.0 ? -fi : fi);
  }

  return *b >= 0.0 ? fi : 360.0 - fi;
}
