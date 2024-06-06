#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fft.h"
#include "ifft.h"

#define PI 3.14159265358979323846

unsigned long * do_decimation(
  unsigned long *N,
  unsigned short *l
);
void initialize_xifft(
  double *xkMagnitude,
  double *xkPhase,
  double *xifftReal,
  double *xifftImag,
  unsigned long *N,
  unsigned long *M
);
void initialize_xp(
  double *xifftReal,
  double *xifftImag,
  double *xpReal,
  double *xpImag,
  unsigned long *indexes,
  unsigned long *N
);
void initialize_w(
  double *WNReal,
  double *WNImag,
  unsigned long *N
);
void calculate_inverse_wn(
  double *WNReal,
  double *WNImag,
  unsigned long *N,
  unsigned short *nPf,
  unsigned short fc
);
void do_butterfly_operation(
  double *xpReal,
  double *xpImag,
  double *xifftReal,
  double *xifftImag,
  double *WNReal,
  double *WNImag,
  unsigned short *BFc,
  unsigned short *nFc
);

/**
 * Calculate the Inverse Fast Fourier Transform of a given signal.
 *
 * @param *xkMagnitude Magnitude part of DFT.
 * @param *xkPhase     Phase part of DFT.
 * @param *N           Size of the signal.
 * @return             The recovered signal.
 */
struct IDFT calculate_ifft(double *xkMagnitude, double *xkPhase, unsigned long *M)
{
  struct IDFT ifft;
  double *xpReal, *xpImag, *xifftReal, *xifftImag, *WNReal, *WNImag;
  unsigned long *indexes;
  unsigned long i, N;
  unsigned short l = 0, sc = 1, nPf = 1, BFc = 0;

  l = (unsigned short) ceil(log((double) *M) / log(2.0));
  N = (unsigned long) pow(2, (double) l);

  xpReal = calloc(N, sizeof(double));
  xpImag = calloc(N, sizeof(double));
  xifftReal = calloc(N, sizeof(double));
  xifftImag = calloc(N, sizeof(double));
  WNReal = calloc(N, sizeof(double));
  WNImag = calloc(N, sizeof(double));

  if (
    NULL == xpReal
    || NULL == xpImag
    || NULL == xifftReal
    || NULL == xifftImag
    || NULL == WNReal
    || NULL == WNImag
  ) {
    printf("Unable to allocate a variable for IFFT\n");
		exit(0);
  }

  indexes = do_decimation(&N, &l); //Calculo de indeces en orden binario inverso
  initialize_xifft(xkMagnitude, xkPhase, xifftReal, xifftImag, &N, M); //Inicializacion y Zero-padding
  initialize_xp(xifftReal, xifftImag, xpReal, xpImag, indexes, &N);
  for (sc = 1; sc < l + 1; sc++) {
    nPf = 1 << (sc - 1); // Calcula # de factores de fase diferentes por etapa
    BFc = 0; //Contador de operaciones mariposa/Indice para resultados (BFc)
    initialize_w(WNReal, WNImag, &N);
    calculate_inverse_wn(WNReal, WNImag, &N, &nPf, (l - sc));

    for (i = 0; i< (unsigned long) (N / 2); i++) {
      do_butterfly_operation(xpReal, xpImag, xifftReal, xifftImag, WNReal, WNImag, &BFc, &nPf);
      // Si BFc es multiplo de nPF(2^sc-1), avanza nPf elementos en arreglo de
      // operandos para realizar mariposa, sino avanza 1 elemento.
      BFc += 0 == (BFc + 1) % nPf ? (1 << (sc - 1)) + 1 : 1;
    }

    for(i = 0; i < N; i++){ // Inicializa valores para la sigueinte etapa
      xpReal[i] = xifftReal[i];
    	xpImag[i] = xifftImag[i];
    }
  }

  for(i = 0; i < N; i++) {
    //Divide entre numero de muestras
    xifftReal[i] = xifftReal[i] / (double) N;//sqrt(pow(xifftReal[i], 2) + pow(xifftImag[i], 2)) / (double) N;
  }


  free(indexes);
  free(xpReal);
  free(xpImag);
  free(xifftImag);
  free(WNReal);
  free(WNImag);

  ifft.idft = xifftReal;
  ifft.size = N;

  return ifft;
}

/**
 * Initialize Xifft values.
 *
 * @param *xkMagnitude   Magnitude part of Xk.
 * @param *xkPhase       Phase part of Xk.
 * @param *xfftReal      Real part of Xfft.
 * @param *xfftImag      Imaginary part of Xfft.
 * @param *N             Size of Xfft.
 * @param *N             Size of Xn.
 */
void initialize_xifft(
  double *xkMagnitude,
  double *xkPhase,
  double *xifftReal,
  double *xifftImag,
  unsigned long *N,
  unsigned long *M
)
{
    unsigned long i;

    for (i = 0; i < *N; i++) {
      xifftReal[i] = i < *M ? xkMagnitude[i] * cos(xkPhase[i] * PI / 180) : 0.0;
      xifftImag[i] = i < *M ? xkMagnitude[i] * sin(xkPhase[i] * PI / 180) : 0.0;
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
void calculate_inverse_wn(
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
    WNrImag = sin(2 * PI * (r) / *N);
    j = *nPf + i; // Evita factores de fase WN^0 al inicio de cada etapa

    while(j < *N) { // Acomoda factores de fase en arreglo
      WNReal[j] = WNrReal;
      WNImag[j] = WNrImag;
      j += *nPf << 1;
    }

    r += 1 << fc; // Modifica exponente para calcular siguiente factor de fase
  }
}
