// TODO: give credits to https://cb.wunderkis.de/wk-pub/www.iowahills.com/P51RootFinder.html
// contact author? or find alternative.

//---------------------------------------------------------------------------

#ifndef QuadRootsRevJH
#define QuadRootsRevJH
//---------------------------------------------------------------------------

//#include "CplxD.h"

#ifndef LDBL_EPSILON // TODO: ifndef makes sense?
#define LDBL_EPSILON 1.084202172485504434E-19L
#endif
// #define M_SQRT3 1.7320508075688772935274463L   // sqrt(3)
#define M_SQRT3_2 0.8660254037844386467637231L // sqrt(3)/2
// #define DBL_EPSILON  2.2204460492503131E-16    // 2^-52  typically defined in the compiler's float.h
#define ZERO_PLUS 8.88178419700125232E-16 // 2^-50 = 4*DBL_EPSILON
#define ZERO_MINUS -8.88178419700125232E-16
#define TINY_VALUE 1.0E-30 // This is what we use to test for zero. Usually to avoid divide by zero.
int QuadCubicRoots(long double *Coeff, int N, long double *RealRoot, long double *ImagRoot);
void QuadRoots(long double *P, long double *RealPart, long double *ImagPart);
void CubicRoots(long double *P, long double *RealPart, long double *ImagPart);
void BiQuadRoots(long double *P, long double *RealPart, long double *ImagPart);
void ReversePoly(long double *P, int N);
void InvertRoots(int N, long double *RealRoot, long double *ImagRoot);

/*
This is a scaled version of QuadRoots. We use b and c because modifying P isn't allowed for BiQuadRoots
void QuadRoots(long double *P, long double *RealRoot, long double *ImagRoot)
{
 long double D, Scalar, b, c;

 if(P[2] == 0.0)
  {
   RealRoot[0] = -P[1];
   RealRoot[1] = ImagRoot[0] = ImagRoot[1] = 0.0;
   return;
  }

 b = P[1];
 if(P[2] > 0.0)c = 1.0;
 else c = -1.0;

 Scalar = sqrtl(fabsl(P[2]));
 b /= Scalar;
 Scalar *= 0.5;

 D = b*b - 4.0*c;
 if(D >= 0.0)  // 2 real roots
  {
   RealRoot[0] = (-b - sqrtl(D)) * Scalar;
   RealRoot[1] = (-b + sqrtl(D)) * Scalar;
   ImagRoot[0] = ImagRoot[1] = 0.0;
  }
 else // 2 complex roots
  {
   RealRoot[0] = RealRoot[1] = -b * Scalar;
   ImagRoot[0] = sqrtl(-D) * Scalar;
   ImagRoot[1] = -ImagRoot[0];
  }
}

*/

//---------------------------------------------------------------------------
int QuadCubicRoots(long double *Coeff, int N, long double *RealRoot, long double *ImagRoot) {
	if (N < 1 || N > 4) {
		//ShowMessage("Invalid Poly Order in QuadCubicRoots()");
		return (0);
	}

	int j;
	long double P[5];

	// Must init to zero, in case N is reduced.
	for (j = 0; j < N; j++)
		RealRoot[j] = ImagRoot[j] = 0.0;
	for (j = 0; j < 5; j++)
		P[j] = 0.0;

	// The functions below modify the coeff array, so we pass P instead of Coeff.
	for (j = 0; j <= N; j++)
		P[j] = (long double)Coeff[j];

	// Remove trailing zeros. A tiny P[N] relative to P[N-1] is not allowed.
	while (N > 0 && fabsl(P[N]) <= TINY_VALUE * fabsl(P[N - 1])) {
		N--;
	}

	// Remove leading zeros.
	while (N > 0 && P[0] == 0.0) {
		for (j = 0; j < N; j++)
			P[j] = P[j + 1];
		N--;
	}

	// P[0] must = 1
	if (P[0] != 1.0) {
		for (j = 1; j <= N; j++)
			P[j] /= P[0];
		P[0] = 1.0;
	}

	// Calculate the roots.
	if (N == 4)
		BiQuadRoots(P, RealRoot, ImagRoot);
	else if (N == 3)
		CubicRoots(P, RealRoot, ImagRoot);
	else if (N == 2)
		QuadRoots(P, RealRoot, ImagRoot);
	else if (N == 1) {
		RealRoot[0] = -P[1] / P[0];
		ImagRoot[0] = 0.0;
	}

	return (N);
}

//---------------------------------------------------------------------------

// This function is the quadratic formula with P[0] = 1
// y = P[0]*x^2 + P[1]*x + P[2]
// Normally, we don't allow P[2] = 0, but this can happen in a call from BiQuadRoots.
// If P[2] = 0, the zero is returned in RealRoot[0].
void QuadRoots(long double *P, long double *RealRoot, long double *ImagRoot) {
	long double D;
	D = P[1] * P[1] - 4.0 * P[2];
	if (D >= 0.0) // 1 or 2 real roots
	{
		RealRoot[0] = (-P[1] - sqrtl(D)) * 0.5; // = -P[1] if P[2] = 0
		RealRoot[1] = (-P[1] + sqrtl(D)) * 0.5; // = 0 if P[2] = 0
		ImagRoot[0] = ImagRoot[1] = 0.0;
	} else // 2 complex roots
	{
		RealRoot[0] = RealRoot[1] = -P[1] * 0.5;
		ImagRoot[0] = sqrtl(-D) * 0.5;
		ImagRoot[1] = -ImagRoot[0];
	}
}

//---------------------------------------------------------------------------
// This finds the roots of y = P0x^3 + P1x^2 + P2x+ P3   P[0] = 1
void CubicRoots(long double *P, long double *RealRoot, long double *ImagRoot) {
	int j;
	long double s, t, b, c, d, Scalar;
	bool CubicPolyReversed = false;

	// Scale the polynomial so that P[N] = +/-1. This moves the roots toward unit circle.
	Scalar = powl(fabsl(P[3]), 1.0 / 3.0);
	for (j = 1; j <= 3; j++)
		P[j] /= powl(Scalar, (long double)j);

	if (fabsl(P[3]) < fabsl(P[2]) && P[2] > 0.0) {
		ReversePoly(P, 3);
		CubicPolyReversed = true;
	}

	s = P[1] / 3.0;
	b = (6.0 * P[1] * P[1] * P[1] - 27.0 * P[1] * P[2] + 81.0 * P[3]) / 162.0;
	t = (P[1] * P[1] - 3.0 * P[2]) / 9.0;
	c = t * t * t;
	d = 2.0 * P[1] * P[1] * P[1] - 9.0 * P[1] * P[2] + 27.0 * P[3];
	d = d * d / 2916.0 - c;

	// if(d > 0) 1 complex and 1 real root. We use LDBL_EPSILON to account for round off err.
	if (d > LDBL_EPSILON) {
		d = powl((sqrtl(d) + fabsl(b)), 1.0 / 3.0);
		if (d != 0.0) {
			if (b > 0)
				b = -d;
			else
				b = d;
			c = t / b;
		}
		d = M_SQRT3_2 * (b - c);
		b = b + c;
		c = -b / 2.0 - s;

		RealRoot[0] = (b - s);
		ImagRoot[0] = 0.0;
		RealRoot[1] = RealRoot[2] = c;
		ImagRoot[1] = d;
		ImagRoot[2] = -ImagRoot[1];
	}

	else // d < 0.0  3 real roots
	{
		if (b == 0.0)
			d = M_PI_2 / 3.0; //  b can be as small as 1.0E-25
		else
			d = atanl(sqrtl(fabsl(d)) / fabsl(b)) / 3.0;

		if (b < 0.0)
			b = 2.0 * sqrtl(fabsl(t));
		else
			b = -2.0 * sqrtl(fabsl(t));

		c = cosl(d) * b;
		t = -M_SQRT3_2 * sinl(d) * b - 0.5 * c;

		RealRoot[0] = (t - s);
		RealRoot[1] = -(t + c + s);
		RealRoot[2] = (c - s);
		ImagRoot[0] = 0.0;
		ImagRoot[1] = 0.0;
		ImagRoot[2] = 0.0;
	}

	// If we reversed the poly, the roots need to be inverted.
	if (CubicPolyReversed)
		InvertRoots(3, RealRoot, ImagRoot);

	// Apply the Scalar to the roots.
	for (j = 0; j < 3; j++)
		RealRoot[j] *= Scalar;
	for (j = 0; j < 3; j++)
		ImagRoot[j] *= Scalar;
}

//---------------------------------------------------------------------------

// This finds the roots of  y = P0x^4 + P1x^3 + P2x^2 + P3x + P4    P[0] = 1
// This function calls CubicRoots and QuadRoots
void BiQuadRoots(long double *P, long double *RealRoot, long double *ImagRoot) {
	int j;
	long double a, b, c, d, e, Q3Limit, Scalar, Q[4], MinRoot;
	bool QuadPolyReversed = false;

	// Scale the polynomial so that P[N] = +/- 1. This moves the roots toward unit circle.
	Scalar = powl(fabsl(P[4]), 0.25);
	for (j = 1; j <= 4; j++)
		P[j] /= powl(Scalar, (long double)j);

	// Having P[1] < P[3] helps with the Q[3] calculation and test.
	if (fabsl(P[1]) > fabsl(P[3])) {
		ReversePoly(P, 4);
		QuadPolyReversed = true;
	}

	a = P[2] - P[1] * P[1] * 0.375;
	b = P[3] + P[1] * P[1] * P[1] * 0.125 - P[1] * P[2] * 0.5;
	c = P[4] + 0.0625 * P[1] * P[1] * P[2] - 0.01171875 * P[1] * P[1] * P[1] * P[1] - 0.25 * P[1] * P[3];
	e = P[1] * 0.25;

	Q[0] = 1.0;
	Q[1] = P[2] * 0.5 - P[1] * P[1] * 0.1875;
	Q[2] = (P[2] * P[2] - P[1] * P[1] * P[2] + 0.1875 * P[1] * P[1] * P[1] * P[1] - 4.0 * P[4] + P[1] * P[3]) * 0.0625;
	Q[3] = -b * b * 0.015625;

	/* The value of Q[3] can cause problems when it should have calculated to zero (just above) but
	 is instead ~ -1.0E-17 because of round off errors. Consequently, we need to determine whether
	 a tiny Q[3] is due to roundoff, or if it is legitimately small. It can legitimately have values
	 of ~ -1E-28. When this happens, we assume Q[2] should also be small. Q[3] can also be tiny with
	 2 sets of equal real roots. Then P[1] and P[3], are approx equal. */

	Q3Limit = ZERO_MINUS;
	if (fabsl(fabsl(P[1]) - fabsl(P[3])) >= ZERO_PLUS &&
			Q[3] > ZERO_MINUS && fabsl(Q[2]) < 1.0E-5)
		Q3Limit = 0.0;

	if (Q[3] < Q3Limit && fabsl(Q[2]) < 1.0E20 * fabsl(Q[3])) {
		CubicRoots(Q, RealRoot, ImagRoot);

		// Find the smallest positive real root. One of the real roots is always positive.
		MinRoot = 1.0E100;
		for (j = 0; j < 3; j++) {
			if (ImagRoot[j] == 0.0 && RealRoot[j] > 0 && RealRoot[j] < MinRoot)
				MinRoot = RealRoot[j];
		}

		d = 4.0 * MinRoot;
		a += d;
		if (a * b < 0.0)
			Q[1] = -sqrtl(d);
		else
			Q[1] = sqrtl(d);
		b = 0.5 * (a + b / Q[1]);
	} else {
		if (Q[2] < 0.0) // 2 sets of equal imag roots
		{
			b = sqrtl(fabsl(c));
			d = b + b - a;
			if (d > 0.0)
				Q[1] = sqrtl(fabsl(d));
			else
				Q[1] = 0.0;
		} else {
			if (Q[1] > 0.0)
				b = 2.0 * sqrtl(fabsl(Q[2])) + Q[1];
			else
				b = -2.0 * sqrtl(fabsl(Q[2])) + Q[1];
			Q[1] = 0.0;
		}
	}

	// Calc the roots from two 2nd order polys and subtract e from the real part.
	if (fabsl(b) > 1.0E-8) {
		Q[2] = c / b;
		QuadRoots(Q, RealRoot, ImagRoot);

		Q[1] = -Q[1];
		Q[2] = b;
		QuadRoots(Q, RealRoot + 2, ImagRoot + 2);

		for (j = 0; j < 4; j++)
			RealRoot[j] -= e;
	} else // b==0 with 4 equal real roots
	{
		for (j = 0; j < 4; j++)
			RealRoot[j] = -e;
		for (j = 0; j < 4; j++)
			ImagRoot[j] = 0.0;
	}

	// If we reversed the poly, the roots need to be inverted.
	if (QuadPolyReversed)
		InvertRoots(4, RealRoot, ImagRoot);

	// Apply the Scalar to the roots.
	for (j = 0; j < 4; j++)
		RealRoot[j] *= Scalar;
	for (j = 0; j < 4; j++)
		ImagRoot[j] *= Scalar;
}

//---------------------------------------------------------------------------

// A reversed polynomial has its roots at the same angle, but reflected about the unit circle.
void ReversePoly(long double *P, int N) {
	int j;
	long double Temp;
	for (j = 0; j <= N / 2; j++) {
		Temp = P[j];
		P[j] = P[N - j];
		P[N - j] = Temp;
	}
	if (P[0] != 0.0) {
		for (j = N; j >= 0; j--)
			P[j] /= P[0];
	}
}

//---------------------------------------------------------------------------
// This is used in conjunction with ReversePoly
void InvertRoots(int N, long double *RealRoot, long double *ImagRoot) {
	int j;
	long double Mag;
	for (j = 0; j < N; j++) {
		// complex math for 1/x
		Mag = RealRoot[j] * RealRoot[j] + ImagRoot[j] * ImagRoot[j];
		if (Mag != 0.0) {
			RealRoot[j] /= Mag;
			ImagRoot[j] /= -Mag;
		}
	}
}
//---------------------------------------------------------------------------

Vector<real_t> solve_quartic(const real_t h0, const real_t h1, const real_t h2, const real_t h3, const real_t h4) {
	Vector<real_t> roots;
	long double coeffs[5] = { h4, h3, h2, h1, h0 }; // Note: reverse convention.
	int N = 4;
	long double real_parts[4] = { 0.0, 0.0, 0.0, 0.0 };
	long double imag_parts[4] = { 0.0, 0.0, 0.0, 0.0 };
	N = QuadCubicRoots(coeffs, N, real_parts, imag_parts);
	for (int i = 0; i < 4; i++) {
		if (Math::is_zero_approx((double)imag_parts[i])) {
			roots.push_back(real_parts[i]);
		}
	}
	return roots;
}

#endif