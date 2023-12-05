package io.github.mtrevisan.astro.core;

import io.github.mtrevisan.astro.helpers.MathHelper;


public class NutationCorrections{

	private static final double[][] NUTATION_COEFFS = {
		{297.85036, 445267.111480, -0.0019142, 1. / 189474},
		{357.52772, 35999.050340, -0.0001603, -1. / 300000},
		{134.96298, 477198.867398, 0.0086972, 1. / 56250},
		{93.27191, 483202.017538, -0.0036825, 1. / 327270},
		{125.04452, -1934.136261, 0.0020708, 1. / 450000}
	};
	private static final double[][] TERMS_Y = {
		{0., 0., 0., 0., 1.},
		{-2., 0., 0., 2., 2.},
		{0., 0., 0., 2., 2.},
		{0., 0., 0., 0., 2.},
		{0., 1., 0., 0., 0.},
		{0., 0., 1., 0., 0.},
		{-2., 1., 0., 2., 2.},
		{0., 0., 0., 2., 1.},
		{0., 0., 1., 2., 2.},
		{-2., -1., 0., 2., 2.},
		{-2., 0., 1., 0., 0.},
		{-2., 0., 0., 2., 1.},
		{0., 0., -1., 2., 2.},
		{2., 0., 0., 0., 0.},
		{0., 0., 1., 0., 1.},
		{2., 0., -1., 2., 2.},
		{0., 0., -1., 0., 1.},
		{0., 0., 1., 2., 1.},
		{-2., 0., 2., 0., 0.},
		{0., 0., -2., 2., 1.},
		{2., 0., 0., 2., 2.},
		{0., 0., 2., 2., 2.},
		{0., 0., 2., 0., 0.},
		{-2., 0., 1., 2., 2.},
		{0., 0., 0., 2., 0.},
		{-2., 0., 0., 2., 0.},
		{0., 0., -1., 2., 1.},
		{0., 2., 0., 0., 0.},
		{2., 0., -1., 0., 1.},
		{-2., 2., 0., 2., 2.},
		{0., 1., 0., 0., 1.},
		{-2., 0., 1., 0., 1.},
		{0., -1., 0., 0., 1.},
		{0., 0., 2., -2., 0.},
		{2., 0., -1., 2., 1.},
		{2., 0., 1., 2., 2.},
		{0., 1., 0., 2., 2.},
		{-2., 1., 1., 0., 0.},
		{0., -1., 0., 2., 2.},
		{2., 0., 0., 2., 1.},
		{2., 0., 1., 0., 0.},
		{-2., 0., 2., 2., 2.},
		{-2., 0., 1., 2., 1.},
		{2., 0., -2., 0., 1.},
		{2., 0., 0., 0., 1.},
		{0., -1., 1., 0., 0.},
		{-2., -1., 0., 2., 1.},
		{-2., 0., 0., 0., 1.},
		{0., 0., 2., 2., 1.},
		{-2., 0., 2., 0., 1.},
		{-2., 1., 0., 2., 1.},
		{0., 0., 1., -2., 0.},
		{-1., 0., 1., 0., 0.},
		{-2., 1., 0., 0., 0.},
		{1., 0., 0., 0., 0.},
		{0., 0., 1., 2., 0.},
		{0., 0., -2., 2., 2.},
		{-1., -1., 1., 0., 0.},
		{0., 1., 1., 0., 0.},
		{0., -1., 1., 2., 2.},
		{2., -1., -1., 2., 2.},
		{0., 0., 3., 2., 2.},
		{2., -1., 0., 2., 2.}
	};
	private static final double[][] TERMS_PE = {
		{-171996., -174.2, 92025., 8.9},
		{-13187., -1.6, 5736., -3.1},
		{-2274., -0.2, 977., -0.5},
		{2062., 0.2, -895., 0.5},
		{1426., -3.4, 54., -0.1},
		{712., 0.1, -7., 0.},
		{-517., 1.2, 224., -0.6},
		{-386., -0.4, 200., 0},
		{-301., 0., 129., -0.1},
		{217., -0.5, -95., 0.3},
		{-158., 0., 0., 0.},
		{129., 0.1, -70., 0.},
		{123., 0., -53., 0.},
		{63., 0., 0., 0.},
		{63., 0.1, -33., 0.},
		{-59., 0., 26., 0.},
		{-58., -0.1, 32., 0.},
		{-51., 0., 27., 0.},
		{48., 0., 0., 0.},
		{46., 0., -24., 0.},
		{-38., 0., 16., 0.},
		{-31., 0., 13., 0.},
		{29., 0., 0., 0.},
		{29., 0., -12., 0.},
		{26., 0., 0., 0.},
		{-22., 0., 0., 0.},
		{21., 0., -10., 0.},
		{17., -0.1, 0., 0.},
		{16., 0., -8., 0.},
		{-16., 0.1, 7., 0.},
		{-15., 0., 9., 0.},
		{-13., 0., 7., 0.},
		{-12., 0., 6., 0.},
		{11., 0., 0., 0.},
		{-10., 0., 5., 0.},
		{-8., 0., 3., 0.},
		{7., 0., -3., 0.},
		{-7., 0., 0., 0.},
		{-7., 0., 3., 0.},
		{-7., 0., 3., 0.},
		{6., 0., 0., 0.},
		{6., 0., -3., 0.},
		{6., 0., -3., 0.},
		{-6., 0., 3., 0.},
		{-6., 0., 3., 0.},
		{5., 0., 0., 0.},
		{-5., 0., 3., 0.},
		{-5., 0., 3., 0.},
		{-5., 0., 3., 0.},
		{4., 0., 0., 0.},
		{4., 0., 0., 0.},
		{4., 0., 0., 0.},
		{-4., 0., 0., 0.},
		{-4., 0., 0., 0.},
		{-4., 0., 0., 0.},
		{3., 0., 0., 0.},
		{-3., 0., 0., 0.},
		{-3., 0., 0., 0.},
		{-3., 0., 0., 0.},
		{-3., 0., 0., 0.},
		{-3., 0., 0., 0.},
		{-3., 0., 0., 0.},
		{-3., 0., 0., 0.}
	};


	/** Calculate <code>Δψ</code> [rad]. */
	private final double deltaPsi;
	/** Calculate <code>Δε</code> [rad]. */
	private final double deltaEpsilon;


	public static NutationCorrections calculate(final double ut){
		return new NutationCorrections(ut);
	}


	/**
	 * @param jce	Julian Century of Terrestrial Time from J2000.0.
	 */
	private NutationCorrections(final double jce){
		//calculate nutation corrections
		final double[] nutationTerms = nutationTerms(jce);
		final double[] deltaPsiCoeffs = deltaPsiCoeffs(jce, nutationTerms);
		final double[] deltaEpsilonCoeffs = deltaEpsilonCoeffs(jce, nutationTerms);
		deltaPsi = deltaPsiEpsilon(deltaPsiCoeffs);
		deltaEpsilon = deltaPsiEpsilon(deltaEpsilonCoeffs);
	}

	private static double[] nutationTerms(final double jce){
		final double[] x = new double[NUTATION_COEFFS.length];
		for(int i = 0; i < x.length; i ++)
			x[i] = MathHelper.polynomial(jce, NUTATION_COEFFS[i]);
		return x;
	}

	private static double[] deltaPsiCoeffs(final double jce, final double[] x){
		final double[] deltaPsiI = new double[TERMS_PE.length];
		for(int i = 0; i < TERMS_PE.length; i ++){
			final double a = TERMS_PE[i][0];
			final double b = TERMS_PE[i][1];
			deltaPsiI[i] = (a + b * jce) * StrictMath.sin(calculateXJYTermSum(i, x));
		}
		return deltaPsiI;
	}

	private static double[] deltaEpsilonCoeffs(final double jce, final double[] x){
		final double[] deltaEpsilonI = new double[TERMS_PE.length];
		for(int i = 0; i < TERMS_PE.length; i ++){
			final double c = TERMS_PE[i][2];
			final double d = TERMS_PE[i][3];
			deltaEpsilonI[i] = (c + d * jce) * StrictMath.cos(calculateXJYTermSum(i, x));
		}
		return deltaEpsilonI;
	}

	private static double calculateXJYTermSum(final int i, final double[] x){
		double sum = 0.;
		for(int j = 0; j < x.length; j ++)
			sum += x[j] * TERMS_Y[i][j];
		return StrictMath.toRadians(sum);
	}

	private static double deltaPsiEpsilon(final double[] deltaPsiOrEpsilonI){
		double sum = 0.;
		for(int i = 0; i < deltaPsiOrEpsilonI.length; i ++)
			sum += deltaPsiOrEpsilonI[i];
		return StrictMath.toRadians(sum / 36_000_000.);
	}


	public double getDeltaPsi(){
		return deltaPsi;
	}

	public double getDeltaEpsilon(){
		return deltaEpsilon;
	}

}
