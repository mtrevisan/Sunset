package io.github.mtrevisan.astro.core;

import io.github.mtrevisan.astro.helpers.MathHelper;
import io.github.mtrevisan.astro.helpers.ResourceReader;

import java.io.IOException;
import java.util.List;
import java.util.Map;


/**
 * @see <a href="https://www.research.unipd.it/retrieve/e14fb26f-e069-3de1-e053-1705fe0ac030/PhD_Thesis_Zoccarato.pdf">Determinazione Orbitale Precisa di satelliti in orbita LEO per la Radio Occultazione attraverso sistemi GNSS</a>
 */
public class NutationCorrections{

	private static Map<String, List<double[]>> NUTATION_DATA;
	static{
		try{
			NUTATION_DATA = ResourceReader.read("nutation.dat");
		}
		catch(final IOException ignored){}
	}

	/**
	 * IAU 2010 theory, <code>D</code> [deg]
	 *
	 * <a href="https://iers-conventions.obspm.fr/content/tn36.pdf">IERS Conventions 2010</a>
	 */
	private static final double[] MOON_MEAN_ELONGATION_COEFFS = {297.850_195_47, 1_602_961_601.209_0/3600., -6.370_6/3600., 0.006_593/3600.,
		-0.000_031_69/3600.};
	/**
	 * IAU 2010 theory, <code>l'</code> [deg]
	 *
	 * <a href="https://iers-conventions.obspm.fr/content/tn36.pdf">IERS Conventions 2010</a>
	 */
	private static final double[] SUN_GEOCENTRIC_MEAN_ANOMALY_COEFFS = {357.529_109_18, 129_596_581.048_1/3600., -0.553_2/3600.,
		0.000_136/3600., -0.000_011_49/3600.};
	/**
	 * IAU 2010 theory, <code>l</code> [deg]
	 *
	 * <a href="https://iers-conventions.obspm.fr/content/tn36.pdf">IERS Conventions 2010</a>
	 */
	private static final double[] MOON_MEAN_ANOMALY_COEFFS = {134.963_402_51, 1_717_915_923.217_8/3600., 31.879_2/3600., 0.051_635/3600.,
		-0.000_244_70/3600.};
	/**
	 * IAU 2010 theory, <code>F</code> [deg]
	 *
	 * <a href="https://iers-conventions.obspm.fr/content/tn36.pdf">IERS Conventions 2010</a>
	 */
	private static final double[] MOON_ARGUMENT_OF_LATITUDE_COEFFS = {93.272_090_62, 1_739_527_262.847_8/3600., -12.751_2/3600.,
		-0.001_037/3600., 0.000_004_17/3600.};
	/**
	 * IAU 2010 theory, <code>Ω</code> [deg]
	 *
	 * <a href="https://iers-conventions.obspm.fr/content/tn36.pdf">IERS Conventions 2010</a>
	 */
	private static final double[] MOON_LONGITUDE_ASCENDING_NODE_COEFFS = {125.044_555_01, -6_962_890.543_1/3600., 7.472_2/3600.,
		0.007_702/3600., -0.000_059_39/3600.};
	private static final double[][] NUTATION_COEFFS = {
		MOON_MEAN_ELONGATION_COEFFS,
		SUN_GEOCENTRIC_MEAN_ANOMALY_COEFFS,
		MOON_MEAN_ANOMALY_COEFFS,
		MOON_ARGUMENT_OF_LATITUDE_COEFFS,
		MOON_LONGITUDE_ASCENDING_NODE_COEFFS
	};


	/** Calculate <code>Δψ</code> [rad]. */
	private final double deltaPsi;
	/** Calculate <code>Δε</code> [rad]. */
	private final double deltaEpsilon;


	public static NutationCorrections calculate(final double ut){
		return new NutationCorrections(ut);
	}


	/**
	 * Calculate nutation corrections following IAU 1980.
	 *
	 * @param jce	Julian Century of Terrestrial Time from J2000.0.
	 */
	private NutationCorrections(final double jce){
		//calculate nutation corrections
		final double[] nutationTerms = nutationTerms(jce);
		final List<double[]> elements = NUTATION_DATA.get("coeffs");
		final double[] deltaPsiCoeffs = deltaPsiCoeffs(jce, nutationTerms, elements);
		final double[] deltaEpsilonCoeffs = deltaEpsilonCoeffs(jce, nutationTerms, elements);
		deltaPsi = deltaPsiEpsilon(deltaPsiCoeffs);
		deltaEpsilon = deltaPsiEpsilon(deltaEpsilonCoeffs);
	}

	private static double[] nutationTerms(final double jce){
		final double[] x = new double[NUTATION_COEFFS.length];
		for(int i = 0; i < NUTATION_COEFFS.length; i ++)
			x[i] = MathHelper.polynomial(jce, NUTATION_COEFFS[i]);
		return x;
	}

	private static double[] deltaPsiCoeffs(final double jce, final double[] nutationTerms, final List<double[]> elements){
		final double[] deltaPsiI = new double[elements.size()];
		for(int i = 0; i < deltaPsiI.length; i ++){
			final double[] params = elements.get(i);
			final double a = params[5];
			final double b = params[6];
			deltaPsiI[i] = (a + b * jce) * StrictMath.sin(calculateXYTermSum(i, nutationTerms, params));
		}
		return deltaPsiI;
	}

	private static double[] deltaEpsilonCoeffs(final double jce, final double[] nutationTerms, final List<double[]> elements){
		final double[] deltaEpsilonI = new double[elements.size()];
		for(int i = 0; i < deltaEpsilonI.length; i ++){
			final double[] params = elements.get(i);
			final double c = params[7];
			final double d = params[8];
			deltaEpsilonI[i] = (c + d * jce) * StrictMath.cos(calculateXYTermSum(i, nutationTerms, params));
		}
		return deltaEpsilonI;
	}

	private static double calculateXYTermSum(final int i, final double[] x, final double[] params){
		double result = 0.;
		for(int j = 0; j < x.length; j ++)
			result += x[j] * params[j];
		return StrictMath.toRadians(result);
	}

	private static double deltaPsiEpsilon(final double[] deltaPsiOrEpsilonI){
		double result = 0.;
		for(int i = 0; i < deltaPsiOrEpsilonI.length; i ++)
			result += deltaPsiOrEpsilonI[i];
		return StrictMath.toRadians(result / 36_000_000.);
	}

	/**
	 * @param jce	Julian Century of Terrestrial Time from J2000.0.
	 * @return	Longitude of the ascending node of the Moon [rad].
	 */
	public static double moonLongitudeAscendingNode(final double jce){
		return StrictMath.toRadians(MathHelper.polynomial(jce, MOON_LONGITUDE_ASCENDING_NODE_COEFFS));
	}


	public double getDeltaPsi(){
		return deltaPsi;
	}

	public double getDeltaEpsilon(){
		return deltaEpsilon;
	}

}
