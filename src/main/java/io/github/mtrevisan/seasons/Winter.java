package io.github.mtrevisan.seasons;

import io.github.mtrevisan.sunset.JulianDay;
import io.github.mtrevisan.sunset.TimeHelper;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.Locale;


public class Winter{

	public static void main(final String[] args){
		final int year = 2022;
		final double longitude = 12.19415;

		//https://www.agopax.it/Libri_astronomia/pdf/Astronomical%20Algorithms.pdf

		//calculate the approximate date
		final double winterSolsticeTDB = winterSolstice(year);

		//TDB ~ TDT (the difference is at most 0.0017 s)
		//https://gssc.esa.int/navipedia/index.php/Transformations_between_Time_Systems#TAI_-_TDT.2C_TCG.2C_TT
		final double tt = JulianDay.centuryJ2000Of(winterSolsticeTDB);
		final double g = Math.toRadians(357.528 + 35999.05 * tt);
		final double winterSolsticeTDT = winterSolsticeTDB - 0.001658 * StrictMath.sin(g + 0.0167 * StrictMath.sin(g)) / 86400.;

		//TODO https://www.agopax.it/Libri_astronomia/pdf/Astronomical%20Algorithms.pdf
		//calculate the Sun's apparent longitude lambda (chap 25 + chap 22 + formula 25.9 + formula 25.10), including the corrections for reduction to the FK5 system, for aberration and for nutation
		//the correction is then 58. * StrictMath.sin(Math.toRadians(270. - lambda[°])) [day]
		//repeat until the correction is small

		//transform TDT into UT
		final double deltaT = TimeHelper.deltaT(year);
		final DecimalFormat decimalFormatter = (DecimalFormat)NumberFormat.getNumberInstance(Locale.US);
		decimalFormatter.applyPattern("0.#");
		System.out.println("∆T = " + decimalFormatter.format(deltaT) + " s");
		final double winterSolsticeUT = TimeHelper.terrestrialTimeToUniversalTime(winterSolsticeTDT, deltaT);
		LocalDateTime winterSolsticeDateTime = JulianDay.dateTimeOf(winterSolsticeUT);
		System.out.println(DateTimeFormatter.ISO_LOCAL_DATE_TIME.format(winterSolsticeDateTime) + " UTC");

		//convert UT into local time
		final double winterSolsticeLocalTime = winterSolsticeUT + longitude / (24. * 15.);
		winterSolsticeDateTime = JulianDay.dateTimeOf(winterSolsticeLocalTime);
		decimalFormatter.applyPattern("0.######");
		System.out.println(DateTimeFormatter.ISO_LOCAL_DATE_TIME.format(winterSolsticeDateTime) + " @ " + decimalFormatter.format(longitude) + "°"
			+ (longitude > 0.? " E": (longitude < 0.? " W": "")));
	}

	/**
	 * @return	JD of Winter Solstice in Barycentric Dynamical Time.
	 */
	private static double winterSolstice(final int year){
		//instant of the "mean" solstice for the years [1000, 3000]
		final double y = (year - 2000.) /1000.;
		final double jde0 = 2451900.05952 + (365242.74049 + (0.06223 - (0.00823 - 0.00032 * y) * y) * y) * y;

		final double tt = JulianDay.centuryJ2000Of(jde0);
		final double w = Math.toRadians(35999.373 * tt - 2.47);
		final double deltaLambda = 1. + 0.0334 * StrictMath.cos(w) + 0.0007 * StrictMath.cos(2. * w);
		return jde0 + (0.00001 * periodicTerms(tt)) / deltaLambda
			- (66. + (year - 2000.)) / 86400.;
	}

	private static double periodicTerms(final double t){
		double x = 485. * StrictMath.cos(Math.toRadians(324.96 + 1934.136 * t));
		x += 203. * StrictMath.cos(Math.toRadians(337.23 + 32964.467 * t));
		x += 199. * StrictMath.cos(Math.toRadians(342.08 + 20.186 * t));
		x += 182. * StrictMath.cos(Math.toRadians(27.85 + 445267.112 * t));
		x += 156. * StrictMath.cos(Math.toRadians(73.14 + 45036.886 * t));
		x += 136. * StrictMath.cos(Math.toRadians(171.52 + 22518.443 * t));
		x += 77. * StrictMath.cos(Math.toRadians(222.54 + 65928.934 * t));
		x += 74. * StrictMath.cos(Math.toRadians(296.72 + 3034.906 * t));
		x += 70. * StrictMath.cos(Math.toRadians(243.58 + 9037.513 * t));
		x += 58. * StrictMath.cos(Math.toRadians(119.81 + 33718.147 * t));
		x += 52. * StrictMath.cos(Math.toRadians(297.17 + 150.678 * t));
		x += 50. * StrictMath.cos(Math.toRadians(21.02 + 2281.226 * t));

		x += 45. * StrictMath.cos(Math.toRadians(247.54 + 29929.562 * t));
		x += 44. * StrictMath.cos(Math.toRadians(325.15 + 31555.956 * t));
		x += 29. * StrictMath.cos(Math.toRadians(60.93 + 4443.417 * t));
		x += 18. * StrictMath.cos(Math.toRadians(155.12 + 67555.328 * t));

		x += 17. * StrictMath.cos(Math.toRadians(288.79 + 4562.452 * t));
		x += 16. * StrictMath.cos(Math.toRadians(198.04 + 62894.029 * t));
		x += 14. * StrictMath.cos(Math.toRadians(199.76 + 31436.921 * t));
		x += 12. * StrictMath.cos(Math.toRadians(95.39 + 14577.848 * t));
		x += 12. * StrictMath.cos(Math.toRadians(287.11 + 31931.756 * t));
		x += 12. * StrictMath.cos(Math.toRadians(320.81 + 34777.259 * t));
		x += 9. * StrictMath.cos(Math.toRadians(227.73 + 1222.114 * t));
		x += 8. * StrictMath.cos(Math.toRadians(15.45 + 16859.074 * t));
		return x;
	}

}
