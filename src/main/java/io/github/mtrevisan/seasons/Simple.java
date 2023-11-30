/**
 * Copyright (c) 2021 Mauro Trevisan
 * <p>
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 * <p>
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * <p>
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */
package io.github.mtrevisan.seasons;

import io.github.mtrevisan.sunset.JulianDay;
import io.github.mtrevisan.sunset.MathHelper;
import io.github.mtrevisan.sunset.TimeHelper;
import io.github.mtrevisan.sunset.core.SunPosition;

import java.time.LocalDateTime;


public class Simple{

	private static final double[] EARTH_WINTER_SOLSTICE_MINUS_1000_1000 = {1_721_414.399_87, 365_242.88257, -0.007_69, -0.009_33,
		-0.000_06};
	private static final double[] EARTH_WINTER_SOLSTICE_1000_3000 = {2_451_900.059_52, 365_242.740_49, -0.062_23, -0.008_23, 0.000_32};

	private static final double[] PERIOD_TERMS_A = new double[]{485., 203., 199., 182., 156., 136., 77., 74., 70., 58., 52., 50., 45., 44.,
		29., 18., 17., 16., 14., 12., 12., 12., 9., 8.};
	private static final double[] PERIOD_TERMS_B = new double[]{324.96, 337.23, 342.08, 27.85, 73.14, 171.52, 222.54, 296.72, 243.58,
		119.81, 297.17, 21.02, 247.54, 325.15, 60.93, 155.12, 288.79, 198.04, 199.76, 95.39, 287.11, 320.81, 227.73, 15.45};
	private static final double[] PERIOD_TERMS_C = new double[]{1934.136, 32964.467, 20.186, 445267.112, 45036.886, 22518.443, 65928.934,
		3034.906, 9037.513, 33718.147, 150.678, 2281.226, 29929.562, 31555.956, 4443.417, 67555.328, 4562.452, 62894.029, 31436.921,
		14577.848, 31931.756, 34777.259, 1222.114, 16859.074};

	//range of years in lookup table:
	private static final int DELTA_T_TABLE_FIRST_YEAR = 1620;
	private static final int DELTA_T_TABLE_LAST_YEAR = 2002;
	//corrections [s]
	private static final double[] DELTA_T_TABLE = {
		/*1620*/
		121, 112, 103, 95, 88, 82, 77, 72, 68, 63, 60, 56, 53, 51, 48, 46, 44, 42, 40, 38,
		/*1660*/
		35, 33, 31, 29, 26, 24, 22, 20, 18, 16, 14, 12, 11, 10, 9, 8, 7, 7, 7, 7,
		/*1700*/
		7, 7, 8, 8, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11,
		/*1740*/
		11, 11, 12, 12, 12, 12, 13, 13, 13, 14, 14, 14, 14, 15, 15, 15, 15, 15, 16, 16,
		/*1780*/
		16, 16, 16, 16, 16, 16, 15, 15, 14, 13,
		/*1800*/
		13.1, 12.5, 12.2, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 11.9, 11.6, 11.0, 10.2, 9.2, 8.2,
		/*1830*/
		7.1, 6.2, 5.6, 5.4, 5.3, 5.4, 5.6, 5.9, 6.2, 6.5, 6.8, 7.1, 7.3, 7.5, 7.6,
		/*1860*/
		7.7, 7.3, 6.2, 5.2, 2.7, 1.4, -1.2, -2.8, -3.8, -4.8, -5.5, -5.3, -5.6, -5.7, -5.9,
		/*1890*/
		-6.0, -6.3, -6.5, -6.2, -4.7, -2.8, -0.1, 2.6, 5.3, 7.7, 10.4, 13.3, 16.0, 18.2, 20.2,
		/*1920*/
		21.1, 22.4, 23.5, 23.8, 24.3, 24.0, 23.9, 23.9, 23.7, 24.0, 24.3, 25.3, 26.2, 27.3, 28.2,
		/*1950*/
		29.1, 30.0, 30.7, 31.4, 32.2, 33.1, 34.0, 35.0, 36.5, 38.3, 40.2, 42.2, 44.5, 46.5, 48.5,
		/*1980*/
		//values for ΔT for 2000 thru 2002 from NASA
		50.5, 52.5, 53.8, 54.9, 55.8, 56.9, 58.3, 60.0, 61.6, 63.0, 63.8, 64.3};


	public static void main(String[] args){
		int year = 2022;

		simple(year);
	}

	private static void simple(int year){
		//calculate the approximate date
		final double jed = winterSolstice(year);

		//just an innocuous approximation
		final double tdt = jed;

		//correct TDT to UTC
		final double deltaT = deltaT(year);
		//NOTE: the error is less than 41 seconds for the years 1951-2050
		double utc = TimeHelper.terrestrialTimeToUniversalTime(tdt, deltaT);
System.out.println(utc);

		final LocalDateTime dateTime = JulianDay.dateTimeOf(utc);
System.out.println(dateTime);

		while(true){
			//calculate the ecliptical mean longitude of the Sun, referred to the mean equinox of the date
			final double jme = JulianDay.millenniumJ2000Of(utc);
			final double meanEclipticalLongitude = SunPosition.meanEclipticalLongitude(jme);
			//calculate the distance between the center of the Sun and the center of the Earth, R
			final double radiusVector = SunPosition.radiusVector(jme);

			//calculate corrections of nutation in longitude, ∆ψ
			final double sunMeanAnomaly = SunPosition.geocentricMeanAnomaly(jme * 10.);
			final double[] nutation = SunPosition.nutationCorrection(sunMeanAnomaly, jme);

			//calculate the correction for aberration, ∆τ
			final double aberration = SunPosition.aberrationCorrection(radiusVector);

			final double correctionFK5 = Math.toRadians(-0.090_33 / 3600.);

			//calculate the apparent longitude of the Sun, Lapp = λ
			final double apparentGeocentricLongitude = SunPosition.apparentGeocentricLongitude(meanEclipticalLongitude, nutation[0],
				aberration) + correctionFK5;

			final double deltaJDE = 58. * StrictMath.sin(3. * Math.PI / 2. - apparentGeocentricLongitude);
			utc += deltaJDE;

			if(Math.abs(deltaJDE) * JulianDay.SECONDS_IN_DAY < 42.)
				break;
		}
System.out.println(utc);

		final LocalDateTime dateTime2 = JulianDay.dateTimeOf(utc);
System.out.println(dateTime2);
System.out.println("2022-12-21T21:49:22");

		//TODO https://www.agopax.it/Libri_astronomia/pdf/Astronomical%20Algorithms.pdf, chap 25
	}

	/**
	 * Instant of the "mean" solstice for the years [-1000, 3000].
	 *
	 * @return	JD of Winter Solstice in Julian Ephemeris Days (hence, in Dynamical Time).
	 */
	private static double winterSolstice(final int year){
		if(year < -1000 || year > 3000)
			throw new IllegalArgumentException("The year must be in [-1000, 3000]");

		final double deltaYear2000 = year - 2000.;
		final double jde0;
		if(year < 1000)
			jde0 = MathHelper.eval(year /1000., EARTH_WINTER_SOLSTICE_MINUS_1000_1000);
		else
			jde0 = MathHelper.eval(deltaYear2000 /1000., EARTH_WINTER_SOLSTICE_1000_3000);

		final double jce = JulianDay.centuryJ2000Of(jde0);
		final double w = Math.toRadians(35_999.373 * jce - 2.47);
		final double deltaLambda = 1. + 0.033_4 * StrictMath.cos(w) + 0.000_7 * StrictMath.cos(2. * w);
		return jde0 + periodicTerms(jce) / deltaLambda;
	}

	private static double periodicTerms(final double t){
		double sum = 0.;
		for(int i = 0; i < PERIOD_TERMS_A.length; i ++)
			sum += PERIOD_TERMS_A[i] * StrictMath.cos(Math.toRadians(PERIOD_TERMS_B[i] + PERIOD_TERMS_C[i] * t));
		return 0.000_01 * sum;
	}

	private static double deltaT(final int year){
		//deltaT = TDT - UTC [s]
		double deltaT;
		//centuries from the epoch J2000.0
		final double t = (year - 2000.) / 100.;

		//find correction in table
		if(year >= DELTA_T_TABLE_FIRST_YEAR && year <= DELTA_T_TABLE_LAST_YEAR){
			if((year % 2) != 0)
				//odd year - interpolate
				deltaT = (DELTA_T_TABLE[(year - DELTA_T_TABLE_FIRST_YEAR - 1) / 2] + DELTA_T_TABLE[(year - DELTA_T_TABLE_FIRST_YEAR + 1) / 2])
					/ 2.;
			else
				//even year - direct table lookup
				deltaT = DELTA_T_TABLE[(year - DELTA_T_TABLE_FIRST_YEAR) / 2];
		}
		else if(year < 948)
			deltaT = MathHelper.eval(t, new double[]{2177., 497., 44.1});
		else{
			deltaT = MathHelper.eval(t, new double[]{102., 102., 25.3});
			//special correction to avoid discontinuity in 2000
			if(year >= 2000 && year <= 2100)
				deltaT += 0.37 * (year - 2100.);
		}
		return deltaT;
	}

}
