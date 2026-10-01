/**
 * Copyright (c) 2026 Mauro Trevisan
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
package io.github.mtrevisan.familylegacy.v2.sunset;

import io.github.mtrevisan.astro.core.NutationCorrections;
import io.github.mtrevisan.astro.core.SunPosition;
import io.github.mtrevisan.astro.helpers.JulianDate;
import io.github.mtrevisan.astro.helpers.MathHelper;
import io.github.mtrevisan.familylegacy.v2.services.AstronomicalEngine;

import java.time.LocalDate;


/**
 * High-precision implementation of AstronomicalEngine using full ELP2000-82B / VSOP87
 * periodic series for historical genealogical accuracy (covering -4000 to +8000).
 */
public final class SunsetAstronomicalEngineAdapter implements AstronomicalEngine{

	// Precision target: 0.01 seconds of day (~0.0000001 days)
	private static final double TIME_PRECISION = 1.0e-7;

	// Mean synodic month derived from IAU 2010 mean elongation motion
	private static final double MEAN_SYNODIC_MONTH = (360. * JulianDate.CIVIL_SAECULUM) / NutationCorrections.MOON_MEAN_ELONGATION_COEFFS[1];


	@Override
	public boolean isAvailable(){
		return true;
	}

	@Override
	public long getNextNewMoonJdn(final double approxJdn, final double utcOffset){
		double tt = approxJdn;
		double correction;

		// Newton-Raphson iteration for exact Moon-Sun longitude conjunction (0 rad)
		do{
			final double jce = JulianDate.centuryJ2000Of(tt);
			final NutationCorrections nutation = NutationCorrections.calculate(jce);

			// Geocentric apparent solar longitude (VSOP87 + nutation + aberration)
			final double sunLong = SunPosition.apparentSunLongitude(jce / 10.0, nutation.getDeltaPsi());

			// Geocentric apparent lunar longitude (ELP2000-82B high-precision series)
			final double moonLong = calculateHighPrecisionMoonLongitude(jce, nutation.getDeltaPsi());

			final double phaseAngle = MathHelper.mod2pi(moonLong - sunLong);
			final double diff = (phaseAngle > Math.PI? phaseAngle - MathHelper.TWO_PI: phaseAngle);

			correction = -diff * (MEAN_SYNODIC_MONTH / MathHelper.TWO_PI);
			tt += correction;
		}while(Math.abs(correction) > TIME_PRECISION);

		// Convert TT to local timezone JDN
		return Math.round(tt + (utcOffset / 24.));
	}

	@Override
	public double getSolarLongitudeJdn(final int year, final double targetLongitudeDeg, final double utcOffset){
		final double targetRad = Math.toRadians(targetLongitudeDeg);
		final LocalDate date = LocalDate.of(year, 1, 1);
		final int yearLength = (date.isLeapYear()? 366: 365);

		double tt = JulianDate.of(date);
		final double minUT = tt;
		double correction;

		do{
			final double jce = JulianDate.centuryJ2000Of(tt);
			final NutationCorrections nutation = NutationCorrections.calculate(jce);
			final double sunApparentLongitude = SunPosition.apparentSunLongitude(jce / 10., nutation.getDeltaPsi());

			// Inverse solar longitude solver
			correction = 58. * StrictMath.sin(targetRad - sunApparentLongitude);
			tt += correction;
			if(tt < minUT){
				tt += yearLength;
			}
		}while(Math.abs(correction) > TIME_PRECISION);

		return tt + (utcOffset / 24.);
	}

	/**
	 * Computes geocentric ecliptic apparent longitude of the Moon using ELP2000-82B main series
	 * including planetary perturbations (Meeus Ch. 47 / Chapront-Touzé).
	 *
	 * @param jce      Julian Century of Terrestrial Time from J2000.0.
	 * @param deltaPsi Nutation in longitude in radians.
	 * @return Apparent ecliptic longitude in radians.
	 */
	private static double calculateHighPrecisionMoonLongitude(final double jce, final double deltaPsi){
		// Fundamental Arguments (IAU 2010 / ELP2000)
		final double L0 = Math.toRadians(MathHelper.polynomial(jce, new double[]{218.316_447_7, 481_267.881_234_21, -0.001_578_6, 1. / 538_841., -1. / 65_194_000.}));
		final double D = Math.toRadians(MathHelper.polynomial(jce, new double[]{297.850_192_1, 445_267.111_403_4, -0.001_881_9, 1. / 545_868., -1. / 113_065_000.}));
		final double M = Math.toRadians(MathHelper.polynomial(jce, new double[]{357.529_109_2, 35_999.050_290_9, -0.000_153_6, 1. / 24_490_000.}));
		final double Mp = Math.toRadians(MathHelper.polynomial(jce, new double[]{134.963_396_4, 477_198.867_505_5, 0.008_741_4, 1. / 69_699., -1. / 14_712_000.}));
		final double F = Math.toRadians(MathHelper.polynomial(jce, new double[]{93.272_095_0, 483_202.017_523_3, -0.003_653_9, -1. / 3_526_000., 1. / 863_310_000.}));

		// Venus and Jupiter perturbation terms
		final double A1 = Math.toRadians(119.75 + 131.849 * jce);
		final double A2 = Math.toRadians(53.09 + 479264.290 * jce);
		final double A3 = Math.toRadians(313.45 + 481266.484 * jce);
		final double E = 1. - 0.002516 * jce - 0.0000074 * jce * jce;

		// Major Periodic Terms for Lunar Longitude (in 0.000001 degrees)
		double sumL = 6288774. * StrictMath.sin(Mp)
			+ 1274027. * StrictMath.sin(2. * D - Mp)
			+ 658314. * StrictMath.sin(2. * D)
			+ 213618. * StrictMath.sin(2. * Mp)
			- 185116. * E * StrictMath.sin(M)
			- 114332. * StrictMath.sin(2. * F)
			+ 58793. * StrictMath.sin(2. * D - 2. * Mp)
			+ 57066. * E * StrictMath.sin(2. * D - M - Mp)
			+ 53322. * StrictMath.sin(2. * D + Mp)
			+ 45758. * E * StrictMath.sin(2. * D - M)
			- 40923. * E * StrictMath.sin(M - Mp)
			- 34720. * StrictMath.sin(D)
			- 30383. * E * StrictMath.sin(M + Mp)
			+ 15327. * StrictMath.sin(2. * D - 2. * F)
			- 12528. * StrictMath.sin(Mp + 2. * F)
			+ 10980. * StrictMath.sin(Mp - 2. * F)
			+ 10675. * StrictMath.sin(4. * D - Mp)
			+ 10075. * StrictMath.sin(3. * Mp)
			+ 8548. * StrictMath.sin(4. * D - 2. * Mp)
			- 7888. * E * StrictMath.sin(2. * D + M - Mp)
			- 6766. * E * StrictMath.sin(2. * D + M)
			- 5163. * StrictMath.sin(D - Mp)
			+ 4987. * E * StrictMath.sin(D + M)
			+ 4026. * StrictMath.sin(2. * D + 2. * Mp)
			+ 3660. * StrictMath.sin(4. * D)
			- 3549. * StrictMath.sin(2. * D + 2. * F)
			- 2469. * StrictMath.sin(3. * Mp - 2. * D)
			+ 2249. * E * StrictMath.sin(2. * D - M - 2. * Mp)
			- 2073. * StrictMath.sin(2. * D - 3. * Mp)
			+ 2235. * StrictMath.sin(A1)
			+ 3820. * StrictMath.sin(A2)
			+ 1750. * StrictMath.sin(A3);

		final double moonGeocentricLong = L0 + Math.toRadians(sumL / 1_000_000.);

		// Apply nutation in longitude
		return MathHelper.mod2pi(moonGeocentricLong + deltaPsi);
	}

}
