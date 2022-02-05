/*
 * Copyright (c) 2020-2022. * mauro Trevisan
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */
package io.github.mtrevisan.sunset;

import io.github.mtrevisan.sunset.coordinates.EquatorialCoordinate;

import java.io.IOException;
import java.util.Collection;
import java.util.Map;


/*
https://github.com/mikereedell/sunrisesunsetlib-java

https://web.archive.org/web/20161202180207/http://williams.best.vwh.net/sunrise_sunset_algorithm.htm
https://ebvalaim.pl/en/2015/12/22/calculating-sunrise-and-sunset-times/
https://en.wikipedia.org/wiki/Position_of_the_Sun
https://edwilliams.org/sunrise_sunset_algorithm.htm

others
https://github.com/MenoData/Time4J/blob/master/base/src/main/java/net/time4j/calendar/astro/SolarTime.java
https://github.com/shred/commons-suncalc/blob/master/src/main/java/org/shredzone/commons/suncalc/SunTimes.java
https://github.com/buelowp/sunset/blob/master/src/sunset.cpp

https://frinklang.org/frinksamp/sun.frink
https://www.astrouw.edu.pl/~jskowron/pracownia/praca/sunspot_answerbook_expl/expl-5.html
http://co2.aos.wisc.edu/data/code/idl-lib/util/sunrise.pro
https://ebvalaim.pl/en/2015/12/22/calculating-sunrise-and-sunset-times/

https://www.iers.org/SharedDocs/Publikationen/EN/IERS/Publications/tn/TechnNote13/tn13.pdf?__blob=publicationFile&v=1
https://squarewidget.com/solar-coordinates/
*/
public final class SunPosition{

	private static Map<String, Collection<Double[]>> EARTH_HELIOCENTRIC_DATA;
	private static Map<String, Collection<Double[]>> EARTH_RADIUS_VECTOR_DATA;
	private static Map<String, Collection<Double[]>> NUTATION_DATA;
	static{
		try{
			EARTH_HELIOCENTRIC_DATA = ResourceReader.read("earthHeliocentric.dat");
			EARTH_RADIUS_VECTOR_DATA = ResourceReader.read("earthRadiusVector.dat");
			NUTATION_DATA = ResourceReader.read("nutation.dat");
		}
		catch(final IOException ignored){}
	}

	static final double EARTH_FLATTENING = 1. / 298.25642;
	//[m]
	static final double EARTH_EQUATORIAL_RADIUS = 6378140.;

	private static final double[] MOON_MEAN_ELONGATION_PARAMETERS = {297.8503631, 445267.1114800, - 0.0019142, 1. / 189474.};
	private static final double[] SUN_MEAN_ANOMALY_PARAMETERS = {357.5277233, 35999.0503400, - 0.0001603, - 1. / 300000.};
	private static final double[] MOON_MEAN_ANOMALY_PARAMETERS = {134.9629814, 477198.8673981, 0.0086972, 1. / 56250.};
	private static final double[] MOON_ARGUMENT_OF_LATITUDE = {93.2719103, 483202.0175381, - 0.00368250, 1. / 327270.};
	private static final double[] MOON_LONGITUDE_ASCENDING_NODE = {125.0445222, - 1934.1362608, 0.002070833, 1. / 450000.};
	private static final double[] MEAN_ECLIPTIC_OBLIQUITY_PARAMETERS = {21.448, - 4680.93, - 1.55, 1999.25, - 51.38, - 249.67, - 39.05, 7.12,
		27.87, 5.79, 2.45};


	private SunPosition(){}


	/**
	 * Calculated the Sun position.
	 *
	 * @param tt	Julian Century of Terrestrial Time from J2000.0.
	 * @return	The Sun position.
	 *
	 * @see <a href="https://squarewidget.com/solar-coordinates/">Solar coordinates</>
	 */
	public static EquatorialCoordinate sunPosition(final double tt){
		//calculate the geometric mean longitude L0 of the Sun referred to the mean equinox of the time T: L0
		final double geometricMeanLongitude = geometricMeanLongitude(tt);
		//calculate the nutation in longitude and obliquity
		final double[] nutation = nutationCorrection(tt);
		final double radiusVector = radiusVector(tt);
		//calculate the aberration correction
		final double aberration = aberrationCorrection(radiusVector);
		//calculate the apparent Sun longitude: Ltrue = L0 + C
		final double apparentGeometricLongitude = apparentGeometricLongitude(geometricMeanLongitude, nutation[0], aberration);
		//calculate the obliquity of the ecliptic (the inclination of the Earth’s equator with respect to the plane at which the Sun
		//and planets appear to move across the sky): ɛ0
		final double meanEclipticObliquity = meanEclipticObliquity(tt);
		//calculate the true obliquity of the ecliptic
		final double trueEclipticObliquity = trueEclipticObliquity(meanEclipticObliquity, nutation[1]);
		final double geometricMeanLatitude = geometricMeanLatitude(tt);
		return EquatorialCoordinate.createFromEcliptical(geometricMeanLatitude, apparentGeometricLongitude, trueEclipticObliquity);
	}

	/**
	 * Calculate the geometric mean longitude of the Sun, referred to the mean equinox of the date, L.
	 *
	 * @param tt	Julian Century of Terrestrial Time from J2000.0.
	 * @return	The geometric mean longitude of the Sun [°].
	 *
	 * @see <a href="https://squarewidget.com/solar-coordinates/">Solar coordinates</>
	 */
	static double geometricMeanLongitude(final double tt){
		final double jme = tt / 10.;
		final double[] parameters = new double[6];
		for(int i = 0; i < parameters.length; i ++){
			double parameter = 0.;
			final Collection<Double[]> elements = EARTH_HELIOCENTRIC_DATA.get("L" + i);
			for(final Double[] element : elements)
				parameter += element[0] * StrictMath.cos(element[1] + element[2] * jme);
			parameters[i] = parameter;
		}
		final double longitude = MathHelper.eval(jme, parameters) / 100_000_000.;
		return MathHelper.correctRangeDegree(StrictMath.toDegrees(longitude) + 180.);
	}

	/**
	 * Calculate corrections of nutation in longitude (∆ψ) and obliquity (∆ε).
	 *
	 * @param tt	Julian Century of Terrestrial Time from J2000.0.
	 * @return	An array where the first element is ∆ψ [°], and the second ∆ε [°].
	 *
	 * http://physics.uwyo.edu/~wiro/planet/nutate.c
	 * https://iopscience.iop.org/article/10.1086/375641/pdf
	 */
	static double[] nutationCorrection(final double tt){
		//mean elongation of the Moon from the Sun [°]
		final double d = StrictMath.toRadians(MathHelper.correctRangeDegree(
			MathHelper.eval(tt, MOON_MEAN_ELONGATION_PARAMETERS))
		);
		//mean anomaly of the Sun [°]
		final double m = StrictMath.toRadians(MathHelper.correctRangeDegree(
			MathHelper.eval(tt, SUN_MEAN_ANOMALY_PARAMETERS))
		);
		//mean anomaly of the Moon [°]
		final double mp = StrictMath.toRadians(MathHelper.correctRangeDegree(
			MathHelper.eval(tt, MOON_MEAN_ANOMALY_PARAMETERS))
		);
		//Moon's argument of Latitude [°]
		final double f = StrictMath.toRadians(MathHelper.correctRangeDegree(
			MathHelper.eval(tt, MOON_ARGUMENT_OF_LATITUDE))
		);
		//Longitude of the ascending node of the Moon's mean orbit on the ecliptic measured from the mean equinox of the date [rad]
		final double omega = StrictMath.toRadians(MathHelper.correctRangeDegree(
			MathHelper.eval(tt, MOON_LONGITUDE_ASCENDING_NODE))
		);

		final Collection<Double[]> elements = NUTATION_DATA.get("coeffs");
		final double[] x = {d, m, mp, f, omega};
		double deltaPsi = 0.;
		double deltaEpsilon = 0.;
		for(int i = 0; i < elements.size(); i ++)
			for(final Double[] element : elements){
				double parameter = 0.;
				for(int j = 0; j < x.length; j ++)
					parameter += x[j] * element[j];
				deltaPsi += (element[x.length] + element[x.length + 1] * tt) * StrictMath.sin(parameter);
				deltaEpsilon += (element[x.length + 2] + element[x.length + 3] * tt) * StrictMath.cos(parameter);
			}
		//FIXME /63 ?!?!?!
		//[°]
		deltaPsi = MathHelper.toDegrees(0, 0, deltaPsi / 10_000.) / 63.;
		//[°]
		deltaEpsilon = MathHelper.toDegrees(0, 0, deltaEpsilon / 10_000.) / 63.;
		return new double[]{deltaPsi, deltaEpsilon};
	}

	/**
	 * Calculate the distance between the center of the Sun and the center of the Earth, R.
	 * <p>U.S. Naval Observatory function.</p>
	 *
	 * @param tt	Julian Century of Terrestrial Time from J2000.0.
	 * @return	Distance between the center of the Sun and the center of the Earth [AU].
	 */
	static double radiusVector(final double tt){
		final double jme = tt / 10.;
		final double[] parameters = new double[5];
		for(int i = 0; i < parameters.length; i ++){
			double parameter = 0.;
			final Collection<Double[]> elements = EARTH_RADIUS_VECTOR_DATA.get("R" + i);
			for(final Double[] element : elements)
				parameter += element[0] * StrictMath.cos(element[1] + element[2] * jme);
			parameters[i] = parameter;
		}
		return MathHelper.eval(jme, parameters) / 100_000_000.;
	}

	/**
	 * Calculate the correction for aberration (∆τ).
	 *
	 * @param earthRadiusVector	Earth radius vector [AU].
	 * @return	The correction for aberration [°].
	 */
	static double aberrationCorrection(final double earthRadiusVector){
		return -20.4898 / (JulianDay.SECONDS_IN_HOUR * earthRadiusVector);
	}

	/**
	 * Calculate the apparent longitude of the Sun, Lapp = λ.
	 *
	 * @param geometricMeanLongitude	Sun's true geometric longitude [°].
	 * @param deltaPsi	Nutation in longitude [°].
	 * @param deltaAberration	Aberration [°].
	 * @return	Apparent longitude of the Sun [°].
	 */
	static double apparentGeometricLongitude(final double geometricMeanLongitude, final double deltaPsi, final double deltaAberration){
		return geometricMeanLongitude + deltaPsi + deltaAberration;
	}

	/**
	 * Calculate the mean obliquity of the ecliptic, ɛ0.
	 *
	 * @param tt	Julian Century of Terrestrial Time from J2000.0.
	 * @return	Apparent longitude of the Sun [°].
	 */
	static double meanEclipticObliquity(final double tt){
		final double u = tt / 100.;
		final double seconds = MathHelper.eval(u, MEAN_ECLIPTIC_OBLIQUITY_PARAMETERS);
		return MathHelper.toDegrees(23, 26, seconds);
	}

	/**
	 * True obliquity of the ecliptic corrected for nutation, ɛ.
	 *
	 * @param meanEclipticObliquity	Obliquity of the ecliptic, corrected for parallax [°].
	 * @param deltaEpsilon	Nutation in obliquity [°].
	 * @return	Apparent longitude of the Sun [°].
	 */
	static double trueEclipticObliquity(final double meanEclipticObliquity, final double deltaEpsilon){
		return meanEclipticObliquity + deltaEpsilon;
	}

	/**
	 * Calculate the geometric mean latitude of the Sun, referred to the mean equinox of the date, β.
	 *
	 * @param tt	Julian Century of Terrestrial Time from J2000.0.
	 * @return	The geometric mean latitude of the Sun [°].
	 *
	 * @see <a href="https://squarewidget.com/solar-coordinates/">Solar coordinates</>
	 */
	static double geometricMeanLatitude(final double tt){
		final double jme = tt / 10.;
		final double[] parameters = new double[2];
		for(int i = 0; i < parameters.length; i ++){
			double parameter = 0.;
			final Collection<Double[]> elements = EARTH_HELIOCENTRIC_DATA.get("B" + i);
			for(final Double[] element : elements)
				parameter += element[0] * StrictMath.cos(element[1] + element[2] * jme);
			parameters[i] = parameter;
		}
		final double latitude = MathHelper.eval(jme, parameters) / 100_000_000.;
		return MathHelper.correctRangeDegree(-StrictMath.toDegrees(latitude));
	}

}
