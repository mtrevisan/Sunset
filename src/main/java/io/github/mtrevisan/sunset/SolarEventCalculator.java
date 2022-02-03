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
import io.github.mtrevisan.sunset.coordinates.GNSSLocation;

import java.time.LocalDate;
import java.time.LocalDateTime;
import java.time.LocalTime;
import java.time.temporal.ChronoUnit;


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
*/
public class SolarEventCalculator{

	// Calculates the approximate set time of a body which has the specified right ascension and declination.
	// The resultant value will be close to the specified date.
	// Return values are undef if the object is circumpolar for that date.
	//	private double approxRiseSet(date, lat, long, ra, decl, h0 = -0.5667 degrees){
	//		transit = approxTransit[date, long, ra, decl];
	//		H0 = calcHourAngle[lat, decl, h0];
	//
	//		if(H0 == undef)
	//			return [undef, undef]
	//
	//		Htime = H0 / (360 degrees/day);
	//		//   println["transit is $transit"]
	//		//   println["H0 is $H0"]
	//		//   println["Htime is " + (Htime -> "hours")]
	//		set = transit + Htime;
	//		return set;
	//	}

	//https://frinklang.org/frinksamp/sun.frink
	//https://www.astrouw.edu.pl/~jskowron/pracownia/praca/sunspot_answerbook_expl/expl-5.html
	//http://co2.aos.wisc.edu/data/code/idl-lib/util/sunrise.pro
	//https://ebvalaim.pl/en/2015/12/22/calculating-sunrise-and-sunset-times/

	//---

	private final GNSSLocation location;


	/**
	 * Constructs a new instance using the given parameters.
	 *
	 * @param location	Location of the place.
	 */
	public static SolarEventCalculator create(final GNSSLocation location){
		return new SolarEventCalculator(location);
	}


	private SolarEventCalculator(final GNSSLocation location){
		this.location = location;
	}


//---
	//https://squarewidget.com/solar-coordinates/

	/**
	 * Calculated the Sun position.
	 *
	 * @param jd	The Julian Day [day].
	 * @return	The Sun position.
	 */
	public static EquatorialCoordinate sunPosition(final double jd){
		final double t = JulianDay.centuryJ2000Of(jd);
		final double geometricMeanLongitude = geometricMeanLongitude(t);
		final double meanAnomaly = geometricMeanAnomaly(t);
		final double equationOfCenter = equationOfCenter(meanAnomaly, t);
		final double trueGeometricLongitude = trueGeometricLongitude(geometricMeanLongitude, equationOfCenter);
		final double apparentLongitude = apparentGeometricLongitude(trueGeometricLongitude, t);
		final double meanEclipticObliquity = meanEclipticObliquity(t);
		final double apparentEclipticObliquity = apparentEclipticObliquity(meanEclipticObliquity, t);
		final double apparentRightAscension = apparentRightAscension(apparentEclipticObliquity, apparentLongitude);
		final double apparentDeclination = apparentDeclination(apparentEclipticObliquity, apparentLongitude);
		return EquatorialCoordinate.create(apparentRightAscension, apparentDeclination);
	}

	//TODO
	//https://github.com/MenoData/Time4J/blob/master/base/src/main/java/net/time4j/calendar/astro/SunPosition.java
	/**
	 * Computes the sunset time.
	 *
	 * @param solarZenith	Enumeration corresponding to the type of sunset to compute.
	 * @param date	Date to compute the sunset for.
	 * @return	The sunset time or {@code null} for no sunset.
	 * @throws SolarEventException	Whenever the Sun never rises or sets.
	 */
	public final LocalDateTime sunset(final LocalDate date, final Zenith solarZenith) throws SolarEventException{
		final double jd = JulianDay.of(date);
		final double t = JulianDay.centuryJ2000Of(jd);
		final EquatorialCoordinate coord = sunPosition(jd);
		//[hrs]
		final double ra = degToHrs(coord.getRightAscension());
		final double decl = coord.getLongitude();

/*
1/1/2022
jd = 2459580,5
n = 8036
longitude lw = 12
mean solar time J* = 8035,9666666666666666666666666667
solar mean anomaly M = 357,78009673733333333333333333333
equation of the center C = -0,07575266651750332420396923728956
ecliptic longitude lambda = 280,64154407081583000912936409604
solar transit Jtransit = 2459580,9714704432394846054500724
sinDecl = -0,39094722575833922205307251034751
declination of the Sun decl = -23,013451399537735607652540761379
latitute phi = 45
cosw0 = 0,40249461588178552731510832196043
hour angle w0 = 66,265778126410893799009823280512
sunset Jset = 2459581.1555420491461815326695441 = 15:43:59
*/
		final double geometricMeanLongitude = geometricMeanLongitude(t);
		final double meanAnomaly = geometricMeanAnomaly(t);
		final double eccentricity = earthOrbitEccentricity(t);
		final double equationOfCenter = equationOfCenter(meanAnomaly, t);
		final double trueGeometricLongitude = trueGeometricLongitude(geometricMeanLongitude, equationOfCenter);
		final double trueAnomaly = trueAnomaly(meanAnomaly, equationOfCenter);
		final double radiusVector = radiusVector(meanAnomaly);
		final double radiusVector2 = radiusVector(eccentricity, trueAnomaly);
		final double apparentLongitude = apparentGeometricLongitude(trueGeometricLongitude, t);
		final double meanEclipticObliquity = meanEclipticObliquity(t);
		final double apparentEclipticObliquity = apparentEclipticObliquity(meanEclipticObliquity, t);
		final double longitudeOfEarthPerihelion = longitudeOfEarthPerihelion(t);

		final double greenwichMeanSiderealTime = greenwichMeanSiderealTime(t);
		final double greenwichApparentSiderealTime = greenwichApparentSiderealTime(greenwichMeanSiderealTime, apparentEclipticObliquity, t);
		final double apparentLocalSiderealTime = apparentLocalSiderealTime(greenwichApparentSiderealTime);
		final double localHourAngle = localHourAngle(apparentLocalSiderealTime, coord.getRightAscension());

		LocalDateTime sunsetLocalDateTime = LocalDateTime.of(date, LocalTime.MIDNIGHT);
		LocalDateTime localDateTime;
		do{
			localDateTime = sunsetLocalDateTime;

			//[hrs]
			final LocalTime localTime = computeSolarEventTime(localDateTime, solarZenith, false);
			sunsetLocalDateTime = LocalDateTime.of(date, localTime);
		}while(!localDateTime.equals(sunsetLocalDateTime));
		return sunsetLocalDateTime;
	}

	/**
	 * Calculate the geometric mean longitude of the Sun, referred to the mean equinox of the date, L0.
	 *
	 * @param t	Julian Century in J2000.0 epoch.
	 * @return	The geometric mean longitude of the Sun [°].
	 */
	private static double geometricMeanLongitude(final double t){
		return correctRangeDegree(eval(t, new double[]{toDegrees(280, 27, 59.26), 36000.76983, 0.0003032}));
	}

	/**
	 * Calculate the geometric mean anomaly of the Sun, M.
	 *
	 * @param t	Julian Century in J2000.0 epoch.
	 * @return	The geometric mean anomaly of the Sun [°].
	 */
	private static double geometricMeanAnomaly(final double t){
		return correctRangeDegree(eval(t, new double[]{357.52911, 35999.05029, -0.0001537}));
	}

	/**
	 * Calculate the eccentricity of Earth's orbit, e.
	 *
	 * @param t	Julian Century in J2000.0 epoch.
	 * @return	The eccentricity of Earth's orbit.
	 */
	private static double earthOrbitEccentricity(final double t){
		return eval(t, new double[]{0.016708634, -0.000042037, -0.0000001267});
	}

	/**
	 * Calculate the Sun's equation of center, C.
	 *
	 * @param geometricMeanAnomaly	The mean anomaly of the Sun [°].
	 * @param t	Julian Century in J2000.0 epoch.
	 * @return	The Sun's equation of center [°].
	 */
	private static double equationOfCenter(double geometricMeanAnomaly, final double t){
		geometricMeanAnomaly = degToRad(geometricMeanAnomaly);
		return eval(t, new double[]{1.914602, -0.004817, -0.000014}) * StrictMath.sin(geometricMeanAnomaly)
			+ eval(t, new double[]{0.019993, -0.000101}) * StrictMath.sin(2. * geometricMeanAnomaly)
			+ 0.000289 * StrictMath.sin(3. * geometricMeanAnomaly);
	}

	/**
	 * Calculate the Sun's true geometric longitude, Ltrue = L0 + C.
	 *
	 * @param geometricMeanLongitude	The geometric mean longitude of the Sun, referred to the mean equinox of the date [°].
	 * @param equationOfCenter	The Sun's equation of center [°].
	 * @return	Sun's true geometric longitude [°].
	 */
	private static double trueGeometricLongitude(final double geometricMeanLongitude, final double equationOfCenter){
		return correctRangeDegree(geometricMeanLongitude + equationOfCenter);
	}

	/**
	 * Calculate the Sun's true anomaly longitude, ν = M + C.
	 *
	 * @param meanAnomaly	The mean anomaly of the Sun [°].
	 * @param equationOfCenter	The Sun's equation of center [°].
	 * @return	Sun's true geometric longitude [°].
	 */
	private static double trueAnomaly(final double meanAnomaly, final double equationOfCenter){
		return correctRangeDegree(meanAnomaly + equationOfCenter);
	}

	/**
	 * Calculate the distance between the center of the Sun and the center of the Earth, R.
	 * <p>U.S. Naval Observatory function.</p>
	 *
	 * @param meanAnomaly	The mean anomaly of the Sun [°].
	 * @return	Distance between the center of the Sun and the center of the Earth [AU].
	 */
	private static double radiusVector(double meanAnomaly){
		meanAnomaly = degToRad(meanAnomaly);
		return 1.00014
			- 0.01671 * StrictMath.cos(meanAnomaly)
			- 0.00014 * StrictMath.cos(meanAnomaly * 2.);
	}

	/**
	 * Calculate the distance between the center of the Sun and the center of the Earth, R.
	 *
	 * @param eccentricity	The eccentricity of Earth's orbit.
	 * @param trueAnomaly	The true anomaly of the Sun [°].
	 * @return	Distance between the center of the Sun and the center of the Earth [AU].
	 */
	private static double radiusVector(final double eccentricity, final double trueAnomaly){
		return 1.000001018 * (1. - eccentricity * eccentricity) / (1. + eccentricity * StrictMath.cos(degToRad(trueAnomaly)));
	}

	/**
	 * Calculate the apparent longitude of the Sun, Lapp = Ω.
	 *
	 * @param trueGeometricLongitude	Sun's true geometric longitude [°].
	 * @param t	Julian Century in J2000.0 epoch.
	 * @return	Apparent longitude of the Sun [°].
	 */
	private static double apparentGeometricLongitude(final double trueGeometricLongitude, final double t){
		final double correction = nutationAndAberrationCorrection(t);
		return trueGeometricLongitude - correction;
	}

	/**
	 * Calculate the mean obliquity of the ecliptic, ɛ0.
	 *
	 * @param t	Julian Century in J2000.0 epoch.
	 * @return	Apparent longitude of the Sun [°].
	 */
	private static double meanEclipticObliquity(final double t){
		final double u = t / 100.;
		final double seconds = eval(u, new double[]{21.448, -4680.93, -1.55, 1999.25, -51.38, -249.67, -39.05, 7.12, 27.87, 5.79, 2.45});
		return toDegrees(23, 26, seconds);
	}

	/**
	 * Calculate the mean obliquity of the ecliptic, corrected for parallax, ɛ'.
	 *
	 * @param meanEclipticObliquity	Mean obliquity of the ecliptic [°].
	 * @param t	Julian Century in J2000.0 epoch.
	 * @return	Apparent longitude of the Sun [°].
	 */
	private static double apparentEclipticObliquity(final double meanEclipticObliquity, final double t){
		final double correction = nutationAndAberrationCorrection(t);
		return meanEclipticObliquity + correction;
	}

	/**
	 * True obliquity of the ecliptic corrected for nutation, ɛ.
	 *
	 * @param meanEclipticObliquity	Obliquity of the ecliptic, corrected for parallax [°].
	 * @param t	Julian Century in J2000.0 epoch.
	 * @return	Apparent longitude of the Sun [°].
	 */
	private static double trueEclipticObliquity(final double meanEclipticObliquity, final double t){
		final double[] deltaPsiEpsilon = highAccuracyNutation(t);
		return meanEclipticObliquity + toDegrees(0, 0, deltaPsiEpsilon[1]);
	}

	/**
	 * Calculate Nutation in longitude (delta psi) and obliquity (delta epsilon).
	 *
	 * @param t	Julian Century in J2000.0 epoch.
	 * @return	An array where the first element is delta psi, and the second delta epsilon ["].
	 */
	private static double[] highAccuracyNutation(final double t){
		//mean elongation of the Moon from the Sun [°]
		final double d = eval(t, new double[]{297.85036, 445267.111480, -0.0019142, 1. / 189474.});
		//mean anomaly of the Sun [°]
		final double m = eval(t, new double[]{357.52772, 35999.050340, -0.0001603, - 1. / 300000.});
		//mean anomaly of the Moon [°]
		final double mp = eval(t, new double[]{134.96298, 477198.867398, 0.0086972, 1. / 56250.});
		//Moon's argument of Latitude [°]
		final double f = eval(t, new double[]{93.27191, 483202.017538, -0.0036825, 1. / 327270.});
		//Longitude of the ascending node of the Moon's mean orbit on the ecliptic measured from the mean equinox of the date [rad]
		final double omega = degToRad(eval(t, new double[]{125.04452, -1934.136261, 0.0020708, 1. / 450000.}));

		//these lines generated by iau1980.frink and pasted in here ["]
		@SuppressWarnings("OverlyComplexArithmeticExpression")
		final double deltaPsi = 0.0001 * ((-171996 - 174.2 * t) * StrictMath.sin(omega)
			+ (-13187. - 1.6 * t) * StrictMath.sin(-2. * d + 2. * f + 2. * omega)
			+ (-2274. - 0.2 * t) * StrictMath.sin(2. * f + 2. * omega)
			+ (2062. + 0.2 * t) * StrictMath.sin(2. * omega)
			+ (1426. + -3.4 * t) * StrictMath.sin(m)
			+ (712. + 0.1 * t) * StrictMath.sin(mp)
			+ (-517. + 1.2 * t) * StrictMath.sin(-2. * d + m + 2. * f + 2. * omega)
			+ (-386. - 0.4*  t) * StrictMath.sin(2. * f + omega)
			- 301. * StrictMath.sin(mp + 2. * f + 2. * omega)
			+ (217. - 0.5 * t) * StrictMath.sin(-2. * d - m + 2. * f + 2. * omega)
			- 158. * StrictMath.sin(-2. * d + mp)
			+ (129. + 0.1 * t) * StrictMath.sin(-2. * d + 2. * f + omega)
			+ 123. * StrictMath.sin(-mp + 2. * f + 2. * omega)
			+ 63. * StrictMath.sin(2. * d)
			+ (63. + 0.1 * t) * StrictMath.sin(mp + omega)
			- 59. * StrictMath.sin(2. * d - mp + 2. * f + 2. * omega)
			+ (-58. - 0.1 * t) * StrictMath.sin(-mp + omega)
			- 51. * StrictMath.sin(mp + 2. * f + omega)
			+ 48. * StrictMath.sin(-2. * d + 2. * mp)
			+ 46. * StrictMath.sin(-2. * mp + 2. * f + omega)
			- 38. * StrictMath.sin(2. * d + 2. * f + 2. * omega)
			- 31. * StrictMath.sin(2. * mp + 2. * f + 2. * omega)
			+ 29. * StrictMath.sin(2. * mp)
			+ 29. * StrictMath.sin(-2. * d + mp + 2. * f + 2. * omega)
			+ 26. * StrictMath.sin(2. * f)
			- 22. * StrictMath.sin(-2. * d + 2. * f)
			+ 21. * StrictMath.sin(-mp + 2. * f + omega)
			+ (17. - 0.1 * t) * StrictMath.sin(2. * m)
			+ 16. * StrictMath.sin(2. * d - mp + omega)
			+ (-16. + 0.1 * t) * StrictMath.sin(-2. * d + 2. * m + 2. * f + 2. * omega)
			- 15. * StrictMath.sin(m + omega)
			- 13. * StrictMath.sin(-2. * d + mp + omega)
			- 12. * StrictMath.sin(-m + omega)
			+ 11. * StrictMath.sin(2. * mp - 2. * f)
			- 10. * StrictMath.sin(2. * d - mp + 2. * f + omega)
			-8. * StrictMath.sin(2. * d + mp + 2. * f + 2. * omega)
			+ 7. * StrictMath.sin(m + 2. * f + 2. * omega)
			- 7. * StrictMath.sin(-2. * d + m + mp)
			- 7. * StrictMath.sin(-m + 2. * f + 2. * omega)
			- 8. * StrictMath.sin(2. * d + 2. * f + omega)
			+ 6. * StrictMath.sin(2. * d + mp)
			+ 6. * StrictMath.sin(-2. * d + 2. * mp + 2. * f + 2. * omega)
			+ 6. * StrictMath.sin(-2. * d + mp + 2. * f + omega)
			- 6. * StrictMath.sin(2. * d - 2. * mp + omega)
			- 6. * StrictMath.sin(2. * d + omega)
			+ 5. * StrictMath.sin(-m + mp)
			- 5. * StrictMath.sin(-2. * d - m + 2. * f + omega)
			- 5. * StrictMath.sin(-2. * d + omega)
			- 5. * StrictMath.sin(2. * mp + 2. * f + omega)
			+ 4. * StrictMath.sin(-2. * d + 2. * mp + omega)
			+ 4. * StrictMath.sin(-2. * d + m + 2. * f + omega)
			+ 4. * StrictMath.sin(mp - 2. * f)
			- 4. * StrictMath.sin(-d + mp)
			- 4. * StrictMath.sin(-2. * d + m)
			- 4. * StrictMath.sin(d)
			+ 3. * StrictMath.sin(mp + 2. * f)
			- 3. * StrictMath.sin(-2. * mp + 2. * f + 2. * omega)
			- 3. * StrictMath.sin(-d - m + mp)
			- 3. * StrictMath.sin(m + mp)
			- 3. * StrictMath.sin(-m + mp + 2. * f + 2. * omega)
			- 3. * StrictMath.sin(2. * d - m - mp + 2. * f + 2. * omega)
			- 3. * StrictMath.sin(3. * mp + 2. * f + 2. * omega)
			- 3. * StrictMath.sin(2. * d - m + 2. * f + 2. * omega));
		@SuppressWarnings("OverlyComplexArithmeticExpression")
		final double deltaEpsilon = 0.0001 * ( (92025 + 8.9 * t) * StrictMath.cos(omega)
			+ (5736. - 3.1 * t) * StrictMath.cos(-2. * d + 2. * f + 2. * omega)
			+ (977. - 0.5 * t) * StrictMath.cos(2. * f + 2. * omega)
			+ (-895. + 0.5 * t) * StrictMath.cos(2. * omega)
			+ (54. - 0.1 * t) * StrictMath.cos(m)
			- 7. * StrictMath.cos(mp)
			+ (224. - 0.6 * t) * StrictMath.cos(-2. * d + m + 2. * f + 2. * omega)
			+ 200. * StrictMath.cos(2. * f + omega)
			+ (129. - 0.1 * t) * StrictMath.cos(mp + 2. * f + 2. * omega)
			+ (-95. + 0.3 * t) * StrictMath.cos(-2. * d - m + 2. * f + 2. * omega)
			- 70. * StrictMath.cos(-2. * d + 2. * f + omega)
			- 53. * StrictMath.cos(-mp + 2. * f + 2. * omega)
			- 33. * StrictMath.cos(mp + omega)
			+ 26. * StrictMath.cos(2. * d - mp + 2. * f + 2. * omega)
			+ 32. * StrictMath.cos(-mp + omega)
			+ 27. * StrictMath.cos(mp + 2. * f + omega)
			- 24. * StrictMath.cos(-2. * mp + 2. * f + omega)
			+ 16. * StrictMath.cos(2. * d + 2. * f + 2. * omega)
			+ 13. * StrictMath.cos(2. * mp + 2. * f + 2. * omega)
			- 12. * StrictMath.cos(-2. * d + mp + 2. * f + 2. * omega)
			- 10. * StrictMath.cos(-mp + 2. * f + omega)
			- 8. * StrictMath.cos(2. * d - mp + omega)
			+ 7. * StrictMath.cos(-2. * d + 2. * m + 2. * f + 2. * omega)
			+ 9. * StrictMath.cos(m + omega)
			+ 7. * StrictMath.cos(-2. * d + mp + omega)
			+ 6. * StrictMath.cos(-m + omega)
			+ 5. * StrictMath.cos(2. * d - mp + 2. * f + omega)
			+ 3. * StrictMath.cos(2. * d + mp + 2. * f + 2. * omega)
			- 3. * StrictMath.cos(m + 2. * f + 2. * omega)
			+ 3. * StrictMath.cos(-m + 2. * f + 2. * omega)
			+ 3. * StrictMath.cos(2. * d + 2. * f + omega)
			- 3. * StrictMath.cos(-2. * d + 2. * mp + 2. * f + 2. * omega)
			- 3. * StrictMath.cos(-2. * d + mp + 2. * f + omega)
			+ 3. * StrictMath.cos(2. * d -2. * mp + omega)
			+ 3. * StrictMath.cos(2. * d + omega)
			+ 3. * StrictMath.cos(-2. * d - m + 2. * f + omega)
			+ 3. * StrictMath.cos(-2. * d + omega)
			+ 3. * StrictMath.cos(2. * mp + 2. * f + omega));

		return new double[]{deltaPsi, deltaEpsilon};
	}

	/**
	 * Calculate the Sun's apparent right ascension, AR.
	 *
	 * @param apparentEclipticObliquity	Apparent obliquity of the ecliptic [°].
	 * @param apparentLongitude	Apparent longitude of the Sun [°].
	 * @return	Sun's right ascension [°].
	 */
	private static double apparentRightAscension(final double apparentEclipticObliquity, double apparentLongitude){
		apparentLongitude = degToRad(apparentLongitude);
		return correctRangeDegree(radToDeg(StrictMath.atan2(
			StrictMath.cos(degToRad(apparentEclipticObliquity)) * StrictMath.sin(apparentLongitude),
			StrictMath.cos(apparentLongitude))));
	}

	/**
	 * Calculate the Sun's apparent declination, δ.
	 *
	 * @param apparentEclipticObliquity	Obliquity of the ecliptic, corrected for parallax [°].
	 * @param apparentLongitude	Apparent longitude of the Sun [°].
	 * @return	Sun's declination [°].
	 */
	private static double apparentDeclination(final double apparentEclipticObliquity, final double apparentLongitude){
		return radToDeg(StrictMath.asin(
			StrictMath.sin(degToRad(apparentEclipticObliquity)) * StrictMath.sin(degToRad(apparentLongitude))
		));
	}

	/**
	 * Calculate the longitude of the perihelion of the orbit, π.
	 *
	 * @param t	Julian Century in J2000.0 epoch.
	 * @return	The longitude of the perihelion of the orbit [°].
	 */
	private static double longitudeOfEarthPerihelion(final double t){
		return eval(t, new double[]{102.93735, 1.71946, 0.00046});
	}

	/**
	 * Calculate the Sun's radius angle, corrected for distance, but not for refraction.
	 *
	 * @param radius	Sun radius [^].
	 * @param distance	Sun distance [AU].
	 * @return	The longitude of the perihelion of the orbit [°].
	 */
	private static double radiusAngle(final double radius, final double distance){
		return StrictMath.asin(radius / (radius + distance));
	}

	private static double nutationAndAberrationCorrection(final double t){
		return 0.00569 + 0.00478 * StrictMath.sin(degToRad(125.04 - 1934.136 * t));
	}

	@SuppressWarnings("OverlyComplexArithmeticExpression")
	private static double equationOfTime(final double t){
		final double e0 = meanEclipticObliquity(t);
		final double epsilon = apparentEclipticObliquity(e0, t);
		final double l0 = degToRad(geometricMeanLongitude(t));
		final double e = earthOrbitEccentricity(t);
		final double m = degToRad(geometricMeanAnomaly(t));
		double y = StrictMath.tan(degToRad(epsilon) / 2.);
		y *= y;

		//[hrs]
		return radToDeg(y * StrictMath.sin(2. * l0)
			- 2. * e * StrictMath.sin(m)
			+ 4. * e * y * StrictMath.sin(m) * StrictMath.cos(2. * l0)
			- 0.5 * y * y * StrictMath.sin(4. * l0)
			- 1.25 * e * e * StrictMath.sin(2. * m)) * 4. / 60.;
	}

	/**
	 * Calculate mean Sidereal time at Greenwich, θ0.
	 *
	 * @param t	Julian Century in J2000.0 epoch.
	 * @return	mean Sidereal time at Greenwich [°].
	 */
	private static double greenwichMeanSiderealTime(final double t){
		//FIXME
//		return 18.697374558 + 24.06570982441908 * t;
		return correctRangeDegree(eval(t, new double[]{280.46061837, 360.98564736629 * 36525., 0.000387933, -1. / 38710000.}));
	}

	/**
	 * Calculate apparent Sidereal time at Greenwich, θ0.
	 *
	 * @param greenwichMeanSiderealTime	Apparent Sidereal time at Greenwich [°].
	 * @param apparentEclipticObliquity	Obliquity of the ecliptic, corrected for parallax [°].
	 * @param t	Julian Century in J2000.0 epoch.
	 * @return	apparent Sidereal time at Greenwich [°].
	 */
	private static double greenwichApparentSiderealTime(final double greenwichMeanSiderealTime, final double apparentEclipticObliquity,
			final double t){
		final double[] deltaPsiEpsilon = highAccuracyNutation(t);
		final double correction = toDegrees(0, 0, deltaPsiEpsilon[0]) * StrictMath.cos(degToRad(apparentEclipticObliquity));
		return correctRangeDegree(greenwichMeanSiderealTime + correction);
	}

	/**
	 * Calculate apparent local Sidereal time at Greenwich.
	 *
	 * @param greenwichApparentSiderealTime	Apparent Sidereal time at Greenwich [°].
	 * @return	The apparent local Sidereal time at Greenwich [°].
	 */
	private double apparentLocalSiderealTime(final double greenwichApparentSiderealTime){
		return greenwichApparentSiderealTime - location.getLongitude();
	}

	/**
	 * Calculate the hour angle of a body given its right ascension, H.
	 *
	 * @param localSiderealTime	Local Sidereal time at Greenwich [°].
	 * @param rightAscension	Right ascension [°].
	 * @return	The hour angle [°].
	 */
	private static double localHourAngle(final double localSiderealTime, final double rightAscension){
		return localSiderealTime - rightAscension;
	}


	private LocalTime computeSolarEventTime(final LocalDateTime date, final Zenith solarZenith, final boolean sunrise)
			throws SolarEventException{
		final double longitudeHour = getLongitudeHour(date, sunrise);

		final double jd = JulianDay.of(date);
		final double t = JulianDay.centuryJ2000Of(jd);
		final double equationOfTime = equationOfTime(t);
		final double geometricMeanLongitude = geometricMeanLongitude(t);
		final double meanAnomaly = geometricMeanAnomaly(t);
		final double equationOfCenter = equationOfCenter(meanAnomaly, t);
		final double trueGeometricLongitude = trueGeometricLongitude(geometricMeanLongitude, equationOfCenter);
		final double apparentGeometricLongitude = apparentGeometricLongitude(trueGeometricLongitude, t);

//		delta = m_longitude + radToDeg(hourAngle);
//		timeDiff = 4 * delta;
//		timeUTC = 720 - timeDiff - eqTime; // in minutes

final double meanEclipticObliquity = meanEclipticObliquity(t);
final double apparentEclipticObliquity = apparentEclipticObliquity(meanEclipticObliquity, t);
final double apparentDeclination = apparentDeclination(apparentEclipticObliquity, apparentGeometricLongitude);
		final double localHourAngle = localHourAngle(trueGeometricLongitude, solarZenith, sunrise, apparentDeclination);
		final double localMeanTime = getLocalMeanTime(trueGeometricLongitude, longitudeHour, localHourAngle);
		return getLocalTime(localMeanTime - equationOfTime);
	}

	/**
	 * Computes the longitude time.
	 *
	 * @return	Longitudinal time [day].
	 */
	private double getLongitudeHour(final LocalDateTime date, final Boolean sunrise){
		final double dividend = (sunrise? 6: 18) - getBaseLongitudeHour();
		return getDayOfYear(date) + dividend / 24.;
	}

	/**
	 * Computes the base longitude hour.
	 *
	 * @return	The longitude of the location of the solar event divided by 15 (deg/hour).
	 */
	private double getBaseLongitudeHour(){
		return degToHrs(location.getLongitude());
	}

	private static double getDayOfYear(final LocalDateTime date){
		return date.getDayOfYear() + JulianDay.timeOf(date);
	}

	/**
	 * @param trueLong	Sun's true longitude [°].
	 * @param zenith	Zenith enumeration.
	 * @param sunrise	Whether it's sunrise or sunset.
	 * @return	The Sun local hour [hrs].
	 * @throws SolarEventException	Whenever the Sun never rises or sets.
	 */
	private double localHourAngle(final double trueLong, final Zenith zenith, final Boolean sunrise,
			double apparentDeclination) throws SolarEventException{
		final double cosLocalHour = cosineLocalHour(trueLong, zenith);

double latRad = degToRad(location.getLatitude());
double sdRad  = degToRad(apparentDeclination);
double cosHA = StrictMath.sin(zenith.getRadians())
	/ (StrictMath.cos(latRad) * StrictMath.cos(sdRad)) - StrictMath.tan(latRad) * StrictMath.tan(sdRad);

		if(cosLocalHour < -1.)
			//the sun never sets on this location on the specified date
			throw SolarEventException.create(SolarEventError.NEVER_SETS);
		if(cosLocalHour > 1.)
			//the sun never rises on this location on the specified date
			throw SolarEventException.create(SolarEventError.NEVER_RISES);

		final double localHour = radToDeg(StrictMath.acos(cosLocalHour));
//		final double localHour = radToDeg(StrictMath.acos(cosHA));
		return degToHrs(sunrise? 360. - localHour: localHour);
	}

	private double cosineLocalHour(final double trueLong, final Zenith zenith){
		final double sinDeclination = 0.39782 * StrictMath.sin(degToRad(trueLong));
		final double cosDeclination = StrictMath.cos(StrictMath.asin(sinDeclination));

		final double sinZenith = StrictMath.sin(zenith.getRadians());
		final double sinLatitude = StrictMath.sin(degToRad(location.getLatitude()));
		final double cosLatitude = StrictMath.cos(degToRad(location.getLongitude()));

		return (sinZenith - sinLatitude * sinDeclination) / (cosLatitude * cosDeclination);
	}

	/**
	 * Computes the local mean time.
	 *
	 * @param trueGeometricLongitude	Suns true longitude [°].
	 * @param longitudeHour	Longitude hour [hrs].
	 * @param localHourAngle	Local hour of the Sun [hrs].
	 * @return	The local mean time [hrs].
	 */
	private static double getLocalMeanTime(final double trueGeometricLongitude, final double longitudeHour, final double localHourAngle){
		final double rightAscension = getRightAscension(trueGeometricLongitude);
		return correctRangeHour(localHourAngle + rightAscension - 0.06571 * longitudeHour - 6.622);
	}

	/**
	 * Computes the Suns right ascension, adjusting for the quadrant of the true longitude of the Sun and turning it
	 * into degree-hours.
	 *
	 * @param sunTrueLongitude	Suns true longitude [°].
	 * @return	Suns right ascension in degree-hours.
	 */
	@SuppressWarnings("NumericCastThatLosesPrecision")
	private static double getRightAscension(final double sunTrueLongitude){
		final double tanL = StrictMath.tan(degToRad(sunTrueLongitude));
		final double rightAscension = correctRangeDegree(radToDeg(StrictMath.atan(0.91764 * tanL)));

		final double longitudeQuadrant = 90. * (int)(sunTrueLongitude / 90.);
		final double rightAscensionQuadrant = 90. * (int)(rightAscension / 90.);
		return (rightAscension + longitudeQuadrant - rightAscensionQuadrant) / 15.;
	}

	@SuppressWarnings("NumericCastThatLosesPrecision")
	private LocalTime getLocalTime(final double localMeanTime){
		//adjust back to UTC
		final double utcTime = correctRangeHour(localMeanTime - getBaseLongitudeHour());
		return LocalTime.of(0, 0)
			.plus((long)(utcTime * 60. * 60.), ChronoUnit.SECONDS);
	}



	public static double toDegrees(final int degree, final int minute, final double second){
		return degree + (minute + second / 60.) / 60.;
	}

	private static double correctRangeDegree(double degree){
		degree %= 360.;
		return (degree < 0.? degree + 360.: degree);
	}

	private static double correctRangeHour(double degree){
		degree %= 24;
		return (degree < 0.? degree + 24.: degree);
	}

	private static double degToRad(final double degrees){
		return degrees * (StrictMath.PI / 180.);
	}

	private static double radToDeg(final double radians){
		return radians * (180. / StrictMath.PI);
	}

	private static double degToHrs(final double degrees){
		return degrees / 15.;
	}

	// use Horner's method to compute and return the polynomial evaluated at x
	// p[0) + p[1] x^1 + p[2] x^2 + ... + p[n-1] x^n-1
	private static double eval(final double x, final double[] p){
		double result = 0.;
		for(int i = p.length - 1; i >= 0; i --)
			result = p[i] + x * result;
		return result;
	}

}
