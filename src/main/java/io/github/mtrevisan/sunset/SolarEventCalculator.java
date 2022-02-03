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

import java.time.LocalDate;
import java.time.LocalDateTime;
import java.time.LocalTime;
import java.time.temporal.ChronoUnit;


/*
https://github.com/mikereedell/sunrisesunsetlib-java

https://web.archive.org/web/20161202180207/http://williams.best.vwh.net/sunrise_sunset_algorithm.htm
https://ebvalaim.pl/en/2015/12/22/calculating-sunrise-and-sunset-times/
https://en.wikipedia.org/wiki/Position_of_the_Sun

others
https://github.com/caarmen/SunriseSunset
https://github.com/MenoData/Time4J/blob/master/base/src/main/java/net/time4j/calendar/astro/SolarTime.java
https://github.com/shred/commons-suncalc/blob/master/src/main/java/org/shredzone/commons/suncalc/SunTimes.java
*/
public class SolarEventCalculator{

	private final Location location;


	/**
	 * Constructs a new instance using the given parameters.
	 *
	 * @param location	Location of the place.
	 */
	public static SolarEventCalculator create(final Location location){
		return new SolarEventCalculator(location);
	}


	private SolarEventCalculator(final Location location){
		this.location = location;
	}


//---
	//https://squarewidget.com/solar-coordinates/

	/**
	 * Calculate Julian Day at 0 UTC.
	 *
	 * @param date	The date.
	 * @return	The Julian Day [day].
	 */
	public static double julianDay(final LocalDate date){
		return julianDay(date.getYear(), date.getMonthValue(), date.getDayOfMonth());
	}

	/**
	 * Calculate Julian Day at 0 UTC.
	 *
	 * @param year	The year.
	 * @param month	The month (1 is January).
	 * @param day	The day.
	 * @return	The Julian Day [day].
	 */
	static double julianDay(int year, int month, final int day){
		if(month <= 2){
			year --;
			month += 12;
		}

		final double a = StrictMath.floor(year / 100.);
		final double b = 2. - a + StrictMath.floor(a / 4.);
		return StrictMath.floor(365.25 * (year + 4716)) + StrictMath.floor(30.6001 * (month + 1)) + day + b - 1524.5;
	}

	/**
	 * Calculated the Sun position.
	 *
	 * @param jd	The Julian Day [day].
	 * @return	The Sun position.
	 */
	public final EquatorialCoordinate sunPosition(final double jd){
		final double t = julianCentury(jd);
		final double geometricMeanLongitude = geometricMeanLongitude(t);
		final double meanAnomaly = meanAnomaly(t);
		final double equationOfCenter = equationOfCenter(meanAnomaly, t);
		final double trueGeometricLongitude = trueGeometricLongitude(geometricMeanLongitude, equationOfCenter);
		final double apparentLongitude = apparentGeometricLongitude(trueGeometricLongitude, t);
		final double meanEclipticObliquity = meanEclipticObliquity(t);
		final double apparentEclipticObliquity = apparentEclipticObliquity(meanEclipticObliquity, apparentLongitude);
		final double apparentRightAscension = apparentRightAscension(apparentEclipticObliquity, apparentLongitude);
		final double apparentDeclination = apparentDeclination(apparentEclipticObliquity, apparentLongitude);
		return EquatorialCoordinate.create(apparentRightAscension, apparentDeclination);
	}

	//TODO
	/**
	 * Computes the sunset time.
	 *
	 * @param solarZenith	Enumeration corresponding to the type of sunset to compute.
	 * @param date	Date to compute the sunset for.
	 * @return	The sunset time or {@code null} for no sunset.
	 * @throws SolarEventException	Whenever the Sun never rises or sets.
	 */
	public final LocalDateTime sunset(final LocalDate date, final Zenith solarZenith) throws SolarEventException{
/*
// Calculate the time when the upper limb of the sun just crosses the
// horizon.
// Returns undef if the sun does not rise or set on that date at that latitude.
sunset[date, lat, long, temp = 283 K, pressure=1010 millibars ] :=
{
   edate = approxMidnight[date,long] + 17.5 hours
   [ra, decl] = sunApparentRADecl[edate]
   [rise, set] = approxRiseSet[edate, lat, long, ra, decl]
   if set == undef
      return undef
   else
      return sunSecantAltitude[set, lat, long, -16 arcmin, temp, pressure]
}*/

		final double jd = julianDay(date);
		final double t = julianCentury(jd);
		final EquatorialCoordinate coord = sunPosition(jd);
		//[hrs]
		final double ra = coord.getRightAscension() / 15.;
		final double decl = coord.getLongitude();

		final double geometricMeanLongitude = geometricMeanLongitude(t);
		final double meanAnomaly = meanAnomaly(t);
		final double eccentricity = earthOrbitEccentricity(t);
		final double equationOfCenter = equationOfCenter(meanAnomaly, t);
		final double trueGeometricLongitude = trueGeometricLongitude(geometricMeanLongitude, equationOfCenter);
		final double trueAnomaly = trueAnomaly(meanAnomaly, equationOfCenter);
		final double radiusVector = radiusVector(meanAnomaly);
		final double radiusVector2 = radiusVector(eccentricity, trueAnomaly);
		final double apparentLongitude = apparentGeometricLongitude(trueGeometricLongitude, t);
		final double meanEclipticObliquity = meanEclipticObliquity(t);
		final double apparentEclipticObliquity = apparentEclipticObliquity(meanEclipticObliquity, apparentLongitude);
		final double longitudeOfEarthPerihelion = longitudeOfEarthPerihelion(t);

		System.out.println(coord);

		final double greenwichMeanSiderealTime = greenwichMeanSiderealTime(t);
		final double greenwichApparentSiderealTime = greenwichApparentSiderealTime(greenwichMeanSiderealTime, apparentEclipticObliquity, t);
		final double apparentLocalSiderealTime = apparentLocalSiderealTime(greenwichApparentSiderealTime, location);
		final double localHourAngle = localHourAngle(apparentLocalSiderealTime, coord.getRightAscension());

		//[hrs]
		final LocalTime localTime = computeSolarEventTime(solarZenith, date, false);
		return LocalDateTime.of(date, localTime);
	}

	/**
	 * Calculated the Julian Century since JDE2451545, that is J2000.0.
	 *
	 * @param jd	The Julian Day [day].
	 * @return	The Julian Century.
	 */
	private double julianCentury(final double jd){
		return (jd - 2451545.) / 36525.;
	}

	/**
	 * Calculate the geometric mean longitude of the Sun, referred to the mean equinox of the date, L0.
	 *
	 * @param t	Julian Century in J2000.0 epoch.
	 * @return	The geometric mean longitude of the Sun [°].
	 */
	private double geometricMeanLongitude(final double t){
		return correctRangeDegree(eval(t, new double[]{toDegrees(280, 27, 59.26), 36000.76983, 0.0003032}));
	}

	/**
	 * Calculate the mean anomaly of the Sun, M.
	 *
	 * @param t	Julian Century in J2000.0 epoch.
	 * @return	The mean anomaly of the Sun [°].
	 */
	private double meanAnomaly(final double t){
		return correctRangeDegree(eval(t, new double[]{357.52911, 35999.05029, -0.0001537}));
	}

	/**
	 * Calculate the eccentricity of Earth's orbit, e.
	 *
	 * @param t	Julian Century in J2000.0 epoch.
	 * @return	The eccentricity of Earth's orbit.
	 */
	private double earthOrbitEccentricity(final double t){
		return eval(t, new double[]{0.016708634, -0.000042037, -0.0000001267});
	}

	/**
	 * Calculate the Sun's equation of center, C.
	 *
	 * @param meanAnomaly	The mean anomaly of the Sun [°].
	 * @param t	Julian Century in J2000.0 epoch.
	 * @return	The Sun's equation of center [°].
	 */
	private double equationOfCenter(double meanAnomaly, final double t){
		meanAnomaly = convertDegreesToRadians(meanAnomaly);
		return eval(t, new double[]{1.914602, -0.004817, -0.000014}) * StrictMath.sin(meanAnomaly)
			+ eval(t, new double[]{0.019993, -0.000101}) * StrictMath.sin(meanAnomaly * 2.)
			+ 0.000289 * StrictMath.sin(meanAnomaly * 3.);
	}

	/**
	 * Calculate the Sun's true geometric longitude, Ltrue = L0 + C.
	 *
	 * @param geometricMeanLongitude	The geometric mean longitude of the Sun, referred to the mean equinox of the date [°].
	 * @param equationOfCenter	The Sun's equation of center [°].
	 * @return	Sun's true geometric longitude [°].
	 */
	private double trueGeometricLongitude(final double geometricMeanLongitude, final double equationOfCenter){
		return correctRangeDegree(geometricMeanLongitude + equationOfCenter);
	}

	/**
	 * Calculate the Sun's true anomaly longitude, ν = M + C.
	 *
	 * @param meanAnomaly	The mean anomaly of the Sun [°].
	 * @param equationOfCenter	The Sun's equation of center [°].
	 * @return	Sun's true geometric longitude [°].
	 */
	private double trueAnomaly(final double meanAnomaly, final double equationOfCenter){
		return correctRangeDegree(meanAnomaly + equationOfCenter);
	}

	/**
	 * Calculate the distance between the center of the Sun and the center of the Earth, R.
	 * <p>U.S. Naval Observatory function.</p>
	 *
	 * @param meanAnomaly	The mean anomaly of the Sun [°].
	 * @return	Distance between the center of the Sun and the center of the Earth [AU].
	 */
	private double radiusVector(double meanAnomaly){
		meanAnomaly = convertDegreesToRadians(meanAnomaly);
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
	private double radiusVector(final double eccentricity, final double trueAnomaly){
		return 1.000001018 * (1. - eccentricity * eccentricity) / (1. + eccentricity * StrictMath.cos(convertDegreesToRadians(trueAnomaly)));
	}

	/**
	 * Calculate the apparent longitude of the Sun, Lapp = Ω.
	 *
	 * @param trueGeometricLongitude	Sun's true geometric longitude [°].
	 * @param t	Julian Century in J2000.0 epoch.
	 * @return	Apparent longitude of the Sun [°].
	 */
	private double apparentGeometricLongitude(final double trueGeometricLongitude, final double t){
		//correction for nutation and aberration
		//[rad]
		final double omega = convertDegreesToRadians(125.04 - 1934.136 * t);
		return trueGeometricLongitude - 0.00569 - 0.00478 * StrictMath.sin(omega);
	}

	/**
	 * Calculate the mean obliquity of the ecliptic, ɛ0.
	 *
	 * @param t	Julian Century in J2000.0 epoch.
	 * @return	Apparent longitude of the Sun [°].
	 */
	private double meanEclipticObliquity(final double t){
		final double u = t / 100.;
		return toDegrees(23, 26, 21.448)
			+ toDegrees(0, 0, eval(u, new double[]{0., -4680.93, -1.55, 1999.25, -51.38, -249.67, -39.05, 7.12, 27.87, 5.79,
			2.45}));
	}

	/**
	 * Calculate the mean obliquity of the ecliptic, corrected for parallax, ɛ'.
	 *
	 * @param meanEclipticObliquity	Mean obliquity of the ecliptic [°].
	 * @param apparentLongitude	Apparent longitude of the Sun [°].
	 * @return	Apparent longitude of the Sun [°].
	 */
	private double apparentEclipticObliquity(final double meanEclipticObliquity, final double apparentLongitude){
		return meanEclipticObliquity + 0.00256 * StrictMath.cos(convertDegreesToRadians(apparentLongitude));
	}

	/**
	 * True obliquity of the ecliptic corrected for nutation, ɛ.
	 *
	 * @param meanEclipticObliquity	Obliquity of the ecliptic, corrected for parallax [°].
	 * @param t	Julian Century in J2000.0 epoch.
	 * @return	Apparent longitude of the Sun [°].
	 */
	private double trueEclipticObliquity(final double meanEclipticObliquity, final double t){
		final double[] deltaPsiEpsilon = highAccuracyNutation(t);
		return meanEclipticObliquity + toDegrees(0, 0, deltaPsiEpsilon[1]);
	}

	/**
	 * Calculate Nutation in longitude (delta psi) and obliquity (delta epsilon).
	 *
	 * @param t	Julian Century in J2000.0 epoch.
	 * @return	An array where the first element is delta psi, and the second delta epsilon ["].
	 */
	private double[] highAccuracyNutation(final double t){
		//mean elongation of the Moon from the Sun [°]
		final double d = eval(t, new double[]{297.85036, 445267.111480, -0.0019142, 1. / 189474.});
		//mean anomaly of the Sun [°]
		final double m = eval(t, new double[]{357.52772, 35999.050340, -0.0001603, - 1. / 300000.});
		//mean anomaly of the Moon [°]
		final double mp = eval(t, new double[]{134.96298, 477198.867398, 0.0086972, 1. / 56250.});
		//Moon's argument of Latitude [°]
		final double f = eval(t, new double[]{93.27191, 483202.017538, -0.0036825, 1. / 327270.});
		//Longitude of the ascending node of the Moon's mean orbit on the ecliptic measured from the mean equinox of the date [rad]
		final double omega = convertDegreesToRadians(eval(t, new double[]{125.04452, -1934.136261, 0.0020708, 1. / 450000.}));

		//these lines generated by iau1980.frink and pasted in here ["]
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
	private double apparentRightAscension(final double apparentEclipticObliquity, double apparentLongitude){
		apparentLongitude = convertDegreesToRadians(apparentLongitude);
		return correctRangeDegree(convertRadiansToDegrees(StrictMath.atan2(
			StrictMath.cos(convertDegreesToRadians(apparentEclipticObliquity)) * StrictMath.sin(apparentLongitude),
			StrictMath.cos(apparentLongitude))));
	}

	/**
	 * Calculate the Sun's apparent declination, δ.
	 *
	 * @param apparentEclipticObliquity	Obliquity of the ecliptic, corrected for parallax [°].
	 * @param apparentLongitude	Apparent longitude of the Sun [°].
	 * @return	Sun's declination [°].
	 */
	private double apparentDeclination(final double apparentEclipticObliquity, final double apparentLongitude){
		return convertRadiansToDegrees(StrictMath.asin(
			StrictMath.sin(convertDegreesToRadians(apparentEclipticObliquity)) * StrictMath.sin(convertDegreesToRadians(apparentLongitude))
		));
	}

	/**
	 * Calculate the longitude of the perihelion of the orbit, π.
	 *
	 * @param t	Julian Century in J2000.0 epoch.
	 * @return	The longitude of the perihelion of the orbit [°].
	 */
	private double longitudeOfEarthPerihelion(final double t){
		return eval(t, new double[]{102.93735, 1.71946, 0.00046});
	}

	/**
	 * Calculate the Sun's radius angle, corrected for distance, but not for refraction.
	 *
	 * @param radius	Sun radius [^].
	 * @param distance	Sun distance [AU].
	 * @return	The longitude of the perihelion of the orbit [°].
	 */
	private double radiusAngle(final double radius, final double distance){
		return StrictMath.asin(radius / (radius + distance));
	}

	/**
	 * Calculate mean Sidereal time at Greenwich, θ0.
	 *
	 * @param t	Julian Century in J2000.0 epoch.
	 * @return	mean Sidereal time at Greenwich [°].
	 */
	private double greenwichMeanSiderealTime(final double t){
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
	private double greenwichApparentSiderealTime(final double greenwichMeanSiderealTime, final double apparentEclipticObliquity,
			final double t){
		final double[] deltaPsiEpsilon = highAccuracyNutation(t);
		final double correction = toDegrees(0, 0, deltaPsiEpsilon[0]) * StrictMath.cos(convertDegreesToRadians(apparentEclipticObliquity));
		return correctRangeDegree(greenwichMeanSiderealTime + correction);
	}

	/**
	 * Calculate apparent local Sidereal time at Greenwich.
	 *
	 * @param greenwichApparentSiderealTime	Apparent Sidereal time at Greenwich [°].
	 * @param location	Location of the place.
	 * @return	The apparent local Sidereal time at Greenwich [°].
	 */
	private double apparentLocalSiderealTime(final double greenwichApparentSiderealTime, final Location location){
		return greenwichApparentSiderealTime - location.getLongitude();
	}

	/**
	 * Calculate the hour angle of a body given its right ascension, H.
	 *
	 * @param localSiderealTime	Local Sidereal time at Greenwich [°].
	 * @param rightAscension	Right ascension [°].
	 * @return	The hour angle [°].
	 */
	private double localHourAngle(final double localSiderealTime, final double rightAscension){
		return localSiderealTime - rightAscension;
	}


	private LocalTime computeSolarEventTime(final Zenith solarZenith, final LocalDate date, final boolean sunrise)
			throws SolarEventException{
		final double longitudeHour = getLongitudeHour(date, sunrise);

		final double jd = julianDay(date);
		final double t = julianCentury(jd);
		final double geometricMeanLongitude = geometricMeanLongitude(t);
		final double meanAnomaly = meanAnomaly(t);
		final double equationOfCenter = equationOfCenter(meanAnomaly, t);
		final double trueGeometricLongitude = trueGeometricLongitude(geometricMeanLongitude, equationOfCenter);
		final double apparentGeometricLongitude = apparentGeometricLongitude(trueGeometricLongitude, t);

		final double cosLocalHour = getCosineSunLocalHour(trueGeometricLongitude, solarZenith);
		if(cosLocalHour < -1.)
			//the sun never sets on this location on the specified date
			throw SolarEventException.create(SolarEventError.NEVER_SETS);
		if(cosLocalHour > 1.)
			//the sun never rises on this location on the specified date
			throw SolarEventException.create(SolarEventError.NEVER_RISES);

		final double sunLocalHour = getSunLocalHour(cosLocalHour, sunrise);
		final double localMeanTime = getLocalMeanTime(trueGeometricLongitude, longitudeHour, sunLocalHour);
		return getLocalTime(localMeanTime);
	}

	/**
	 * Calculate the parallactic angle for a given body from a point on the Earth, q.
	 *
	 * @param hourAngle	Hour angle [°].
	 * @param declination	Declination [°].
	 * @param location	Location of the place.
	 * @param t	Julian Century in J2000.0 epoch.
	 * @return	The hour angle [°].
	 */
	private double parallatticAngle(double hourAngle, final double declination, final Location location, final double t){
		hourAngle = convertDegreesToRadians(hourAngle);
		final double latitude = convertDegreesToRadians(location.getLatitude());
		final double decl = convertDegreesToRadians(declination);
		return StrictMath.atan2(StrictMath.sin(hourAngle), StrictMath.tan(latitude) * StrictMath.cos(decl) - StrictMath.sin(decl) * StrictMath.cos(hourAngle));
	}

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

	/**
	 * Computes the longitude time.
	 *
	 * @return	Longitudinal time [day].
	 */
	private double getLongitudeHour(final LocalDate date, final Boolean sunrise){
		final double dividend = (sunrise? 6: 18) - getBaseLongitudeHour();
		return getDayOfYear(date) + dividend / 24.;
	}

	/**
	 * Computes the base longitude hour.
	 *
	 * @return	The longitude of the location of the solar event divided by 15 (deg/hour).
	 */
	private double getBaseLongitudeHour(){
		return location.getLongitude() / 15.;
	}

	private int getDayOfYear(final LocalDate date){
		return date.getDayOfYear();
	}

	private double getCosineSunLocalHour(final double sunTrueLong, final Zenith zenith){
		final double sinDeclination = 0.39782 * StrictMath.sin(convertDegreesToRadians(sunTrueLong));
		final double cosDeclination = StrictMath.cos(StrictMath.asin(sinDeclination));

		final double cosZenith = StrictMath.cos(zenith.getRadians());
		final double sinLatitude = StrictMath.sin(convertDegreesToRadians(location.getLatitude()));
		final double cosLatitude = StrictMath.cos(convertDegreesToRadians(location.getLongitude()));

		return (cosZenith - sinDeclination * sinLatitude) / (cosDeclination * cosLatitude);
	}

	/**
	 * @param cosSunLocalHour	Cosine of Sun's local hour.
	 * @param sunrise	Whether it's sunrise or sunset.
	 * @return	The Sun local hour [hrs].
	 */
	private double getSunLocalHour(final double cosSunLocalHour, final Boolean sunrise){
		double localHour = convertRadiansToDegrees(StrictMath.acos(cosSunLocalHour));
		if(sunrise)
			localHour = 360. - localHour;

		return localHour / 15.;
	}

	/**
	 *
	 * @param sunTrueLongitude
	 * @param longitudeHour
	 * @param sunLocalHour
	 * @return	The local mean time [hrs].
	 */
	private double getLocalMeanTime(final double sunTrueLongitude, final double longitudeHour, final double sunLocalHour){
		final double rightAscension = getRightAscension(sunTrueLongitude);
		double localMeanTime = sunLocalHour + rightAscension - 0.06571 * longitudeHour - 6.622;
		if(localMeanTime < 0.)
			localMeanTime += 24.;
		else if(localMeanTime > 24.)
			localMeanTime -= 24.;
		return localMeanTime;
	}

	/**
	 * Computes the Suns right ascension, adjusting for the quadrant of the true longitude of the Sun and turning it
	 * into degree-hours.
	 *
	 * @param sunTrueLong	Suns true longitude [°].
	 * @return	Suns right ascension in degree-hours.
	 */
	private double getRightAscension(final double sunTrueLong){
		final double tanL = StrictMath.tan(convertDegreesToRadians(sunTrueLong));
		double rightAscension = StrictMath.atan(convertDegreesToRadians(0.91764 * convertRadiansToDegrees(tanL)));
		rightAscension = convertRadiansToDegrees(rightAscension);
		if(rightAscension < 0.)
			rightAscension += 360.;
		else if(rightAscension > 360.)
			rightAscension -= 360.;

		final double longitudeQuadrant = 90. * (int)(sunTrueLong / 90.);
		final double rightAscensionQuadrant = 90. * (int)(rightAscension / 90.);
		return (rightAscension + longitudeQuadrant - rightAscensionQuadrant) / 15.;
	}

	private LocalTime getLocalTime(final double localMeanTime){
		double utcTime = localMeanTime - getBaseLongitudeHour();
		utcTime = adjustForDST(utcTime);
		return LocalTime.of(0, 0)
			.plus((long)(utcTime * 24. * 60. * 60.), ChronoUnit.SECONDS);
	}

	private double adjustForDST(final double localMeanTime){
		double localTime = localMeanTime;
		if(localTime > 24.)
			localTime -= 24.;
		return localTime;
	}



	private static double toDegrees(final int degree, final int minute, final double second){
		return degree + (minute + second / 60.) / 60.;
	}

	private static double correctRangeDegree(double degree){
		degree %= 360.;
		return (degree < 0.? degree + 360.: degree);
	}

	private static double correctRangeHour(double degree){
		degree %= 360.;
		return (degree < 0.? degree + 360.: degree);
	}

	private static double convertDegreesToRadians(final double degrees){
		return degrees * (StrictMath.PI / 180.);
	}

	private static double convertRadiansToDegrees(final double radians){
		return radians * (180. / StrictMath.PI);
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
