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

import java.io.IOException;
import java.time.LocalDate;
import java.time.LocalDateTime;
import java.time.LocalTime;
import java.time.temporal.ChronoUnit;
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

https://www.iers.org/SharedDocs/Publikationen/EN/IERS/Publications/tn/TechnNote13/tn13.pdf?__blob=publicationFile&v=1
*/
public class SolarEventCalculator{

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

	private static final double EARTH_FLATTENING = 1. / 298.25642;
	//[m]
	private static final double EARTH_EQUATORIAL_RADIUS = 6378140.;


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
	 *
	 * @see <a href="https://squarewidget.com/solar-coordinates/">Solar coordinates</>
	 */
	public static EquatorialCoordinate sunPosition(final double jd){
		final double t = JulianDay.centuryJ2000Of(jd);

		//calculate the geometric mean longitude L0 of the Sun referred to the mean equinox of the time T: L0
		final double geometricMeanLongitude = geometricMeanLongitude(t);
		//calculate the mean anomaly of the Sun at time T: M
		final double meanAnomaly = geometricMeanAnomaly(t);
		//calculate the Sun’s equation of the center at time T: C
		final double equationOfCenter = equationOfCenter(meanAnomaly, t);
		//calculate the Sun’s true longitude: Ltrue = L0 + C
		final double trueGeometricLongitude = correctRangeDegree(geometricMeanLongitude + equationOfCenter);
		//calculate the Sun’s true anomaly: ν = M + C
		final double trueAnomaly = correctRangeDegree(meanAnomaly + equationOfCenter);
		//correct for nutation and aberration in order to get the Sun’s apparent longitude referred to the true equinox of time T: Lapp
		final double[] nutationInLongitudeAndObliquity = correctionNutationInLongitudeAndObliquity(t);
		final double aberration = correctionAberration(radiusVector(t));
		final double apparentGeometricLatitude = apparentGeometricLongitude(geometricMeanLongitude, nutationInLongitudeAndObliquity[0],
			aberration);
		final double apparentGeometricLongitude = apparentGeometricLongitude(geometricMeanLongitude, nutationInLongitudeAndObliquity[0],
			aberration);
		//calculate the obliquity of the ecliptic (the inclination of the Earth’s equator with respect to the plane at which the Sun
		//and planets appear to move across the sky): ɛ0
		final double meanEclipticObliquity = meanEclipticObliquity(t);

		//calculate the apparent position of the Sun on the celestial sphere at time T:
		final double apparentEclipticObliquity = apparentEclipticObliquity(meanEclipticObliquity, t);
		final double apparentRightAscension = rightAscension(apparentGeometricLatitude, apparentEclipticObliquity, apparentGeometricLongitude);
		final double apparentDeclination = declination(apparentGeometricLatitude, apparentGeometricLongitude, apparentEclipticObliquity);
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
		//[h]
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
		//Ltrue = L0 + C
		final double trueGeometricLongitude = correctRangeDegree(geometricMeanLongitude + equationOfCenter);
		//ν = M + C
		final double trueAnomaly = correctRangeDegree(meanAnomaly + equationOfCenter);
		final double radiusVector = radiusVector(t);
		final double radiusVector2 = radiusVector(eccentricity, trueAnomaly);
		final double[] nutationInLongitudeAndObliquity = correctionNutationInLongitudeAndObliquity(t);
		final double aberration = correctionAberration(radiusVector(t));
		final double apparentGeometricLongitude = apparentGeometricLongitude(geometricMeanLongitude, nutationInLongitudeAndObliquity[0],
			aberration);
		final double meanEclipticObliquity = meanEclipticObliquity(t);
		final double apparentEclipticObliquity = apparentEclipticObliquity(meanEclipticObliquity, t);
		final double longitudeOfEarthPerihelion = longitudeOfEarthPerihelion(t);

		final double meanSiderealTime = meanSiderealTime(t);
		final double apparentSiderealTime = apparentSiderealTime(meanSiderealTime, apparentEclipticObliquity, t);
		final double apparentLocalSiderealTime = localMeanSiderealTime(apparentSiderealTime, location);
		final double localHourAngle = localHourAngle(apparentLocalSiderealTime, coord.getRightAscension());

		LocalDateTime sunsetLocalDateTime = LocalDateTime.of(date, LocalTime.MIDNIGHT);
		LocalDateTime localDateTime;
		do{
			localDateTime = sunsetLocalDateTime;

			//[h]
			final LocalTime localTime = computeSolarEventTime(localDateTime, solarZenith, false);
			sunsetLocalDateTime = LocalDateTime.of(date, localTime);
		}while(!localDateTime.equals(sunsetLocalDateTime));
		return sunsetLocalDateTime;
	}

	/**
	 * Calculate the geometric mean longitude of the Sun, referred to the mean equinox of the date, L.
	 *
	 * @param tt	Julian Century of Terrestrial Time from J2000.0.
	 * @return	The geometric mean longitude of the Sun [°].
	 *
	 * @see <a href="https://squarewidget.com/solar-coordinates/">Solar coordinates</>
	 */
	private static double geometricMeanLongitude(final double tt){
		final double jme = tt / 10.;
		final double[] parameters = new double[6];
		for(int i = 0; i < parameters.length; i ++){
			double parameter = 0.;
			final Collection<Double[]> elements = EARTH_HELIOCENTRIC_DATA.get("L" + i);
			for(final Double[] element : elements)
				parameter += element[0] * StrictMath.cos(element[1] + element[2] * jme);
			parameters[i] = parameter;
		}
		final double longitude = eval(jme, parameters) / 100_000_000.;
		return correctRangeDegree(StrictMath.toDegrees(longitude) + 180.);
	}

	/**
	 * Calculate the geometric mean latitude of the Sun, referred to the mean equinox of the date, β.
	 *
	 * @param tt	Julian Century of Terrestrial Time from J2000.0.
	 * @return	The geometric mean latitude of the Sun [°].
	 *
	 * @see <a href="https://squarewidget.com/solar-coordinates/">Solar coordinates</>
	 */
	private static double geometricMeanLatitude(final double tt){
		final double jme = tt / 10.;
		final double[] parameters = new double[2];
		for(int i = 0; i < parameters.length; i ++){
			double parameter = 0.;
			final Collection<Double[]> elements = EARTH_HELIOCENTRIC_DATA.get("B" + i);
			for(final Double[] element : elements)
				parameter += element[0] * StrictMath.cos(element[1] + element[2] * jme);
			parameters[i] = parameter;
		}
		final double latitude = eval(jme, parameters) / 100_000_000.;
		return correctRangeDegree(-StrictMath.toDegrees(latitude));
	}

	/**
	 * Calculate the geometric mean anomaly of the Sun, M.
	 *
	 * @param tt	Julian Century of Terrestrial Time from J2000.0.
	 * @return	The geometric mean anomaly of the Sun [°].
	 */
	private static double geometricMeanAnomaly(final double tt){
		return correctRangeDegree(eval(tt, new double[]{357.52911, 35999.05029, -0.0001537}));
	}

	/**
	 * Calculate the Sun's equation of center, C.
	 *
	 * @param geometricMeanAnomaly	The mean anomaly of the Sun [°].
	 * @param tt	Julian Century of Terrestrial Time from J2000.0.
	 * @return	The Sun's equation of center [°].
	 */
	private static double equationOfCenter(double geometricMeanAnomaly, final double tt){
		geometricMeanAnomaly = StrictMath.toRadians(geometricMeanAnomaly);
		return eval(tt, new double[]{1.914602, -0.004817, -0.000014}) * StrictMath.sin(geometricMeanAnomaly)
			+ eval(tt, new double[]{0.019993, -0.000101}) * StrictMath.sin(2. * geometricMeanAnomaly)
			+ 0.000289 * StrictMath.sin(3. * geometricMeanAnomaly);
	}

	/**
	 * Calculate the mean obliquity of the ecliptic, ɛ0.
	 *
	 * @param tt	Julian Century of Terrestrial Time from J2000.0.
	 * @return	Apparent longitude of the Sun [°].
	 */
	private static double meanEclipticObliquity(final double tt){
		final double u = tt / 100.;
		final double seconds = eval(u, new double[]{21.448, -4680.93, -1.55, 1999.25, -51.38, -249.67, -39.05, 7.12, 27.87, 5.79, 2.45});
		return toDegrees(23, 26, seconds);
	}

	/**
	 * Calculate the apparent longitude of the Sun, Lapp = λ.
	 *
	 * @param geometricMeanLongitude	Sun's true geometric longitude [°].
	 * @param deltaPsi	Nutation in longitude [°].
	 * @param deltaAberration	Aberration [°].
	 * @return	Apparent longitude of the Sun [°].
	 */
	private static double apparentGeometricLongitude(final double geometricMeanLongitude, final double deltaPsi,
			final double deltaAberration){
		return geometricMeanLongitude + deltaPsi + deltaAberration;
	}

	/**
	 * Calculate the mean obliquity of the ecliptic, corrected for parallax, ɛ'.
	 *
	 * @param meanEclipticObliquity	Mean obliquity of the ecliptic [°].
	 * @param tt	Julian Century of Terrestrial Time from J2000.0.
	 * @return	Apparent longitude of the Sun [°].
	 */
	private static double apparentEclipticObliquity(final double meanEclipticObliquity, final double tt){
		return meanEclipticObliquity + 0.00256 * StrictMath.cos(StrictMath.toRadians(omega(tt)));
	}

	private static double nutationAndAberrationCorrection(final double tt){
		return 0.00569 + 0.00478 * StrictMath.sin(StrictMath.toRadians(omega(tt)));
	}

	private static double omega(final double tt){
		return 125.04 - 1934.136 * tt;
	}

	/**
	 * Calculate the Sun's right ascension, AR.
	 *
	 * @param geometricMeanLatitude   Geometric mean latitude of the Sun [°].
	 * @param eclipticObliquity	Obliquity of the ecliptic [°].
	 * @param longitude	Longitude of the Sun [°].
	 * @return	Sun's right ascension [°].
	 */
	private static double rightAscension(final double geometricMeanLatitude, double eclipticObliquity, double longitude){
		longitude = StrictMath.toRadians(longitude);
		eclipticObliquity = StrictMath.toRadians(eclipticObliquity);
		return correctRangeDegree(StrictMath.toDegrees(StrictMath.atan2(
			StrictMath.sin(longitude) * StrictMath.cos(eclipticObliquity)
				- StrictMath.tan(StrictMath.toRadians(geometricMeanLatitude)) * StrictMath.sin(eclipticObliquity),
			StrictMath.cos(longitude))));
	}

	/**
	 * Calculate the Sun's declination, δ.
	 *
	 * @param geometricMeanLatitude   Geometric mean latitude of the Sun [°].
	 * @param apparentGeometricLongitude   Apparent longitude of the Sun [°].
	 * @param trueEclipticObliquity   True obliquity of the ecliptic [°].
	 * @return	Sun's declination [°].
	 */
	private static double declination(double geometricMeanLatitude, final double apparentGeometricLongitude, double trueEclipticObliquity){
		geometricMeanLatitude = StrictMath.toRadians(geometricMeanLatitude);
		trueEclipticObliquity = StrictMath.toRadians(trueEclipticObliquity);
		return StrictMath.toDegrees(StrictMath.asin(
			StrictMath.sin(geometricMeanLatitude) * StrictMath.cos(trueEclipticObliquity)
			+ StrictMath.cos(geometricMeanLatitude) * StrictMath.sin(trueEclipticObliquity) * StrictMath.sin(StrictMath.toRadians(apparentGeometricLongitude))
		));
	}

	//https://www.nrel.gov/docs/fy08osti/34302.pdf
	public static void main(String[] args){
//		final double jd = JulianDay.of(1997, 8, 7) + JulianDay.timeOf(LocalTime.of(11, 0));
final GNSSLocation location = GNSSLocation.create(39.742476, -105.1786);
final double observerElevation = 1830.14;
final double ut = JulianDay.of(2003, 10, 17)
	+ JulianDay.timeOf(LocalTime.of(19, 30, 30));
final double jd = TimeHelper.universalTimeToTerrestrialTime(ut, 67.);
final double tt = JulianDay.centuryJ2000Of(jd);

		//calculate the geometric mean longitude L0 of the Sun referred to the mean equinox of the time T: L0
		final double geometricMeanLongitude = geometricMeanLongitude(tt);
		if(Math.abs(geometricMeanLongitude - 204.0182616917) > 0.0000000001)
			throw new IllegalArgumentException("geometricMeanLongitude: " + (geometricMeanLongitude - 204.0182616917));
		//calculate the nutation in longitude and obliquity
		final double[] nutationInLongitudeAndObliquity = correctionNutationInLongitudeAndObliquity(tt);
		if(Math.abs(nutationInLongitudeAndObliquity[0] - -0.00399840) > 0.00000001)
			throw new IllegalArgumentException("nutationInLongitude: " + (nutationInLongitudeAndObliquity[0] - -0.00399840));
		if(Math.abs(nutationInLongitudeAndObliquity[1] - 0.00166657) > 0.00000001)
			throw new IllegalArgumentException("nutationInObliquity: " + (nutationInLongitudeAndObliquity[1] - 0.00166657));
		final double radiusVector = radiusVector(tt);
		if(Math.abs(radiusVector - 0.9965422974) > 0.0000000001)
			throw new IllegalArgumentException("radiusVector: " + (radiusVector - 0.9965422974));
		//calculate the aberration correction
		final double aberration = correctionAberration(radiusVector);
		//calculate the apparent Sun longitude: Ltrue = L0 + C
		final double apparentGeometricLongitude = apparentGeometricLongitude(geometricMeanLongitude, nutationInLongitudeAndObliquity[0],
			aberration);
		if(Math.abs(apparentGeometricLongitude - 204.0085519281) > 0.0000000002)
			throw new IllegalArgumentException("apparentGeometricLongitude: " + (apparentGeometricLongitude - 204.0085519281));
		//calculate the obliquity of the ecliptic (the inclination of the Earth’s equator with respect to the plane at which the Sun
		//and planets appear to move across the sky): ɛ0
		final double meanEclipticObliquity = meanEclipticObliquity(tt);
		//calculate the true obliquity of the ecliptic
		final double trueEclipticObliquity = trueEclipticObliquity(meanEclipticObliquity, nutationInLongitudeAndObliquity[1]);
		if(Math.abs(trueEclipticObliquity - 23.440465) > 0.000001)
			throw new IllegalArgumentException("trueEclipticObliquity: " + (trueEclipticObliquity - 23.440465));
		final double geometricMeanLatitude = geometricMeanLatitude(tt);
		if(Math.abs(geometricMeanLatitude - 0.0001011219) > 0.0000000001)
			throw new IllegalArgumentException("geometricMeanLatitude: " + (geometricMeanLatitude - 0.0001011219));
		//calculate the geocentric Sun Right Ascension
		final double rightAscension = rightAscension(geometricMeanLatitude, trueEclipticObliquity, apparentGeometricLongitude);
		if(Math.abs(rightAscension - 202.22741) > 0.00001)
			throw new IllegalArgumentException("rightAscension: " + (rightAscension - 202.22741));
		//calculate the geocentric Sun declination
		final double declination = declination(geometricMeanLatitude, apparentGeometricLongitude, trueEclipticObliquity);
		if(Math.abs(declination - -9.31434) > 0.00001)
			throw new IllegalArgumentException("declination: " + (declination - -9.31434));
		EquatorialCoordinate coord = EquatorialCoordinate.create(rightAscension, declination);
		System.out.println(coord);

		final double meanSiderealTime = meanSiderealTime(ut);
		//calculate the apparent sidereal time at Greenwich at any given time
		final double apparentSiderealTime = apparentSiderealTime(meanSiderealTime, trueEclipticObliquity, nutationInLongitudeAndObliquity[0]);
		final double localMeanSiderealTime = localMeanSiderealTime(apparentSiderealTime, location);
		//calculate the observer local hour angle: H
		final double localHourAngle = localHourAngle(localMeanSiderealTime, rightAscension);
		if(Math.abs(localHourAngle - 11.105902) > 0.000001)
			throw new IllegalArgumentException("localHourAngle: " + (localHourAngle - 11.105902));
		final double equatorialHorizontalParallax = equatorialHorizontalParallax(radiusVector);
		final double latitude = StrictMath.toRadians(location.getLatitude());
		final double u = StrictMath.atan((1. - EARTH_FLATTENING) * StrictMath.tan(latitude));
		final double chi = StrictMath.cos(u) + (observerElevation / EARTH_EQUATORIAL_RADIUS) * StrictMath.cos(latitude);
		final double y = 0.99664719 * StrictMath.sin(u) + (observerElevation / EARTH_EQUATORIAL_RADIUS) * StrictMath.sin(latitude);
		//calculate the parallax in the sun right ascension: Δα
		final double horizontalParallax = StrictMath.toRadians(equatorialHorizontalParallax);
		final double lha = StrictMath.toRadians(localHourAngle);
		final double decl = StrictMath.toRadians(declination);
		final double deltaRightAscension = StrictMath.toDegrees(StrictMath.atan2(
			-chi * StrictMath.sin(horizontalParallax) * StrictMath.sin(lha),
			StrictMath.cos(decl) - chi * StrictMath.sin(horizontalParallax) * StrictMath.cos(lha)
		));
		//calculate the topocentric Sun Right Ascension: α'
		final double rightAscensionTopocentric = rightAscension + deltaRightAscension;
		if(Math.abs(rightAscensionTopocentric - 202.22704) > 0.00001)
			throw new IllegalArgumentException("rightAscensionTopocentric: " + (rightAscensionTopocentric - 202.22704));
		final double declinationTopocentric = StrictMath.toDegrees(StrictMath.atan2(
			((StrictMath.sin(decl) - y * StrictMath.sin(horizontalParallax)) * StrictMath.cos(StrictMath.toRadians(deltaRightAscension))),
			StrictMath.cos(decl) - chi * StrictMath.sin(horizontalParallax) * StrictMath.cos(lha)
		));
		if(Math.abs(declinationTopocentric - -9.316179) > 0.000001)
			throw new IllegalArgumentException("declinationTopocentric: " + (declinationTopocentric - -9.316179));
		//calculate the topocentric local hour angle: H’
		final double localHourAngleTopocentric = localHourAngle - deltaRightAscension;
		if(Math.abs(localHourAngleTopocentric - 11.10629) > 0.00002)
			throw new IllegalArgumentException("localHourAngleTopocentric: " + (localHourAngleTopocentric - 11.10629));
		//TODO
		//calculate the topocentric zenith angle: θ
		final double zenithTopocentric = 0;



//		EquatorialCoordinate coord2 = sunPosition(jd);

		System.out.println(coord);
	}

	/**
	 * Calculate the eccentricity of Earth's orbit, e.
	 *
	 * @param tt	Julian Century of Terrestrial Time from J2000.0.
	 * @return	The eccentricity of Earth's orbit.
	 */
	private static double earthOrbitEccentricity(final double tt){
		return eval(tt, new double[]{0.016708634, -0.000042037, -0.0000001267});
	}

	/**
	 * Calculate the distance between the center of the Sun and the center of the Earth, R.
	 * <p>U.S. Naval Observatory function.</p>
	 *
	 * @param tt	Julian Century of Terrestrial Time from J2000.0.
	 * @return	Distance between the center of the Sun and the center of the Earth [AU].
	 */
	private static double radiusVector(final double tt){
		final double jme = tt / 10.;
		final double[] parameters = new double[5];
		for(int i = 0; i < parameters.length; i ++){
			double parameter = 0.;
			final Collection<Double[]> elements = EARTH_RADIUS_VECTOR_DATA.get("R" + i);
			for(final Double[] element : elements)
				parameter += element[0] * StrictMath.cos(element[1] + element[2] * jme);
			parameters[i] = parameter;
		}
		return eval(jme, parameters) / 100_000_000.;
	}

	/**
	 * Calculate the distance between the center of the Sun and the center of the Earth, R.
	 *
	 * @param eccentricity	The eccentricity of Earth's orbit.
	 * @param trueAnomaly	The true anomaly of the Sun [°].
	 * @return	Distance between the center of the Sun and the center of the Earth [AU].
	 */
	private static double radiusVector(final double eccentricity, final double trueAnomaly){
		return 1.000001018 * (1. - eccentricity * eccentricity) / (1. + eccentricity * StrictMath.cos(StrictMath.toRadians(trueAnomaly)));
	}

	/**
	 * Calculate the equatorial horizontal parallax of the Sun, ξ.
	 *
	 * @param radiusVector	Radius vector of the Earth [AU].
	 * @return	The equatorial horizontal parallax of the Sun [°].
	 */
	private static double equatorialHorizontalParallax(double radiusVector){
		return 8.794 / (JulianDay.SECONDS_IN_HOUR * radiusVector);
	}

	/**
	 * True obliquity of the ecliptic corrected for nutation, ɛ.
	 *
	 * @param meanEclipticObliquity	Obliquity of the ecliptic, corrected for parallax [°].
	 * @param deltaEpsilon	Nutation in obliquity [°].
	 * @return	Apparent longitude of the Sun [°].
	 */
	private static double trueEclipticObliquity(final double meanEclipticObliquity, final double deltaEpsilon){
		return meanEclipticObliquity + deltaEpsilon;
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
	private static double[] correctionNutationInLongitudeAndObliquity(final double tt){
		//mean elongation of the Moon from the Sun [°]
		final double d = StrictMath.toRadians(correctRangeDegree(eval(tt, new double[]{297.8503631, 445267.1114800, -0.0019142, 1. / 189474.})));
		//mean anomaly of the Sun [°]
		final double m = StrictMath.toRadians(correctRangeDegree(eval(tt, new double[]{357.5277233, 35999.0503400, -0.0001603, - 1. / 300000.})));
		//mean anomaly of the Moon [°]
		final double mp = StrictMath.toRadians(correctRangeDegree(eval(tt, new double[]{134.9629814, 477198.8673981, 0.0086972, 1. / 56250.})));
		//Moon's argument of Latitude [°]
		final double f = StrictMath.toRadians(correctRangeDegree(eval(tt, new double[]{93.2719103, 483202.0175381, -0.00368250, 1. / 327270.})));
		//Longitude of the ascending node of the Moon's mean orbit on the ecliptic measured from the mean equinox of the date [rad]
		final double omega = StrictMath.toRadians(correctRangeDegree(eval(tt, new double[]{125.0445222, -1934.1362608, 0.002070833, 1. / 450000.})));

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
		deltaPsi = toDegrees(0, 0, deltaPsi / 10_000.) / 63.;
		//[°]
		deltaEpsilon = toDegrees(0, 0, deltaEpsilon / 10_000.) / 63.;
		return new double[]{deltaPsi, deltaEpsilon};
	}

	/**
	 * Calculate corrections of aberration (∆τ).
	 *
	 * @param earthRadiusVector	Earth radius vector [AU].
	 * @return	∆τ [°].
	 */
	private static double correctionAberration(final double earthRadiusVector){
		return -20.4898 / (JulianDay.SECONDS_IN_HOUR * earthRadiusVector);
	}

	/**
	 * Calculate the longitude of the perihelion of the orbit, π.
	 *
	 * @param tt	Julian Century of Terrestrial Time from J2000.0.
	 * @return	The longitude of the perihelion of the orbit [°].
	 */
	private static double longitudeOfEarthPerihelion(final double tt){
		return eval(tt, new double[]{102.93735, 1.71946, 0.00046});
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

	@SuppressWarnings("OverlyComplexArithmeticExpression")
	private static double equationOfTime(final double tt){
		final double e0 = meanEclipticObliquity(tt);
		final double epsilon = apparentEclipticObliquity(e0, tt);
		final double l0 = StrictMath.toRadians(geometricMeanLongitude(tt));
		final double e = earthOrbitEccentricity(tt);
		final double m = StrictMath.toRadians(geometricMeanAnomaly(tt));
		double y = StrictMath.tan(StrictMath.toRadians(epsilon) / 2.);
		y *= y;

		//[h]
		return StrictMath.toDegrees(y * StrictMath.sin(2. * l0)
			- 2. * e * StrictMath.sin(m)
			+ 4. * e * y * StrictMath.sin(m) * StrictMath.cos(2. * l0)
			- 0.5 * y * y * StrictMath.sin(4. * l0)
			- 1.25 * e * e * StrictMath.sin(2. * m)) * 4. / 60.;
	}

	/**
	 * Calculate Greenwich Mean Sidereal Time, ΘGMST.
	 *
	 * @param ut	Julian Day of Universal Time from J2000.0.
	 * @return	mean Sidereal time at Greenwich [°].
	 */
	private static double meanSiderealTime(final double ut){
		final double[] dateAndTime = JulianDay.extractDateAndTime(ut);
		final double ut0 = JulianDay.centuryJ2000Of(dateAndTime[0]);
		//[s]
		final double t = dateAndTime[1] * JulianDay.SECONDS_IN_DAY;

		//Greenwich Sidereal Time at midnight [day]
		final double h0 = eval(ut0, new double[]{24110.54841, 8640184.812866, 0.093104, -6.2e-6}) / JulianDay.SECONDS_IN_DAY;
		final double earthSiderealRotationRate = earthSiderealRotationRate(ut0);
		/*
		This is the difference between UT1 (time using the mean rotating Earth as a clock) and UTC (time that runs at the same rate as
		an atomic clock, but with leap seconds occasionally inserted to keep UT1 and UTC in sync).
		Since dUT1 = |UT1 − UTC| < 0.9 s, it is ignored.
		See <a href="https://en.wikipedia.org/wiki/DUT1">DUT1</a>.
		*/
		//[s]
		final double dUT1 = 0.;
		//[day]
		final double h = frac(h0 + earthSiderealRotationRate * (t - dUT1));
		return h * JulianDay.HOURS_IN_DAY * JulianDay.DEGREES_PER_HOUR;
		//alternative:
		//return correctRangeDegree(eval(JulianDay.centuryJ2000Of(ut), new double[]{280.46061837, 360.98564736629 * JulianDay.CIVIL_SAECULUM, 0.000387933, -1. / 38710000.}));
	}

	private static double frac(final double x){
		return x - (int)x;
	}

	/**
	 * Calculate the Earth's sidereal rotation rate.
	 *
	 * @param ut	Julian Century of Universal Time from J2000.0.
	 * @return	Earth's sidereal rotation rate [sidereal day/UT s].
	 */
	private static double earthSiderealRotationRate(final double ut){
		//[rad/s]
		final double rate = 7.2921158553e-5 + 4.3e-15 * ut;
		return rate / (2. * StrictMath.PI);
	}

	/**
	 * Calculate Greenwich Apparent Sidereal Time, ΘGAST.
	 *
	 * @param meanSiderealTime	Greenwich Mean Sidereal Time [°].
	 * @param trueEclipticObliquity	Obliquity of the ecliptic, corrected for nutation [°].
	 * @param deltaPsi	Nutation in longitude [°].
	 * @return	apparent Sidereal time at Greenwich [°].
	 */
	private static double apparentSiderealTime(final double meanSiderealTime, final double trueEclipticObliquity, final double deltaPsi){
		final double equationOfTheEquinoxes = deltaPsi * StrictMath.cos(StrictMath.toRadians(trueEclipticObliquity));
		return correctRangeDegree(meanSiderealTime + equationOfTheEquinoxes);
	}

	/**
	 * Calculate Local Mean Sidereal Time, ΘLMST.
	 *
	 * @param meanSiderealTime	Greenwich Mean Sidereal Time [°].
	 * @return	The apparent local Sidereal time at Greenwich [°].
	 */
	private static double localMeanSiderealTime(final double meanSiderealTime, final GNSSLocation location){
		return meanSiderealTime + location.getLongitude();
	}

	/**
	 * Calculate the hour angle of a body, H.
	 *
	 * @param localSiderealTime	Local Sidereal Time [°].
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
		//Ltrue = L0 + C
		final double trueGeometricLongitude = correctRangeDegree(geometricMeanLongitude + equationOfCenter);
		final double[] nutationInLongitudeAndObliquity = correctionNutationInLongitudeAndObliquity(t);
		final double aberration = correctionAberration(radiusVector(t));
		final double apparentGeometricLatitude = apparentGeometricLongitude(geometricMeanLongitude, nutationInLongitudeAndObliquity[0],
			aberration);
		final double apparentGeometricLongitude = apparentGeometricLongitude(geometricMeanLongitude, nutationInLongitudeAndObliquity[0],
			aberration);


		//		delta = m_longitude + StrictMath.toDegrees(hourAngle);
//		timeDiff = 4 * delta;
//		timeUTC = 720 - timeDiff - eqTime; // in minutes

final double meanEclipticObliquity = meanEclipticObliquity(t);
final double apparentEclipticObliquity = apparentEclipticObliquity(meanEclipticObliquity, t);
final double apparentDeclination = declination(apparentGeometricLatitude, apparentGeometricLongitude, apparentEclipticObliquity);
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
		return date.getDayOfYear() + JulianDay.timeOf(date) + dividend / 24.;
	}

	/**
	 * Computes the base longitude hour.
	 *
	 * @return	The longitude of the location of the solar event divided by 15 (deg/hour).
	 */
	private double getBaseLongitudeHour(){
		return degToHrs(location.getLongitude());
	}

	/**
	 * @param trueLong	Sun's true longitude [°].
	 * @param zenith	Zenith enumeration.
	 * @param sunrise	Whether it's sunrise or sunset.
	 * @return	The Sun local hour [h].
	 * @throws SolarEventException	Whenever the Sun never rises or sets.
	 */
	private double localHourAngle(final double trueLong, final Zenith zenith, final Boolean sunrise,
			double apparentDeclination) throws SolarEventException{
		final double cosLocalHour = cosineLocalHour(trueLong, zenith);

double latRad = StrictMath.toRadians(location.getLatitude());
double sdRad  = StrictMath.toRadians(apparentDeclination);
double cosHA = StrictMath.sin(zenith.getRadians())
	/ (StrictMath.cos(latRad) * StrictMath.cos(sdRad)) - StrictMath.tan(latRad) * StrictMath.tan(sdRad);

		if(cosLocalHour < -1.)
			//the sun never sets on this location on the specified date
			throw SolarEventException.create(SolarEventError.NEVER_SETS);
		if(cosLocalHour > 1.)
			//the sun never rises on this location on the specified date
			throw SolarEventException.create(SolarEventError.NEVER_RISES);

		final double localHour = StrictMath.toDegrees(StrictMath.acos(cosLocalHour));
//		final double localHour = StrictMath.toDegrees(StrictMath.acos(cosHA));
		return degToHrs(sunrise? 360. - localHour: localHour);
	}

	private double cosineLocalHour(final double trueLong, final Zenith zenith){
		final double sinDeclination = 0.39782 * StrictMath.sin(StrictMath.toRadians(trueLong));
		final double cosDeclination = StrictMath.cos(StrictMath.asin(sinDeclination));

		final double sinZenith = StrictMath.sin(zenith.getRadians());
		final double sinLatitude = StrictMath.sin(StrictMath.toRadians(location.getLatitude()));
		final double cosLatitude = StrictMath.cos(StrictMath.toRadians(location.getLongitude()));

		return (sinZenith - sinLatitude * sinDeclination) / (cosLatitude * cosDeclination);
	}

	/**
	 * Computes the local mean time.
	 *
	 * @param trueGeometricLongitude	Suns true longitude [°].
	 * @param longitudeHour	Longitude hour [h].
	 * @param localHourAngle	Local hour of the Sun [h].
	 * @return	The local mean time [h].
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
		final double tanL = StrictMath.tan(StrictMath.toRadians(sunTrueLongitude));
		final double rightAscension = correctRangeDegree(StrictMath.toDegrees(StrictMath.atan(0.91764 * tanL)));

		final double longitudeQuadrant = 90. * (int)(sunTrueLongitude / 90.);
		final double rightAscensionQuadrant = 90. * (int)(rightAscension / 90.);
		return (rightAscension + longitudeQuadrant - rightAscensionQuadrant) / JulianDay.DEGREES_PER_HOUR;
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

	private static double degToHrs(final double degrees){
		return degrees / JulianDay.DEGREES_PER_HOUR;
	}

	/**
	 * Use Horner's method to compute and return the polynomial evaluated at {@code x}:<br/>
	 * {@code p[0] + p[1] * x^1 + p[2] * x^2 + ... + p[n-1] * x^n-1}
	 *
	 * @param x	The value at which to calculate the polynomial.
	 * @param p	The polynomial coefficients.
	 * @return	The value of the polynomial.
	 */
	private static double eval(final double x, final double[] p){
		double result = 0.;
		for(int i = p.length - 1; i >= 0; i --)
			result = p[i] + result * x;
		return result;
	}

}
