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

	static final double EARTH_FLATTENING = 1. / 298.25642;
	//[m]
	static final double EARTH_EQUATORIAL_RADIUS = 6378140.;

	/** [°C] */
	static final double ABSOLUTE_ZERO = 273.15;


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
	//[m]
	private double observerElevation;
	//[hPa]
	private double pressure;
	//[°C]
	private double temperature;
	//[°]
	private double surfaceSlope;
	//[°]
	private double surfaceAzimuthRotation;


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

	public double getObserverElevation(){
		return observerElevation;
	}

	public SolarEventCalculator withObserverElevation(final double observerElevation){
		this.observerElevation = observerElevation;

		return this;
	}

	public double getPressure(){
		return pressure;
	}

	public SolarEventCalculator withPressure(final double pressure){
		this.pressure = pressure;

		return this;
	}

	public double getTemperature(){
		return temperature;
	}

	public SolarEventCalculator withTemperature(final double temperature){
		this.temperature = temperature;

		return this;
	}

	public double getSurfaceSlope(){
		return surfaceSlope;
	}

	/**
	 * Slope of the surface measured from the horizontal plane.
	 *
	 * @param surfaceSlope	The surface slope [°].
	 * @return	This instance.
	 */
	public SolarEventCalculator withSurfaceSlope(final double surfaceSlope){
		this.surfaceSlope = surfaceSlope;

		return this;
	}

	public double getSurfaceAzimuthRotation(){
		return surfaceAzimuthRotation;
	}

	/**
	 * Surface azimuth rotation angle, measured from south to the projection of the surface normal on the horizontal plane, positive if
	 * oriented west from south.
	 *
	 * @param surfaceAzimuthRotation	The surface azimuth rotation angle [°].
	 * @return	This instance.
	 */
	public SolarEventCalculator withSurfaceAzimuthRotation(final double surfaceAzimuthRotation){
		this.surfaceAzimuthRotation = surfaceAzimuthRotation;

		return this;
	}


	//---
	//https://squarewidget.com/solar-coordinates/

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
		final EquatorialCoordinate coord = SunPosition.sunPosition(jd);
		//[h]
		final double ra = MathHelper.degToHrs(coord.getRightAscension());
		final double decl = coord.getDeclination();

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
		final double geometricMeanLongitude = SunPosition.geometricMeanLongitude(t);
		final double meanAnomaly = geometricMeanAnomaly(t);
		final double eccentricity = earthOrbitEccentricity(t);
		final double equationOfCenter = equationOfCenter(meanAnomaly, t);
		//Ltrue = L0 + C
		final double trueGeometricLongitude = MathHelper.correctRangeDegree(geometricMeanLongitude + equationOfCenter);
		//ν = M + C
		final double trueAnomaly = MathHelper.correctRangeDegree(meanAnomaly + equationOfCenter);
		final double radiusVector = SunPosition.radiusVector(t);
		final double radiusVector2 = radiusVector(eccentricity, trueAnomaly);
		final double[] nutation = SunPosition.nutationCorrection(t);
		final double aberration = SunPosition.aberrationCorrection(SunPosition.radiusVector(t));
		final double apparentGeometricLongitude = SunPosition.apparentGeometricLongitude(geometricMeanLongitude, nutation[0], aberration);
		final double meanEclipticObliquity = SunPosition.meanEclipticObliquity(t);
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
	 * Calculate the geometric mean anomaly of the Sun, M.
	 *
	 * @param tt	Julian Century of Terrestrial Time from J2000.0.
	 * @return	The geometric mean anomaly of the Sun [°].
	 */
	private static double geometricMeanAnomaly(final double tt){
		return MathHelper.correctRangeDegree(MathHelper.eval(tt, new double[]{357.52911, 35999.05029, -0.0001537}));
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
		return MathHelper.eval(tt, new double[]{1.914602, -0.004817, -0.000014}) * StrictMath.sin(geometricMeanAnomaly)
			+ MathHelper.eval(tt, new double[]{0.019993, -0.000101}) * StrictMath.sin(2. * geometricMeanAnomaly)
			+ 0.000289 * StrictMath.sin(3. * geometricMeanAnomaly);
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

	//https://www.nrel.gov/docs/fy08osti/34302.pdf
	public static void main(String[] args){
		final GNSSLocation location = GNSSLocation.create(39.742476, -105.1786);
		//[m]
		final double observerElevation = 1830.14;
		//[hPa]
		final double pressure = 820.;
		//[°C]
		final double temperature = 11.;
		//slope of the surface measured from the horizontal plane [°]
		final double surfaceSlope = 30.;
		//surface azimuth rotation angle, measured from south to the projection of the surface normal on the horizontal plane, positive if
		//oriented west from south [°]
		final double surfaceAzimuthRotation = -10.;

		final double ut = JulianDay.of(2003, 10, 17)
			+ JulianDay.timeOf(LocalTime.of(19, 30, 30));
		TimeHelper.deltaT(2003);
		final double jd = TimeHelper.universalTimeToTerrestrialTime(ut, 67.);
		final double tt = JulianDay.centuryJ2000Of(jd);

		EquatorialCoordinate coord = SunPosition.sunPosition(tt);
		final double[] nutation = SunPosition.nutationCorrection(tt);
		final double radiusVector = SunPosition.radiusVector(tt);
		//calculate the obliquity of the ecliptic (the inclination of the Earth’s equator with respect to the plane at which the Sun
		//and planets appear to move across the sky): ɛ0
		final double meanEclipticObliquity = SunPosition.meanEclipticObliquity(tt);
		//calculate the true obliquity of the ecliptic
		final double trueEclipticObliquity = SunPosition.trueEclipticObliquity(meanEclipticObliquity, nutation[1]);

		//calculate the mean sidereal time at Greenwich at any given time: ΘGMST
		final double meanSiderealTime = meanSiderealTime(ut);
		//calculate the apparent sidereal time at Greenwich at any given time: ΘGAST
		final double apparentSiderealTime = apparentSiderealTime(meanSiderealTime, trueEclipticObliquity, nutation[0]);
		//calculate the local mean sidereal time at Greenwich at any given time: ΘLMST
		final double localMeanSiderealTime = localMeanSiderealTime(apparentSiderealTime, location);
		//calculate the observer local hour angle: H
		final double localHourAngle = localHourAngle(localMeanSiderealTime, coord.getRightAscension());
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
		final double decl = StrictMath.toRadians(coord.getDeclination());
		final double deltaRightAscension = StrictMath.toDegrees(StrictMath.atan2(
			-chi * StrictMath.sin(horizontalParallax) * StrictMath.sin(lha),
			StrictMath.cos(decl) - chi * StrictMath.sin(horizontalParallax) * StrictMath.cos(lha)
		));
		//calculate the topocentric Sun Right Ascension: α'
		final double rightAscensionTopocentric = coord.getRightAscension() + deltaRightAscension;
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
		//calculate the topocentric elevation angle without atmospheric refraction correction: e0
		final double lhaTopocentric = StrictMath.toRadians(localHourAngleTopocentric);
		final double e0 = StrictMath.toDegrees(StrictMath.asin(
			StrictMath.sin(latitude) * StrictMath.sin(StrictMath.toRadians(declinationTopocentric))
				+ StrictMath.cos(latitude) * StrictMath.cos(StrictMath.toRadians(declinationTopocentric)) * StrictMath.cos(lhaTopocentric)
		));
		final double deltaE = atmosphericRefractionCorrection(pressure, temperature, e0);
		//calculate the topocentric elevation angle: e
		final double elevationTopocentric = e0 + deltaE;
		//calculate the topocentric zenith angle: θ
		final double zenithTopocentric = 90. - elevationTopocentric;
		if(Math.abs(zenithTopocentric - 50.11162) > 0.00001)
			throw new IllegalArgumentException("zenithTopocentric: " + (zenithTopocentric - 50.11162));
		//calculate the topocentric astronomers azimuth angle (measured westward from south): Γ
		final double azimuthTopocentric = MathHelper.correctRangeDegree(StrictMath.toDegrees(StrictMath.atan2(
			StrictMath.sin(lhaTopocentric),
			StrictMath.cos(lhaTopocentric) * StrictMath.sin(latitude)
				- StrictMath.tan(StrictMath.toRadians(declinationTopocentric)) * StrictMath.cos(latitude)
		)));
		//calculate the topocentric azimuth angle (measured westward from north): M
		final double azimuthTopocentricNavigators = MathHelper.correctRangeDegree(azimuthTopocentric + 180.);
		if(Math.abs(azimuthTopocentricNavigators - 194.34024) > 0.00001)
			throw new IllegalArgumentException("azimuthTopocentricNavigators: " + (azimuthTopocentricNavigators - 194.34024));
		//calculate the incidence angle for a surface oriented in any direction, I
		final double zenith = StrictMath.toRadians(zenithTopocentric);
		final double slope = StrictMath.toRadians(surfaceSlope);
		final double incidence = StrictMath.toDegrees(StrictMath.acos(
			StrictMath.cos(zenith) * StrictMath.cos(slope)
				+ StrictMath.sin(zenith) * StrictMath.sin(slope)
				* StrictMath.cos(StrictMath.toRadians(azimuthTopocentric - surfaceAzimuthRotation))
		));
		if(Math.abs(incidence - 25.18700) > 0.00001)
			throw new IllegalArgumentException("incidence: " + (incidence - 25.18700));


		//		EquatorialCoordinate coord2 = sunPosition(jd);
		//		System.out.println(coord2);
	}

	/**
	 * Calculate the atmospheric refraction correction, Δe.
	 *
	 * @param pressure	The pressure [hPa].
	 * @param temperature	The temperature [°C].
	 * @param e0	The topocentric elevation angle without atmospheric refraction correction [°].
	 * @return	The correction [°].
	 */
	static double atmosphericRefractionCorrection(final double pressure, final double temperature, final double e0){
		return (pressure / 1010.)
			* ((ABSOLUTE_ZERO + 10.) / (ABSOLUTE_ZERO + temperature))
			* (1.02 / (60. * StrictMath.tan(StrictMath.toRadians(e0 + 10.3 / (e0 + 5.11)))));
	}

	/**
	 * Calculate the eccentricity of Earth's orbit, e.
	 *
	 * @param tt	Julian Century of Terrestrial Time from J2000.0.
	 * @return	The eccentricity of Earth's orbit.
	 */
	private static double earthOrbitEccentricity(final double tt){
		return MathHelper.eval(tt, new double[]{0.016708634, -0.000042037, -0.0000001267});
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
	static double equatorialHorizontalParallax(double radiusVector){
		return 8.794 / (JulianDay.SECONDS_IN_HOUR * radiusVector);
	}

	/**
	 * Calculate the longitude of the perihelion of the orbit, π.
	 *
	 * @param tt	Julian Century of Terrestrial Time from J2000.0.
	 * @return	The longitude of the perihelion of the orbit [°].
	 */
	private static double longitudeOfEarthPerihelion(final double tt){
		return MathHelper.eval(tt, new double[]{102.93735, 1.71946, 0.00046});
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
		final double e0 = SunPosition.meanEclipticObliquity(tt);
		final double epsilon = apparentEclipticObliquity(e0, tt);
		final double l0 = StrictMath.toRadians(SunPosition.geometricMeanLongitude(tt));
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
	static double meanSiderealTime(final double ut){
		final double[] dateAndTime = JulianDay.extractDateAndTime(ut);
		final double ut0 = JulianDay.centuryJ2000Of(dateAndTime[0]);
		//[s]
		final double t = dateAndTime[1] * JulianDay.SECONDS_IN_DAY;

		//Greenwich Sidereal Time at midnight [day]
		final double h0 = MathHelper.eval(ut0, new double[]{24110.54841, 8640184.812866, 0.093104, -6.2e-6}) / JulianDay.SECONDS_IN_DAY;
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
		final double h = MathHelper.frac(h0 + earthSiderealRotationRate * (t - dUT1));
		return h * JulianDay.HOURS_IN_DAY * JulianDay.DEGREES_PER_HOUR;
		//alternative:
		//return correctRangeDegree(eval(JulianDay.centuryJ2000Of(ut), new double[]{280.46061837, 360.98564736629 * JulianDay.CIVIL_SAECULUM, 0.000387933, -1. / 38710000.}));
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
	static double apparentSiderealTime(final double meanSiderealTime, final double trueEclipticObliquity, final double deltaPsi){
		final double equationOfTheEquinoxes = deltaPsi * StrictMath.cos(StrictMath.toRadians(trueEclipticObliquity));
		return MathHelper.correctRangeDegree(meanSiderealTime + equationOfTheEquinoxes);
	}

	/**
	 * Calculate Local Mean Sidereal Time, ΘLMST.
	 *
	 * @param meanSiderealTime	Greenwich Mean Sidereal Time [°].
	 * @return	The apparent local Sidereal time at Greenwich [°].
	 */
	static double localMeanSiderealTime(final double meanSiderealTime, final GNSSLocation location){
		return meanSiderealTime + location.getLongitude();
	}

	/**
	 * Calculate the hour angle of a body, H.
	 *
	 * @param localSiderealTime	Local Sidereal Time [°].
	 * @param rightAscension	Right ascension [°].
	 * @return	The hour angle [°].
	 */
	static double localHourAngle(final double localSiderealTime, final double rightAscension){
		return localSiderealTime - rightAscension;
	}


	private LocalTime computeSolarEventTime(final LocalDateTime date, final Zenith solarZenith, final boolean sunrise)
		throws SolarEventException{
		final double longitudeHour = getLongitudeHour(date, sunrise);

		final double jd = JulianDay.of(date);
		final double t = JulianDay.centuryJ2000Of(jd);
		final double equationOfTime = equationOfTime(t);
		final double geometricMeanLongitude = SunPosition.geometricMeanLongitude(t);
		final double meanAnomaly = geometricMeanAnomaly(t);
		final double equationOfCenter = equationOfCenter(meanAnomaly, t);
		//Ltrue = L0 + C
		final double trueGeometricLongitude = MathHelper.correctRangeDegree(geometricMeanLongitude + equationOfCenter);
		final double[] nutation = SunPosition.nutationCorrection(t);
		final double aberration = SunPosition.aberrationCorrection(SunPosition.radiusVector(t));
		final double apparentGeometricLatitude = SunPosition.apparentGeometricLongitude(geometricMeanLongitude, nutation[0],
			aberration);
		final double apparentGeometricLongitude = SunPosition.apparentGeometricLongitude(geometricMeanLongitude, nutation[0],
			aberration);


//		delta = m_longitude + StrictMath.toDegrees(hourAngle);
//		timeDiff = 4 * delta;
//		timeUTC = 720 - timeDiff - eqTime; // in minutes

		final double meanEclipticObliquity = SunPosition.meanEclipticObliquity(t);
		final double apparentEclipticObliquity = apparentEclipticObliquity(meanEclipticObliquity, t);
		EquatorialCoordinate coord = SunPosition.toEquatorialCoordinate(apparentGeometricLatitude, apparentGeometricLongitude,
			apparentEclipticObliquity);
		final double apparentDeclination = coord.getDeclination();
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
		return MathHelper.degToHrs(location.getLongitude());
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
		return MathHelper.degToHrs(sunrise? 360. - localHour: localHour);
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
		return MathHelper.correctRangeHour(localHourAngle + rightAscension - 0.06571 * longitudeHour - 6.622);
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
		final double rightAscension = MathHelper.correctRangeDegree(StrictMath.toDegrees(StrictMath.atan(0.91764 * tanL)));

		final double longitudeQuadrant = 90. * (int)(sunTrueLongitude / 90.);
		final double rightAscensionQuadrant = 90. * (int)(rightAscension / 90.);
		return (rightAscension + longitudeQuadrant - rightAscensionQuadrant) / JulianDay.DEGREES_PER_HOUR;
	}

	@SuppressWarnings("NumericCastThatLosesPrecision")
	private LocalTime getLocalTime(final double localMeanTime){
		//adjust back to UTC
		final double utcTime = MathHelper.correctRangeHour(localMeanTime - getBaseLongitudeHour());
		return LocalTime.of(0, 0)
			.plus((long)(utcTime * 60. * 60.), ChronoUnit.SECONDS);
	}

}
