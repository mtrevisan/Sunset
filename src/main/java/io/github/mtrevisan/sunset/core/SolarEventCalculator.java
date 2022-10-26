/*
 * Copyright (c) 2022-2022 Mauro Trevisan
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
package io.github.mtrevisan.sunset.core;

import io.github.mtrevisan.sunset.AtmosphericModel;
import io.github.mtrevisan.sunset.JulianDay;
import io.github.mtrevisan.sunset.MathHelper;
import io.github.mtrevisan.sunset.SolarEventError;
import io.github.mtrevisan.sunset.SolarEventException;
import io.github.mtrevisan.sunset.StringHelper;
import io.github.mtrevisan.sunset.TimeHelper;
import io.github.mtrevisan.sunset.Zenith;
import io.github.mtrevisan.sunset.coordinates.EclipticCoordinate;
import io.github.mtrevisan.sunset.coordinates.EquatorialCoordinate;
import io.github.mtrevisan.sunset.coordinates.GNSSLocation;

import java.time.LocalDate;
import java.time.LocalDateTime;
import java.time.LocalTime;
import java.time.temporal.ChronoUnit;


public class SolarEventCalculator{

	private static final double EARTH_FLATTENING = 1. / 298.25642;
	//[m]
	private static final double EARTH_EQUATORIAL_RADIUS = 6378140.;

	private static final double[] SUN_GEOCENTRIC_MEAN_ANOMALY = {357.52911, 35_999.050_29, -0.000_1537};


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

	//TODO
	//https://github.com/MenoData/Time4J/blob/master/base/src/main/java/net/time4j/calendar/astro/SunPosition.java
	/**
	 * Computes the sunset time.
	 *
	 * @param solarZenith	Enumeration corresponding to the type of sunset to compute.
	 * @param date	Date to compute the sunset for.
	 * @return	The sunset time or {@code null} for no sunset.
	 * @throws SolarEventException   Whenever the Sun never rises or sets.
	 */
	public final LocalDateTime sunset(final LocalDate date, final Zenith solarZenith) throws SolarEventException{
		final double ut = JulianDay.of(date);
		final double dt = TimeHelper.deltaT(date.getYear());
		double jd = TimeHelper.universalTimeToTerrestrialTime(ut, dt);
		double tt = JulianDay.centuryJ2000Of(jd);

//		//calculate geocentric mean ecliptic coordinate
//		final EclipticCoordinate eclipticCoord = SunPosition.sunEclipticPosition(jd);
//		final EquatorialCoordinate equatorialCoord = SunPosition.sunEquatorialPosition(eclipticCoord, jd);
//		final double[] nutation = SunPosition.nutationCorrection(tt);
//		final double meanEclipticObliquity = SunPosition.meanEclipticObliquity(tt);
//		final double trueEclipticObliquity = SunPosition.trueEclipticObliquity(meanEclipticObliquity, nutation[1]);
//
//		final double meanSiderealTime = TimeHelper.meanSiderealTime(ut);
//		final double apparentSiderealTime = TimeHelper.apparentSiderealTime(meanSiderealTime, trueEclipticObliquity, nutation[0]);
//		final double localMeanSiderealTime = TimeHelper.localMeanSiderealTime(apparentSiderealTime, location);
//		final double rightAscension = equatorialCoord.getRightAscension();
//		final double localHourAngle = TimeHelper.localHourAngle(localMeanSiderealTime, rightAscension);
//
//		//compute the sun position (right ascension and declination) with respect to the observer local position at the Earth surface:
//		final double equatorialHorizontalParallax = equatorialHorizontalParallax(eclipticCoord.getDistance());
//		final double latitude = location.getLatitude();
//		final double u = StrictMath.atan((1. - SunPosition.EARTH_FLATTENING) * StrictMath.tan(latitude));
//		final double height = location.getAltitude() / SunPosition.EARTH_EQUATORIAL_RADIUS;
//		final double x = StrictMath.cos(u) + height * StrictMath.cos(latitude);
//		final double y = 0.99664719 * StrictMath.sin(u) + height * StrictMath.sin(latitude);
//		final double declination = equatorialCoord.getDeclination();
//		final double deltaRightAscension = StrictMath.atan2(
//			-x * StrictMath.sin(equatorialHorizontalParallax) * StrictMath.sin(localHourAngle),
//			StrictMath.cos(declination) - x * StrictMath.sin(equatorialHorizontalParallax)
//				* StrictMath.cos(localHourAngle)
//		);
//		//calculate the topocentric Sun Right Ascension: α'
//		final double topocentricRightAscension = rightAscension + deltaRightAscension;
//		//calculate the topocentric Sun declination: δ'
//		final double topocentricDeclination = StrictMath.atan2(
//			(StrictMath.sin(declination) - y * StrictMath.sin(equatorialHorizontalParallax)) * StrictMath.cos(deltaRightAscension),
//			StrictMath.cos(declination) - y * StrictMath.sin(equatorialHorizontalParallax) * StrictMath.cos(localHourAngle)
//		);
//		//calculate the topocentric local hour angle: H’
//		final double topocentricLocalHourAngle = localHourAngle - deltaRightAscension;
//		//calculate the true elevation angle without atmospheric refraction correction: e0
//		final double trueElevation = StrictMath.asin(
//			StrictMath.sin(latitude) * StrictMath.sin(topocentricDeclination)
//				+ StrictMath.cos(latitude) * StrictMath.cos(topocentricDeclination) * StrictMath.cos(topocentricLocalHourAngle)
//		);
//		//calculate the atmospheric refraction correction: Δe
//		final double deltaElevation = AtmosphereHelper.atmosphericRefractionCorrection(pressure, temperature, trueElevation);
//		//calculate the topocentric elevation angle: e
//		final double topocentricElevation = trueElevation + deltaElevation;
//		//calculate the topocentric zenith angle: θ
//		final double topocentricZenith = StrictMath.PI / 2. - topocentricElevation;
//		//calculate the topocentric astronomers azimuth angle (measured westward from south): Γ
//		final double topocentricAzimuth = StrictMath.atan2(
//			StrictMath.sin(topocentricLocalHourAngle),
//			StrictMath.cos(topocentricLocalHourAngle) * StrictMath.sin(latitude)
//				- StrictMath.tan(topocentricDeclination) * StrictMath.cos(latitude)
//		);
//		//calculate the (navigators) topocentric azimuth angle (measured westward from north): M
//		final double topocentricAzimuthNavigators = MathHelper.mod2pi(topocentricAzimuth + StrictMath.PI);
//		HorizontalCoordinate.create(topocentricAzimuth, topocentricElevation, eclipticCoord.getDistance());

//		//[h]
//		final double ra = MathHelper.degToHrs(coord.getRightAscension());
//		final double decl = coord.getDeclination();

//		final double geocentricMeanLongitude = SunPosition.geocentricMeanLongitude(t);
//		final double meanAnomaly = geocentricMeanAnomaly(t);
//		final double eccentricity = earthOrbitEccentricity(t);
//		final double equationOfCenter = equationOfCenter(meanAnomaly, t);
//		//Ltrue = L + C
//		final double trueGeocentricLongitude = MathHelper.limitRangeDegree(geocentricMeanLongitude + equationOfCenter);
//		//ν = M + C
//		final double trueAnomaly = MathHelper.limitRangeDegree(meanAnomaly + equationOfCenter);
//		final double radiusVector = SunPosition.radiusVector(t);
//		final double radiusVector2 = radiusVector(eccentricity, trueAnomaly);
//		final double aberration = SunPosition.aberrationCorrection(SunPosition.radiusVector(t));
//		final double apparentGeocentricLongitude = SunPosition.apparentGeocentricLongitude(geocentricMeanLongitude, nutation[0], aberration);
//		final double longitudeOfEarthPerihelion = longitudeOfEarthPerihelion(t);

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
	 * Calculate the geocentric mean anomaly of the Sun, M.
	 *
	 * @param tt	Julian Century of Terrestrial Time from J2000.0.
	 * @return	The geocentric mean anomaly of the Sun [rad].
	 */
	private static double geocentricMeanAnomaly(final double tt){
		return MathHelper.mod2pi(StrictMath.toRadians(MathHelper.eval(tt, SUN_GEOCENTRIC_MEAN_ANOMALY)));
	}

	/**
	 * Calculate the Sun's equation of center, C.
	 *
	 * @param geocentricMeanAnomaly	The mean anomaly of the Sun [rad].
	 * @param tt	Julian Century of Terrestrial Time from J2000.0.
	 * @return	The Sun's equation of center [rad].
	 */
	private static double equationOfCenter(double geocentricMeanAnomaly, final double tt){
		return MathHelper.eval(tt, new double[]{1.914602, -0.004817, -0.000014}) * StrictMath.sin(geocentricMeanAnomaly)
			+ MathHelper.eval(tt, new double[]{0.019993, -0.000101}) * StrictMath.sin(2. * geocentricMeanAnomaly)
			+ 0.000289 * StrictMath.sin(3. * geocentricMeanAnomaly);
	}

	/**
	 * Calculate the mean obliquity of the ecliptic, corrected for parallax, ɛ'.
	 *
	 * @param meanEclipticObliquity	Mean obliquity of the ecliptic [rad].
	 * @param tt	Julian Century of Terrestrial Time from J2000.0.
	 * @return	Apparent longitude of the Sun [rad].
	 */
	private static double apparentEclipticObliquity(final double meanEclipticObliquity, final double tt){
		return meanEclipticObliquity + 0.00256 * StrictMath.cos(omega(tt));
	}

	private static double nutationAndAberrationCorrection(final double tt){
		return 0.00569 + 0.00478 * StrictMath.sin(omega(tt));
	}

	private static double omega(final double tt){
		return 125.04 - 1934.136 * tt;
	}

	//WORKS
	//https://stellafane.org/misc/equinox.html
	//<a href="https://www.agopax.it/Libri_astronomia/pdf/Astronomical%20Algorithms.pdf">Meeus, Jean. Astronomical algorithms. 2nd ed. 1998.</a>
	private static void meeus(){
		//2022: 21 Dec 2022 21:47:33 GMT
		//1487: 13 Dec 1487 00:20:31 GMT
		final int year = 2022;

		//calculate an initial guess
		final double jde0 = MathHelper.eval((year - 2000) / 1000., new double[]{2451900.05952, 365242.74049, -0.06223, -0.00823, 0.00032});
		final double t = JulianDay.centuryJ2000Of(jde0);
		final double w = StrictMath.toRadians(35_999.373 * t - 2.47);
		final double deltaLambda = 1. + 0.033_4 * StrictMath.cos(w) + 0.000_7 * StrictMath.cos(2. * w);

		final double[] a = {0.00485, 0.00203, 0.00199, 0.00182, 0.00156, 0.00136, 0.00077, 0.00074, 0.00070, 0.00058, 0.00052, 0.00050, 0.00045, 0.00044, 0.00029, 0.00018, 0.00017, 0.00016, 0.00014,
			0.00012, 0.00012, 0.00012, 0.00009, 0.00008};
		final double[] b = {324.96, 337.23, 342.08, 27.85, 73.14, 171.52, 222.54, 296.72, 243.58, 119.81, 297.17, 21.02,
			247.54, 325.15, 60.93, 155.12, 288.79, 198.04, 199.76, 95.39, 287.11, 320.81, 227.73, 15.45};
		final double[] c = {1934.136, 32964.467, 20.186, 445267.112, 45036.886, 22518.443,
			65928.934, 3034.906, 9037.513, 33718.147, 150.678, 2281.226,
			29929.562, 31555.956, 4443.417, 67555.328, 4562.452, 62894.029,
			31436.921, 14577.848, 31931.756, 34777.259, 1222.114, 16859.074};
		double s = 0.;
		for(int i = 0; i < a.length; i ++)
			s += a[i] * StrictMath.cos(StrictMath.toRadians(b[i] + c[i] * t));

		final double tdt = jde0 + s / deltaLambda;

		final double deltaT = deltaT(year);
		//NOTE: the error is less than one minute for the years 1951-2050
		final double utc = TimeHelper.terrestrialTimeToUniversalTime(tdt, deltaT);

		System.out.println("deltaT: " + deltaT);
		System.out.println("tdt: " + JulianDay.dateTimeOf(tdt));
		System.out.println("utc: " + JulianDay.dateTimeOf(utc));
	}

	private static double deltaT(final int year){
		//from Meeus Astronomical Algorithms Chapter 10
		//Correction lookup table has entry for every even year between tableFirstYear and tableLastYear

		//range of years in lookup table:
		final int tableFirstYear = 1620;
		final int tableLastYear = 2002;
		//corrections [s]
		final double[] table = {
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
			50.5, 52.5, 53.8, 54.9, 55.8, 56.9, 58.3, 60.0, 61.6, 63.0, 63.8, 64.3};
		//values for Delta T for 2000 thru 2002 from NASA
		//deltaT = TDT - UTC (in Seconds)
		double deltaT;
		//centuries from the epoch 2000.0
		final double t = (year - 2000) / 100.;

		//find correction in table
		if(year >= tableFirstYear && year <= tableLastYear){
			//odd year - interpolate
			if((year % 2) != 0)
				deltaT = (table[(year - tableFirstYear - 1) / 2] + table[(year - tableFirstYear + 1) / 2]) / 2.;
			//even year - direct table lookup
			else
				deltaT = table[(year - tableFirstYear) / 2];
		}
		else if(year < 948)
			deltaT = 2177. + 497. * t + 44.1 * StrictMath.pow(t, 2.);
		else{
			deltaT = 102. + 102. * t + 25.3 * StrictMath.pow(t, 2);
			//special correction to avoid discontinuity in 2000
			if(year >= 2000 && year <= 2100)
				deltaT += 0.37 * (year - 2100);
		}
		return deltaT;
	}


	//https://www.nrel.gov/docs/fy08osti/34302.pdf
	//https://midcdmz.nrel.gov/solpos/spa.html
	//http://phpsciencelabs.us/wiki_programs/Sidereal_Time_Calculator.php
	//https://lweb.cfa.harvard.edu/~jzhao/times.html
	//https://www.meteopiateda.it/busteggia/pages/astronomy/equisol.php
	//https://www.researchgate.net/publication/262200491_Calcolo_analitico_della_posizione_del_sole_per_l%27allineamento_di_impianti_solari_ed_altre_applicazioni
	//https://www.suncalc.org/#/45.7149,12.1941,17/2022.06.27/14:23/1/3
	public static void main(String[] args) throws SolarEventException{
		meeus();
		final GNSSLocation location = GNSSLocation.create(
			MathHelper.toDegrees(45, 42, 54.),
			MathHelper.toDegrees(12, 11, 37.),
			40.
		);
		AtmosphericModel atmosphericModel = AtmosphericModel.create(1013.25, 25.);
//		final GNSSLocation location = GNSSLocation.create(39.742476, -105.1786, 1830.14);
//		AtmosphericModel atmosphericModel = AtmosphericModel.create(1021.5, 8.);

		final double ut = JulianDay.of(2022, 12, 22)
			+ JulianDay.timeOf(LocalTime.of(19, 30, 30));
		//[s]
		final double deltaT = deltaT(2022);
		final double jd = TimeHelper.universalTimeToTerrestrialTime(ut, deltaT);
		final double tt = JulianDay.centuryJ2000Of(jd);

		final EclipticCoordinate eclipticCoordBefore = SunPosition.sunEclipticPosition(jd - 1.);
		EquatorialCoordinate coordBefore = SunPosition.sunEquatorialPosition(eclipticCoordBefore, jd - 1.);
		final EclipticCoordinate eclipticCoord = SunPosition.sunEclipticPosition(jd);
		EquatorialCoordinate coord = SunPosition.sunEquatorialPosition(eclipticCoord, jd);
		final EclipticCoordinate eclipticCoordAfter = SunPosition.sunEclipticPosition(jd + 1.);
		EquatorialCoordinate coordAfter = SunPosition.sunEquatorialPosition(eclipticCoordAfter, jd + 1.);
		final double[] nutation = SunPosition.nutationCorrection(tt);
		//calculate the obliquity of the ecliptic (the inclination of the Earth’s equator with respect to the plane at which the Sun
		//and planets appear to move across the sky): ɛ0
		final double meanEclipticObliquity = SunPosition.meanEclipticObliquity(tt);
		//calculate the true obliquity of the ecliptic
		final double trueEclipticObliquity = SunPosition.trueEclipticObliquity(meanEclipticObliquity, nutation[1]);

		//calculate the mean sidereal time at Greenwich at any given time: ΘGMST
		final double meanSiderealTime = meanSiderealTime(ut);
		//calculate the apparent sidereal time at Greenwich at any given time: ΘGAST
		final double apparentSiderealTime = apparentSiderealTime(meanSiderealTime, trueEclipticObliquity, nutation[0]);

		//calculate the approximate sun transit time: m0 [day]
		double m0 = (coord.getRightAscension() - location.getLongitude() - apparentSiderealTime) / (2. * StrictMath.PI);
		//calculate the local hour angle corresponding to the sun elevation equals -0.8333°: H0
		final double latitude = StrictMath.toRadians(location.getLatitude());
		final double declination = StrictMath.toRadians(coord.getDeclination());
		final double cosH0 = (StrictMath.sin(StrictMath.toRadians(-0.8333)) - StrictMath.sin(latitude) * StrictMath.sin(declination))
			/ (StrictMath.cos(latitude) * StrictMath.cos(declination));
		if(cosH0 < -1.)
			//the sun never sets on this location on the specified date
			throw SolarEventException.create(SolarEventError.NEVER_SETS);
		if(cosH0 > 1.)
			//the sun never rises on this location on the specified date
			throw SolarEventException.create(SolarEventError.NEVER_RISES);
		final double h0 = MathHelper.modpi(StrictMath.acos(cosH0));
		//calculate the approximate sunrise time: m1
		final double m1 = MathHelper.limitRangeDay(m0 - h0 / (2. * StrictMath.PI));
		//calculate the approximate sunset time: m2
		final double m2 = MathHelper.limitRangeDay(m0 + h0 / (2. * StrictMath.PI));
		m0 = MathHelper.limitRangeDay(m0);

		//calculate the sidereal time at Greenwich for the sun transit, sunrise, and sunset
		final double v0 = apparentSiderealTime + 360.985647 * m0;
		final double v1 = apparentSiderealTime + 360.985647 * m1;
		final double v2 = apparentSiderealTime + 360.985647 * m2;
		final double n0 = m0 + deltaT / 86400.;
		final double n1 = m1 + deltaT / 86400.;
		final double n2 = m2 + deltaT / 86400.;
		//calculate the Right Ascension and declination
		double a = coord.getRightAscension() - coordBefore.getRightAscension();
		double a_prime = coord.getDeclination() - coordBefore.getDeclination();
		double b = coordAfter.getRightAscension() - coord.getRightAscension();
		double b_prime = coordAfter.getDeclination() - coord.getDeclination();
		double c = b - a;
		double c_prime = b_prime - a_prime;
		if(StrictMath.abs(a) > 2.)
			a = MathHelper.limitRangeDay(a);
		if(StrictMath.abs(a_prime) > 2.)
			a_prime = MathHelper.limitRangeDay(a_prime);
		if(StrictMath.abs(b) > 2.)
			b = MathHelper.limitRangeDay(b);
		if(StrictMath.abs(b_prime) > 2.)
			b_prime = MathHelper.limitRangeDay(b_prime);
		final double alpha0 = coord.getRightAscension() + n0 * (a + b + c * n0) / 2.;
		final double alpha1 = coord.getRightAscension() + n1 * (a + b + c * n1) / 2.;
		final double alpha2 = coord.getRightAscension() + n1 * (a + b + c * n2) / 2.;
		final double delta0 = StrictMath.toRadians(coord.getDeclination() + n0 * (a_prime + b_prime + c_prime * n0) / 2.);
		final double delta1 = StrictMath.toRadians(coord.getDeclination() + n1 * (a_prime + b_prime + c_prime * n1) / 2.);
		final double delta2 = StrictMath.toRadians(coord.getDeclination() + n1 * (a_prime + b_prime + c_prime * n2) / 2.);
		//calculate the local hour angle (measured as positive westward from the meridian): H’
		double h_prime0 = MathHelper.modpi(StrictMath.toRadians(MathHelper.frac(v0 + location.getLongitude() - alpha0) * 2. * StrictMath.PI));
		double h_prime1 = MathHelper.modpi(StrictMath.toRadians(MathHelper.frac(v1 + location.getLongitude() - alpha1) * 2. * StrictMath.PI));
		double h_prime2 = MathHelper.modpi(StrictMath.toRadians(MathHelper.frac(v2 + location.getLongitude() - alpha2) * 2. * StrictMath.PI));
		//calculate the sun altitude
		final double hh0 = StrictMath.asin(latitude) * StrictMath.sin(delta0)
			+ StrictMath.cos(latitude) * StrictMath.cos(delta0) * StrictMath.cos(h_prime0);
		final double hh1 = StrictMath.asin(latitude) * StrictMath.sin(delta1)
			+ StrictMath.cos(latitude) * StrictMath.cos(delta1) * StrictMath.cos(h_prime1);
		final double hh2 = StrictMath.asin(latitude) * StrictMath.sin(delta2)
			+ StrictMath.cos(latitude) * StrictMath.cos(delta2) * StrictMath.cos(h_prime2);
		//calculate the sun transit [UT day]
		final double transit = m0 - h_prime0 / (2. * StrictMath.PI);
		//calculate the sunrise [UT day]
		final double sunrise = m1 + (hh1 - h_prime0) / (2. * StrictMath.PI * StrictMath.cos(delta1) * StrictMath.cos(latitude) * StrictMath.sin(h_prime1));
		//calculate the sunset [UT day]
		final double sunset = m2 + (hh2 - h_prime1) / (2. * StrictMath.PI * StrictMath.cos(delta2) * StrictMath.cos(latitude) * StrictMath.sin(h_prime2));
		System.out.println(StringHelper.degreeToHMSString(sunrise * 2. * StrictMath.PI, 0));
		System.out.println(StringHelper.degreeToHMSString(transit * 2. * StrictMath.PI, 0));
		System.out.println(StringHelper.degreeToHMSString(sunset * 2. * StrictMath.PI, 0));

/*
Topocentric zenith angle	116.01759
Top. azimuth angle (eastward from N)	279.725959
Top. azimuth angle (westward from S)	99.725959
Surface incidence angle	116.01759
Local sunrise time	6.212067
Local sun transit time	11.768045
Local sunset time	17.338667
Earth heliocentric longitude	24.307715
Earth heliocentric latitude	-0.000106
Earth radius vector	0.996462
Geocentric longitude	204.307715
Geocentric latitude	0.000106
Mean elongation (moon-sun)	17189.41681
Mean anomaly (sun)	1723.180685
Mean anomaly (moon)	18237.88633
Argument latitude (moon)	18423.92957
Ascending longitude (moon)	51.671506
Aberration correction	-0.005712
Apparent sun longitude	204.29801
Observer hour angle	116.120588
Sun equatorial horizontal parallax	0.002451
Sun right ascension parallax	-0.001718
Topocentric sun declination	-9.422223
Topocentric sun right ascension	202.498489
Topocentric local hour angle	116.122306
Top. elevation angle (uncorrected)	-26.01759
Atmospheric refraction correction	0
Top. elevation angle (corrected)	-26.01759
Equation of time	14.700254
Sunrise hour angle	-83.496338
Sunset hour angle	83.524274
*/
	}

	public static void main3(String[] args){
		final GNSSLocation location = GNSSLocation.create(39.742476, -105.1786, 1830.14);
		final AtmosphericModel atmosphericModel = AtmosphericModel.create(820., 11.);

		final double ut = JulianDay.of(2003, 10, 17)
			+ JulianDay.timeOf(LocalTime.of(19, 30, 30));
		final double jd = TimeHelper.universalTimeToTerrestrialTime(ut, 67.);
		final double tt = JulianDay.centuryJ2000Of(jd);

		final EclipticCoordinate eclipticCoord = SunPosition.sunEclipticPosition(jd);
		EquatorialCoordinate coord = SunPosition.sunEquatorialPosition(eclipticCoord, jd);
		final double[] nutation = SunPosition.nutationCorrection(tt);
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

		//calculate the parallax in the sun right ascension: Δα
		final double equatorialHorizontalParallax = equatorialHorizontalParallax(eclipticCoord.getDistance());
		final double latitude = StrictMath.toRadians(location.getLatitude());
		final double u = StrictMath.atan((1. - EARTH_FLATTENING) * StrictMath.tan(latitude));
		final double chi = StrictMath.cos(u) + (location.getAltitude() / EARTH_EQUATORIAL_RADIUS) * StrictMath.cos(latitude);
		final double y = 0.99664719 * StrictMath.sin(u) + (location.getAltitude() / EARTH_EQUATORIAL_RADIUS) * StrictMath.sin(latitude);
		final double deltaRightAscension = StrictMath.atan2(
			-chi * StrictMath.sin(equatorialHorizontalParallax) * StrictMath.sin(localHourAngle),
			StrictMath.cos(coord.getDeclination()) - chi * StrictMath.sin(equatorialHorizontalParallax) * StrictMath.cos(localHourAngle)
		);
		//calculate the topocentric Sun Right Ascension: α'
		final double rightAscensionTopocentric = coord.getRightAscension() + deltaRightAscension;
		if(Math.abs(rightAscensionTopocentric - StrictMath.toRadians(202.22704)) > 0.00001)
			throw new IllegalArgumentException("rightAscensionTopocentric: " + rightAscensionTopocentric);
		final double declinationTopocentric = StrictMath.atan2(
			((StrictMath.sin(coord.getDeclination()) - y * StrictMath.sin(equatorialHorizontalParallax))
				* StrictMath.cos(deltaRightAscension)),
			StrictMath.cos(coord.getDeclination()) - chi * StrictMath.sin(equatorialHorizontalParallax) * StrictMath.cos(localHourAngle)
		);
		if(Math.abs(declinationTopocentric - StrictMath.toRadians(-9.316179)) > 0.000001)
			throw new IllegalArgumentException("declinationTopocentric: " + declinationTopocentric);
		//calculate the topocentric local hour angle: H’
		final double localHourAngleTopocentric = localHourAngle - deltaRightAscension;
		if(Math.abs(localHourAngleTopocentric - 11.10629) > 0.00002)
			throw new IllegalArgumentException("localHourAngleTopocentric: " + (localHourAngleTopocentric - 11.10629));
		//calculate the topocentric elevation angle without atmospheric refraction correction: e0
		final double e0 = StrictMath.asin(
			StrictMath.sin(latitude) * StrictMath.sin(declinationTopocentric)
				+ StrictMath.cos(latitude) * StrictMath.cos(declinationTopocentric) * StrictMath.cos(localHourAngleTopocentric)
		);
		final double deltaE = atmosphericModel.atmosphericRefractionCorrection(e0);
		//calculate the topocentric elevation angle: e
		final double elevationTopocentric = e0 + deltaE;
		//calculate the topocentric zenith angle: θ
		final double zenithTopocentric = 90. - elevationTopocentric;
		if(Math.abs(zenithTopocentric - 116.01759) > 0.00001)
			throw new IllegalArgumentException("zenithTopocentric: " + (zenithTopocentric - 116.01759));
		//calculate the topocentric astronomers azimuth angle (measured westward from south): Γ
		final double azimuthTopocentric = MathHelper.mod2pi(StrictMath.atan2(
			StrictMath.sin(localHourAngleTopocentric),
			StrictMath.cos(localHourAngleTopocentric) * StrictMath.sin(latitude)
				- StrictMath.tan(declinationTopocentric) * StrictMath.cos(latitude)
		));
		//calculate the topocentric azimuth angle (measured westward from north): M
		final double azimuthTopocentricNavigators = MathHelper.mod2pi(azimuthTopocentric + StrictMath.PI);
		if(Math.abs(azimuthTopocentricNavigators - 279.725959) > 0.000001)
			throw new IllegalArgumentException("azimuthTopocentricNavigators: " + (azimuthTopocentricNavigators - 279.725959));


//		EquatorialCoordinate coord2 = sunPosition(jd);
//		System.out.println(coord2);
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
	 * Calculate the equatorial horizontal parallax of the Sun, ξ.
	 *
	 * @param radiusVector	Radius vector of the Earth [AU].
	 * @return	The equatorial horizontal parallax of the Sun [rad].
	 */
	static double equatorialHorizontalParallax(double radiusVector){
		return StrictMath.toRadians(8.794 / (JulianDay.SECONDS_IN_HOUR * radiusVector));
	}

	/**
	 * Calculate the longitude of the perihelion of the orbit, π.
	 *
	 * @param tt	Julian Century of Terrestrial Time from J2000.0.
	 * @return	The longitude of the perihelion of the orbit [rad].
	 */
	private static double longitudeOfEarthPerihelion(final double tt){
		return MathHelper.eval(tt, new double[]{102.93735, 1.71946, 0.00046});
	}

	/**
	 * Calculate the Sun's radius angle, corrected for distance, but not for refraction.
	 *
	 * @param radius	Sun radius [^].
	 * @param distance	Sun distance [AU].
	 * @return	The longitude of the perihelion of the orbit [rad].
	 */
	private static double radiusAngle(final double radius, final double distance){
		return StrictMath.asin(radius / (radius + distance));
	}

	/**
	 * Calculate the Equation of Time.
	 *
	 * @param eclipticCoord	Mean ecliptic coordinate of the Sun.
	 * @param tt	Julian Century of Terrestrial Time from J2000.0.
	 * @return	The Equation of Time [h].
	 */
	@SuppressWarnings("OverlyComplexArithmeticExpression")
	private static double equationOfTime(final EclipticCoordinate eclipticCoord, final double tt){
		final double e0 = SunPosition.meanEclipticObliquity(tt);
		final double epsilon = apparentEclipticObliquity(e0, tt);
		final double l0 = eclipticCoord.getLongitude();
		final double e = earthOrbitEccentricity(tt);
		final double m = geocentricMeanAnomaly(tt);
		double y = StrictMath.tan(epsilon / 2.);
		y *= y;

		return y * StrictMath.sin(2. * l0)
			- 2. * e * StrictMath.sin(m)
			+ 4. * e * y * StrictMath.sin(m) * StrictMath.cos(2. * l0)
			- 0.5 * y * y * StrictMath.sin(4. * l0)
			- 1.25 * e * e * StrictMath.sin(2. * m) * 4. / 60.;
	}

	/**
	 * Calculate Greenwich Mean Sidereal Time, ΘGMST.
	 *
	 * @param ut	Julian Day of Universal Time from J2000.0.
	 * @return	mean Sidereal time at Greenwich [rad].
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
		return StrictMath.toRadians(h * JulianDay.HOURS_IN_DAY * JulianDay.DEGREES_PER_HOUR);
		//alternative:
		//return MathHelper.mod2pi(StrictMath.toRadians(
		//	MathHelper.eval(JulianDay.centuryJ2000Of(ut), new double[]{280.46061837, 360.98564736629 * JulianDay.CIVIL_SAECULUM, 0.000387933, -1. / 38710000.})
		//));
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
	 * @param meanSiderealTime	Greenwich Mean Sidereal Time [rad].
	 * @param trueEclipticObliquity	Obliquity of the ecliptic, corrected for nutation [rad].
	 * @param deltaPsi	Nutation in longitude [rad].
	 * @return	apparent Sidereal time at Greenwich [rad].
	 */
	static double apparentSiderealTime(final double meanSiderealTime, final double trueEclipticObliquity, final double deltaPsi){
		final double equationOfTheEquinoxes = deltaPsi * StrictMath.cos(trueEclipticObliquity);
		return MathHelper.mod2pi(meanSiderealTime + equationOfTheEquinoxes);
	}

	/**
	 * Calculate Local Mean Sidereal Time, ΘLMST.
	 *
	 * @param meanSiderealTime	Greenwich Mean Sidereal Time [rad].
	 * @return	The apparent local Sidereal time at Greenwich [rad].
	 */
	static double localMeanSiderealTime(final double meanSiderealTime, final GNSSLocation location){
		return meanSiderealTime + location.getLongitude();
	}

	/**
	 * Calculate the hour angle of a body, H.
	 *
	 * @param localSiderealTime	Local Sidereal Time [rad].
	 * @param rightAscension	Right ascension [rad].
	 * @return	The hour angle [rad].
	 */
	static double localHourAngle(final double localSiderealTime, final double rightAscension){
		return localSiderealTime - rightAscension;
	}


	private LocalTime computeSolarEventTime(final LocalDateTime date, final Zenith solarZenith, final boolean sunrise)
			throws SolarEventException{
		final double longitudeHour = getLongitudeHour(date, sunrise);

		final double jd = JulianDay.of(date);
		final double t = JulianDay.centuryJ2000Of(jd);
		final EclipticCoordinate eclipticCoord = SunPosition.sunEclipticPosition(jd);
		final double equationOfTime = equationOfTime(eclipticCoord, t);
		final double geocentricMeanLongitude = eclipticCoord.getLongitude();
		final double meanAnomaly = geocentricMeanAnomaly(t);
		final double equationOfCenter = equationOfCenter(meanAnomaly, t);
		//Ltrue = L + C
		final double trueGeocentricLongitude = MathHelper.mod2pi(geocentricMeanLongitude + equationOfCenter);
		final double[] nutation = SunPosition.nutationCorrection(t);
		final double aberration = SunPosition.aberrationCorrection(eclipticCoord.getDistance());
		final double apparentGeocentricLatitude = SunPosition.apparentGeocentricLongitude(geocentricMeanLongitude, nutation[0],
			aberration);
		final double apparentGeocentricLongitude = SunPosition.apparentGeocentricLongitude(geocentricMeanLongitude, nutation[0],
			aberration);


//		delta = m_longitude + hourAngle;
//		timeDiff = 4 * delta;
//		timeUTC = 720 - timeDiff - eqTime; // in minutes

		final double meanEclipticObliquity = SunPosition.meanEclipticObliquity(t);
		final double apparentEclipticObliquity = apparentEclipticObliquity(meanEclipticObliquity, t);
		EquatorialCoordinate coord = EquatorialCoordinate.createFromEcliptical(apparentGeocentricLatitude, apparentGeocentricLongitude,
			apparentEclipticObliquity);
		final double apparentDeclination = coord.getDeclination();
		final double localHourAngle = localHourAngle(trueGeocentricLongitude, solarZenith, sunrise, apparentDeclination);
		final double localMeanTime = getLocalMeanTime(trueGeocentricLongitude, longitudeHour, localHourAngle);
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
	 * @param trueLong	Sun's true longitude [rad.
	 * @param zenith	Zenith enumeration.
	 * @param sunrise	Whether it's sunrise or sunset.
	 * @return	The Sun local hour [h].
	 * @throws SolarEventException	Whenever the Sun never rises or sets.
	 */
	private double localHourAngle(final double trueLong, final Zenith zenith, final Boolean sunrise,
		double apparentDeclination) throws SolarEventException{
		final double cosLocalHour = cosineLocalHour(trueLong, zenith);

		double cosHA = StrictMath.sin(zenith.getAngle())
			/ (StrictMath.cos(location.getLatitude()) * StrictMath.cos(apparentDeclination))
			- StrictMath.tan(location.getLatitude()) * StrictMath.tan(apparentDeclination);

		if(cosLocalHour < -1.)
			//the sun never sets on this location on the specified date
			throw SolarEventException.create(SolarEventError.NEVER_SETS);
		if(cosLocalHour > 1.)
			//the sun never rises on this location on the specified date
			throw SolarEventException.create(SolarEventError.NEVER_RISES);

		final double localHour = StrictMath.acos(cosLocalHour);
//		final double localHour = StrictMath.acos(cosHA);
		return MathHelper.degToHrs(sunrise? 2. * StrictMath.PI - localHour: localHour);
	}

	private double cosineLocalHour(final double trueLong, final Zenith zenith){
		final double sinDeclination = 0.39782 * StrictMath.sin(trueLong);
		final double cosDeclination = StrictMath.cos(StrictMath.asin(sinDeclination));

		final double sinZenith = StrictMath.sin(zenith.getAngle());
		final double sinLatitude = StrictMath.sin(location.getLatitude());
		final double cosLatitude = StrictMath.cos(location.getLongitude());

		return (sinZenith - sinLatitude * sinDeclination) / (cosLatitude * cosDeclination);
	}

	/**
	 * Computes the local mean time.
	 *
	 * @param trueGeocentricLongitude	Suns true longitude [rad].
	 * @param longitudeHour	Longitude hour [h].
	 * @param localHourAngle	Local hour of the Sun [h].
	 * @return	The local mean time [h].
	 */
	private static double getLocalMeanTime(final double trueGeocentricLongitude, final double longitudeHour, final double localHourAngle){
		final double rightAscension = getRightAscension(trueGeocentricLongitude);
		return MathHelper.limitRangeHour(localHourAngle + rightAscension - 0.06571 * longitudeHour - 6.622);
	}

	/**
	 * Computes the Suns right ascension, adjusting for the quadrant of the true longitude of the Sun and turning it
	 * into degree-hours.
	 *
	 * @param sunTrueLongitude	Suns true longitude [rad].
	 * @return	Suns right ascension in degree-hours.
	 */
	@SuppressWarnings("NumericCastThatLosesPrecision")
	private static double getRightAscension(final double sunTrueLongitude){
		final double tanL = StrictMath.tan(sunTrueLongitude);
		final double rightAscension = MathHelper.mod2pi(StrictMath.atan(0.91764 * tanL));

		final double longitudeQuadrant = 90. * (int)(sunTrueLongitude / 90.);
		final double rightAscensionQuadrant = 90. * (int)(rightAscension / 90.);
		return (rightAscension + longitudeQuadrant - rightAscensionQuadrant) / JulianDay.DEGREES_PER_HOUR;
	}

	@SuppressWarnings("NumericCastThatLosesPrecision")
	private LocalTime getLocalTime(final double localMeanTime){
		//adjust back to UTC
		final double utcTime = MathHelper.limitRangeHour(localMeanTime - getBaseLongitudeHour());
		return LocalTime.of(0, 0)
			.plus((long)(utcTime * 60. * 60.), ChronoUnit.SECONDS);
	}

}
