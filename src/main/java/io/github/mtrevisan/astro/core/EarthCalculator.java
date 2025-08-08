/*
 * Copyright (c) 2023 Mauro Trevisan
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
package io.github.mtrevisan.astro.core;

import io.github.mtrevisan.astro.helpers.JulianDate;
import io.github.mtrevisan.astro.helpers.MathHelper;
import io.github.mtrevisan.astro.helpers.TimeHelper;
import io.github.mtrevisan.astro.coordinates.EquatorialCoordinate;
import io.github.mtrevisan.astro.coordinates.GeographicLocation;

import java.time.LocalDate;
import java.time.ZoneId;
import java.time.ZonedDateTime;
import java.time.temporal.ChronoUnit;


public class EarthCalculator{

	//	private static final double EARTH_FLATTENING = 1. / 298.25642;
//	//[m]
//	private static final double EARTH_EQUATORIAL_RADIUS = 6378140.;
//
//	private static final double[] SUN_GEOCENTRIC_MEAN_ANOMALY = {357.52911, 35_999.050_29, -0.000_1537};
//
//
//	//https://frinklang.org/frinksamp/sun.frink
//	//https://www.astrouw.edu.pl/~jskowron/pracownia/praca/sunspot_answerbook_expl/expl-5.html
//	//http://co2.aos.wisc.edu/data/code/idl-lib/util/sunrise.pro
//	//https://ebvalaim.pl/en/2015/12/22/calculating-sunrise-and-sunset-times/
//
//	//---

	//1 s
	private static final double TIME_PRECISION = 1. / JulianDate.SECONDS_PER_DAY;
	/**
	 * <a href="https://en.wikipedia.org/wiki/Sidereal_time#Relationship_between_solar_time_and_sidereal_time_intervals">Sidereal time</a>
	 * <a href="https://iers-conventions.obspm.fr/content/tn36.pdf">IERS Conventions 2010</a>
	 */
	//[deg/day]
	private static final double[] EARTH_SIDEREAL_ROTATION_RATE_COEFF = {0., 1.002_737_811_911_354_48 * 360., 5.900_6e-11 * 360., -5.9e-15 * 360.};


	private enum SunVisibility{
		NORMAL, ALWAYS_DAY, ALWAYS_NIGHT
	}


	private final GeographicLocation location;


	/**
	 * Constructs a new instance using the given parameters.
	 *
	 * @param location	Location of the place.
	 */
	public static EarthCalculator create(final GeographicLocation location){
		return new EarthCalculator(location);
	}


	private EarthCalculator(final GeographicLocation location){
		this.location = location;
	}


	/**
	 * Calculate the times of sunlight phases (sunrise, sun transit or solar noon, and sunset) for a given day.
	 * <p>
	 * The definition of sunrise or sunset can be chosen based on a horizon type (defined via its elevation angle).
	 * </p>
	 *
	 * @param date	Date for which sunrise/transit/sunset are to be calculated.
	 * @param solarZenith	Solar zenith (basically, elevation angle) to use as the sunrise/sunset definition.
	 * 	This can be used to calculate twilight times.
	 * @return	The required time instants in Terrestrial Time from J2000.0.
	 *
	 * @see <a href="https://www.nrel.gov/docs/fy08osti/34302.pdf">Solar Position Algorithm for Solar Radiation Applications</a>
	 * @see <a href="https://midcdmz.nrel.gov/spa/">NREL's Solar Position Algorithm (SPA)</a>
	 */
	public final SunlightPhase sunlightPhase(final ZonedDateTime date, final Zenith solarZenith){
		return sunlightPhase(date, TimeHelper.deltaT(date), solarZenith);
	}

	/**
	 * Calculate the times of sunlight phases (sunrise, sun transit or solar noon, and sunset) for a given day.
	 * <p>
	 * The definition of sunrise or sunset can be chosen based on a horizon type (defined via its elevation angle).
	 * </p>
	 *
	 * @param date	Date for which sunrise/transit/sunset are to be calculated.
	 * @param deltaT	Difference between earth rotation time and terrestrial time (or Universal Time and Terrestrial Time) [s].
	 * @param solarZenith	Solar zenith (basically, elevation angle) to use as the sunrise/sunset definition.
	 * 	This can be used to calculate twilight times.
	 * @return	The required time instants in Terrestrial Time from J2000.0.
	 *
	 * @see <a href="https://www.nrel.gov/docs/fy08osti/34302.pdf">Solar Position Algorithm for Solar Radiation Applications</a>
	 * @see <a href="https://midcdmz.nrel.gov/spa/">NREL's Solar Position Algorithm (SPA)</a>
	 */
	public final SunlightPhase sunlightPhase(final ZonedDateTime date, final double deltaT, final Zenith solarZenith){
		final double ut = JulianDate.of(date.truncatedTo(ChronoUnit.DAYS));
		final double jce = JulianDate.centuryJ2000Of(ut);

		//A.2.1. Calculate the apparent sidereal time at Greenwich at 0 UT
		final NutationCorrections nutationCorrections = NutationCorrections.calculate(jce);
		final double meanEclipticObliquity = SunPosition.meanEclipticObliquity(jce);
		final double trueEclipticObliquity = SunPosition.trueEclipticObliquity(meanEclipticObliquity,
			nutationCorrections.getDeltaEpsilon());
		final double moonLongitudeAscendingNode = NutationCorrections.moonLongitudeAscendingNode(jce);
		final double greenwichMeanSiderealTime = TimeHelper.greenwichMeanSiderealTime(jce);
		final double greenwichApparentSiderealTime = TimeHelper.greenwichApparentSiderealTime(greenwichMeanSiderealTime,
			nutationCorrections.getDeltaPsi(), trueEclipticObliquity, moonLongitudeAscendingNode);


		//A.2.2. Calculate the geocentric right ascension and declination at 0 TT for day before, same day, and next day
		final EquatorialCoordinate[] equatorialCoords = new EquatorialCoordinate[3];
		for(int i = 0; i < equatorialCoords.length; i ++){
			final double jme0 = JulianDate.millenniumJ2000Of(ut + i - 1);
			equatorialCoords[i] = SunPosition.sunEquatorialPosition(jme0, nutationCorrections.getDeltaPsi(), trueEclipticObliquity);
		}


		//A.2.3. Calculate the approximate sun transit time, <code>m0</code> [day]
		final double[] m = new double[3];
		m[0] = (StrictMath.toDegrees(equatorialCoords[1].getRightAscension()) - location.getLongitude()
			- greenwichApparentSiderealTime) / 360.;


		//A.2.4. Calculate the local hour angle of the Sun, <code>H0</code>
		final double phi = StrictMath.toRadians(location.getLatitude());
		//the solar zenith angle is the correction for:
		// - observer altitude
		// - atmospheric refraction at sunrise/sunset
		// - the size of the solar disk
		final double correctionAltitude = StrictMath.toRadians(-1.75 * StrictMath.sqrt(location.getAltitude())
			/ JulianDate.MINUTES_PER_HOUR);
		double trueElevation = solarZenith.getElevation() + correctionAltitude;
		if(solarZenith == Zenith.OFFICIAL){
			//Bennett's equation (with a geometric elevation of 0°)
			final double atmosphericRefractionCorrection = location.getAtmosphere()
				.atmosphericRefractionCorrection(0.);
			//[rad]
			final double sunSolarDiskRadius = solarDiskRadius(jce);
			trueElevation -= atmosphericRefractionCorrection + sunSolarDiskRadius;
		}
		final double cosSunLocalHour = (StrictMath.sin(trueElevation)
			- StrictMath.sin(phi) * StrictMath.sin(equatorialCoords[1].getDeclination()))
			/ (StrictMath.cos(phi) * StrictMath.cos(equatorialCoords[1].getDeclination()));

		SunVisibility type = SunVisibility.NORMAL;
		if(cosSunLocalHour < -1.)
			type = SunVisibility.ALWAYS_DAY;
		else if(cosSunLocalHour > 1.)
			type = SunVisibility.ALWAYS_NIGHT;

		//[rad]
		final double sunLocalHour = MathHelper.modpipi(StrictMath.acos(cosSunLocalHour));


		//A.2.5. Calculate the approximate sunrise time, <code>m1</code>, in fraction of day
		m[1] = MathHelper.mod(m[0] - sunLocalHour / MathHelper.TWO_PI, 1.);


		//A.2.6. Calculate the approximate sunset time, <code>m2</code>, in fraction of day
		m[2] = MathHelper.mod(m[0] + sunLocalHour / MathHelper.TWO_PI, 1.);

		m[0] = MathHelper.mod(m[0], 1.);


		//A.2.8. Calculate the mean sidereal time at Greenwich [deg] for the sun transit, sunrise, and sunset
		final double[] greenwichSiderealTime = new double[3];
		for(int i = 0; i < m.length; i ++)
			greenwichSiderealTime[i] = greenwichApparentSiderealTime + MathHelper.polynomial(m[i], EARTH_SIDEREAL_ROTATION_RATE_COEFF);


		//A.2.9. Calculate the terms <code>n_i</code>
		final double[] n = new double[3];
		for(int i = 0; i < m.length; i ++)
			n[i] = m[i] + deltaT / JulianDate.SECONDS_PER_DAY;


		//A.2.10. Calculate the values alpha'i and delta'i [deg]
		final double a = limitIfNecessary(
			StrictMath.toDegrees(equatorialCoords[1].getRightAscension() - equatorialCoords[0].getRightAscension())
		);
		final double aPrime = limitIfNecessary(
			StrictMath.toDegrees(equatorialCoords[1].getDeclination() - equatorialCoords[0].getDeclination())
		);

		final double b = limitIfNecessary(
			StrictMath.toDegrees(equatorialCoords[2].getRightAscension() - equatorialCoords[1].getRightAscension())
		);
		final double bPrime = limitIfNecessary(
			StrictMath.toDegrees(equatorialCoords[2].getDeclination() - equatorialCoords[1].getDeclination())
		);

		final double c = b - a;
		final double cPrime = bPrime - aPrime;

		//apply interpolation
		final EquatorialCoordinate[] localEquatorialCoords = new EquatorialCoordinate[3];
		for(int i = 0; i < localEquatorialCoords.length; i ++)
			if(!Double.isNaN(n[i])){
				final double localRightAscension = StrictMath.toDegrees(equatorialCoords[1].getRightAscension())
					+ (n[i] * (a + b + c * n[i])) / 2.;
				final double localDeclination = StrictMath.toDegrees(equatorialCoords[1].getDeclination())
					+ (n[i] * (aPrime + bPrime + cPrime * n[i])) / 2.;
				localEquatorialCoords[i] = EquatorialCoordinate.create(localRightAscension, localDeclination);
			}


		//A.2.11. Calculate the local hour angle for the sun transit, sunrise, and sunset
		final double[] localHourAngle = new double[3];
		for(int i = 0; i < localHourAngle.length; i ++)
			if(localEquatorialCoords[i] != null){
				final double localSiderealTime = TimeHelper.localSiderealTime(greenwichSiderealTime[i], location.getLongitude());
				localHourAngle[i] = StrictMath.toRadians(
					TimeHelper.localHourAngle(localSiderealTime, localEquatorialCoords[i].getRightAscension())
				);
			}


		//A.2.12. Calculate the sun altitude for the sun transit, sunrise, and sunset, <code>h_i</code> [rad]
		final double[] sunAltitude = new double[3];
		for(int i = 0; i < sunAltitude.length; i ++)
			if(localEquatorialCoords[i] != null){
				final double declination = StrictMath.toRadians(localEquatorialCoords[i].getDeclination());
				sunAltitude[i] = StrictMath.asin(StrictMath.sin(phi) * StrictMath.sin(declination)
					+ StrictMath.cos(phi) * StrictMath.cos(declination) * StrictMath.cos(localHourAngle[i]));
			}


		//A.2.13. Calculate the sun transit, <code>T</code> [day]
		final double ttTransit = m[0] - localHourAngle[0] / MathHelper.TWO_PI;


		//A.2.14. Calculate the sunrise, <code>R</code> [day]
		double ttSunrise = Double.NaN;
		if(localEquatorialCoords[1] != null)
			ttSunrise = m[1] + (sunAltitude[1] - trueElevation)
				/ (MathHelper.TWO_PI
				* StrictMath.cos(StrictMath.toRadians(localEquatorialCoords[1].getDeclination()))
				* StrictMath.cos(phi)
				* StrictMath.sin(localHourAngle[1]));

		//A.2.15. Calculate the sunset, <code>S</code> [day]
		double ttSunset = Double.NaN;
		if(localEquatorialCoords[1] != null)
			ttSunset = m[2] + (sunAltitude[2] - trueElevation)
				/ (MathHelper.TWO_PI
				* StrictMath.cos(StrictMath.toRadians(localEquatorialCoords[2].getDeclination()))
				* StrictMath.cos(phi)
				* StrictMath.sin(localHourAngle[2]));


		return switch(type){
			case NORMAL -> new SunlightPhase.RegularDay(addFractionOfDay(date, ttSunrise), addFractionOfDay(date, ttTransit),
				addFractionOfDay(date, ttSunset));
			case ALWAYS_DAY -> new SunlightPhase.AlwaysDay(addFractionOfDay(date, ttTransit));
			case ALWAYS_NIGHT -> new SunlightPhase.AlwaysNight(addFractionOfDay(date, ttTransit));
		};
	}

	/**
	 * Calculate the apparent Sun's disk radius.
	 *
	 * @param jce	Julian Century of Terrestrial Time from J2000.0.
	 * @return	The Sun disk radius [rad].
	 */
	private static double solarDiskRadius(final double jce){
		final double sunRadiusVector = SunPosition.radiusVector(jce / 10.);
		return StrictMath.atan2(SunPosition.SUN_EQUATORIAL_RADIUS / SunPosition.ASTRONOMICAL_UNIT, sunRadiusVector);
	}

	/**
	 * Limit to 0..1 if absolute value > 2.
	 * <p>
	 * Refer to A.2.10 in NREL report.
	 * </p>
	 */
	private static double limitIfNecessary(final double val){
		return (StrictMath.abs(val) > 2.? MathHelper.mod(val, 1.): val);
	}

	private static ZonedDateTime addFractionOfDay(final ZonedDateTime date, final double fraction){
		return date.truncatedTo(ChronoUnit.DAYS)
			.plus((int)(fraction * JulianDate.MILLISECONDS_PER_DAY), ChronoUnit.MILLIS);
	}


	/**
	 * Calculate the times of season for a given year.
	 *
	 * @param year	Year for which sunrise/transit/sunset are to be calculated.
	 * @param zoneId	The time zone.
	 * @param northernHemisphereSeason	Season for the northern hemisphere to be calculated.
	 * @return	The required time instant in Terrestrial Time from J2000.0.
	 */
	public static ZonedDateTime season(final int year, final ZoneId zoneId, final Season northernHemisphereSeason){
		final double targetSunApparentLongitude = northernHemisphereSeason.ordinal() * Math.PI / 2.;
		final LocalDate date = LocalDate.of(year, 1, 1);
		//[day]
		final int yearLength = (date.isLeapYear()? 366: 365);
		double tt = JulianDate.of(date);
		final double minUT = tt;
		double correction;
		do{
			final double jce = JulianDate.centuryJ2000Of(tt);

			final NutationCorrections nutationCorrections = NutationCorrections.calculate(jce);
			final double sunApparentLongitude = SunPosition.apparentSunLongitude(jce / 10., nutationCorrections.getDeltaPsi());

			//[day]
			correction = 58. * StrictMath.sin(targetSunApparentLongitude - sunApparentLongitude);
			tt += correction;
			if(tt < minUT)
				tt += yearLength;
		}while(Math.abs(correction) > TIME_PRECISION);

		return JulianDate.dateTimeOf(tt)
			.atZone(zoneId);
	}


//	//---
//	//https://squarewidget.com/solar-coordinates/
//
//	//TODO
//	//https://github.com/MenoData/Time4J/blob/master/base/src/main/java/net/time4j/calendar/astro/SunPosition.java
////		//calculate the topocentric Sun Right Ascension: α'
////		final double topocentricRightAscension = rightAscension + deltaRightAscension;
////		//calculate the topocentric Sun declination: δ'
////		final double topocentricDeclination = StrictMath.atan2(
////			(StrictMath.sin(declination) - y * StrictMath.sin(equatorialHorizontalParallax)) * StrictMath.cos(deltaRightAscension),
////			StrictMath.cos(declination) - y * StrictMath.sin(equatorialHorizontalParallax) * StrictMath.cos(localHourAngle)
////		);
////		//calculate the topocentric local hour angle: H’
////		final double topocentricLocalHourAngle = localHourAngle - deltaRightAscension;
////		//calculate the true elevation angle without atmospheric refraction correction: e0
////		final double trueElevation = StrictMath.asin(
////			StrictMath.sin(latitude) * StrictMath.sin(topocentricDeclination)
////				+ StrictMath.cos(latitude) * StrictMath.cos(topocentricDeclination) * StrictMath.cos(topocentricLocalHourAngle)
////		);
////		//calculate the atmospheric refraction correction: Δe
////		final double deltaElevation = AtmosphereHelper.atmosphericRefractionCorrection(pressure, temperature, trueElevation);
////		//calculate the topocentric elevation angle: e
////		final double topocentricElevation = trueElevation + deltaElevation;
////		//calculate the topocentric zenith angle: θ
////		final double topocentricZenith = StrictMath.PI / 2. - topocentricElevation;
////		//calculate the topocentric astronomers azimuth angle (measured westward from south): Γ
////		final double topocentricAzimuth = StrictMath.atan2(
////			StrictMath.sin(topocentricLocalHourAngle),
////			StrictMath.cos(topocentricLocalHourAngle) * StrictMath.sin(latitude)
////				- StrictMath.tan(topocentricDeclination) * StrictMath.cos(latitude)
////		);
////		//calculate the (navigators) topocentric azimuth angle (measured westward from north): M
////		final double topocentricAzimuthNavigators = MathHelper.mod2pi(topocentricAzimuth + StrictMath.PI);
////		HorizontalCoordinate.create(topocentricAzimuth, topocentricElevation, eclipticCoord.getDistance());
//
//	//WORKS
//	//https://stellafane.org/misc/equinox.html
//	//<a href="https://www.agopax.it/Libri_astronomia/pdf/Astronomical%20Algorithms.pdf">Meeus, Jean. Astronomical algorithms. 2nd ed. 1998.</a>
//
//	//https://www.nrel.gov/docs/fy08osti/34302.pdf
//	//https://midcdmz.nrel.gov/solpos/spa.html
//	//http://phpsciencelabs.us/wiki_programs/Sidereal_Time_Calculator.php
//	//https://lweb.cfa.harvard.edu/~jzhao/times.html
//	//https://www.meteopiateda.it/busteggia/pages/astronomy/equisol.php
//	//https://www.researchgate.net/publication/262200491_Calcolo_analitico_della_posizione_del_sole_per_l%27allineamento_di_impianti_solari_ed_altre_applicazioni
//	//https://www.suncalc.org/#/45.7149,12.194179,17/2022.06.27/14:23/1/3
//
//		//calculate the parallax in the sun right ascension: Δα
//		final double equatorialHorizontalParallax = SunPosition.equatorialHorizontalParallax(eclipticCoord.getDistance());
//		final double phi = StrictMath.toRadians(location.getLatitude());
//		final double u = StrictMath.atan((1. - EARTH_FLATTENING) * StrictMath.tan(phi));
//		final double chi = StrictMath.cos(u) + (location.getAltitude() / EARTH_EQUATORIAL_RADIUS) * StrictMath.cos(phi);
//		final double y = 0.99664719 * StrictMath.sin(u) + (location.getAltitude() / EARTH_EQUATORIAL_RADIUS) * StrictMath.sin(phi);
//		final double deltaRightAscension = StrictMath.atan2(
//			-chi * StrictMath.sin(equatorialHorizontalParallax) * StrictMath.sin(localHourAngle),
//			StrictMath.cos(coord.getDeclination()) - chi * StrictMath.sin(equatorialHorizontalParallax) * StrictMath.cos(localHourAngle)
//		);
//		//calculate the topocentric Sun Right Ascension: α'
//		final double rightAscensionTopocentric = coord.getRightAscension() + deltaRightAscension;
//		if(Math.abs(rightAscensionTopocentric - StrictMath.toRadians(202.22704)) > 0.00001)
//			throw new IllegalArgumentException("rightAscensionTopocentric: " + rightAscensionTopocentric);
//		final double declinationTopocentric = StrictMath.atan2(
//			((StrictMath.sin(coord.getDeclination()) - y * StrictMath.sin(equatorialHorizontalParallax))
//				* StrictMath.cos(deltaRightAscension)),
//			StrictMath.cos(coord.getDeclination()) - chi * StrictMath.sin(equatorialHorizontalParallax) * StrictMath.cos(localHourAngle)
//		);
//		if(Math.abs(declinationTopocentric - StrictMath.toRadians(-9.316179)) > 0.000001)
//			throw new IllegalArgumentException("declinationTopocentric: " + declinationTopocentric);
//		//calculate the topocentric local hour angle: H’
//		final double localHourAngleTopocentric = localHourAngle - deltaRightAscension;
//		if(Math.abs(localHourAngleTopocentric - 11.10629) > 0.00002)
//			throw new IllegalArgumentException("localHourAngleTopocentric: " + (localHourAngleTopocentric - 11.10629));
//		//calculate the topocentric elevation angle without atmospheric refraction correction: e0
//		final double e0 = StrictMath.asin(
//			StrictMath.sin(phi) * StrictMath.sin(declinationTopocentric)
//				+ StrictMath.cos(phi) * StrictMath.cos(declinationTopocentric) * StrictMath.cos(localHourAngleTopocentric)
//		);
//		final double deltaE = atmosphericModel.atmosphericRefractionCorrection(e0);
//		//calculate the topocentric elevation angle: e
//		final double elevationTopocentric = e0 + deltaE;
//		//calculate the topocentric zenith angle: θ
//		final double zenithTopocentric = 90. - elevationTopocentric;
//		if(Math.abs(zenithTopocentric - 116.01759) > 0.00001)
//			throw new IllegalArgumentException("zenithTopocentric: " + (zenithTopocentric - 116.01759));
//		//calculate the topocentric astronomers azimuth angle (measured westward from south): Γ
//		final double azimuthTopocentric = MathHelper.mod2pi(StrictMath.atan2(
//			StrictMath.sin(localHourAngleTopocentric),
//			StrictMath.cos(localHourAngleTopocentric) * StrictMath.sin(phi)
//				- StrictMath.tan(declinationTopocentric) * StrictMath.cos(phi)
//		));
//		//calculate the topocentric azimuth angle (measured westward from north): M
//		final double azimuthTopocentricNavigators = MathHelper.mod2pi(azimuthTopocentric + StrictMath.PI);
//		if(Math.abs(azimuthTopocentricNavigators - 279.725959) > 0.000001)
//			throw new IllegalArgumentException("azimuthTopocentricNavigators: " + (azimuthTopocentricNavigators - 279.725959));
//
//	/**
//	 * Calculate the eccentricity of Earth's orbit, e.
//	 *
//	 * @param tt	Julian Century of Terrestrial Time from J2000.0.
//	 * @return	The eccentricity of Earth's orbit.
//	 */
//	private static double earthOrbitEccentricity(final double tt){
//		return MathHelper.polynomial(tt, new double[]{0.016708634, -0.000042037, -0.0000001267});
//	}
//
//	/**
//	 * Calculate the longitude of the perihelion of the orbit, π.
//	 *
//	 * @param tt	Julian Century of Terrestrial Time from J2000.0.
//	 * @return	The longitude of the perihelion of the orbit [rad].
//	 */
//	private static double longitudeOfEarthPerihelion(final double tt){
//		return MathHelper.polynomial(tt, new double[]{102.93735, 1.71946, 0.00046});
//	}
//
//	/**
//	 * Calculate the Sun's radius angle, corrected for distance, but not for refraction.
//	 *
//	 * @param radius	Sun radius [^].
//	 * @param distance	Sun distance [AU].
//	 * @return	The longitude of the perihelion of the orbit [rad].
//	 */
//	private static double radiusAngle(final double radius, final double distance){
//		return StrictMath.asin(radius / (radius + distance));
//	}
//
//	/**
//	 * Calculate the Equation of Time.
//	 *
//	 * @param eclipticCoord	Mean ecliptic coordinate of the Sun.
//	 * @param tt	Julian Century of Terrestrial Time from J2000.0.
//	 * @return	The Equation of Time [h].
//	 */
//	@SuppressWarnings("OverlyComplexArithmeticExpression")
//	public static double equationOfTime(final EclipticCoordinate eclipticCoord, final double tt){
//		final double e0 = SunPosition.meanEclipticObliquity(tt);
//		final double epsilon = apparentEclipticObliquity(e0, tt);
//		final double l0 = eclipticCoord.getLongitude();
//		final double e = earthOrbitEccentricity(tt);
//		final double m = geocentricMeanAnomaly(tt);
//		double y = StrictMath.tan(epsilon / 2.);
//		y *= y;
//
//		return y * StrictMath.sin(2. * l0)
//			- 2. * e * StrictMath.sin(m)
//			+ 4. * e * y * StrictMath.sin(m) * StrictMath.cos(2. * l0)
//			- 0.5 * y * y * StrictMath.sin(4. * l0)
//			- 1.25 * e * e * StrictMath.sin(2. * m) * 4. / JulianDate.MINUTES_PER_HOUR;
//	}
//
//	/**
//	 * Calculate the Equation of Time.
//	 *
//	 * @param sunMeanLongitude	Geocentric mean longitude of the Sun, <code>L0</code>.
//	 * @param sunMeanAnomaly	Geocentric mean anomaly of the Sun, <code>M</code>.
//	 * @param apparentEclipticObliquity	Mean obliquity of the ecliptic, corrected for parallax, <code>ɛ'</code>.
//	 * @param earthEccentricity	Eccentricity of Earth's orbit, <code>e</code>.
//	 * @return	The Equation of Time [min].
//	 */
//	public static double equationOfTime(final double sunMeanLongitude, final double sunMeanAnomaly, final double apparentEclipticObliquity,
//			final double earthEccentricity){
//		double y = StrictMath.tan(apparentEclipticObliquity / 2.);
//		y *= y;
//		final double sin2L0 = StrictMath.sin(2. * sunMeanLongitude);
//		final double cos2L0 = StrictMath.cos(2. * sunMeanLongitude);
//		final double sin4L0 = StrictMath.sin(4. * sunMeanLongitude);
//		final double sinM = StrictMath.sin(sunMeanAnomaly);
//		final double sin2M = StrictMath.sin(2. * sunMeanAnomaly);
//		return 4. * StrictMath.toDegrees(y * sin2L0
//			- 2. * earthEccentricity * sinM
//			+ 4. * earthEccentricity * y * sinM * cos2L0
//			- 0.5 * y * y * sin4L0
//			- 1.25 * earthEccentricity * earthEccentricity * sin2M);
//	}
//
//
//	/**
//	 * Computes the longitude time.
//	 *
//	 * @return	Longitudinal time [day].
//	 */
//	private double getLongitudeHour(final LocalDateTime date, final Boolean sunrise){
//		final double dividend = (sunrise? 6: 18) - getBaseLongitudeHour();
//		return date.getDayOfYear() + JulianDate.timeOf(date) + dividend / JulianDate.HOURS_PER_DAY;
//	}
//
//	/**
//	 * Computes the base longitude hour.
//	 *
//	 * @return	The longitude of the location of the solar event divided by 15 (deg/hour).
//	 */
//	private double getBaseLongitudeHour(){
//		return MathHelper.degToHrs(location.getLongitude());
//	}
//
//	/**
//	 * Computes the Suns right ascension, adjusting for the quadrant of the true longitude of the Sun and turning it
//	 * into degree-hours.
//	 *
//	 * @param sunTrueLongitude	Suns true longitude [rad].
//	 * @return	Suns right ascension in degree-hours.
//	 */
//	@SuppressWarnings("NumericCastThatLosesPrecision")
//	private static double getRightAscension(final double sunTrueLongitude){
//		final double tanL = StrictMath.tan(sunTrueLongitude);
//		final double rightAscension = MathHelper.mod2pi(StrictMath.atan(0.91764 * tanL));
//
//		final double longitudeQuadrant = 90. * (int)(sunTrueLongitude / 90.);
//		final double rightAscensionQuadrant = 90. * (int)(rightAscension / 90.);
//		return (rightAscension + longitudeQuadrant - rightAscensionQuadrant) / JulianDate.DEGREES_PER_HOUR;
//	}
//
//	@SuppressWarnings("NumericCastThatLosesPrecision")
//	private LocalTime getLocalTime(final double localMeanTime){
//		//adjust back to UTC
//		final double utcTime = MathHelper.limitRangeHour(localMeanTime - getBaseLongitudeHour());
//		return LocalTime.of(0, 0)
//			.plus((long)(utcTime * JulianDate.SECONDS_PER_HOUR), ChronoUnit.SECONDS);
//	}

}
