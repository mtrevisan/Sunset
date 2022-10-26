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
import io.github.mtrevisan.sunset.ResourceReader;
import io.github.mtrevisan.sunset.TimeHelper;
import io.github.mtrevisan.sunset.coordinates.EclipticCoordinate;
import io.github.mtrevisan.sunset.coordinates.EquatorialCoordinate;
import io.github.mtrevisan.sunset.coordinates.GNSSLocation;
import io.github.mtrevisan.sunset.coordinates.HorizontalCoordinate;

import java.io.IOException;
import java.time.LocalDateTime;
import java.util.Collection;
import java.util.List;
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

https://dspace.mit.edu/bitstream/handle/1721.1/112396/1013461937-MIT.pdf
https://github.com/asif31iqbal/adhanalarm/blob/master/AdhanAlarm/src/net/sourceforge/jitl/astro/AstroLib.java
https://www.scielo.br/j/jatm/a/ZZ6L3gWKF9LpvW8tLcQPx8p/?lang=en&format=pdf
https://ftp.imcce.fr/pub/ephem/planets/vsop2013/solution/
*/
public final class SunPosition{

	private static Map<String, Collection<Double[]>> EARTH_HELIOCENTRIC_DATA;
	private static Map<ResourceReader.VariableIndex, List<ResourceReader.VSOP2013Data>> EARTH_HELIOCENTRIC_DATA2;
	private static Map<String, Collection<Double[]>> NUTATION_DATA;
	static{
		try{
			EARTH_HELIOCENTRIC_DATA = ResourceReader.read("earthHeliocentric.dat");
			EARTH_HELIOCENTRIC_DATA2 = ResourceReader.readData("VSOP2013p3.dat");
			NUTATION_DATA = ResourceReader.read("nutation.dat");
		}
		catch(final IOException ignored){}
	}

	static final double EARTH_FLATTENING = 1. / 298.25642;
	//[m]
	static final double EARTH_EQUATORIAL_RADIUS = 6378140.;
	private static final double[] EARTH_ORBIT_ECCENTRICITY = {0.016_708_634, -0.000_042_037, -0.000_000_1267};

	private static final double[] SUN_GEOCENTRIC_MEAN_ANOMALY_PARAMETERS = {357.5277233, 35999.0503400, -0.0001603, -1. / 300000.};
	private static final double[] SUN_EQUATION_OF_CENTER_1 = {1.914_602, -0.004_817, -0.000_014};
	private static final double[] SUN_EQUATION_OF_CENTER_2 = {0.019_993, -0.000_101};

	private static final double[] MOON_MEAN_ELONGATION_PARAMETERS = {297.8503631, 445267.1114800, -0.0019142, 1. / 189474.};
	private static final double[] MOON_MEAN_ANOMALY_PARAMETERS = {134.9629814, 477198.8673981, 0.0086972, 1. / 56250.};
	private static final double[] MOON_ARGUMENT_OF_LATITUDE = {93.2719103, 483202.0175381, -0.00368250, 1. / 327270.};
	private static final double[] MOON_LONGITUDE_ASCENDING_NODE = {125.0445222, -1934.1362608, 0.002070833, 1. / 450000.};
	private static final double[] MEAN_ECLIPTIC_OBLIQUITY_PARAMETERS = {84381.448, -4680.93, -1.55, 1999.25, -51.38, -249.67, -39.05, 7.12,
		27.87, 5.79, 2.45};


	private SunPosition(){}


	//https://hal-mines-paristech.archives-ouvertes.fr/hal-00725987/document


	/**
	 * Calculated the Sun ecliptic position.
	 *
	 * @param jd	Julian Day of Terrestrial Time from J2000.0.
	 * @return	The Sun ecliptic position.
	 */
	public static EclipticCoordinate sunEclipticPosition(final double jd){
		final double jme = JulianDay.millenniumJ2000Of(jd);

		final double geocentricMeanLatitude = geocentricMeanLatitude(jme);
		final double geocentricMeanLongitude = geocentricMeanLongitude(jme);
		final double radiusVector = radiusVector(jme);
final double radiusVectorApprox = 0.016704 * StrictMath.cos(2. * StrictMath.PI * jd / 365.254902 + 3.091159) + 1.000140;
		return EclipticCoordinate.create(geocentricMeanLatitude, geocentricMeanLongitude, radiusVector);
	}

	/**
	 * Calculated the Sun equatorial position.
	 *
	 * @param eclipticCoord	Mean ecliptic coordinate of the Sun.
	 * @param jd	Julian Day of Terrestrial Time from J2000.0.
	 * @return	The Sun equatorial position.
	 *
	 * @see <a href="https://squarewidget.com/solar-coordinates/">Solar coordinates</a>
	 */
	public static EquatorialCoordinate sunEquatorialPosition(final EclipticCoordinate eclipticCoord, final double jd){
		final double tt = JulianDay.centuryJ2000Of(jd);

		//calculate the nutation in longitude and obliquity
		final double[] nutation = nutationCorrection(tt);
		//calculate the aberration correction
		final double aberration = aberrationCorrection(eclipticCoord.getDistance());
		//calculate the apparent Sun longitude: Ltrue = L + C
		final double apparentGeocentricLongitude = apparentGeocentricLongitude(eclipticCoord.getLongitude(), nutation[0], aberration);

		//calculate the obliquity of the ecliptic (the inclination of the Earth’s equator with respect to the plane at which the Sun
		//and planets appear to move across the sky): ɛ0
		final double meanEclipticObliquity = meanEclipticObliquity(tt);
		//calculate the true obliquity of the ecliptic
		final double trueEclipticObliquity = trueEclipticObliquity(meanEclipticObliquity, nutation[1]);

		return EquatorialCoordinate.createFromEcliptical(eclipticCoord.getLatitude(), apparentGeocentricLongitude, trueEclipticObliquity);
	}

	/**
	 * Calculated the Sun topocentric position.
	 *
	 * @param location	Location of the place.
	 * @param atmosphericModel	Atmospheric model of the place.
	 * @param eclipticCoord	Mean ecliptic coordinate of the Sun.
	 * @param jd	Julian Day of Terrestrial Time from J2000.0.
	 * @return	The Sun topocentric position.
	 */
	public static HorizontalCoordinate sunTopocentricPosition(final GNSSLocation location, final AtmosphericModel atmosphericModel,
			final EclipticCoordinate eclipticCoord, final double jd){
		final LocalDateTime date = JulianDay.dateTimeOf(jd);
		final double dt = TimeHelper.deltaT(date.getYear());
		final double ut = TimeHelper.terrestrialTimeToUniversalTime(jd, dt);
		final double tt = JulianDay.centuryJ2000Of(jd);

		final EquatorialCoordinate equatorialCoord = sunEquatorialPosition(eclipticCoord, jd);
		final double[] nutation = nutationCorrection(tt);
		final double meanEclipticObliquity = meanEclipticObliquity(tt);
		final double trueEclipticObliquity = trueEclipticObliquity(meanEclipticObliquity, nutation[1]);

		final double meanSiderealTime = TimeHelper.meanSiderealTime(ut);
		final double apparentSiderealTime = TimeHelper.apparentSiderealTime(meanSiderealTime, trueEclipticObliquity, nutation[0]);
		final double localMeanSiderealTime = TimeHelper.localMeanSiderealTime(apparentSiderealTime, location);
		final double rightAscension = equatorialCoord.getRightAscension();
		final double localHourAngle = TimeHelper.localHourAngle(localMeanSiderealTime, rightAscension);

		//compute the sun position (right ascension and declination) with respect to the observer local position at the Earth surface:
		final double equatorialHorizontalParallax = equatorialHorizontalParallax(eclipticCoord.getDistance());
		final double latitude = location.getLatitude();
		final double u = StrictMath.atan((1. - EARTH_FLATTENING) * StrictMath.tan(latitude));
		final double height = location.getAltitude() / EARTH_EQUATORIAL_RADIUS;
		final double x = StrictMath.cos(u) + height * StrictMath.cos(latitude);
		final double y = 0.99664719 * StrictMath.sin(u) + height * StrictMath.sin(latitude);
		final double declination = equatorialCoord.getDeclination();
		final double deltaRightAscension = StrictMath.atan2(
			-x * StrictMath.sin(equatorialHorizontalParallax) * StrictMath.sin(localHourAngle),
			StrictMath.cos(declination) - x * StrictMath.sin(equatorialHorizontalParallax)
				* StrictMath.cos(localHourAngle)
		);
		//calculate the topocentric Sun Right Ascension: α'
		//final double topocentricRightAscension = rightAscension + deltaRightAscension;
		//calculate the topocentric Sun declination: δ'
		final double topocentricDeclination = StrictMath.atan2(
			(StrictMath.sin(declination) - y * StrictMath.sin(equatorialHorizontalParallax)) * StrictMath.cos(deltaRightAscension),
			StrictMath.cos(declination) - y * StrictMath.sin(equatorialHorizontalParallax) * StrictMath.cos(localHourAngle)
		);
		//calculate the topocentric local hour angle: H’
		final double topocentricLocalHourAngle = localHourAngle - deltaRightAscension;
		//calculate the true elevation angle without atmospheric refraction correction: e0
		final double trueElevation = StrictMath.asin(
			StrictMath.sin(latitude) * StrictMath.sin(topocentricDeclination)
				+ StrictMath.cos(latitude) * StrictMath.cos(topocentricDeclination) * StrictMath.cos(topocentricLocalHourAngle)
		);
		//calculate the atmospheric refraction correction: Δe
		final double deltaElevation = atmosphericModel.atmosphericRefractionCorrection(trueElevation);
		//calculate the topocentric elevation angle: e
		final double topocentricElevation = MathHelper.modpi(trueElevation + deltaElevation);
		//calculate the topocentric zenith angle: θ
		//final double topocentricZenith = StrictMath.PI / 2. - topocentricElevation;
		//calculate the topocentric astronomers azimuth angle (measured westward from south): Γ
		final double topocentricAzimuth = StrictMath.atan2(
			StrictMath.sin(topocentricLocalHourAngle),
			StrictMath.cos(topocentricLocalHourAngle) * StrictMath.sin(latitude)
				- StrictMath.tan(topocentricDeclination) * StrictMath.cos(latitude)
		);
		//calculate the (navigators) topocentric azimuth angle (measured westward from north): M
		final double topocentricAzimuthNavigators = MathHelper.mod2pi(topocentricAzimuth + StrictMath.PI);
		return HorizontalCoordinate.create(topocentricAzimuthNavigators, topocentricElevation, eclipticCoord.getDistance());
	}


	/**
	 * Calculate the geocentric mean latitude of the Sun, referred to the mean equinox of the date, β.
	 *
	 * @param jme	Julian Ephemeris Millennium of Terrestrial Time from J2000.0.
	 * @return	The geocentric mean latitude of the Sun [rad].
	 *
	 * @see <a href="https://squarewidget.com/solar-coordinates/">Solar coordinates</a>
	 */
	private static double geocentricMeanLatitude(final double jme){
		final double[] parameters = new double[6];
		for(int i = 0; i < parameters.length; i ++){
			final Collection<Double[]> elements = EARTH_HELIOCENTRIC_DATA.get("B" + i);
			if(elements != null){
				double parameter = 0.;
				for(final Double[] element : elements)
					parameter += element[0] * StrictMath.cos(element[1] + element[2] * jme);
				parameters[i] = parameter;
			}
		}
		//transform heliocentric to geocentric
		return -MathHelper.modpi(
			MathHelper.eval(jme, parameters) / 100_000_000.
		);
	}

	/**
	 * Calculate the geocentric mean longitude of the Sun, referred to the mean equinox of the date, L.
	 *
	 * @param jme	Julian Ephemeris Millennium of Terrestrial Time from J2000.0.
	 * @return	The geocentric mean longitude of the Sun [rad].
	 *
	 * @see <a href="https://squarewidget.com/solar-coordinates/">Solar coordinates</a>
	 */
	public static double geocentricMeanLongitude(final double jme){
		final double[] parameters = new double[6];
		for(int i = 0; i < parameters.length; i ++){
			final Collection<Double[]> elements = EARTH_HELIOCENTRIC_DATA.get("L" + i);
			if(elements != null){
				double parameter = 0.;
				for(final Double[] element : elements)
					parameter += element[0] * StrictMath.cos(element[1] + element[2] * jme);
				parameters[i] = parameter;
			}
		}
		//transform heliocentric to geocentric
		return MathHelper.mod2pi(
			MathHelper.eval(jme, parameters) / 100_000_000. + StrictMath.PI
		);
	}

	/**
	 * Calculate the distance between the center of the Sun and the center of the Earth, R.
	 * <p>U.S. Naval Observatory function.</p>
	 *
	 * @param jme	Julian Ephemeris Millennium of Terrestrial Time from J2000.0.
	 * @return	Distance between the center of the Sun and the center of the Earth [AU].
	 */
	public static double radiusVector(final double jme){
		final double[] parameters = new double[6];
		for(int i = 0; i < parameters.length; i ++){
			final Collection<Double[]> elements = EARTH_HELIOCENTRIC_DATA.get("R" + i);
			if(elements != null){
				double parameter = 0.;
				for(final Double[] element : elements)
					parameter += element[0] * StrictMath.cos(element[1] + element[2] * jme);
				parameters[i] = parameter;
			}
		}
		return MathHelper.eval(jme, parameters) / 100_000_000.;
	}

	/**
	 * Calculate the Sun's equation of center, C.
	 *
	 * @param geocentricMeanAnomaly	The mean anomaly of the Sun [rad].
	 * @param tt	Julian Century of Terrestrial Time from J2000.0.
	 * @return	The Sun's equation of center [rad].
	 */
	public static double equationOfCenter(final double geocentricMeanAnomaly, final double tt){
		return Math.toRadians(
			MathHelper.eval(tt, SUN_EQUATION_OF_CENTER_1) * StrictMath.sin(geocentricMeanAnomaly)
				+ MathHelper.eval(tt, SUN_EQUATION_OF_CENTER_2) * StrictMath.sin(2. * geocentricMeanAnomaly)
				+ 0.000_289 * StrictMath.sin(3. * geocentricMeanAnomaly)
		);
	}

	/**
	 * Calculate the eccentricity of Earth's orbit, e.
	 *
	 * @param tt	Julian Century of Terrestrial Time from J2000.0.
	 * @return	The eccentricity of Earth's orbit.
	 */
	public static double earthOrbitEccentricity(final double tt){
		return MathHelper.eval(tt, EARTH_ORBIT_ECCENTRICITY);
	}


	/**
	 * Calculate corrections of nutation in longitude (∆ψ) and obliquity (∆ε).
	 *
	 * @param tt	Julian Century of Terrestrial Time from J2000.0.
	 * @return	An array where the first element is ∆ψ [rad], and the second ∆ε [rad].
	 *
	 * http://physics.uwyo.edu/~wiro/planet/nutate.c
	 * https://iopscience.iop.org/article/10.1086/375641/pdf
	 */
	public static double[] nutationCorrection(final double tt){
		final double d = meanElongationMoonSun(tt);
		final double m = meanAnomalySun(tt);
		final double mp = meanAnomalyMoon(tt);
		final double f = argumentLatitudeMoon(tt);
		final double omega = ascendingLongitudeMoon(tt);

		final Collection<Double[]> elements = NUTATION_DATA.get("coeffs");
		final double[] x = {d, m, mp, f, omega};
		double deltaPsi = 0.;
		double deltaEpsilon = 0.;
		for(final Double[] element : elements){
			final double parameter = xyTermSummation(x, element);
			deltaPsi += (element[x.length] + element[x.length + 1] * tt) * StrictMath.sin(parameter);
			deltaEpsilon += (element[x.length + 2] + element[x.length + 3] * tt) * StrictMath.cos(parameter);
		}
		//[°]
		deltaPsi /= JulianDay.SECONDS_IN_HOUR * 10_000.;
		//[°]
		deltaEpsilon /= JulianDay.SECONDS_IN_HOUR * 10_000.;
		return new double[]{StrictMath.toRadians(deltaPsi), StrictMath.toRadians(deltaEpsilon)};
	}

	//mean elongation of the Moon from the Sun [rad]
	private static double meanElongationMoonSun(final double tt){
		return MathHelper.mod2pi(StrictMath.toRadians(
			MathHelper.eval(tt, MOON_MEAN_ELONGATION_PARAMETERS)
		));
	}

	/**
	 * Calculate the geocentric mean anomaly of the Sun, M.
	 *
	 * @param tdb	Julian Century of Terrestrial Time from J2000.0.
	 * @return	The geocentric mean anomaly of the Sun [rad].
	 */
	public static double meanAnomalySun(final double tdb){
		return MathHelper.mod2pi(StrictMath.toRadians(
			MathHelper.eval(tdb, SUN_GEOCENTRIC_MEAN_ANOMALY_PARAMETERS)
		));
	}

	//mean anomaly of the Moon [rad]
	private static double meanAnomalyMoon(final double tt){
		return MathHelper.mod2pi(StrictMath.toRadians(
			MathHelper.eval(tt, MOON_MEAN_ANOMALY_PARAMETERS)
		));
	}

	//Moon's argument of Latitude [rad]
	private static double argumentLatitudeMoon(final double tt){
		return MathHelper.mod2pi(StrictMath.toRadians(
			MathHelper.eval(tt, MOON_ARGUMENT_OF_LATITUDE)
		));
	}

	//Longitude of the ascending node of the Moon's mean orbit on the ecliptic measured from the mean equinox of the date [rad]
	private static double ascendingLongitudeMoon(final double tt){
		return MathHelper.mod2pi(StrictMath.toRadians(
			MathHelper.eval(tt, MOON_LONGITUDE_ASCENDING_NODE)
		));
	}

	private static double xyTermSummation(final double[] x, final Double[] terms){
		double result = 0.;
		for(int j = 0; j < x.length; j ++)
			result += x[j] * terms[j];
		return result;
	}

	/**
	 * Calculate the correction for aberration (∆τ).
	 *
	 * @param earthRadiusVector	Earth radius vector [AU].
	 * @return	The correction for aberration [rad].
	 */
	static double aberrationCorrection(final double earthRadiusVector){
		return StrictMath.toRadians(
			-20.4898 / (JulianDay.SECONDS_IN_HOUR * earthRadiusVector)
		);
	}

	/**
	 * Calculate the apparent longitude of the Sun, Lapp = λ.
	 *
	 * @param geocentricMeanLongitude	Sun's true geocentric longitude [rad].
	 * @param deltaPsi	Nutation in longitude [rad].
	 * @param deltaAberration	Aberration [rad].
	 * @return	Apparent longitude of the Sun [rad].
	 */
	static double apparentGeocentricLongitude(final double geocentricMeanLongitude, final double deltaPsi, final double deltaAberration){
		return geocentricMeanLongitude + deltaPsi + deltaAberration;
	}

	/**
	 * Calculate the mean obliquity of the ecliptic, ɛ0.
	 *
	 * @param tt	Julian Century of Terrestrial Time from J2000.0.
	 * @return	Apparent longitude of the Sun [rad].
	 */
	public static double meanEclipticObliquity(final double tt){
		final double u = tt / 100.;
		return StrictMath.toRadians(
			MathHelper.eval(u, MEAN_ECLIPTIC_OBLIQUITY_PARAMETERS) / JulianDay.SECONDS_IN_HOUR
		);
	}

	/**
	 * True obliquity of the ecliptic corrected for nutation, ɛ.
	 *
	 * @param meanEclipticObliquity	Obliquity of the ecliptic, corrected for parallax [rad].
	 * @param deltaEpsilon	Nutation in obliquity [rad.
	 * @return	Apparent longitude of the Sun [rad].
	 */
	public static double trueEclipticObliquity(final double meanEclipticObliquity, final double deltaEpsilon){
		return meanEclipticObliquity + deltaEpsilon;
	}


	/**
	 * Calculate the equatorial horizontal parallax of the Sun, ξ.
	 *
	 * @param radiusVector	Radius vector of the Earth [AU].
	 * @return	The equatorial horizontal parallax of the Sun [rad].
	 */
	private static double equatorialHorizontalParallax(double radiusVector){
		return StrictMath.toRadians(8.794 / (JulianDay.SECONDS_IN_HOUR * radiusVector));
	}

}
