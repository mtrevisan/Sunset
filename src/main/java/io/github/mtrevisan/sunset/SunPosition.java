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
package io.github.mtrevisan.sunset;

import io.github.mtrevisan.sunset.coordinates.EclipticCoordinate;
import io.github.mtrevisan.sunset.coordinates.EquatorialCoordinate;
import io.github.mtrevisan.sunset.coordinates.OrbitalElements;

import java.io.IOException;
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

	private static final double[] MOON_MEAN_ELONGATION_PARAMETERS = {297.8503631, 445267.1114800, -0.0019142, 1. / 189474.};
	private static final double[] SUN_MEAN_ANOMALY_PARAMETERS = {357.5277233, 35999.0503400, -0.0001603, -1. / 300000.};
	private static final double[] MOON_MEAN_ANOMALY_PARAMETERS = {134.9629814, 477198.8673981, 0.0086972, 1. / 56250.};
	private static final double[] MOON_ARGUMENT_OF_LATITUDE = {93.2719103, 483202.0175381, -0.00368250, 1. / 327270.};
	private static final double[] MOON_LONGITUDE_ASCENDING_NODE = {125.0445222, -1934.1362608, 0.002070833, 1. / 450000.};
	private static final double[] MEAN_ECLIPTIC_OBLIQUITY_PARAMETERS = {84381.448, -4680.93, -1.55, 1999.25, -51.38, -249.67, -39.05, 7.12,
		27.87, 5.79, 2.45};


	private SunPosition(){}


	/**
	 * Calculated the Sun ecliptic position.
	 *
	 * @param jd	Julian Day of Terrestrial Time from J2000.0.
	 * @return	The Sun ecliptic position.
	 */
	public static EclipticCoordinate sunEclipticPosition(final double jd){
		final double tt = JulianDay.centuryJ2000Of(jd);
		final double geometricMeanLatitude = geometricMeanLatitude(tt);
geometricMeanLatitude2(jd);
		final double geometricMeanLongitude = geometricMeanLongitude(tt);
		final double radiusVector = radiusVector(tt);
		final EclipticCoordinate coord = EclipticCoordinate.create(geometricMeanLatitude, geometricMeanLongitude, radiusVector);

		//conversion to the FK5 System:
//		final double lambdaPrime = geometricMeanLongitude + MathHelper.eval(tt, new double[]{0., -1.397, -0.00031});
//		final double deltaLatitude = (0.03916 / JulianDay.SECONDS_IN_HOUR) * (StrictMath.cos(lambdaPrime) - StrictMath.sin(lambdaPrime));
//		final double deltaLongitude = -0.09033 / JulianDay.SECONDS_IN_HOUR;
//		geometricMeanLatitude += deltaLatitude;
//		geometricMeanLongitude += deltaLongitude;

		return coord;
	}

	/**
	 * Calculated the Sun equatorial position.
	 *
	 * @param eclipticCoord	Geometric mean ecliptic coordinate of the Sun.
	 * @param jd	Julian Day of Terrestrial Time from J2000.0.
	 * @return	The Sun equatorial position.
	 *
	 * @see <a href="https://squarewidget.com/solar-coordinates/">Solar coordinates</>
	 */
	public static EquatorialCoordinate sunEquatorialPosition(final EclipticCoordinate eclipticCoord, final double jd){
		final double tt = JulianDay.centuryJ2000Of(jd);
		//calculate the nutation in longitude and obliquity
		final double[] nutation = nutationCorrection(tt);
		//calculate the aberration correction
		final double aberration = aberrationCorrection(eclipticCoord.getDistance());
		//calculate the apparent Sun longitude: Ltrue = L + C
		final double apparentGeometricLongitude = apparentGeometricLongitude(eclipticCoord.getLongitude(), nutation[0], aberration);

		//calculate the obliquity of the ecliptic (the inclination of the Earth’s equator with respect to the plane at which the Sun
		//and planets appear to move across the sky): ɛ0
		final double meanEclipticObliquity = meanEclipticObliquity(tt);
		//calculate the true obliquity of the ecliptic
		final double trueEclipticObliquity = trueEclipticObliquity(meanEclipticObliquity, nutation[1]);

		return EquatorialCoordinate.createFromEcliptical(eclipticCoord.getLatitude(), apparentGeometricLongitude, trueEclipticObliquity);
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
		final double[] parameters = new double[6];
		for(int i = 0; i < parameters.length; i ++){
			double parameter = 0.;
			final Collection<Double[]> elements = EARTH_HELIOCENTRIC_DATA.get("B" + i);
			if(elements != null){
				for(final Double[] element : elements)
					parameter += element[0] * StrictMath.cos(element[1] + element[2] * jme);
				parameters[i] = parameter;
			}
		}
		return -StrictMath.toDegrees(MathHelper.eval(jme, parameters) / 100_000_000.);
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
			if(elements != null){
				for(final Double[] element : elements)
					parameter += element[0] * StrictMath.cos(element[1] + element[2] * jme);
				parameters[i] = parameter;
			}
		}
		final double longitude = StrictMath.toDegrees(MathHelper.eval(jme, parameters) / 100_000_000.);
		return MathHelper.limitRangeDegree(longitude + 180.);
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
		final double[] parameters = new double[6];
		for(int i = 0; i < parameters.length; i ++){
			double parameter = 0.;
			final Collection<Double[]> elements = EARTH_HELIOCENTRIC_DATA.get("R" + i);
			if(elements != null){
				for(final Double[] element : elements)
					parameter += element[0] * StrictMath.cos(element[1] + element[2] * jme);
				parameters[i] = parameter;
			}
		}
		return MathHelper.eval(jme, parameters) / 100_000_000.;
	}


	/**
	 * Calculate the geometric mean latitude of the Sun, referred to the inertial frame defined by the dynamical equinox and ecliptic J2000,
	 * β.
	 *
	 * @param jd	Julian Day of Terrestrial Time from J2000.0.
	 * @return	The geometric mean latitude of the Sun [°].
	 *
	 * @see <a href="https://squarewidget.com/solar-coordinates/">Solar coordinates</>
	 * https://github.com/timmyd7777/SSCore/tree/a915f57f8754b206614d3f0b5055d1a7b56e9e70/SSCode/VSOP2013
	 */
	private static double geometricMeanLatitude2(final double jd){
		final double tt = JulianDay.centuryJ2000Of(jd);
		final double jme = tt / 10.;

		final double[] planetLongitudes = new double[17];
		planetLongitudes[0] =  4.402608631669 + 26087.90314068555 * jme;     // Mercury
		planetLongitudes[1] =  3.176134461576 + 10213.28554743445 * jme;     // Venus
		planetLongitudes[2] =  1.753470369433 +  6283.075850353215 * jme;    // Earth-Moon Barycenter
		planetLongitudes[3] =  6.203500014141 +  3340.612434145457 * jme;    // Mars
		planetLongitudes[4] =  4.091360003050 +  1731.170452721855 * jme;    // Vesta
		planetLongitudes[5] =  1.713740719173 +  1704.450855027201 * jme;    // Iris
		planetLongitudes[6] =  5.598641292287 +  1428.948917844273 * jme;    // Bamberga
		planetLongitudes[7] =  2.805136360408 +  1364.756513629990 * jme;    // Ceres
		planetLongitudes[8] =  2.326989734620 +  1361.923207632842 * jme;    // Pallas
		planetLongitudes[9]  = 0.599546107035 +   529.6909615623250 * jme;   // Jupiter
		planetLongitudes[10] = 0.874018510107 +   213.2990861084880 * jme;   // Saturn
		planetLongitudes[11] = 5.481225395663 +    74.78165903077800 * jme;   // Uranus
		planetLongitudes[12] = 5.311897933164 +    38.13297222612500 * jme;   // Neptune
		planetLongitudes[13] =                      0.3595362285049309 * jme; // Pluto (mu)
		planetLongitudes[14] = 5.198466400630 + 77713.7714481804 * jme;       // Moon (D)
		planetLongitudes[15] = 1.627905136020 + 84334.6615717837 * jme;       // Moon (F)
		planetLongitudes[16] = 2.355555638750 + 83286.9142477147 * jme;       // Moon (l)

		final double a = evalSeries(EARTH_HELIOCENTRIC_DATA2.get(ResourceReader.VariableIndex.A), jme, planetLongitudes);
		final double l = evalSeries(EARTH_HELIOCENTRIC_DATA2.get(ResourceReader.VariableIndex.L), jme, planetLongitudes);
		final double k = evalSeries(EARTH_HELIOCENTRIC_DATA2.get(ResourceReader.VariableIndex.K), jme, planetLongitudes);
		final double h = evalSeries(EARTH_HELIOCENTRIC_DATA2.get(ResourceReader.VariableIndex.H), jme, planetLongitudes);
		final double q = evalSeries(EARTH_HELIOCENTRIC_DATA2.get(ResourceReader.VariableIndex.Q), jme, planetLongitudes);
		final double p = evalSeries(EARTH_HELIOCENTRIC_DATA2.get(ResourceReader.VariableIndex.P), jme, planetLongitudes);

		final double e = StrictMath.sqrt(k * k + h * h);
		final double w = StrictMath.atan2(h, k);
		final double n = StrictMath.atan2(p, q);
		final double i = 2. * StrictMath.asin(StrictMath.sqrt(q * q + p * p));
		final double mm = StrictMath.sqrt(8.9970116036316091182e-10 + 2.9591220836841438269e-04) / StrictMath.pow(a, 1.5);

		final OrbitalElements orbit = new OrbitalElements();
		orbit.t = tt * JulianDay.CIVIL_SAECULUM;
		orbit.eccentricity = e;
		orbit.semimajorAxis = a;
		orbit.inclination = i;
		orbit.longitudeAscendingNode = MathHelper.mod2pi(n);
		orbit.argumentOfPerihelion = MathHelper.mod2pi(w - n);
		orbit.meanAnomaly = MathHelper.mod2pi(l - w);
		orbit.meanMotion = mm;
		//1.5149480825765068E-4
		return -StrictMath.toDegrees(a);
	}

	private static double evalSeries(final List<ResourceReader.VSOP2013Data> elements, final double jme, final double[] planetLongitudes){
		double sum = 0.;
		final double[] parameters = new double[elements.size()];
		for(int i = 0; i < elements.size(); i ++){
			final ResourceReader.VSOP2013Data data = elements.get(i);
			for(final ResourceReader.VSOP2013Coeffs coeff : data.coeffs){
				double phi = 0.;
				for(int j = 0; j < 17; j ++)
					phi += coeff.iphi[j] * planetLongitudes[j];
				sum += coeff.sine * StrictMath.sin(phi) + coeff.cosine * StrictMath.cos(phi);
			}
			parameters[data.timePower] = sum;
		}
		return MathHelper.eval(jme, parameters);
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
		final double d = StrictMath.toRadians(meanElongationMoonSun(tt));
		final double m = StrictMath.toRadians(meanAnomalySun(tt));
		final double mp = StrictMath.toRadians(meanAnomalyMoon(tt));
		final double f = StrictMath.toRadians(argumentLatitudeMoon(tt));
		final double omega = StrictMath.toRadians(ascendingLongitudeMoon(tt));

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
		return new double[]{deltaPsi, deltaEpsilon};
	}

	//mean elongation of the Moon from the Sun [°]
	private static double meanElongationMoonSun(final double tt){
		return MathHelper.eval(tt, MOON_MEAN_ELONGATION_PARAMETERS);
	}

	//mean anomaly of the Sun [°]
	private static double meanAnomalySun(final double tt){
		return MathHelper.eval(tt, SUN_MEAN_ANOMALY_PARAMETERS);
	}

	//mean anomaly of the Moon [°]
	private static double meanAnomalyMoon(final double tt){
		return MathHelper.eval(tt, MOON_MEAN_ANOMALY_PARAMETERS);
	}

	//Moon's argument of Latitude [°]
	private static double argumentLatitudeMoon(final double tt){
		return MathHelper.eval(tt, MOON_ARGUMENT_OF_LATITUDE);
	}

	//Longitude of the ascending node of the Moon's mean orbit on the ecliptic measured from the mean equinox of the date [rad]
	private static double ascendingLongitudeMoon(final double tt){
		return MathHelper.eval(tt, MOON_LONGITUDE_ASCENDING_NODE);
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
		return MathHelper.eval(u, MEAN_ECLIPTIC_OBLIQUITY_PARAMETERS) / JulianDay.SECONDS_IN_HOUR;
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

}
