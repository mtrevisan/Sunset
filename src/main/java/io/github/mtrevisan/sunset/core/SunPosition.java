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

import io.github.mtrevisan.sunset.coordinates.AtmosphericModel;
import io.github.mtrevisan.sunset.JulianDate;
import io.github.mtrevisan.sunset.MathHelper;
import io.github.mtrevisan.sunset.ResourceReader;
import io.github.mtrevisan.sunset.TimeHelper;
import io.github.mtrevisan.sunset.coordinates.EclipticCoordinate;
import io.github.mtrevisan.sunset.coordinates.EquatorialCoordinate;
import io.github.mtrevisan.sunset.coordinates.GeographicLocation;
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

	private static final double[] OBLIQUITY_COEFFS = {
		84381.448, -4680.93, -1.55, 1999.25, 51.38, -249.67, -39.05, 7.12, 27.87, 5.79, 2.45
	};

	static final double EARTH_FLATTENING = 1. / 298.25642;
	//[m]
	static final double EARTH_EQUATORIAL_RADIUS = 6378136.6;
	private static final double[] EARTH_ORBIT_ECCENTRICITY = {0.016_708_634, -0.000_042_037, -0.000_000_126_7};

	//[m]
	public static final double SUN_EQUATORIAL_RADIUS = 6.95700e8;
	//[m]
	public static final double ASTRONOMICAL_UNIT = 1.495978707e11;

	private static final double[][][] TERMS_L = {
		{
			{175347046., 0., 0.},
			{3341656., 4.6692568, 6283.07585},
			{34894., 4.6261, 12566.1517},
			{3497., 2.7441, 5753.3849},
			{3418., 2.8289, 3.5231},
			{3136., 3.6277, 77713.7715},
			{2676., 4.4181, 7860.4194},
			{2343., 6.1352, 3930.2097},
			{1324., 0.7425, 11506.7698},
			{1273., 2.0371, 529.691},
			{1199., 1.1096, 1577.3435},
			{990., 5.233, 5884.927},
			{902., 2.045, 26.298},
			{857., 3.508, 398.149},
			{780., 1.179, 5223.694},
			{753., 2.533, 5507.553},
			{505., 4.583, 18849.228},
			{492., 4.205, 775.523},
			{357., 2.92, 0.067},
			{317., 5.849, 11790.629},
			{284., 1.899, 796.298},
			{271., 0.315, 10977.079},
			{243., 0.345, 5486.778},
			{206., 4.806, 2544.314},
			{205., 1.869, 5573.143},
			{202., 2.458, 6069.777},
			{156., 0.833, 213.299},
			{132., 3.411, 2942.463},
			{126., 1.083, 20.775},
			{115., 0.645, 0.98},
			{103., 0.636, 4694.003},
			{102., 0.976, 15720.839},
			{102., 4.267, 7.114},
			{99., 6.21, 2146.17},
			{98., 0.68, 155.42},
			{86., 5.98, 161000.69},
			{85., 1.3, 6275.96},
			{85., 3.67, 71430.7},
			{80., 1.81, 17260.15},
			{79., 3.04, 12036.46},
			{75., 1.76, 5088.63},
			{74., 3.5, 3154.69},
			{74., 4.68, 801.82},
			{70., 0.83, 9437.76},
			{62., 3.98, 8827.39},
			{61., 1.82, 7084.9},
			{57., 2.78, 6286.6},
			{56., 4.39, 14143.5},
			{56., 3.47, 6279.55},
			{52., 0.19, 12139.55},
			{52., 1.33, 1748.02},
			{51., 0.28, 5856.48},
			{49., 0.49, 1194.45},
			{41., 5.37, 8429.24},
			{41., 2.4, 19651.05},
			{39., 6.17, 10447.39},
			{37., 6.04, 10213.29},
			{37., 2.57, 1059.38},
			{36., 1.71, 2352.87},
			{36., 1.78, 6812.77},
			{33., 0.59, 17789.85},
			{30., 0.44, 83996.85},
			{30., 2.74, 1349.87},
			{25., 3.16, 4690.48}
		},
		{
			{628331966747., 0., 0.},
			{206059., 2.678235, 6283.07585},
			{4303., 2.6351, 12566.1517},
			{425., 1.59, 3.523},
			{119., 5.796, 26.298},
			{109., 2.966, 1577.344},
			{93., 2.59, 18849.23},
			{72., 1.14, 529.69},
			{68., 1.87, 398.15},
			{67., 4.41, 5507.55},
			{59., 2.89, 5223.69},
			{56., 2.17, 155.42},
			{45., 0.4, 796.3},
			{36., 0.47, 775.52},
			{29., 2.65, 7.11},
			{21., 5.34, 0.98},
			{19., 1.85, 5486.78},
			{19., 4.97, 213.3},
			{17., 2.99, 6275.96},
			{16., 0.03, 2544.31},
			{16., 1.43, 2146.17},
			{15., 1.21, 10977.08},
			{12., 2.83, 1748.02},
			{12., 3.26, 5088.63},
			{12., 5.27, 1194.45},
			{12., 2.08, 4694},
			{11., 0.77, 553.57},
			{10., 1.3, 6286.6},
			{10., 4.24, 1349.87},
			{9., 2.7, 242.73},
			{9., 5.64, 951.72},
			{8., 5.3, 2352.87},
			{6., 2.65, 9437.76},
			{6., 4.67, 4690.48}
		},
		{
			{52919., 0., 0.},
			{8720., 1.0721, 6283.0758},
			{309., 0.867, 12566.152},
			{27., 0.05, 3.52},
			{16., 5.19, 26.3},
			{16., 3.68, 155.42},
			{10., 0.76, 18849.23},
			{9., 2.06, 77713.77},
			{7., 0.83, 775.52},
			{5., 4.66, 1577.34},
			{4., 1.03, 7.11},
			{4., 3.44, 5573.14},
			{3., 5.14, 796.3},
			{3., 6.05, 5507.55},
			{3., 1.19, 242.73},
			{3., 6.12, 529.69},
			{3., 0.31, 398.15},
			{3., 2.28, 553.57},
			{2., 4.38, 5223.69},
			{2., 3.75, 0.98}
		},
		{
			{289., 5.844, 6283.076},
			{35., 0., 0.},
			{17., 5.49, 12566.15},
			{3., 5.2, 155.42},
			{1., 4.72, 3.52},
			{1., 5.3, 18849.23},
			{1., 5.97, 242.73}
		},
		{
			{114., 3.142, 0.},
			{8., 4.13, 6283.08},
			{1., 3.84, 12566.15}
		},
		{
			{1., 3.14, 0.}
		}
	};
	private static final double[][][] TERMS_B = {
		{
			{280., 3.199, 84334.662},
			{102., 5.422, 5507.553},
			{80., 3.88, 5223.69},
			{44., 3.7, 2352.87},
			{32., 4., 1577.34}
		},
		{
			{9., 3.9, 5507.55},
			{6., 1.73, 5223.69}
		}
	};
	private static final double[][][] TERMS_R = {
		{
			{100013989., 0., 0.},
			{1670700., 3.0984635, 6283.07585},
			{13956., 3.05525, 12566.1517},
			{3084., 5.1985, 77713.7715},
			{1628., 1.1739, 5753.3849},
			{1576., 2.8469, 7860.4194},
			{925., 5.453, 11506.77},
			{542., 4.564, 3930.21},
			{472., 3.661, 5884.927},
			{346., 0.964, 5507.553},
			{329., 5.9, 5223.694},
			{307., 0.299, 5573.143},
			{243., 4.273, 11790.629},
			{212., 5.847, 1577.344},
			{186., 5.022, 10977.079},
			{175., 3.012, 18849.228},
			{110., 5.055, 5486.778},
			{98., 0.89, 6069.78},
			{86., 5.69, 15720.84},
			{86., 1.27, 161000.69},
			{65., 0.27, 17260.15},
			{63., 0.92, 529.69},
			{57., 2.01, 83996.85},
			{56., 5.24, 71430.7},
			{49., 3.25, 2544.31},
			{47., 2.58, 775.52},
			{45., 5.54, 9437.76},
			{43., 6.01, 6275.96},
			{39., 5.36, 4694},
			{38., 2.39, 8827.39},
			{37., 0.83, 19651.05},
			{37., 4.9, 12139.55},
			{36., 1.67, 12036.46},
			{35., 1.84, 2942.46},
			{33., 0.24, 7084.9},
			{32., 0.18, 5088.63},
			{32., 1.78, 398.15},
			{28., 1.21, 6286.6},
			{28., 1.9, 6279.55},
			{26., 4.59, 10447.39}
		},
		{
			{103019., 1.10749, 6283.07585},
			{1721., 1.0644, 12566.1517},
			{702., 3.142, 0.},
			{32., 1.02, 18849.23},
			{31., 2.84, 5507.55},
			{25., 1.32, 5223.69},
			{18., 1.42, 1577.34},
			{10., 5.91, 10977.08},
			{9., 1.42, 6275.96},
			{9., 0.27, 5486.78}
		},
		{
			{4359., 5.7846, 6283.0758},
			{124., 5.579, 12566.152},
			{12., 3.14, 0.},
			{9., 3.63, 77713.77},
			{6., 1.87, 5573.14},
			{3., 5.47, 18849.23}
		},
		{
			{145., 4.273, 6283.076},
			{7., 3.92, 12566.15}
		},
		{
			{4., 2.56, 6283.08}
		}
	};

	private static final double[] SUN_GEOCENTRIC_MEAN_LONGITUDE_PARAMETERS = {280.466_46, 36_000.769_83, 0.000_303_2};
	private static final double[] SUN_GEOCENTRIC_MEAN_ANOMALY_PARAMETERS = {357.527_723_33, 35_999.050_34, -0.000_160_28, -0.000_003_33};
	private static final double[] SUN_EQUATION_OF_CENTER_1 = {1.914_602, -0.004_817, -0.000_014};
	private static final double[] SUN_EQUATION_OF_CENTER_2 = {0.019_993, -0.000_101};

	//http://vadimchazov.narod.ru/besmfore.htm
	private static final double[] MOON_MEAN_ELONGATION_PARAMETERS = {297.850_363_06, 445_267.111_48, -0.001_914_17, 0.000_005_28};
	private static final double[] MOON_MEAN_ANOMALY_PARAMETERS = {134.962_981_39, 477_198.867_398_055_6, 0.008_697_22, 0.000_017_78};
	private static final double[] MOON_ARGUMENT_OF_LATITUDE = {93.271_910_28, 483_202.017_538_055_5, -0.003_682_50, 0.000_003_06};
//	private static final double[] MOON_LONGITUDE_ASCENDING_NODE = {125.044_522_22, -1934.1362608, 0.002070833, 1. / 450000.};
	private static final double[] MOON_LONGITUDE_ASCENDING_NODE = {125.04, -1934.136};
	//https://www.aanda.org/articles/aa/pdf/2003/48/aa4068.pdf
	private static final double[] MEAN_ECLIPTIC_OBLIQUITY_PARAMETERS = {84381.448, -4680.93, -1.55, 1999.25, -51.38, -249.67, -39.05, 7.12,
		27.87, 5.79, 2.45};
	private static final double[] MOON_GEOCENTRIC_MEAN_LONGITUDE = {218.3165, 481_267.8813};

	//at J2000
	private static final double ABERRATION_CONSTANT = 20.49552;


	private SunPosition(){}


	//https://hal-mines-paristech.archives-ouvertes.fr/hal-00725987/document


	/**
	 * Calculated the Sun ecliptic position.
	 *
	 * @param jd	Julian Day of Terrestrial Time from J2000.0.
	 * @return	The Sun ecliptic position.
	 */
	public static EclipticCoordinate sunEclipticPosition(final double jd){
		final double jme = JulianDate.millenniumJ2000Of(jd);

		final double meanEclipticalLatitude = meanEclipticalLatitude(jme);
		final double meanEclipticalLongitude = meanEclipticalLongitude(jme);
		final double radiusVector = radiusVector(jme);
final double radiusVectorApprox = 0.016704 * StrictMath.cos(2. * StrictMath.PI * jd / 365.254902 + 3.091159) + 1.000140;
		return EclipticCoordinate.create(meanEclipticalLatitude, meanEclipticalLongitude, radiusVector);
	}

	public static EquatorialCoordinate sunEquatorialPosition(final double jme, final double deltaPsi,
			final double trueEclipticObliquity){
		//calculate Earth heliocentric latitude, <code>B</code> [deg]
		final double[] bTerms = calculateLBRTerms(jme, TERMS_B);
		final double earthHeliocentricLatitude = MathHelper.mod(StrictMath.toDegrees(calculateLBRPolynomial(jme, bTerms)), 360.);

		//calculate Earth radius vector, <code>R</code> [AU]
		final double[] rTerms = calculateLBRTerms(jme, TERMS_R);
		final double earthRadiusVector = calculateLBRPolynomial(jme, rTerms);
		assert earthRadiusVector != 0;

		//calculate Earth heliocentric longitude, <code>L</code> [deg]
		final double[] lTerms = calculateLBRTerms(jme, TERMS_L);
		final double earthHeliocentricLongitude = MathHelper.mod(StrictMath.toDegrees(calculateLBRPolynomial(jme, lTerms)), 360.);

		//calculate geocentric longitude, <code>theta</code> [deg]
		final double sunGeocentricLongitude = MathHelper.mod(earthHeliocentricLongitude + 180., 360.);

		//calculate geocentric latitude, <code>beta</code> [rad]
		final double sunGeocentricLatitude = StrictMath.toRadians(-earthHeliocentricLatitude);

		//calculate aberration correction, <code>delta tau</code>
		final double aberrationCorrection = -20.4898 / (JulianDate.SECONDS_PER_HOUR * earthRadiusVector);

		//calculate the apparent sun longitude, <code>lambda</code>
		//note: should be 0 for Spring equinox, 90 for Summer solstice, 180 for Autumn equinox, 270 for Winter solstice
		final double apparentSunLongitude = StrictMath.toRadians(sunGeocentricLongitude + deltaPsi + aberrationCorrection);

		//calculate the geocentric sun right ascension [deg]
		final double eps = StrictMath.toRadians(trueEclipticObliquity);
		final double rightAscension = calculateGeocentricSunRightAscension(sunGeocentricLatitude, eps, apparentSunLongitude);
		//calculate geocentric sun declination [deg]
		final double declination = StrictMath.toDegrees(calculateGeocentricSunDeclination(sunGeocentricLatitude, eps, apparentSunLongitude));

		return EquatorialCoordinate.create(rightAscension, declination);
	}

	private static double[] calculateLBRTerms(final double jme, final double[][][] termCoeffs){
		final double[] lbrTerms = {0., 0., 0., 0., 0., 0.};
		//L0, L1, ... Ln
		for(int i = 0; i < termCoeffs.length; i ++){
			double lbrSum = 0;
			//rows of each Li
			for(int v = 0; v < termCoeffs[i].length; v ++){
				//coefficients
				final double a = termCoeffs[i][v][0];
				final double b = termCoeffs[i][v][1];
				final double c = termCoeffs[i][v][2];
				lbrSum += a * StrictMath.cos(b + c * jme);
			}
			lbrTerms[i] = lbrSum;
		}
		return lbrTerms;
	}

	private static double calculateLBRPolynomial(final double jme, final double[] terms){
		return MathHelper.polynomial(jme, terms) / 1.e8;
	}

	private static double calculateGeocentricSunRightAscension(final double beta, final double epsilon, final double lambda){
		final double rightAscension = StrictMath.atan2(StrictMath.sin(lambda) * StrictMath.cos(epsilon)
			- StrictMath.tan(beta) * StrictMath.sin(epsilon), StrictMath.cos(lambda));
		return MathHelper.mod(StrictMath.toDegrees(rightAscension), 360.);
	}

	private static double calculateGeocentricSunDeclination(final double beta, final double epsilon, final double lambda){
		return StrictMath.asin(StrictMath.sin(beta) * StrictMath.cos(epsilon)
			+ StrictMath.cos(beta) * StrictMath.sin(epsilon) * StrictMath.sin(lambda));
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
		final double jce = JulianDate.centuryJ2000Of(jd);

		final double sunMeanAnomaly = SunPosition.geocentricMeanAnomaly(jce);
		//calculate the nutation in longitude and obliquity
		final double[] nutation = nutationCorrection(sunMeanAnomaly, jce);
		//calculate the aberration correction
		final double aberration = aberrationCorrection(eclipticCoord.getDistance());
		//calculate the apparent Sun longitude: Ltrue = L + C
		final double apparentGeocentricLongitude = apparentGeocentricLongitude(eclipticCoord.getLongitude(), nutation[0], aberration);

		//calculate the obliquity of the ecliptic (the inclination of the Earth’s equator with respect to the plane at which the Sun
		//and planets appear to move across the sky): ɛ0
		final double meanEclipticObliquity = meanEclipticObliquity(jce);
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
	public static HorizontalCoordinate sunTopocentricPosition(final GeographicLocation location, final AtmosphericModel atmosphericModel,
			final EclipticCoordinate eclipticCoord, final double jd){
		final LocalDateTime date = JulianDate.dateTimeOf(jd);
		final double dt = TimeHelper.deltaT(date.getYear());
		final double ut = TimeHelper.terrestrialTimeToUniversalTime(jd, dt);
		final double jce = JulianDate.centuryJ2000Of(jd);

		final EquatorialCoordinate equatorialCoord = sunEquatorialPosition(eclipticCoord, jd);
		final double sunMeanAnomaly = SunPosition.geocentricMeanAnomaly(jce);
		final double[] nutation = nutationCorrection(sunMeanAnomaly, jce);
		final double meanEclipticObliquity = meanEclipticObliquity(jce);
		final double trueEclipticObliquity = trueEclipticObliquity(meanEclipticObliquity, nutation[1]);

		final double meanSiderealTime = TimeHelper.greenwichMeanSiderealTime(ut);
		final double apparentSiderealTime = TimeHelper.greenwichApparentSiderealTime(meanSiderealTime, trueEclipticObliquity, nutation[0]);
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
	 * Calculate the ecliptical mean latitude of the Sun, referred to the mean equinox of the date, <code>β</code>.
	 *
	 * @param jme	Julian Ephemeris Millennium of Terrestrial Time from J2000.0.
	 * @return	The ecliptical mean latitude of the Sun [rad].
	 *
	 * @see <a href="https://squarewidget.com/solar-coordinates/">Solar coordinates</a>
	 */
	public static double meanEclipticalLatitude(final double jme){
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
			MathHelper.polynomial(jme, parameters) / 100_000_000.
		);
	}

	/**
	 * Calculate the ecliptical mean longitude of the Sun, referred to the mean equinox of the date, <code>L</code>.
	 *
	 * @param jme	Julian Ephemeris Millennium of Terrestrial Time from J2000.0.
	 * @return	The ecliptical mean longitude of the Sun [rad].
	 *
	 * @see <a href="https://squarewidget.com/solar-coordinates/">Solar coordinates</a>
	 */
	public static double meanEclipticalLongitude(final double jme){
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
			MathHelper.polynomial(jme, parameters) / 100_000_000. + StrictMath.PI
		);
	}

	/**
	 * Calculate the distance between the center of the Sun and the center of the Earth, <code>R</code>.
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
		return MathHelper.polynomial(jme, parameters) / 100_000_000.;
	}

	/**
	 * Calculate the Sun's equation of center, <code>C</code>.
	 *
	 * @param geocentricMeanAnomaly	The mean anomaly of the Sun [rad].
	 * @param tt	Julian Century of Terrestrial Time from J2000.0.
	 * @return	The Sun's equation of center [rad].
	 */
	public static double equationOfCenter(final double geocentricMeanAnomaly, final double tt){
		return Math.toRadians(
			MathHelper.polynomial(tt, SUN_EQUATION_OF_CENTER_1) * StrictMath.sin(geocentricMeanAnomaly)
				+ MathHelper.polynomial(tt, SUN_EQUATION_OF_CENTER_2) * StrictMath.sin(2. * geocentricMeanAnomaly)
				+ 0.000_289 * StrictMath.sin(3. * geocentricMeanAnomaly)
		);
	}

	/**
	 * Calculate the Sun's true longitude, <code>☉ = Ltrue</code>.
	 *
	 * @param meanLongitude	The geocentric mean longitude of the Sun [rad].
	 * @param equationOfCenter	The Sun's equation of center [rad].
	 * @return	The Sun's true longitude [rad].
	 */
	public static double trueLongitude(final double meanLongitude, final double equationOfCenter){
		return MathHelper.mod2pi(meanLongitude + equationOfCenter);
	}

	/**
	 * Calculate the Sun's true anomaly, <code>v</code>.
	 *
	 * @param meanAnomaly	The geocentric mean longitude of the Sun [rad].
	 * @param equationOfCenter	The Sun's equation of center [rad].
	 * @return	The Sun's true anomaly [rad].
	 */
	public static double trueAnomaly(final double meanAnomaly, final double equationOfCenter){
		return MathHelper.mod2pi(meanAnomaly + equationOfCenter);
	}

	/**
	 * Calculate the geocentric apparent longitude of the Sun, referred to the true equinox of the date, <code>λ = Lapp</code>.
	 *
	 * @param geocentricTrueLongitude	Suns geocentric true longitude [rad].
	 * @param ascendingLongitudeMoon	Longitude of the ascending node of the Moon's mean orbit on the ecliptic [rad].
	 * @param sunMeanAnomaly	The geocentric mean anomaly of the Sun [rad].
	 * @param jce	Julian Century of Terrestrial Time from J2000.0.
	 * @return	The apparent true longitude of the Sun [rad].
	 */
	public static double apparentLongitude(final double geocentricTrueLongitude, final double ascendingLongitudeMoon, final double sunMeanAnomaly, final double jce){
		final double nutationCorrection = nutationAndAberrationCorrection(ascendingLongitudeMoon, sunMeanAnomaly, jce);
		return geocentricTrueLongitude - nutationCorrection;
	}

	/**
	 * Calculate the Sun's apparent declination.
	 *
	 * @param meanEclipticObliquity	Mean obliquity of the ecliptic [rad].
	 * @param apparentLongitude	The Sun's apparent longitude [rad].
	 * @return	The Sun's apparent declination [rad].
	 */
	public static double apparentDeclination(final double meanEclipticObliquity, final double apparentLongitude){
		return StrictMath.asin(StrictMath.sin(meanEclipticObliquity) * StrictMath.sin(apparentLongitude));
	}

	private static double nutationAndAberrationCorrection(final double moonAscendingLongitude, final double sunMeanAnomaly, final double jce){
		final double deltaPsi = Math.toDegrees(nutationCorrection(sunMeanAnomaly, jce)[0]);

		//Sun's radius vector, r [AU]
		final double earthRadiusVector = SunPosition.radiusVector(jce / 10.);
		//[deg]
		final double aberration = 20.4898 / (earthRadiusVector * JulianDate.SECONDS_PER_HOUR);

		return StrictMath.toRadians(aberration + deltaPsi * StrictMath.sin(moonAscendingLongitude));
	}

	/**
	 * Calculate the geocentric mean longitude of the Moon, referred to the mean equinox of the date, <code>L'</code>.
	 *
	 * @param tt	Julian Century of Terrestrial Time from J2000.0.
	 * @return	The geocentric mean longitude of the Moon [rad].
	 */
	private static double moonGeocentricMeanLongitude(final double tt){
		return MathHelper.mod2pi(StrictMath.toRadians(MathHelper.polynomial(tt, MOON_GEOCENTRIC_MEAN_LONGITUDE)));
	}

	/**
	 * Calculate the eccentricity of Earth's orbit, <code>e</code>.
	 *
	 * @param jce	Julian Century of Terrestrial Time from J2000.0.
	 * @return	The eccentricity of Earth's orbit.
	 */
	public static double earthOrbitEccentricity(final double jce){
		return MathHelper.polynomial(jce, EARTH_ORBIT_ECCENTRICITY);
	}


	/**
	 * Calculate corrections of nutation in longitude (<code>∆ψ</code>) and obliquity (<code>∆ε</code>).
	 *
	 * @param sunMeanAnomaly	The geocentric mean anomaly of the Sun [rad].
	 * @param jce	Julian Century of Terrestrial Time from J2000.0.
	 * @return	An array where the first element is <code>∆ψ</code> [rad], and the second <code>∆ε</code> [rad].
	 *
	 * http://physics.uwyo.edu/~wiro/planet/nutate.c
	 * https://iopscience.iop.org/article/10.1086/375641/pdf
	 */
	public static double[] nutationCorrection(final double sunMeanAnomaly, final double jce){
		final double d = meanElongationMoonSun(jce);
		final double m = sunMeanAnomaly;
		final double mp = meanAnomalyMoon(jce);
		final double f = argumentLatitudeMoon(jce);
		final double omega = ascendingLongitudeMoon(jce);

		final Collection<Double[]> elements = NUTATION_DATA.get("coeffs");
		final double[] x = {d, m, mp, f, omega};
		double deltaPsi = 0.;
		double deltaEpsilon = 0.;
		for(final Double[] element : elements){
			final double parameter = xyTermSummation(x, element);
			deltaPsi += (element[x.length] + element[x.length + 1] * jce) * StrictMath.sin(parameter);
			deltaEpsilon += (element[x.length + 2] + element[x.length + 3] * jce) * StrictMath.cos(parameter);
		}
		return new double[]{
			StrictMath.toRadians(deltaPsi / (JulianDate.SECONDS_PER_HOUR * 10_000.)),
			StrictMath.toRadians(deltaEpsilon / (JulianDate.SECONDS_PER_HOUR * 10_000.))
		};
	}

	//mean elongation of the Moon from the Sun [rad]
	private static double meanElongationMoonSun(final double tt){
		return MathHelper.mod2pi(StrictMath.toRadians(
			MathHelper.polynomial(tt, MOON_MEAN_ELONGATION_PARAMETERS)
		));
	}

	/**
	 * Calculate the geocentric mean longitude of the Sun, <code>L0</code>.
	 *
	 * @param tdb	Julian Century of Terrestrial Time from J2000.0.
	 * @return	The geocentric mean longitude of the Sun [rad].
	 */
	public static double geocentricMeanLongitude(final double tdb){
		return MathHelper.mod2pi(StrictMath.toRadians(
			MathHelper.polynomial(tdb, SUN_GEOCENTRIC_MEAN_LONGITUDE_PARAMETERS)
		));
	}

	/**
	 * Calculate the geocentric mean anomaly of the Sun, <code>M</code>.
	 *
	 * @param tdb	Julian Century of Terrestrial Time from J2000.0.
	 * @return	The geocentric mean anomaly of the Sun [rad].
	 */
	public static double geocentricMeanAnomaly(final double tdb){
		return MathHelper.mod2pi(StrictMath.toRadians(
			MathHelper.polynomial(tdb, SUN_GEOCENTRIC_MEAN_ANOMALY_PARAMETERS)
		));
	}

	//mean anomaly of the Moon [rad]
	private static double meanAnomalyMoon(final double tt){
		return MathHelper.mod2pi(StrictMath.toRadians(
			MathHelper.polynomial(tt, MOON_MEAN_ANOMALY_PARAMETERS)
		));
	}

	//Moon's argument of Latitude [rad]
	private static double argumentLatitudeMoon(final double tt){
		return MathHelper.mod2pi(StrictMath.toRadians(
			MathHelper.polynomial(tt, MOON_ARGUMENT_OF_LATITUDE)
		));
	}

	/**
	 * Longitude of the ascending node of the Moon's mean orbit on the ecliptic measured from the mean equinox of the date, <code>☊</code>.
	 *
	 * @param tt	Julian Century of Terrestrial Time from J2000.0.
	 * @return	Longitude of the ascending node of the Moon's mean orbit on the ecliptic [rad].
	 */
	public static double ascendingLongitudeMoon(final double tt){
		return MathHelper.mod2pi(StrictMath.toRadians(
			MathHelper.polynomial(tt, MOON_LONGITUDE_ASCENDING_NODE)
		));
	}

	private static double xyTermSummation(final double[] x, final Double[] terms){
		double result = 0.;
		for(int j = 0; j < x.length; j ++)
			result += x[j] * terms[j];
		return result;
	}

	/**
	 * Calculate the correction for aberration, the annual aberration (<code>∆τ</code>).
	 * <p>
	 * As the Earth revolves around the Sun, it is moving at a velocity of approximately 29.78 km/s. The speed of light is approximately
	 * 300,000 km/s. In the special case where the Earth is moving perpendicularly to the direction of the star, the angle of displacement,
	 * would therefore be (in radians) the ratio of the two velocities, i.e. vₑ = 2 ⋅ π ⋅ 1 AU / (365.25 ⋅ d) = 29.75, and vₑ / c = 0.00009935
	 * rad, or about 20.5 arcseconds.
	 * <br/>
	 * This quantity is known as the constant of aberration, and is conventionally represented by κ. Its precise accepted value is 20".49552
	 * (at J2000).
	 * </p>
	 *
	 * @param earthRadiusVector	Earth radius vector [AU].
	 * @return	The correction for aberration [rad].
	 */
	public static double aberrationCorrection(final double earthRadiusVector){
		return StrictMath.toRadians(
			-ABERRATION_CONSTANT / (JulianDate.SECONDS_PER_HOUR * earthRadiusVector)
		);
	}

	/**
	 * Calculate the apparent longitude of the Sun, <code>Lapp = λ</code>.
	 *
	 * @param geocentricMeanLongitude	Sun's true geocentric longitude [rad].
	 * @param deltaPsi	Nutation in longitude [rad].
	 * @param deltaAberration	Aberration [rad].
	 * @return	Apparent longitude of the Sun [rad].
	 */
	public static double apparentGeocentricLongitude(final double geocentricMeanLongitude, final double deltaPsi, final double deltaAberration){
		return geocentricMeanLongitude + deltaPsi + deltaAberration;
	}

	/**
	 * Calculate the mean obliquity of the ecliptic, <code>ɛ0</code>.
	 *
	 * @param jce	Julian Century of Terrestrial Time from J2000.0.
	 * @return	Apparent longitude of the Sun [rad].
	 */
	public static double meanEclipticObliquity(final double jce){
		return MathHelper.polynomial(jce / 100., OBLIQUITY_COEFFS) / JulianDate.SECONDS_PER_HOUR;
	}

	/**
	 * Calculate the mean obliquity of the ecliptic, corrected for parallax, <code>ɛ'</code>.
	 *
	 * @param meanEclipticObliquity	Mean obliquity of the ecliptic [rad].
	 * @param tt	Julian Century of Terrestrial Time from J2000.0.
	 * @return	Apparent longitude of the Sun [rad].
	 */
	public static double apparentEclipticObliquity(final double meanEclipticObliquity, final double tt){
		return meanEclipticObliquity + StrictMath.toRadians(0.002_56 * StrictMath.cos(ascendingLongitudeMoon(tt)));
	}

	/**
	 * True obliquity of the ecliptic corrected for nutation, <code>ɛ</code>.
	 *
	 * @param meanEclipticObliquity	Obliquity of the ecliptic, corrected for parallax [deg].
	 * @param obliquityNutation	Corrections of nutation in obliquity (∆ε) [deg].
	 * @return	Apparent longitude of the Sun [deg].
	 */
	public static double trueEclipticObliquity(final double meanEclipticObliquity, final double obliquityNutation){
		return meanEclipticObliquity + obliquityNutation;
	}


	/**
	 * Calculate the equatorial horizontal parallax of the Sun, <code>ξ</code>.
	 *
	 * @param radiusVector	Radius vector of the Earth [AU].
	 * @return	The equatorial horizontal parallax of the Sun [rad].
	 */
	public static double equatorialHorizontalParallax(double radiusVector){
		return StrictMath.toRadians(8.794 / (JulianDate.SECONDS_PER_HOUR * radiusVector));
	}

}
