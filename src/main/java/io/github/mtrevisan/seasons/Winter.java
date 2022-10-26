package io.github.mtrevisan.seasons;

import io.github.mtrevisan.sunset.JulianDay;
import io.github.mtrevisan.sunset.MathHelper;
import io.github.mtrevisan.sunset.TimeHelper;
import io.github.mtrevisan.sunset.coordinates.EclipticCoordinate;
import io.github.mtrevisan.sunset.coordinates.EquatorialCoordinate;
import io.github.mtrevisan.sunset.coordinates.GNSSLocation;
import io.github.mtrevisan.sunset.core.SunPosition;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.time.LocalDateTime;
import java.time.LocalTime;
import java.time.format.DateTimeFormatter;
import java.util.Locale;


public class Winter{

	private static final double[] SUN_GEOCENTRIC_MEAN_LONGITUDE = {280.46646, 36_000.769_83, 0.000_3032};
	private static final double[] SUN_GEOCENTRIC_MEAN_ANOMALY_PARAMETERS = {357.52911, 35_999.050_29, -0.000_1537};
	private static final double[] MOON_GEOCENTRIC_MEAN_LONGITUDE = {218.3165, 481_267.8813};

	private static final double[] EARTH_WINTER_SOLSTICE = {2451900.05952, 365242.74049, 0.06223, -0.00823, 0.00032};

	//1 s
	private static final double TIME_PRECISION = 1. / JulianDay.SECONDS_IN_DAY;

	private static final double ABERRATION_CONSTANT = 20.49552;


	//https://www.agopax.it/Libri_astronomia/pdf/Astronomical%20Algorithms.pdf
	//https://gml.noaa.gov/grad/solcalc/calcdetails.html
	public static void main(final String[] args){
		final int year = 2022;
		final double longitude = 12.19415;

		final DecimalFormat decimalFormatter = (DecimalFormat)NumberFormat.getNumberInstance(Locale.US);

		LocalDateTime winterSolsticeDateTime = winterSolstice(year);
		System.out.println(DateTimeFormatter.ISO_LOCAL_DATE_TIME.format(winterSolsticeDateTime) + " UTC");

		winterSolsticeDateTime = winterSolstice(year, longitude);
		decimalFormatter.applyPattern("0.######");
		System.out.println(DateTimeFormatter.ISO_LOCAL_DATE_TIME.format(winterSolsticeDateTime)
			+ " @ " + decimalFormatter.format(longitude) + "°" + (longitude > 0.? " E": (longitude < 0.? " W": "")));

		sunset(winterSolsticeDateTime);
	}


	private static LocalTime sunset(final LocalDateTime dateTime){
		//pag 109
		final GNSSLocation location = GNSSLocation.create(45.65, 12.19, 100.);

		final double ut = JulianDay.of(dateTime);
		final double tt = JulianDay.centuryJ2000Of(ut);

		final EclipticCoordinate eclipticCoord = SunPosition.sunEclipticPosition(ut);
		//EclipticCoordinate eclipticCoord = SunPosition.sunEclipticPosition(ut);
		final EquatorialCoordinate equatorialCoord = SunPosition.sunEquatorialPosition(eclipticCoord, ut);
		//Sun's geometric mean longitude (referred to the mean equinox of the date)
		final double meanLongitudeSun = sunGeocentricMeanLongitude(tt);
//		final double meanLongitudeSun2 = SunPosition.geocentricMeanLongitude(tt);
		//Sun's mean anomaly
		final double meanAnomaly = sunGeocentricMeanAnomaly(tt);
		//final double meanAnomaly = SunPosition.meanAnomalySun(tt);
		//EquatorialCoordinate equatorialCoord = SunPosition.sunEquatorialPosition(eclipticCoord, ut);
		final double omega = sunGeocentricApparentLongitude(0., new double[]{0., 0.}, 0.);
		final double meanLongitudeMoon = moonGeocentricMeanLongitude(tt);
		final double[] nutation = {Math.toRadians((
			(-17.1996 - 0.01742 * tt) * StrictMath.sin(omega)
				+ (-1.3187 - 0.000_16 * tt) * StrictMath.sin(2. * meanLongitudeSun)
				+ (-0.2274 - 0.000_02 * tt) * StrictMath.sin(2. * meanLongitudeMoon)
				+ (0.2062 + 0.000_02 * tt) * StrictMath.sin(2. * omega)
				+ (0.1426 - 0.000_34 * tt) * StrictMath.sin(meanAnomaly)
		) / 3600.), 0.};
		//double[] nutation = SunPosition.nutationCorrection(tt);
		final double meanEclipticObliquity = SunPosition.meanEclipticObliquity(tt);
		final double trueEclipticObliquity = SunPosition.trueEclipticObliquity(meanEclipticObliquity, nutation[1]);

		final double meanSiderealTime = TimeHelper.meanSiderealTime(ut);
		final double apparentSiderealTime = TimeHelper.apparentSiderealTime(meanSiderealTime, trueEclipticObliquity, nutation[0]);
		final double localMeanSiderealTime = TimeHelper.localMeanSiderealTime(apparentSiderealTime, location);
		final double localHourAngle = TimeHelper.localHourAngle(localMeanSiderealTime, equatorialCoord.getRightAscension());

		return null;
	}/**/


	private static LocalDateTime winterSolstice(final int year){
		return winterSolstice(year, Double.NaN);
	}

	private static LocalDateTime winterSolstice(final int year, final double longitude){
		//calculate the approximate date
		double winterSolsticeTDB = winterSolsticeTDB(year);
		double tt = JulianDay.centuryJ2000Of(winterSolsticeTDB);

		//improve precision:
		double correction;
		do{
			//Sun's geometric mean longitude (referred to the mean equinox of the date)
			final double sunMeanLongitude = sunGeocentricMeanLongitude(tt);
			//FIXME WRONG
//			final double sunMeanLongitude2 = SunPosition.geocentricMeanLongitude(tt);
			//Sun's mean anomaly
			final double sunMeanAnomaly = SunPosition.geocentricMeanAnomaly(tt);
			//Sun's equation of center
			final double equationOfCenter = SunPosition.equationOfCenter(sunMeanAnomaly, tt);
			//eccentricity of the Earth's orbit
			final double eccentricity = SunPosition.earthOrbitEccentricity(tt);
			final double trueAnomaly = sunMeanAnomaly + equationOfCenter;
			//Sun's radius vector [AU]
			final double earthRadiusVector = earthRadiusVector(eccentricity, trueAnomaly);
			//FIXME WRONG
//			final double earthRadiusVector2 = SunPosition.radiusVector(tt);
			//[arcsec]
			final double aberration = ABERRATION_CONSTANT * 1.000_001_018 * (1. - eccentricity * eccentricity);
			final double sunTrueLongitude = sunMeanLongitude + equationOfCenter
				- Math.toRadians(
					//reduction to the FK5 system
					(0.09033
					//correction for aberration
					+ aberration / earthRadiusVector) / 3600.);
			final double meanLongitudeMoon = moonGeocentricMeanLongitude(tt);
			final double ascendingLongitudeMoon = SunPosition.ascendingLongitudeMoon(tt);
			final double[] nutationCorrection = {Math.toRadians((
				(-17.1996 - 0.01742 * tt) * StrictMath.sin(ascendingLongitudeMoon)
				+ (-1.3187 - 0.000_16 * tt) * StrictMath.sin(2. * sunMeanLongitude)
				+ (-0.2274 - 0.000_02 * tt) * StrictMath.sin(2. * meanLongitudeMoon)
				+ (0.2062 + 0.000_02 * tt) * StrictMath.sin(2. * ascendingLongitudeMoon)
				+ (0.1426 - 0.000_34 * tt) * StrictMath.sin(sunMeanAnomaly)
			) / 3600.), 0.};
			//final double[] nutationCorrection = SunPosition.nutationCorrection(tt);
			//Sun's apparent longitude, referred to the true equinox of the date
			final double apparentLongitudeSun = sunGeocentricApparentLongitude(sunTrueLongitude, nutationCorrection, ascendingLongitudeMoon);

			//[day]
			correction = 58. * StrictMath.sin(1.5 * Math.PI - apparentLongitudeSun);

			winterSolsticeTDB += correction;
			tt = JulianDay.centuryJ2000Of(winterSolsticeTDB);
		}while(correction > TIME_PRECISION);

		//TDB ~ TDT (the difference is at most 0.0017 s)
		//https://gssc.esa.int/navipedia/index.php/Transformations_between_Time_Systems#TAI_-_TDT.2C_TCG.2C_TT
		final double g = Math.toRadians(357.528 + 35999.05 * tt);
		final double winterSolsticeTDT = winterSolsticeTDB
			- 0.001658 * StrictMath.sin(g + 0.0167 * StrictMath.sin(g)) / JulianDay.SECONDS_IN_DAY;

		//transform TDT into UT
		//∆T [s]
		final double deltaT = TimeHelper.deltaT(year);
		final double winterSolsticeUT = TimeHelper.terrestrialTimeToUniversalTime(winterSolsticeTDT, deltaT);
		LocalDateTime winterSolsticeDateTime = JulianDay.dateTimeOf(winterSolsticeUT);

		if(!Double.isNaN(longitude)){
			//convert UT into local time
			final double winterSolsticeLocalTime = winterSolsticeUT + longitude / (JulianDay.HOURS_IN_DAY * JulianDay.DEGREES_PER_HOUR);
			winterSolsticeDateTime = JulianDay.dateTimeOf(winterSolsticeLocalTime);
		}

		return winterSolsticeDateTime;
	}

	/**
	 * Calculate the geocentric mean longitude of the Sun, referred to the mean equinox of the date, L.
	 *
	 * @param tt	Julian Century of Terrestrial Time from J2000.0.
	 * @return	The geocentric mean longitude of the Sun [rad].
	 *
	 * @see <a href="https://squarewidget.com/solar-coordinates/">Solar coordinates</a>
	 */
	private static double sunGeocentricMeanLongitude(final double tt){
		return MathHelper.mod2pi(StrictMath.toRadians(MathHelper.eval(tt, SUN_GEOCENTRIC_MEAN_LONGITUDE)));
	}

	/**
	 * Calculate the geocentric mean anomaly of the Sun, M.
	 *
	 * @param tt	Julian Century of Terrestrial Time from J2000.0.
	 * @return	The geocentric mean anomaly of the Sun [rad].
	 */
	private static double sunGeocentricMeanAnomaly(final double tt){
		return MathHelper.mod2pi(StrictMath.toRadians(MathHelper.eval(tt, SUN_GEOCENTRIC_MEAN_ANOMALY_PARAMETERS)));
	}

	/**
	 * Calculate the geocentric apparent longitude of the Sun, referred to the true equinox of the date.
	 *
	 * @param geocentricTrueLongitude	Suns geocentric true longitude [rad].
	 * @param nutationCorrection	Corrections of nutation in longitude (∆ψ) and obliquity (∆ε) [rad].
	 * @param ascendingLongitudeMoon	Longitude of the ascending node of the Moon's mean orbit on the ecliptic [rad].
	 * @return	The apparent true longitude of the Sun [rad].
	 */
	private static double sunGeocentricApparentLongitude(final double geocentricTrueLongitude, final double[] nutationCorrection,
			final double ascendingLongitudeMoon){
		return MathHelper.mod2pi(
			geocentricTrueLongitude
			+ nutationCorrection[0]
			- Math.toRadians(0.005_69 + 0.004_78 * StrictMath.sin(ascendingLongitudeMoon))
		);
	}

	/**
	 * Calculate the distance between the center of the Sun and the center of the Earth, R.
	 * <p>U.S. Naval Observatory function.</p>
	 *
	 * @param eccentricity	Eccentricity of the Earth's orbit.
	 * @param geocentricTrueAnomaly	The true anomaly of the Sun [rad].
	 * @return	Distance between the center of the Sun and the center of the Earth [AU].
	 */
	private static double earthRadiusVector(final double eccentricity, final double geocentricTrueAnomaly){
		return 1.000_001_018 * (1. - eccentricity * eccentricity) / (1. + eccentricity * StrictMath.cos(geocentricTrueAnomaly));
	}

	/**
	 * Calculate the geocentric mean longitude of the Moon, referred to the mean equinox of the date, L'.
	 *
	 * @param tt	Julian Century of Terrestrial Time from J2000.0.
	 * @return	The geocentric mean longitude of the Moon [rad].
	 */
	private static double moonGeocentricMeanLongitude(final double tt){
		return MathHelper.mod2pi(StrictMath.toRadians(MathHelper.eval(tt, MOON_GEOCENTRIC_MEAN_LONGITUDE)));
	}

	/**
	 * Instant of the "mean" solstice for the years [1000, 3000].
	 *
	 * @return	JD of Winter Solstice in Barycentric Dynamical Time.
	 */
	private static double winterSolsticeTDB(final int year){
		final double jde0 = MathHelper.eval((year - 2000.) /1000., EARTH_WINTER_SOLSTICE);

		final double tt = JulianDay.centuryJ2000Of(jde0);
		final double w = Math.toRadians(35999.373 * tt - 2.47);
		final double deltaLambda = 1. + 0.0334 * StrictMath.cos(w) + 0.0007 * StrictMath.cos(2. * w);
		return jde0 + (0.00001 * periodicTerms(tt)) / deltaLambda
			- (66. + (year - 2000.)) / JulianDay.SECONDS_IN_DAY;
	}

	private static double periodicTerms(final double t){
		double x = 485. * StrictMath.cos(Math.toRadians(324.96 + 1934.136 * t));
		x += 203. * StrictMath.cos(Math.toRadians(337.23 + 32964.467 * t));
		x += 199. * StrictMath.cos(Math.toRadians(342.08 + 20.186 * t));
		x += 182. * StrictMath.cos(Math.toRadians(27.85 + 445267.112 * t));
		x += 156. * StrictMath.cos(Math.toRadians(73.14 + 45036.886 * t));
		x += 136. * StrictMath.cos(Math.toRadians(171.52 + 22518.443 * t));
		x += 77. * StrictMath.cos(Math.toRadians(222.54 + 65928.934 * t));
		x += 74. * StrictMath.cos(Math.toRadians(296.72 + 3034.906 * t));
		x += 70. * StrictMath.cos(Math.toRadians(243.58 + 9037.513 * t));
		x += 58. * StrictMath.cos(Math.toRadians(119.81 + 33718.147 * t));
		x += 52. * StrictMath.cos(Math.toRadians(297.17 + 150.678 * t));
		x += 50. * StrictMath.cos(Math.toRadians(21.02 + 2281.226 * t));

		x += 45. * StrictMath.cos(Math.toRadians(247.54 + 29929.562 * t));
		x += 44. * StrictMath.cos(Math.toRadians(325.15 + 31555.956 * t));
		x += 29. * StrictMath.cos(Math.toRadians(60.93 + 4443.417 * t));
		x += 18. * StrictMath.cos(Math.toRadians(155.12 + 67555.328 * t));

		x += 17. * StrictMath.cos(Math.toRadians(288.79 + 4562.452 * t));
		x += 16. * StrictMath.cos(Math.toRadians(198.04 + 62894.029 * t));
		x += 14. * StrictMath.cos(Math.toRadians(199.76 + 31436.921 * t));
		x += 12. * StrictMath.cos(Math.toRadians(95.39 + 14577.848 * t));
		x += 12. * StrictMath.cos(Math.toRadians(287.11 + 31931.756 * t));
		x += 12. * StrictMath.cos(Math.toRadians(320.81 + 34777.259 * t));
		x += 9. * StrictMath.cos(Math.toRadians(227.73 + 1222.114 * t));
		x += 8. * StrictMath.cos(Math.toRadians(15.45 + 16859.074 * t));
		return x;
	}

}
