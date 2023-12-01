package io.github.mtrevisan.seasons;

import io.github.mtrevisan.sunset.JulianDate;
import io.github.mtrevisan.sunset.MathHelper;
import io.github.mtrevisan.sunset.SolarEventException;
import io.github.mtrevisan.sunset.TimeHelper;
import io.github.mtrevisan.sunset.Zenith;
import io.github.mtrevisan.sunset.coordinates.GeographicLocation;
import io.github.mtrevisan.sunset.core.SolarEvent;
import io.github.mtrevisan.sunset.core.SolarEventCalculator;
import io.github.mtrevisan.sunset.core.SunPosition;
import net.e175.klaus.solarpositioning.SPA;
import net.e175.klaus.solarpositioning.SunriseResult;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.time.LocalDate;
import java.time.LocalDateTime;
import java.time.LocalTime;
import java.time.ZoneId;
import java.time.ZonedDateTime;
import java.time.format.DateTimeFormatter;
import java.util.Locale;


public class Winter{

	private static final double[] MOON_GEOCENTRIC_MEAN_LONGITUDE = {218.3165, 481_267.8813};

	private static final double[] EARTH_WINTER_SOLSTICE = {2451900.05952, 365242.74049, 0.06223, -0.00823, 0.00032};

	//1 s
	private static final double TIME_PRECISION = 1. / JulianDate.SECONDS_PER_DAY;

	private static final double ABERRATION_CONSTANT = 20.49552;


	//https://www.agopax.it/Libri_astronomia/pdf/Astronomical%20Algorithms.pdf
	//https://gml.noaa.gov/grad/solcalc/calcdetails.html
	//https://www.nrel.gov/docs/fy08osti/34302.pdf
	//https://www.sunearthtools.com/dp/tools/pos_sun.php
	public static void main(final String[] args) throws SolarEventException{
		final int year = 2023;
		final GeographicLocation location = GeographicLocation.create(45.714920, 12.194179, 16.);

		final DecimalFormat decimalFormatter = (DecimalFormat)NumberFormat.getNumberInstance(Locale.US);


		LocalDateTime winterSolstice = winterSolstice(year);
		System.out.println(DateTimeFormatter.ISO_LOCAL_DATE_TIME.format(winterSolstice) + " UTC");
winterSolstice = winterSolstice2(year);
System.out.println(DateTimeFormatter.ISO_LOCAL_DATE_TIME.format(winterSolstice) + " UTC");

		winterSolstice = winterSolstice(year, location);
		decimalFormatter.applyPattern("0.######");
		System.out.println(DateTimeFormatter.ISO_LOCAL_DATE_TIME.format(winterSolstice) + " local");
winterSolstice = winterSolstice2(year, location);
decimalFormatter.applyPattern("0.######");
System.out.println(DateTimeFormatter.ISO_LOCAL_DATE_TIME.format(winterSolstice) + " local");


		final double deltaT = TimeHelper.deltaT(winterSolstice.toLocalDate());
//ZonedDateTime date = winterSolstice
//	.atZone(ZoneId.of("Etc/UTC"));
//SunriseResult result = SPA.calculateSunriseTransitSet(date, 45.714920, 12.194179, deltaT, SPA.Horizon.CIVIL_TWILIGHT);
//if(result instanceof SunriseResult.RegularDay regular)
//	System.out.println("spa " + regular.sunset());
		final LocalDate winterSolsticeDate = winterSolstice.toLocalDate();
		SolarEventCalculator calculator = SolarEventCalculator.create(location);
		SolarEvent solarEvent = calculator.solarEvent(winterSolsticeDate, deltaT, Zenith.CIVIL);
		if(solarEvent instanceof SolarEvent.RegularDay event){
			System.out.println(DateTimeFormatter.ISO_LOCAL_TIME.format(event.sunset()) + " UTC");
		}
		System.out.println("16:04:23.347");

//		winterSolsticeSunset = sunset(winterSolsticeDate, location);
//		decimalFormatter.applyPattern("0.######");
//		System.out.println(DateTimeFormatter.ISO_LOCAL_TIME.format(winterSolsticeSunset) + " local");
	}



	public static LocalTime sunset(final LocalDate dateTime, final GeographicLocation location, final double deltaT, final Zenith zenith){
		LocalTime time = LocalTime.ofSecondOfDay(0);
		for(int i = 0; i < 2; i ++){
			final double ut = JulianDate.of(dateTime.atTime(time));
			time = sunsetInternal(ut, location, deltaT, zenith);
		}
		return time;
	}

	//https://github.com/buelowp/sunset/blob/master/src/sunset.cpp
	//http://www.internetsv.info/UniClock.html
	//https://www.nrel.gov/docs/fy08osti/34302.pdf
	private static LocalTime sunsetInternal(final double ut, final GeographicLocation location, final double deltaT, final Zenith zenith){
		double jce = JulianDate.centuryJ2000Of(ut);

		//Sun's geometric mean longitude (referred to the mean equinox of the date)
		final double sunMeanLongitude = SunPosition.geocentricMeanLongitude(jce);
		//Sun's mean anomaly
		final double sunMeanAnomaly = SunPosition.geocentricMeanAnomaly(jce);
		//eccentricity of the Earth's orbit
		final double earthEccentricity = SunPosition.earthOrbitEccentricity(jce);
		//Sun's equation of center
		final double equationOfCenter = SunPosition.equationOfCenter(sunMeanAnomaly, jce);
		final double sunTrueLongitude = SunPosition.trueLongitude(sunMeanLongitude, equationOfCenter);
		final double ascendingLongitudeMoon = SunPosition.ascendingLongitudeMoon(jce);
		//note: Sun's apparent longitude should be 0 for Spring equinox, 90 for Summer solstice, 180 for Autumn equinox, 270 for Winter solstice
		final double sunApparentLongitude = SunPosition.apparentLongitude(sunTrueLongitude, ascendingLongitudeMoon, sunMeanAnomaly, jce);

		final double meanEclipticObliquity = SunPosition.meanEclipticObliquity(jce);
		final double apparentEclipticObliquity = SunPosition.apparentEclipticObliquity(meanEclipticObliquity, jce);

		final double sunApparentDeclination = SunPosition.apparentDeclination(meanEclipticObliquity, sunApparentLongitude);
		final double equationOfTime = SolarEventCalculator.equationOfTime(sunMeanLongitude, sunMeanAnomaly, apparentEclipticObliquity, earthEccentricity);
		//[min]
		final double trueSolarTime = MathHelper.mod((ut - Math.floor(ut)) * 1440. + equationOfTime + 4. * location.getLongitude(), 1440.);
		//[deg]
		final double hourAngle = (trueSolarTime < 0.? trueSolarTime / 4. + 180.: trueSolarTime / 4. - 180.);

		//[day]
		final double sunriseHourAngle = SolarEventCalculator.sunriseHourAngle(location.getLatitude(), sunApparentDeclination, hourAngle)
			/ (2. * StrictMath.PI);
		//difference between the local solar time and the mean solar time
		final double timeOffset = equationOfTime + 4. * location.getLongitude();
		//solar noon in Local Sidereal Time [day]
//11:09:24
		final double noon = (60. * 12. - timeOffset) / JulianDate.MINUTES_PER_DAY;
		//sunrise/sunset in Local Sidereal Time:
		final double sunrise = noon - sunriseHourAngle;
//15:30:00
		final double sunset = noon + sunriseHourAngle;
		return LocalTime.ofSecondOfDay(Math.round(sunset * JulianDate.SECONDS_PER_DAY));
	}


	public static LocalDateTime winterSolstice2(final int year){
		return winterSolstice2(year, null);
	}

	public static LocalDateTime winterSolstice2(final int year, final GeographicLocation location){
		//calculate the approximate date
		double winterSolsticeTDB = winterSolsticeTDB(year);

		//improve precision:
		double jce;
		double correction = 0.;
		do{
			winterSolsticeTDB += correction;
			jce = JulianDate.centuryJ2000Of(winterSolsticeTDB);

			//Sun's geometric mean longitude (referred to the mean equinox of the date)
			final double sunMeanLongitude = SunPosition.geocentricMeanLongitude(jce);
			//Sun's mean anomaly
			final double sunMeanAnomaly = SunPosition.geocentricMeanAnomaly(jce);
			//Sun's equation of center
			final double equationOfCenter = SunPosition.equationOfCenter(sunMeanAnomaly, jce);
			final double sunTrueLongitude = SunPosition.trueLongitude(sunMeanLongitude, equationOfCenter);
			final double ascendingLongitudeMoon = SunPosition.ascendingLongitudeMoon(jce);
			//note: Sun's apparent longitude should be 0 for Spring equinox, 90 for Summer solstice, 180 for Autumn equinox, 270 for Winter solstice
			final double sunApparentLongitude = SunPosition.apparentLongitude(sunTrueLongitude, ascendingLongitudeMoon, sunMeanAnomaly, jce);


			//[day]
			correction = 58. * StrictMath.sin(1.5 * Math.PI - sunApparentLongitude);
		}while(Math.abs(correction) > TIME_PRECISION);

		//TDB ~ TDT (the difference is at most 0.0017 s)
		//https://gssc.esa.int/navipedia/index.php/Transformations_between_Time_Systems#TAI_-_TDT.2C_TCG.2C_TT
		final double g = Math.toRadians(357.528 + 35999.05 * jce);
		final double winterSolsticeTDT = winterSolsticeTDB
			- 0.001658 * StrictMath.sin(g + 0.0167 * StrictMath.sin(g)) / JulianDate.SECONDS_PER_DAY;

		//transform TDT into UT
		//∆T [s]
		final double deltaT = TimeHelper.deltaT(year);
		final double winterSolsticeUT = TimeHelper.terrestrialTimeToUniversalTime(winterSolsticeTDT, deltaT);
		LocalDateTime winterSolsticeDateTime = JulianDate.dateTimeOf(winterSolsticeUT);

		if(location != null){
			//convert UT into local time
			final double winterSolsticeLocalTime = winterSolsticeUT + location.getLongitude() / (JulianDate.HOURS_PER_DAY * JulianDate.DEGREES_PER_HOUR);
			winterSolsticeDateTime = JulianDate.dateTimeOf(winterSolsticeLocalTime);
		}

		return winterSolsticeDateTime;
	}

	public static LocalDateTime winterSolstice(final int year){
		return winterSolstice(year, null);
	}

	public static LocalDateTime winterSolstice(final int year, final GeographicLocation location){
		//calculate the approximate date
		double winterSolsticeTDB = winterSolsticeTDB(year);
//winterSolsticeTDB = 2437837.38589;
//winterSolsticeTDB = 2448908.5;

		//improve precision:
		double jce;
		double correction = 0.;
		do{
			winterSolsticeTDB += correction;
			jce = JulianDate.centuryJ2000Of(winterSolsticeTDB);
			final double jme = jce / 10.;

			//Sun's geometric mean longitude (referred to the mean equinox of the date), L0
			final double meanLongitude = SunPosition.geocentricMeanLongitude(jce);
			//etc etc...
			//Sun's mean anomaly, M
			final double meanAnomaly = SunPosition.geocentricMeanAnomaly(jce);
			//eccentricity of the Earth's orbit, e
			final double eccentricity = SunPosition.earthOrbitEccentricity(jce);
			//Sun's equation of center, C
			final double equationOfCenter = SunPosition.equationOfCenter(meanAnomaly, jce);
			//Sun's true geometric longitude, ☉ = Ltrue
			final double trueLongitude = MathHelper.mod2pi(meanLongitude + equationOfCenter);
			//Sun's true geometric anomaly, ν
			final double trueAnomaly = MathHelper.mod2pi(meanAnomaly + equationOfCenter);
			//Sun's radius vector, r
			final double earthRadiusVector = SunPosition.radiusVector(jme);

			//longitude of the ascending node of the Moon's mean orbit on the ecliptic, ☊
			final double moonAscendingLongitude = SunPosition.ascendingLongitudeMoon(jce);
			final double meanLongitudeMoon = moonGeocentricMeanLongitude(jce);
			final double[] nutationCorrection = {Math.toRadians((
				(-17.1996 - 0.01742 * jce) * StrictMath.sin(moonAscendingLongitude)
					+ (-1.3187 - 0.000_16 * jce) * StrictMath.sin(2. * meanLongitude)
					+ (-0.2274 - 0.000_02 * jce) * StrictMath.sin(2. * meanLongitudeMoon)
					+ (0.2062 + 0.000_02 * jce) * StrictMath.sin(2. * moonAscendingLongitude)
					+ (0.1426 - 0.000_34 * jce) * StrictMath.sin(meanAnomaly)
			) / JulianDate.SECONDS_PER_HOUR), 0.};
			//final double[] nutationCorrection = SunPosition.nutationCorrection(tt);

			final double aberration = Math.toRadians(ABERRATION_CONSTANT * 1.000_001_018 * (1. - eccentricity * eccentricity) / JulianDate.SECONDS_PER_HOUR);
			final double sunMeanAnomaly = SunPosition.geocentricMeanAnomaly(jce);
			final double sunTrueAnomaly = sunMeanAnomaly + equationOfCenter;
			//Sun's apparent longitude, referred to the true equinox of the date, λ = Lapp
			final double apparentLongitude = MathHelper.mod2pi(
				SunPosition.apparentLongitude(trueLongitude, moonAscendingLongitude, sunMeanAnomaly, jce)
				//correction for nutation in longitude, ∆ψ
				- nutationCorrection[0]
				//reduction to the FK5 system
				- Math.toRadians(0.090_33 / JulianDate.SECONDS_PER_HOUR)
				//correction for aberration
				- aberration / earthRadiusVector);

			//mean obliquity of the ecliptic, ɛ0
			final double meanEclipticObliquity = SunPosition.meanEclipticObliquity(jce);
			//mean obliquity of the ecliptic, ɛ = ɛ0 + Δɛ
			final double trueEclipticObliquity = meanEclipticObliquity + nutationCorrection[1];

			//[day]
			correction = 58. * StrictMath.sin(1.5 * Math.PI - apparentLongitude);
		}while(Math.abs(correction) > TIME_PRECISION);

		//TDB ~ TDT (the difference is at most 0.0017 s)
		//https://gssc.esa.int/navipedia/index.php/Transformations_between_Time_Systems#TAI_-_TDT.2C_TCG.2C_TT
		final double g = Math.toRadians(357.528 + 35999.05 * jce);
		final double winterSolsticeTDT = winterSolsticeTDB
			- 0.001658 * StrictMath.sin(g + 0.0167 * StrictMath.sin(g)) / JulianDate.SECONDS_PER_DAY;

		//transform TDT into UT
		//∆T [s]
		final double deltaT = TimeHelper.deltaT(year);
		final double winterSolsticeUT = TimeHelper.terrestrialTimeToUniversalTime(winterSolsticeTDT, deltaT);
		LocalDateTime winterSolsticeDateTime = JulianDate.dateTimeOf(winterSolsticeUT);

		if(location != null){
			//convert UT into local time
			final double winterSolsticeLocalTime = winterSolsticeUT + location.getLongitude() / (JulianDate.HOURS_PER_DAY * JulianDate.DEGREES_PER_HOUR);
			winterSolsticeDateTime = JulianDate.dateTimeOf(winterSolsticeLocalTime);
		}

		return winterSolsticeDateTime;
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
		return MathHelper.mod2pi(StrictMath.toRadians(MathHelper.polynomial(tt, MOON_GEOCENTRIC_MEAN_LONGITUDE)));
	}

	/**
	 * Instant of the "mean" solstice for the years [1000, 3000].
	 *
	 * @return	JD of Winter Solstice in Barycentric Dynamical Time.
	 */
	private static double winterSolsticeTDB(final int year){
		final double jde0 = MathHelper.polynomial((year - 2000.) /1000., EARTH_WINTER_SOLSTICE);

		final double jce = JulianDate.centuryJ2000Of(jde0);
		final double w = Math.toRadians(35999.373 * jce - 2.47);
		final double deltaLambda = 1. + 0.0334 * StrictMath.cos(w) + 0.0007 * StrictMath.cos(2. * w);
		return jde0 + (0.00001 * periodicTerms(jce)) / deltaLambda
			- (66. + (year - 2000.)) / JulianDate.SECONDS_PER_DAY;
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
