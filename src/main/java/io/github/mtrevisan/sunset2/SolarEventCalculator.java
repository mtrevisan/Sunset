package io.github.mtrevisan.sunset2;

import io.github.mtrevisan.sunset.JulianDate;
import io.github.mtrevisan.sunset.MathHelper;
import io.github.mtrevisan.sunset2.coordinates.GeographicalLocation;

import java.time.LocalTime;
import java.time.ZoneId;
import java.time.ZonedDateTime;


//Uses the algorithm found at https://web.archive.org/web/20161202180207/http://williams.best.vwh.net/sunrise_sunset_algorithm.htm
//https://github.com/mikereedell/sunrisesunsetlib-java
public class SolarEventCalculator{

	private static final double[] SUN_GEOCENTRIC_MEAN_LONGITUDE = {280.46646, 36_000.769_83, 0.000_3032};
	private static final double[] SUN_GEOCENTRIC_MEAN_ANOMALY_PARAMETERS = {357.52911, 35_999.050_29, -0.000_1537};
	private static final double[] MOON_GEOCENTRIC_MEAN_LONGITUDE = {218.3165, 481_267.8813};


	final private GeographicalLocation location;
	final private ZoneId timeZone;


	public static void main(final String[] args){
		final double latitude = 45.651727;
		final double longitude = 12.212561;
		final double elevation = 16.;

		SolarEventCalculator solarEventCalculator = new SolarEventCalculator(GeographicalLocation.of(latitude, longitude, elevation),
			ZoneId.of("Europe/Rome"));
		ZonedDateTime officialSunset = solarEventCalculator.sunset(Zenith.CIVIL, ZonedDateTime.now());
		//17:03
		System.out.println(officialSunset.toLocalTime());
	}

	/**
	 * Constructs a new <code>SolarEventCalculator</code> using the given parameters.
	 *
	 * @param location	<code>Location</code> of the place where the solar event should be calculated from.
	 * @param timeZone	Time zone of the location parameter.
	 */
	public SolarEventCalculator(final GeographicalLocation location, final ZoneId timeZone){
		this.location = location;
		this.timeZone = timeZone;
	}

	/**
	 * Computes the sunrise time for the given zenith at the given date.
	 *
	 * @param solarZenith	<code>Zenith</code> enum corresponding to the type of sunrise to compute.
	 * @param date	Calendar object representing the date to compute the sunrise for.
	 * @return	The sunrise time as a calendar or null for no sunrise
	 */
	public ZonedDateTime sunrise(final Zenith solarZenith, final ZonedDateTime date){
		final double localTimeParam = computeSolarEventTime(solarZenith, date, true);
		return getLocalTimeAsCalendar(localTimeParam, date);
	}

	/**
	 * Computes the sunset time for the given zenith at the given date.
	 *
	 * @param solarZenith	<code>Zenith</code> enum corresponding to the type of sunset to compute.
	 * @param date	Calendar object representing the date to compute the sunset for.
	 * @return	The sunset time as a Calendar or null for no sunset.
	 */
	public ZonedDateTime sunset(final Zenith solarZenith, final ZonedDateTime date){
		final double localTimeParam = computeSolarEventTime(solarZenith, date, false);
		return getLocalTimeAsCalendar(localTimeParam, date);
	}

	private double computeSolarEventTime(final Zenith solarZenith, ZonedDateTime date, final boolean isSunrise){
		date = date.withZoneSameInstant(timeZone);
		final double longitudeHour = getLongitudeHour(date, isSunrise);
final double jd = JulianDate.of(date.toLocalDateTime());

		final double meanAnomaly = getSunMeanAnomaly(longitudeHour);
final double terrestrialTime = JulianDate.centuryJ2000Of(jd);
final double sunGeocentricMeanAnomaly = StrictMath.toDegrees(sunGeocentricMeanAnomaly(terrestrialTime));
System.out.println(meanAnomaly - sunGeocentricMeanAnomaly);
		final double sunTrueLong = getSunTrueLongitude(meanAnomaly);
		final double cosSunLocalHour = getCosSunLocalHour(sunTrueLong, solarZenith);
		if(cosSunLocalHour < -1. || cosSunLocalHour > 1.)
			return Double.NaN;

		final double sunLocalHour = getSunLocalHour(cosSunLocalHour, isSunrise);
		final double localMeanTime = getLocalMeanTime(sunTrueLong, longitudeHour, sunLocalHour);
		return getLocalTime(localMeanTime, date);
	}

	/**
	 * Computes the longitude time, <code>t</code> in the algorithm.
	 *
	 * @return	Longitude time [day].
	 */
	private double getLongitudeHour(final ZonedDateTime date, final Boolean isSunrise){
		final int offset = (isSunrise? 6: 18);
		return date.getDayOfYear() + (offset - location.getLongitude() / 15.) / 24.;
	}

	/**
	 * Computes the mean anomaly of the Sun, <code>M</code> in the algorithm.
	 *
	 * @return	The suns mean anomaly, <code>M</code> [deg].
	 */
	//https://www.astrouw.edu.pl/~jskowron/pracownia/praca/sunspot_answerbook_expl/expl-5.html
	private double getSunMeanAnomaly(final double longitudeHour){
		return 0.9856 * longitudeHour - 3.289;
	}

	/**
	 * Calculate the geocentric mean anomaly of the Sun, M.
	 *
	 * @param tt	Julian Century of Terrestrial Time from J2000.0.
	 * @return	The geocentric mean anomaly of the Sun [rad].
	 *
	 * https://squarewidget.com/solar-coordinates/
	 */
	private static double sunGeocentricMeanAnomaly(final double tt){
		return MathHelper.mod2pi(StrictMath.toRadians(MathHelper.polynomial(tt, SUN_GEOCENTRIC_MEAN_ANOMALY_PARAMETERS)));
	}

	/**
	 * Computes the true longitude of the sun, L in the algorithm, at the given location, adjusted to fit in
	 * the range [0-360].
	 *
	 * @param meanAnomaly	The suns mean anomaly [deg].
	 * @return	The suns true longitude [deg].
	 */
	private double getSunTrueLongitude(final double meanAnomaly){
		final double sinMeanAnomaly = Math.sin(Math.toRadians(meanAnomaly));
		final double sinDoubleMeanAnomaly = Math.sin(Math.toRadians(meanAnomaly) * 2.);

		double trueLongitude = meanAnomaly + sinMeanAnomaly * 1.916
			+ sinDoubleMeanAnomaly * 0.020 + 282.634;

		if(trueLongitude > 360.)
			trueLongitude = trueLongitude - 360.;
		return trueLongitude;
	}

	/**
	 * Computes the suns right ascension, RA in the algorithm, adjusting for the quadrant of L and turning it
	 * into degree-hours. Will be in the range [0,360].
	 *
	 * @param sunTrueLong	Suns true longitude [deg].
	 * @return	Suns right ascension in degree-hours.
	 */
	private double getRightAscension(final double sunTrueLong){
		final double tanL = Math.tan(Math.toRadians(sunTrueLong));

		final double innerParens = Math.toDegrees(tanL) * 0.91764;
		double rightAscension = Math.atan(Math.toRadians(innerParens));
		rightAscension = Math.toDegrees(rightAscension);

		if(rightAscension < 0.)
			rightAscension = rightAscension + 360.;
		else if(rightAscension > 360.)
			rightAscension = rightAscension - 360.;

		final double ninety = 90.;
		double longitudeQuadrant = Math.floor(sunTrueLong / ninety);
		longitudeQuadrant = longitudeQuadrant * ninety;

		double rightAscensionQuadrant = Math.floor(rightAscension / ninety);
		rightAscensionQuadrant = rightAscensionQuadrant * ninety;

		final double augend = longitudeQuadrant - rightAscensionQuadrant;
		return (rightAscension + augend) / 15.;
	}

	private double getCosSunLocalHour(final double sunTrueLong, final Zenith zenith){
		final double sinSunDeclination = Math.sin(Math.toRadians(sunTrueLong)) * 0.39782;
		final double cosSunDeclination = Math.cos(Math.asin(sinSunDeclination));

		//[rad]
		final double zenithAngle = Math.toRadians(zenith.getAngle());
		final double cosZenith = Math.sin(zenithAngle);
		final double sinLatitude = Math.sin(Math.toRadians(location.getLatitude()));
		final double cosLatitude = Math.cos(Math.toRadians(location.getLatitude()));

		final double sinDeclinationTimesSinLat = sinSunDeclination * sinLatitude;
		final double dividend = cosZenith - sinDeclinationTimesSinLat;
		final double divisor = cosSunDeclination * cosLatitude;

		return dividend / divisor;
	}

	private double getSunLocalHour(final double cosSunLocalHour, final Boolean isSunrise){
		final double arccosOfCosHourAngle = Math.acos(cosSunLocalHour);
		double localHour = Math.toDegrees(arccosOfCosHourAngle);
		if(isSunrise)
			localHour = 360. - localHour;
		return localHour / 15.;
	}

	private double getLocalMeanTime(final double sunTrueLong, final double longitudeHour, final double sunLocalHour){
		final double rightAscension = getRightAscension(sunTrueLong);
		final double innerParens = longitudeHour * 0.06571;
		double localMeanTime = sunLocalHour + rightAscension - innerParens;
		localMeanTime = localMeanTime - 6.622;

		if(localMeanTime < 0.)
			localMeanTime = localMeanTime + 24.;
		else if(localMeanTime > 24.)
			localMeanTime = localMeanTime - 24.;
		return localMeanTime;
	}

	private double getLocalTime(final double localMeanTime, final ZonedDateTime date){
		final double utcTime = localMeanTime - location.getLongitude() / 15.;
		final double utcOffSet = date.getOffset().getTotalSeconds() / JulianDate.SECONDS_PER_HOUR;
		final double utcOffSetTime = utcTime + utcOffSet;
		return adjustForDST(utcOffSetTime, date);
	}

	private double adjustForDST(final double localMeanTime, final ZonedDateTime date){
		double localTime = localMeanTime;
		if(timeZone.getRules().isDaylightSavings(date.toInstant()))
			localTime = localTime + 1.;
		if(localTime > 24.)
			localTime = localTime - 24.;
		return localTime;
	}

	/**
	 * Returns the local rise/set time in the form HH:MM.
	 *
	 * @param localTimeParam	Representation of the local rise/set time.
	 * @return	Calendar representation of the local time as a calendar, or <code>null</code> for none.
	 */
	protected ZonedDateTime getLocalTimeAsCalendar(final double localTimeParam, final ZonedDateTime date){
		if(Double.isNaN(localTimeParam))
			return null;

		//create a clone of the input calendar, so we get locale/timezone information
		ZonedDateTime resultTime = ZonedDateTime.ofInstant(date.toInstant(), date.getZone());

		double localTime = localTimeParam;
		if(localTime < 0.){
			localTime = localTime + 24.;
			resultTime = resultTime.minusDays(1);
		}
		final int intPart = (int)localTime;
		final double decimalPart = localTime - intPart;
		int hour = intPart;
		int minutes = (int)(decimalPart * 60.);
		if(minutes == 60){
			minutes = 0;
			hour ++;
		}
		if(hour == 24){
			hour = 0;
		}

		//set the local time
		return resultTime.with(LocalTime.of(hour, minutes, 0));
	}

}
