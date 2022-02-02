/*
 * Copyright (c) 2020-2022 Mauro Trevisan
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

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.time.LocalDate;
import java.time.LocalDateTime;
import java.time.LocalTime;
import java.time.temporal.ChronoUnit;
import java.time.temporal.TemporalAmount;


/*
https://github.com/mikereedell/sunrisesunsetlib-java
*/
public class SolarEventCalculator{

	private final Location location;


	/**
	 * Constructs a new instance using the given parameters.
	 *
	 * @param location	Location of the place where the solar event should be calculated from.
	 */
	public static SolarEventCalculator create(final Location location){
		return new SolarEventCalculator(location);
	}


	private SolarEventCalculator(final Location location){
		this.location = location;
	}


	/**
	 * Computes the sunset time for the given zenith at the given date.
	 *
	 * @param solarZenith	Enumeration corresponding to the type of sunset to compute.
	 * @param date	Date to compute the sunset for.
	 * @return	The sunset time or {@code null} for no sunset.
	 */
	public LocalDateTime computeSunsetCalendar(final Zenith solarZenith, final LocalDate date){
		//[hrs]
		final LocalTime localTime = computeSolarEventTime(solarZenith, date, false);
		return LocalDateTime.of(date, localTime);
	}

	private LocalTime computeSolarEventTime(final Zenith solarZenith, final LocalDate date, final boolean sunrise){
		final BigDecimal longitudeHour = getLongitudeHour(date, sunrise);

		final BigDecimal meanAnomaly = getMeanAnomaly(longitudeHour);
		final BigDecimal sunTrueLong = getSunTrueLongitude(meanAnomaly);
		final BigDecimal cosineSunLocalHour = getCosineSunLocalHour(sunTrueLong, solarZenith);
		if((cosineSunLocalHour.doubleValue() < -1.0) || (cosineSunLocalHour.doubleValue() > 1.0))
			return null;

		final BigDecimal sunLocalHour = getSunLocalHour(cosineSunLocalHour, sunrise);
		final BigDecimal localMeanTime = getLocalMeanTime(sunTrueLong, longitudeHour, sunLocalHour);
		final LocalTime localTime = getLocalTime(localMeanTime);

		return localTime;
	}

	/**
	 * Computes the longitude time.
	 *
	 * @return	Longitudinal time.
	 */
	private BigDecimal getLongitudeHour(final LocalDate date, final Boolean sunrise){
		int offset = 18;
		if(sunrise)
			offset = 6;

		final BigDecimal dividend = BigDecimal.valueOf(offset).subtract(getBaseLongitudeHour());
		final BigDecimal addend = divideBy(dividend, BigDecimal.valueOf(24));
		final BigDecimal longitudeHour = getDayOfYear(date).add(addend);
		return setScale(longitudeHour);
	}

	/**
	 * Computes the base longitude hour.
	 *
	 * @return	The longitude of the location of the solar event divided by 15 (deg/hour).
	 */
	private BigDecimal getBaseLongitudeHour(){
		return divideBy(location.getLongitude(), BigDecimal.valueOf(15));
	}

	private BigDecimal getDayOfYear(final LocalDate date){
		return new BigDecimal(date.getDayOfYear());
	}

	/**
	 * Computes the mean anomaly of the Sun.
	 *
	 * @return	The Suns mean anomaly [°].
	 */
	private BigDecimal getMeanAnomaly(final BigDecimal longitudeHour){
		final BigDecimal meanAnomaly = multiplyBy(new BigDecimal("0.9856"), longitudeHour).subtract(new BigDecimal("3.289"));
		return setScale(meanAnomaly);
	}

	/**
	 * Computes the true longitude of the Sun at the given location.
	 *
	 * @param meanAnomaly	The Suns mean anomaly [°].
	 * @return	The Suns true longitude [°].
	 */
	private BigDecimal getSunTrueLongitude(BigDecimal meanAnomaly){
		meanAnomaly = convertDegreesToRadians(meanAnomaly);
		final BigDecimal sinMeanAnomaly = new BigDecimal(Math.sin(meanAnomaly.doubleValue()));
		final BigDecimal sinDoubleMeanAnomaly = new BigDecimal(Math.sin(multiplyBy(meanAnomaly, BigDecimal.valueOf(2)).doubleValue()));

		final BigDecimal firstPart = meanAnomaly.add(multiplyBy(sinMeanAnomaly, new BigDecimal("1.916")));
		final BigDecimal secondPart = multiplyBy(sinDoubleMeanAnomaly, new BigDecimal("0.020")).add(new BigDecimal("282.634"));
		BigDecimal trueLongitude = firstPart.add(secondPart);
		if(trueLongitude.doubleValue() > 360.)
			trueLongitude = trueLongitude.subtract(BigDecimal.valueOf(360));

		return setScale(trueLongitude);
	}

	private BigDecimal getCosineSunLocalHour(final BigDecimal sunTrueLong, final Zenith zenith){
		final BigDecimal sinSunDeclination = getSinOfSunDeclination(sunTrueLong);
		final BigDecimal cosineSunDeclination = getCosOfSunDeclination(sinSunDeclination);

		final BigDecimal zenithInRadians = zenith.getRadians();
		final BigDecimal cosineZenith = BigDecimal.valueOf(Math.cos(zenithInRadians.doubleValue()));
		final BigDecimal sinLatitude = BigDecimal.valueOf(Math.sin(convertDegreesToRadians(location.getLatitude()).doubleValue()));
		final BigDecimal cosLatitude = BigDecimal.valueOf(Math.cos(convertDegreesToRadians(location.getLatitude()).doubleValue()));

		final BigDecimal sinDeclinationTimesSinLat = sinSunDeclination.multiply(sinLatitude);
		final BigDecimal dividend = cosineZenith.subtract(sinDeclinationTimesSinLat);
		final BigDecimal divisor = cosineSunDeclination.multiply(cosLatitude);

		return setScale(divideBy(dividend, divisor));
	}

	private BigDecimal getSinOfSunDeclination(final BigDecimal sunTrueLong){
		final BigDecimal sinTrueLongitude = BigDecimal.valueOf(Math.sin(convertDegreesToRadians(sunTrueLong).doubleValue()));
		final BigDecimal sinOfDeclination = sinTrueLongitude.multiply(new BigDecimal("0.39782"));
		return setScale(sinOfDeclination);
	}

	private BigDecimal getCosOfSunDeclination(final BigDecimal sinSunDeclination) {
		final BigDecimal arcSinOfSinDeclination = BigDecimal.valueOf(Math.asin(sinSunDeclination.doubleValue()));
		final BigDecimal cosDeclination = BigDecimal.valueOf(Math.cos(arcSinOfSinDeclination.doubleValue()));
		return setScale(cosDeclination);
	}

	private BigDecimal getSunLocalHour(final BigDecimal cosSunLocalHour, final Boolean sunrise){
		final BigDecimal arcCosineOfCosineHourAngle = getArcCosFor(cosSunLocalHour);
		BigDecimal localHour = convertRadiansToDegrees(arcCosineOfCosineHourAngle);
		if(sunrise)
			localHour = BigDecimal.valueOf(360).subtract(localHour);

		return divideBy(localHour, BigDecimal.valueOf(15));
	}

	private BigDecimal getLocalMeanTime(final BigDecimal sunTrueLong, final BigDecimal longHour, final BigDecimal sunLocalHour){
		final BigDecimal rightAscension = getRightAscension(sunTrueLong);
		final BigDecimal innerParens = longHour.multiply(new BigDecimal("0.06571"));
		BigDecimal localMeanTime = sunLocalHour.add(rightAscension).subtract(innerParens);
		localMeanTime = localMeanTime.subtract(new BigDecimal("6.622"));
		if(localMeanTime.doubleValue() < 0)
			localMeanTime = localMeanTime.add(BigDecimal.valueOf(24));
		else if(localMeanTime.doubleValue() > 24)
			localMeanTime = localMeanTime.subtract(BigDecimal.valueOf(24));
		return setScale(localMeanTime);
	}

	/**
	 * Computes the Suns right ascension, adjusting for the quadrant of the true longitude of the Sun and turning it
	 * into degree-hours.
	 *
	 * @param sunTrueLong	Suns true longitude [°].
	 * @return	Suns right ascension in degree-hours.
	 */
	private BigDecimal getRightAscension(final BigDecimal sunTrueLong){
		final BigDecimal tanL = new BigDecimal(Math.tan(convertDegreesToRadians(sunTrueLong).doubleValue()));

		final BigDecimal innerParens = multiplyBy(convertRadiansToDegrees(tanL), new BigDecimal("0.91764"));
		BigDecimal rightAscension = new BigDecimal(Math.atan(convertDegreesToRadians(innerParens).doubleValue()));
		rightAscension = setScale(convertRadiansToDegrees(rightAscension));
		if(rightAscension.doubleValue() < 0)
			rightAscension = rightAscension.add(BigDecimal.valueOf(360));
		else if(rightAscension.doubleValue() > 360)
			rightAscension = rightAscension.subtract(BigDecimal.valueOf(360));

		final BigDecimal ninety = BigDecimal.valueOf(90);
		BigDecimal longitudeQuadrant = sunTrueLong.divide(ninety, 0, RoundingMode.FLOOR);
		longitudeQuadrant = longitudeQuadrant.multiply(ninety);

		BigDecimal rightAscensionQuadrant = rightAscension.divide(ninety, 0, RoundingMode.FLOOR);
		rightAscensionQuadrant = rightAscensionQuadrant.multiply(ninety);

		final BigDecimal augend = longitudeQuadrant.subtract(rightAscensionQuadrant);
		return divideBy(rightAscension.add(augend), BigDecimal.valueOf(15));
	}

	private LocalTime getLocalTime(final BigDecimal localMeanTime){
		BigDecimal utcTime = localMeanTime.subtract(getBaseLongitudeHour());
		utcTime = adjustForDST(utcTime);
		return LocalTime.of(0, 0)
			.plus(utcTime.longValue() * 24 * 60 * 60, ChronoUnit.SECONDS);
	}

	private BigDecimal adjustForDST(final BigDecimal localMeanTime){
		BigDecimal localTime = localMeanTime;
		if(localTime.doubleValue() > 24.0)
			localTime = localTime.subtract(BigDecimal.valueOf(24));
		return localTime;
	}


	private BigDecimal multiplyBy(final BigDecimal multiplicand, final BigDecimal multiplier){
		return multiplicand.multiply(multiplier);
	}

	private BigDecimal divideBy(final BigDecimal dividend, final BigDecimal divisor){
		return dividend.divide(divisor, 4, RoundingMode.HALF_EVEN);
	}

	private BigDecimal setScale(final BigDecimal number){
		return number.setScale(4, RoundingMode.HALF_EVEN);
	}

	private BigDecimal getArcCosFor(final BigDecimal radians){
		final BigDecimal arcCosine = BigDecimal.valueOf(Math.acos(radians.doubleValue()));
		return setScale(arcCosine);
	}

	private BigDecimal convertDegreesToRadians(final BigDecimal degrees){
		return multiplyBy(degrees, BigDecimal.valueOf(Math.PI / 180.0));
	}

	private BigDecimal convertRadiansToDegrees(final BigDecimal radians){
		return multiplyBy(radians, new BigDecimal(180 / Math.PI));
	}

}
