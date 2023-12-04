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

import io.github.mtrevisan.sunset.Zenith;
import io.github.mtrevisan.sunset.coordinates.AtmosphericModel;
import io.github.mtrevisan.sunset.coordinates.GeographicLocation;
import net.e175.klaus.solarpositioning.SPA;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.time.LocalDate;
import java.time.LocalDateTime;
import java.time.ZoneOffset;
import java.time.ZonedDateTime;
import java.time.format.DateTimeFormatter;


@SuppressWarnings("ALL")
class SolarEventCalculatorTest{

	@Test
	void testExampleSunriseTransitSet(){
		LocalDate time = LocalDate.of(2003, 10, 17);

		AtmosphericModel atmosphere = AtmosphericModel.create(1010., 10.);
		final GeographicLocation location = GeographicLocation.create(39.742476, -105.1786, 0.)
			.withAtmosphere(atmosphere);
		SolarEventCalculator calculator = SolarEventCalculator.create(location);
		SolarEvent solarEvent = calculator.solarEvent(time, 67., Zenith.OFFICIAL);

//		compare(solarEvent, SolarEvent.RegularDay.class, ZoneOffset.ofHours(-7), "2003-10-17T06:12:43-07:00", "2003-10-17T11:46:04-07:00", "2003-10-17T17:18:51-07:00");
		compare(solarEvent, SolarEvent.RegularDay.class, ZoneOffset.ofHours(-7), "2003-10-17T06:17:08.662-07:00", "2003-10-17T11:46:04.958-07:00", "2003-10-17T17:15:54.134-07:00");
	}

	@Test
	void testAllDay(){
		LocalDate time = LocalDate.of(2015, 6, 17);

		AtmosphericModel atmosphere = AtmosphericModel.create(1010., 10.);
		//location is Honningsvåg, Norway (near North Cape)
		final GeographicLocation location = GeographicLocation.create(70.978056, 25.974722, 0.)
			.withAtmosphere(atmosphere);
		SolarEventCalculator calculator = SolarEventCalculator.create(location);
		SolarEvent solarEvent = calculator.solarEvent(time, 0., Zenith.OFFICIAL);

		compare(solarEvent, SolarEvent.AlwaysDay.class, ZoneOffset.ofHours(2), null, "2015-06-17T12:16:55.717+02:00", null);
	}

	@Test
	void testAllNight(){
		LocalDate time = LocalDate.of(2015, 1, 17);

		AtmosphericModel atmosphere = AtmosphericModel.create(1010., 10.);
		//location is Honningsvåg, Norway (near North Cape)
		final GeographicLocation location = GeographicLocation.create(70.978056, 25.974722, 0.)
			.withAtmosphere(atmosphere);
		SolarEventCalculator calculator = SolarEventCalculator.create(location);
		SolarEvent solarEvent = calculator.solarEvent(time, 0., Zenith.OFFICIAL);

		compare(solarEvent, SolarEvent.AlwaysNight.class, ZoneOffset.ofHours(2), null, null, null);
	}

	@Test
	void testNonZeroSunriseTransitSet(){
		LocalDate time = LocalDate.of(2015, 6, 17);

		AtmosphericModel atmosphere = AtmosphericModel.create(1010., 10.);
		final GeographicLocation location = GeographicLocation.create(-36.8406, 174.74, 0.)
			.withAtmosphere(atmosphere);
		SolarEventCalculator calculator = SolarEventCalculator.create(location);
		SolarEvent solarEvent = calculator.solarEvent(time, 0., Zenith.OFFICIAL);

		//NOAA: 7:32, 12:21:41, 17:11
//		compare(solarEvent, SolarEvent.RegularDay.class, ZoneOffset.ofHours(12), "2015-06-17T07:32:26+12:00", "2015-06-17T12:21:46+12:00", "2015-06-17T17:11:03+12:00");
		compare(solarEvent, SolarEvent.RegularDay.class, ZoneOffset.ofHours(12), "2015-06-17T07:37:32.354+12:00", "2015-06-17T12:21:46.655+12:00", "2015-06-17T17:06:17.792+12:00");
	}

	@Test
	void testDSToffDayBerlin(){
		LocalDate time = LocalDate.of(2015, 10, 25);

		AtmosphericModel atmosphere = AtmosphericModel.create(1010., 10.);
		final GeographicLocation location = GeographicLocation.create(52.33, 13.3, 0.)
			.withAtmosphere(atmosphere);
		SolarEventCalculator calculator = SolarEventCalculator.create(location);
		SolarEvent solarEvent = calculator.solarEvent(time, 68., Zenith.OFFICIAL);

		//NOAA: 6:49, 11:50:53, 16:52
//		compare(solarEvent, SolarEvent.RegularDay.class, ZoneOffset.ofHours(1), "2015-10-25T06:49:02+01:00", "2015-10-25T11:50:55+01:00", "2015-10-25T16:51:59+01:00");
		compare(solarEvent, SolarEvent.RegularDay.class, ZoneOffset.ofHours(1), "2015-10-25T06:54:48.063+01:00", "2015-10-25T11:50:55.328+01:00", "2015-10-25T16:46:14.078+01:00");
	}

	@Test
	void testDSTonDayBerlin(){
		LocalDate time = LocalDate.of(2016, 3, 27);

		AtmosphericModel atmosphere = AtmosphericModel.create(1010., 10.);
		final GeographicLocation location = GeographicLocation.create(52.33, 13.3, 0.)
			.withAtmosphere(atmosphere);
		SolarEventCalculator calculator = SolarEventCalculator.create(location);
		SolarEvent solarEvent = calculator.solarEvent(time, 68., Zenith.OFFICIAL);

		//NOAA: 06:52, 13:12:01, 19:33
		compare(solarEvent, SolarEvent.RegularDay.class, ZoneOffset.ofHours(1), "2016-03-27T06:52:19+02:00", "2016-03-27T13:12:02+02:00", "2016-03-27T19:32:49+02:00");
	}

	@Test
	void testDSToffDayAuckland(){
		LocalDate time = LocalDate.of(2016, 4, 3);

		AtmosphericModel atmosphere = AtmosphericModel.create(1010., 10.);
		final GeographicLocation location = GeographicLocation.create(-36.84, 174.74, 0.)
			.withAtmosphere(atmosphere);
		SolarEventCalculator calculator = SolarEventCalculator.create(location);
		SolarEvent solarEvent = calculator.solarEvent(time, 68., Zenith.OFFICIAL);

		//NOAA: 06:36, same, 18:12
//		compare(solarEvent, SolarEvent.RegularDay.class, ZoneOffset.ofHours(12), "2016-04-03T06:36:09+12:00", "2016-04-03T12:24:19+12:00", "2016-04-03T18:11:55+12:00");
		compare(solarEvent, SolarEvent.RegularDay.class, ZoneOffset.ofHours(12), "2016-04-03T06:41:12.981+12:00", "2016-04-03T12:24:19.285+12:00", "2016-04-03T18:07:44.036+12:00");
	}

	@Test
	void testDSTonDayAuckland(){
		LocalDate time = LocalDate.of(2015, 9, 27);

		AtmosphericModel atmosphere = AtmosphericModel.create(1010., 10.);
		final GeographicLocation location = GeographicLocation.create(-36.84, 174.74, 0.)
			.withAtmosphere(atmosphere);
		SolarEventCalculator calculator = SolarEventCalculator.create(location);
		SolarEvent solarEvent = calculator.solarEvent(time, 68., Zenith.OFFICIAL);

		//NOAA: 07:04, 13:12:19, 19:21
		compare(solarEvent, SolarEvent.RegularDay.class, ZoneOffset.ofHours(12), "2015-09-27T07:04:14+13:00", "2015-09-27T13:12:17+13:00", "2015-09-27T19:20:56+13:00");
	}


	private static void compare(SolarEvent result, Class<?> refClass, ZoneOffset zoneOffset, String refSunrise, String refTransit, String refSunset){
		if(refClass != null)
			Assertions.assertEquals(refClass, result.getClass());
		if(result instanceof SolarEvent.RegularDay regularDay){
			Assertions.assertEquals(refSunrise, toZoned(regularDay.sunrise(), zoneOffset).toString());
			Assertions.assertEquals(refSunset, toZoned(regularDay.sunset(), zoneOffset).toString());
		}
		if(refTransit != null)
			Assertions.assertEquals(refTransit, toZoned(result.transit(), zoneOffset).toString());
	}
	private static ZonedDateTime toZoned(LocalDateTime local, ZoneOffset zoneOffset){
		ZonedDateTime date = ZonedDateTime.ofLocal(local.plusSeconds(zoneOffset.getTotalSeconds()), zoneOffset.normalized(), zoneOffset);
		if(date.getDayOfMonth() < local.getDayOfMonth())
			date = date.plusDays(1l);
		else if(date.getDayOfMonth() > local.getDayOfMonth())
			date = date.minusDays(1l);
		return date;
	}


	//	@Test
//	void localSunPosition() throws SolarEventException{
//		GeographicLocation location = GeographicLocation.create(39.742476, -105.1786, 1830.14);
//		SolarEventCalculator calc = SolarEventCalculator.create(location);
//		final AtmosphericModel atmosphericModel = AtmosphericModel.create(820, 11.);
//
//		double ut = JulianDate.of(2003, 10, 17)
//			+ JulianDate.timeOf(LocalTime.of(19, 30, 30));
//		double jd = TimeHelper.universalTimeToTerrestrialTime(ut, 67.);
//		double jce = JulianDate.centuryJ2000Of(jd);
//
//		EclipticCoordinate eclipticCoord = SunPosition.sunEclipticPosition(jd);
//		EquatorialCoordinate equatorialCoord = SunPosition.sunEquatorialPosition(eclipticCoord, jd);
//		double[] nutation = SunPosition.nutationCorrection(jce);
//		double meanEclipticObliquity = SunPosition.meanEclipticObliquity(jce);
//		double trueEclipticObliquity = SunPosition.trueEclipticObliquity(meanEclipticObliquity, nutation[1]);
//
//		double meanSiderealTime = TimeHelper.greenwichMeanSiderealTime(ut);
//		double apparentSiderealTime = TimeHelper.greenwichApparentSiderealTime(meanSiderealTime, trueEclipticObliquity, nutation[0]);
//		double localMeanSiderealTime = TimeHelper.localMeanSiderealTime(apparentSiderealTime, location);
//		double localHourAngle = TimeHelper.localHourAngle(localMeanSiderealTime, equatorialCoord.getRightAscension());
//
//		//compute the sun position (rigth ascension and declination) with respect to the observer local position at the Earth surface:
//		double equatorialHorizontalParallax = SunPosition.equatorialHorizontalParallax(eclipticCoord.getDistance());
//		double u = StrictMath.atan((1. - SunPosition.EARTH_FLATTENING) * StrictMath.tan(location.getLatitude()));
//		double height = location.getAltitude() / SunPosition.EARTH_EQUATORIAL_RADIUS;
//		double x = StrictMath.cos(u) + height * StrictMath.cos(location.getLatitude());
//		double y = 0.99664719 * StrictMath.sin(u) + height * StrictMath.sin(location.getLatitude());
//		double deltaRightAscension = StrictMath.atan2(
//			-x * StrictMath.sin(equatorialHorizontalParallax) * StrictMath.sin(localHourAngle),
//			StrictMath.cos(equatorialCoord.getDeclination()) - x * StrictMath.sin(equatorialHorizontalParallax)
//				* StrictMath.cos(localHourAngle)
//		);
//		//calculate the topocentric Sun Right Ascension: α'
//		double rightAscensionTopocentric = equatorialCoord.getRightAscension() + deltaRightAscension;
//		double declinationTopocentric = StrictMath.atan2(
//			(StrictMath.sin(equatorialCoord.getDeclination()) - y * StrictMath.sin(equatorialHorizontalParallax))
//				* StrictMath.cos(deltaRightAscension),
//			StrictMath.cos(equatorialCoord.getDeclination()) - y * StrictMath.sin(equatorialHorizontalParallax) * StrictMath.cos(localHourAngle)
//		);
//		//calculate the topocentric local hour angle: H’
//		double localHourAngleTopocentric = localHourAngle - deltaRightAscension;
//		//calculate the topocentric elevation angle without atmospheric refraction correction: e0
//		double trueElevation = StrictMath.asin(
//			StrictMath.sin(location.getLatitude()) * StrictMath.sin(declinationTopocentric)
//				+ StrictMath.cos(location.getLatitude()) * StrictMath.cos(declinationTopocentric)
//				* StrictMath.cos(localHourAngleTopocentric)
//		);
//		//calculate the atmospheric refraction correction: Δe
//		final double deltaElevation = atmosphericModel.atmosphericRefractionCorrection(trueElevation);
//		//calculate the topocentric elevation angle: e
//		double elevationTopocentric = trueElevation + deltaElevation;
//		//calculate the topocentric zenith angle: θ
//		double zenithTopocentric = StrictMath.PI / 2. - elevationTopocentric;
//		//calculate the topocentric astronomers azimuth angle (measured westward from south): Γ
//		double azimuthTopocentric = MathHelper.mod2pi(StrictMath.atan2(
//			StrictMath.sin(localHourAngleTopocentric),
//			StrictMath.cos(localHourAngleTopocentric) * StrictMath.sin(location.getLatitude())
//			- StrictMath.tan(declinationTopocentric) * StrictMath.cos(location.getLatitude())
//		));
//		//calculate the topocentric azimuth angle (measured westward from north): M
//		double azimuthTopocentricNavigators = MathHelper.mod2pi(azimuthTopocentric + StrictMath.PI);
//		AstroLib lib = new AstroLib();
//		AstroDay ad = new AstroDay();
//		lib.computeAstroDay(jd, ad);
//	}

//	@Test
//	void astronomicalSunset() throws SolarEventException{
//		GNSSLocation location = GNSSLocation.create(45.714920, 12.194179);
//		SolarEventCalculator calc = SolarEventCalculator.create(location);
//
//		LocalDateTime datetime = calc.sunset(LocalDate.of(2022, 12, 25), Zenith.ASTRONOMICAL);
//
//		assertTimeEquals("17:20:00", datetime);
//	}

//	@Test
//	void civilSunset() throws SolarEventException{
//		GNSSLocation location = GNSSLocation.create(45.714920, 12.194179);
//		SolarEventCalculator calc = SolarEventCalculator.create(location);
//
//		LocalDateTime datetime = calc.sunset(LocalDate.of(2022, 12, 25), Zenith.CIVIL);
//
//		assertTimeEquals("16:06:00", datetime);
//	}

//	@Test
//	void officialSunset() throws SolarEventException{
//		GNSSLocation location = GNSSLocation.create(45.714920, 12.194179);
//		SolarEventCalculator calc = SolarEventCalculator.create(location);
//
//		LocalDateTime datetime = calc.sunset(LocalDate.of(2022, 12, 25), Zenith.OFFICIAL);
//
//		assertTimeEquals("15:32:00", datetime);
//	}


	private static void assertTimeEquals(final String expectedTime, final LocalDateTime actualTime){
		int expectedSeconds = getSeconds(expectedTime);
		int actualSeconds = getSeconds(actualTime.toLocalTime().toString());
		Assertions.assertEquals(expectedSeconds, actualSeconds, 1,
			"Expected " + expectedTime + ", got " + actualTime.format(DateTimeFormatter.ofPattern("HH:MM:ss")));
	}

	private static int getSeconds(final String time){
		String[] timeParts = time.split("\\:");
		return 60 * (60 * Integer.valueOf(timeParts[0]) + Integer.valueOf(timeParts[1])) + Integer.valueOf(timeParts[2]);
	}

}
