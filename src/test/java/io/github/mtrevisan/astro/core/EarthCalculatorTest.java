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

import io.github.mtrevisan.astro.coordinates.AtmosphericModel;
import io.github.mtrevisan.astro.coordinates.GeographicLocation;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.time.LocalDateTime;
import java.time.ZoneOffset;
import java.time.ZonedDateTime;
import java.time.format.DateTimeFormatter;


@SuppressWarnings("ALL")
class EarthCalculatorTest{

	@Test
	void testExampleSunriseTransitSet(){
		ZonedDateTime date = ZonedDateTime.of(2003, 10, 17, 0, 0, 0, 0,
			ZoneOffset.ofHours(-7));

		AtmosphericModel atmosphere = AtmosphericModel.create(1010., 10.);
		final GeographicLocation location = GeographicLocation.create(39.742476, -105.1786, 0.)
			.withAtmosphere(atmosphere);
		EarthCalculator calculator = EarthCalculator.create(location);
		SunlightPhase sunlightPhase = calculator.sunlightPhase(date, 67., Zenith.OFFICIAL);

		compare(sunlightPhase, SunlightPhase.RegularDay.class, date.getOffset(),
			"2003-10-17T06:12:43.931-07:00",
			"2003-10-17T11:46:04.959-07:00",
			"2003-10-17T17:18:51.207-07:00");
	}

	@Test
	void testAllDay(){
		ZonedDateTime date = ZonedDateTime.of(2015, 6, 17, 0, 0, 0, 0,
			ZoneOffset.ofHours(2));

		AtmosphericModel atmosphere = AtmosphericModel.create(1010., 10.);
		//location is Honningsvåg, Norway (near North Cape)
		final GeographicLocation location = GeographicLocation.create(70.978056, 25.974722, 0.)
			.withAtmosphere(atmosphere);
		EarthCalculator calculator = EarthCalculator.create(location);
		SunlightPhase sunlightPhase = calculator.sunlightPhase(date, 0., Zenith.OFFICIAL);

		compare(sunlightPhase, SunlightPhase.AlwaysDay.class, date.getOffset(),
			null,
			"2015-06-17T12:16:55.717+02:00",
			null);
	}

	@Test
	void testAllNight(){
		ZonedDateTime date = ZonedDateTime.of(2015, 1, 17, 0, 0, 0, 0,
			ZoneOffset.ofHours(2));

		AtmosphericModel atmosphere = AtmosphericModel.create(1010., 10.);
		//location is Honningsvåg, Norway (near North Cape)
		final GeographicLocation location = GeographicLocation.create(70.978056, 25.974722, 0.)
			.withAtmosphere(atmosphere);
		EarthCalculator calculator = EarthCalculator.create(location);
		SunlightPhase sunlightPhase = calculator.sunlightPhase(date, 0., Zenith.OFFICIAL);

		compare(sunlightPhase, SunlightPhase.AlwaysNight.class, date.getOffset(),
			null,
			null,
			null);
	}

	@Test
	void testNonZeroSunriseTransitSet(){
		ZonedDateTime date = ZonedDateTime.of(2015, 6, 17, 0, 0, 0, 0,
			ZoneOffset.ofHours(12));

		AtmosphericModel atmosphere = AtmosphericModel.create(1010., 10.);
		final GeographicLocation location = GeographicLocation.create(-36.8406, 174.74, 0.)
			.withAtmosphere(atmosphere);
		EarthCalculator calculator = EarthCalculator.create(location);
		SunlightPhase sunlightPhase = calculator.sunlightPhase(date, 0., Zenith.OFFICIAL);

		//NOAA: 7:32, 12:21:41, 17:11
		compare(sunlightPhase, SunlightPhase.RegularDay.class, date.getOffset(),
			"2015-06-17T07:32:29.115+12:00",
			"2015-06-17T12:21:46.651+12:00",
			"2015-06-17T17:11:01.320+12:00");
	}

	@Test
	void testDSToffDayBerlin(){
		ZonedDateTime date = ZonedDateTime.of(2015, 10, 25, 0, 0, 0, 0,
			ZoneOffset.ofHours(1));

		AtmosphericModel atmosphere = AtmosphericModel.create(1010., 10.);
		final GeographicLocation location = GeographicLocation.create(52.33, 13.3, 0.)
			.withAtmosphere(atmosphere);
		EarthCalculator calculator = EarthCalculator.create(location);
		SunlightPhase sunlightPhase = calculator.sunlightPhase(date, 68., Zenith.OFFICIAL);

		//NOAA: 6:49, 11:50:53, 16:52
		compare(sunlightPhase, SunlightPhase.RegularDay.class, date.getOffset(),
			"2015-10-25T06:49:02.979+01:00",
			"2015-10-25T11:50:55.328+01:00",
			"2015-10-25T16:51:59.301+01:00");
	}

	@Test
	void testDSTonDayBerlin(){
		ZonedDateTime date = ZonedDateTime.of(2016, 3, 27, 0, 0, 0, 0,
			ZoneOffset.ofHours(1));

		AtmosphericModel atmosphere = AtmosphericModel.create(1010., 10.);
		final GeographicLocation location = GeographicLocation.create(52.33, 13.3, 0.)
			.withAtmosphere(atmosphere);
		EarthCalculator calculator = EarthCalculator.create(location);
		SunlightPhase sunlightPhase = calculator.sunlightPhase(date, 68., Zenith.OFFICIAL);

		//NOAA: 06:52, 13:12:01, 19:33
		compare(sunlightPhase, SunlightPhase.RegularDay.class, date.getOffset(),
			"2016-03-27T05:52:19.922+01:00",
			"2016-03-27T12:12:02.218+01:00",
			"2016-03-27T18:32:48.972+01:00");
	}

	@Test
	void testDSToffDayAuckland(){
		ZonedDateTime date = ZonedDateTime.of(2016, 4, 3, 0, 0, 0, 0,
			ZoneOffset.ofHours(12));

		AtmosphericModel atmosphere = AtmosphericModel.create(1010., 10.);
		final GeographicLocation location = GeographicLocation.create(-36.84, 174.74, 0.)
			.withAtmosphere(atmosphere);
		EarthCalculator calculator = EarthCalculator.create(location);
		SunlightPhase sunlightPhase = calculator.sunlightPhase(date, 68., Zenith.OFFICIAL);

		//NOAA: 06:36, same, 18:12
		compare(sunlightPhase, SunlightPhase.RegularDay.class, date.getOffset(),
			"2016-04-03T06:36:10.302+12:00",
			"2016-04-03T12:24:19.284+12:00",
			"2016-04-03T18:11:54.701+12:00");
	}

	@Test
	void testDSTonDayAuckland(){
		ZonedDateTime date = ZonedDateTime.of(2015, 9, 27, 0, 0, 0, 0,
			ZoneOffset.ofHours(12));

		AtmosphericModel atmosphere = AtmosphericModel.create(1010., 10.);
		final GeographicLocation location = GeographicLocation.create(-36.84, 174.74, 0.)
			.withAtmosphere(atmosphere);
		EarthCalculator calculator = EarthCalculator.create(location);
		SunlightPhase sunlightPhase = calculator.sunlightPhase(date, 68., Zenith.OFFICIAL);

		//NOAA: 07:04, 13:12:19, 19:21
		compare(sunlightPhase, SunlightPhase.RegularDay.class, date.getOffset(),
			"2015-09-27T06:04:15.194+12:00",
			"2015-09-27T12:12:17.453+12:00",
			"2015-09-27T18:20:55.543+12:00");
	}


	private static void compare(SunlightPhase result, Class<?> refClass, ZoneOffset zoneOffset, String refSunrise, String refTransit, String refSunset){
		if(refClass != null)
			Assertions.assertEquals(refClass, result.getClass());
		if(result instanceof SunlightPhase.RegularDay regularDay){
			Assertions.assertEquals(refSunrise, regularDay.sunrise().toString());
			Assertions.assertEquals(refSunset, regularDay.sunset().toString());
		}
		if(refTransit != null)
			Assertions.assertEquals(refTransit, result.transit().toString());
	}


	//	@Test
//	void localSunPosition() throws SunlightPhaseException{
//		GeographicLocation location = GeographicLocation.create(39.742476, -105.1786, 1830.14);
//		EarthCalculator calc = EarthCalculator.create(location);
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
//	void astronomicalSunset() throws SunlightPhaseException{
//		GNSSLocation location = GNSSLocation.create(45.714920, 12.194179);
//		EarthCalculator calc = EarthCalculator.create(location);
//
//		LocalDateTime datetime = calc.sunset(LocalDate.of(2022, 12, 25), Zenith.ASTRONOMICAL);
//
//		assertTimeEquals("17:20:00", datetime);
//	}

//	@Test
//	void civilSunset() throws SunlightPhaseException{
//		GNSSLocation location = GNSSLocation.create(45.714920, 12.194179);
//		EarthCalculator calc = EarthCalculator.create(location);
//
//		LocalDateTime datetime = calc.sunset(LocalDate.of(2022, 12, 25), Zenith.CIVIL);
//
//		assertTimeEquals("16:06:00", datetime);
//	}

//	@Test
//	void officialSunset() throws SunlightPhaseException{
//		GNSSLocation location = GNSSLocation.create(45.714920, 12.194179);
//		EarthCalculator calc = EarthCalculator.create(location);
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
