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
import io.github.mtrevisan.sunset.SolarEventException;
import io.github.mtrevisan.sunset.TimeHelper;
import io.github.mtrevisan.sunset.coordinates.EclipticCoordinate;
import io.github.mtrevisan.sunset.coordinates.EquatorialCoordinate;
import io.github.mtrevisan.sunset.coordinates.GeographicLocation;
import io.github.mtrevisan.sunset.test.AstroDay;
import io.github.mtrevisan.sunset.test.AstroLib;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.time.LocalDateTime;
import java.time.LocalTime;
import java.time.format.DateTimeFormatter;


@SuppressWarnings("ALL")
class SolarEventCalculatorTest{

	@Test
	void localSunPosition() throws SolarEventException{
		GeographicLocation location = GeographicLocation.create(39.742476, -105.1786, 1830.14);
		SolarEventCalculator calc = SolarEventCalculator.create(location);
		final AtmosphericModel atmosphericModel = AtmosphericModel.create(820, 11.);

		double ut = JulianDay.of(2003, 10, 17)
			+ JulianDay.timeOf(LocalTime.of(19, 30, 30));
		double jd = TimeHelper.universalTimeToTerrestrialTime(ut, 67.);
		double jce = JulianDay.centuryJ2000Of(jd);

		EclipticCoordinate eclipticCoord = SunPosition.sunEclipticPosition(jd);
		EquatorialCoordinate equatorialCoord = SunPosition.sunEquatorialPosition(eclipticCoord, jd);
		double[] nutation = SunPosition.nutationCorrection(jce);
		double meanEclipticObliquity = SunPosition.meanEclipticObliquity(jce);
		double trueEclipticObliquity = SunPosition.trueEclipticObliquity(meanEclipticObliquity, nutation[1]);

		double meanSiderealTime = TimeHelper.meanSiderealTime(ut);
		double apparentSiderealTime = TimeHelper.apparentSiderealTime(meanSiderealTime, trueEclipticObliquity, nutation[0]);
		double localMeanSiderealTime = TimeHelper.localMeanSiderealTime(apparentSiderealTime, location);
		double localHourAngle = TimeHelper.localHourAngle(localMeanSiderealTime, equatorialCoord.getRightAscension());

		//compute the sun position (rigth ascension and declination) with respect to the observer local position at the Earth surface:
		double equatorialHorizontalParallax = calc.equatorialHorizontalParallax(eclipticCoord.getDistance());
		double u = StrictMath.atan((1. - SunPosition.EARTH_FLATTENING) * StrictMath.tan(location.getLatitude()));
		double height = location.getAltitude() / SunPosition.EARTH_EQUATORIAL_RADIUS;
		double x = StrictMath.cos(u) + height * StrictMath.cos(location.getLatitude());
		double y = 0.99664719 * StrictMath.sin(u) + height * StrictMath.sin(location.getLatitude());
		double deltaRightAscension = StrictMath.atan2(
			-x * StrictMath.sin(equatorialHorizontalParallax) * StrictMath.sin(localHourAngle),
			StrictMath.cos(equatorialCoord.getDeclination()) - x * StrictMath.sin(equatorialHorizontalParallax)
				* StrictMath.cos(localHourAngle)
		);
		//calculate the topocentric Sun Right Ascension: α'
		double rightAscensionTopocentric = equatorialCoord.getRightAscension() + deltaRightAscension;
		double declinationTopocentric = StrictMath.atan2(
			(StrictMath.sin(equatorialCoord.getDeclination()) - y * StrictMath.sin(equatorialHorizontalParallax))
				* StrictMath.cos(deltaRightAscension),
			StrictMath.cos(equatorialCoord.getDeclination()) - y * StrictMath.sin(equatorialHorizontalParallax) * StrictMath.cos(localHourAngle)
		);
		//calculate the topocentric local hour angle: H’
		double localHourAngleTopocentric = localHourAngle - deltaRightAscension;
		//calculate the topocentric elevation angle without atmospheric refraction correction: e0
		double trueElevation = StrictMath.asin(
			StrictMath.sin(location.getLatitude()) * StrictMath.sin(declinationTopocentric)
				+ StrictMath.cos(location.getLatitude()) * StrictMath.cos(declinationTopocentric)
				* StrictMath.cos(localHourAngleTopocentric)
		);
		//calculate the atmospheric refraction correction: Δe
		final double deltaElevation = atmosphericModel.atmosphericRefractionCorrection(trueElevation);
		//calculate the topocentric elevation angle: e
		double elevationTopocentric = trueElevation + deltaElevation;
		//calculate the topocentric zenith angle: θ
		double zenithTopocentric = StrictMath.PI / 2. - elevationTopocentric;
		//calculate the topocentric astronomers azimuth angle (measured westward from south): Γ
		double azimuthTopocentric = MathHelper.mod2pi(StrictMath.atan2(
			StrictMath.sin(localHourAngleTopocentric),
			StrictMath.cos(localHourAngleTopocentric) * StrictMath.sin(location.getLatitude())
			- StrictMath.tan(declinationTopocentric) * StrictMath.cos(location.getLatitude())
		));
		//calculate the topocentric azimuth angle (measured westward from north): M
		double azimuthTopocentricNavigators = MathHelper.mod2pi(azimuthTopocentric + StrictMath.PI);
		AstroLib lib = new AstroLib();
		AstroDay ad = new AstroDay();
		lib.computeAstroDay(jd, ad);
	}

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
