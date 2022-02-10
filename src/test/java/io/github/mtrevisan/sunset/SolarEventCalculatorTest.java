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
import io.github.mtrevisan.sunset.coordinates.GNSSLocation;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.time.LocalDate;
import java.time.LocalDateTime;
import java.time.LocalTime;
import java.time.format.DateTimeFormatter;


@SuppressWarnings("ALL")
class SolarEventCalculatorTest{

	@Test
	void localSunPosition() throws SolarEventException{
		GNSSLocation location = GNSSLocation.create(39.742476, -105.1786);
		SolarEventCalculator calc = SolarEventCalculator.create(location)
			.withObserverElevation(1830.14)
			.withPressure(820.)
			.withTemperature(11.)
			.withSurfaceSlopeAndAzimuthRotation(30., -10.);

		double ut = JulianDay.of(2003, 10, 17)
			+ JulianDay.timeOf(LocalTime.of(19, 30, 30));
		double jd = TimeHelper.universalTimeToTerrestrialTime(ut, 67.);
		double tt = JulianDay.centuryJ2000Of(jd);

		EclipticCoordinate eclipticCoord = SunPosition.sunEclipticPosition(jd);
		EquatorialCoordinate coord = SunPosition.sunEquatorialPosition(eclipticCoord, jd);
		double[] nutation = SunPosition.nutationCorrection(tt);
		double meanEclipticObliquity = SunPosition.meanEclipticObliquity(tt);
		double trueEclipticObliquity = SunPosition.trueEclipticObliquity(meanEclipticObliquity, nutation[1]);

		double meanSiderealTime = TimeHelper.meanSiderealTime(ut);
		double apparentSiderealTime = TimeHelper.apparentSiderealTime(meanSiderealTime, trueEclipticObliquity, nutation[0]);
		double localMeanSiderealTime = TimeHelper.localMeanSiderealTime(apparentSiderealTime, location);
		double localHourAngle = TimeHelper.localHourAngle(localMeanSiderealTime, coord.getRightAscension());
		if(Math.abs(localHourAngle - 11.105902) > 0.000001)
			throw new IllegalArgumentException("localHourAngle: " + (localHourAngle - 11.105902));

		//compute the sun position (rigth ascension and declination) with respect to the observer local position at the Earth surface:
		double equatorialHorizontalParallax = calc.equatorialHorizontalParallax(eclipticCoord.getDistance());
		double u = StrictMath.atan((1. - SunPosition.EARTH_FLATTENING) * StrictMath.tan(location.getLatitude()));
		double x = StrictMath.cos(u)
			+ (calc.getObserverElevation() / SunPosition.EARTH_EQUATORIAL_RADIUS) * StrictMath.cos(location.getLatitude());
		double y = 0.99664719 * StrictMath.sin(u)
			+ (calc.getObserverElevation() / SunPosition.EARTH_EQUATORIAL_RADIUS) * StrictMath.sin(location.getLatitude());
		double deltaRightAscension = StrictMath.atan2(
			-x * StrictMath.sin(equatorialHorizontalParallax) * StrictMath.sin(localHourAngle),
			StrictMath.cos(coord.getDeclination()) - x * StrictMath.sin(equatorialHorizontalParallax) * StrictMath.cos(localHourAngle)
		);
		//calculate the topocentric Sun Right Ascension: α'
		double rightAscensionTopocentric = coord.getRightAscension() + deltaRightAscension;
		if(Math.abs(rightAscensionTopocentric - StrictMath.toRadians(202.22704)) > 0.00001)
			throw new IllegalArgumentException("rightAscensionTopocentric: " + rightAscensionTopocentric);
		double declinationTopocentric = StrictMath.atan2(
			((StrictMath.sin(coord.getDeclination()) - y * StrictMath.sin(equatorialHorizontalParallax))
				* StrictMath.cos(deltaRightAscension)),
			StrictMath.cos(coord.getDeclination()) - y * StrictMath.sin(equatorialHorizontalParallax) * StrictMath.cos(localHourAngle)
		);
		if(Math.abs(declinationTopocentric - StrictMath.toRadians(-9.316179)) > 0.000001)
			throw new IllegalArgumentException("declinationTopocentric: " + declinationTopocentric);
		//calculate the topocentric local hour angle: H’
		double localHourAngleTopocentric = localHourAngle - deltaRightAscension;
		if(Math.abs(localHourAngleTopocentric - StrictMath.toRadians(11.10629)) > 0.00002)
			throw new IllegalArgumentException("localHourAngleTopocentric: " + localHourAngleTopocentric);
		//calculate the topocentric elevation angle without atmospheric refraction correction: e0
		double trueElevation = StrictMath.asin(
			StrictMath.sin(location.getLatitude()) * StrictMath.sin(declinationTopocentric)
				+ StrictMath.cos(location.getLatitude()) * StrictMath.cos(declinationTopocentric)
				* StrictMath.cos(localHourAngleTopocentric)
		);
		//calculate the atmospheric refraction correction: Δe
		double deltaElevation = AtmosphereHelper.atmosphericRefractionCorrection(calc.getPressure(), calc.getTemperature(),
			trueElevation);
		//calculate the topocentric elevation angle: e
		double elevationTopocentric = trueElevation + deltaElevation;
		//calculate the topocentric zenith angle: θ
		double zenithTopocentric = 90. - elevationTopocentric;
		if(Math.abs(zenithTopocentric - 50.11162) > 0.00001)
			throw new IllegalArgumentException("zenithTopocentric: " + (zenithTopocentric - 50.11162));
		//calculate the topocentric astronomers azimuth angle (measured westward from south): Γ
		double azimuthTopocentric = MathHelper.mod2pi(StrictMath.atan2(
			StrictMath.sin(localHourAngleTopocentric),
			StrictMath.cos(localHourAngleTopocentric) * StrictMath.sin(location.getLatitude())
			- StrictMath.tan(declinationTopocentric) * StrictMath.cos(location.getLatitude())
		));
		//calculate the topocentric azimuth angle (measured westward from north): M
		double azimuthTopocentricNavigators = MathHelper.mod2pi(azimuthTopocentric + StrictMath.PI);
		if(Math.abs(azimuthTopocentricNavigators - 194.34024) > 0.00001)
			throw new IllegalArgumentException("azimuthTopocentricNavigators: " + (azimuthTopocentricNavigators - 194.34024));
		//calculate the incidence angle for a surface oriented in any direction, I
		double incidence = StrictMath.acos(
			StrictMath.cos(zenithTopocentric) * StrictMath.cos(calc.getSurfaceSlope())
				+ StrictMath.sin(zenithTopocentric) * StrictMath.sin(calc.getSurfaceSlope())
				* StrictMath.cos(azimuthTopocentric - calc.getSurfaceAzimuthRotation())
		);
		if(Math.abs(incidence - StrictMath.toRadians(25.18700)) > 0.00001)
			throw new IllegalArgumentException("incidence: " + incidence);
	}

	@Test
	void astronomicalSunset() throws SolarEventException{
		GNSSLocation location = GNSSLocation.create(45.65, 12.19);
		SolarEventCalculator calc = SolarEventCalculator.create(location);

		LocalDateTime datetime = calc.sunset(LocalDate.of(2022, 12, 25), Zenith.ASTRONOMICAL);

		assertTimeEquals("17:20:00", datetime);
	}

	@Test
	void civilSunset() throws SolarEventException{
		GNSSLocation location = GNSSLocation.create(45.65, 12.19);
		SolarEventCalculator calc = SolarEventCalculator.create(location);

		LocalDateTime datetime = calc.sunset(LocalDate.of(2022, 12, 25), Zenith.CIVIL);

		assertTimeEquals("16:06:00", datetime);
	}

	@Test
	void officialSunset() throws SolarEventException{
		GNSSLocation location = GNSSLocation.create(45.65, 12.19);
		SolarEventCalculator calc = SolarEventCalculator.create(location);

		LocalDateTime datetime = calc.sunset(LocalDate.of(2022, 12, 25), Zenith.OFFICIAL);

		assertTimeEquals("15:32:00", datetime);
	}


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
