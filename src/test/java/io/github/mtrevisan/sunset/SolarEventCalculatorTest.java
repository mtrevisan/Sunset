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
		final GNSSLocation location = GNSSLocation.create(39.742476, -105.1786);
		SolarEventCalculator calc = SolarEventCalculator.create(location)
			.withObserverElevation(1830.14)
			.withPressure(820.)
			.withTemperature(11.)
			.withSurfaceSlopeAndAzimuthRotation(30., -10.);

		final double ut = JulianDay.of(2003, 10, 17)
			+ JulianDay.timeOf(LocalTime.of(19, 30, 30));
		final double jd = TimeHelper.universalTimeToTerrestrialTime(ut, 67.);
		final double tt = JulianDay.centuryJ2000Of(jd);

		EquatorialCoordinate coord = SunPosition.sunPosition(tt);
		final double[] nutation = SunPosition.nutationCorrection(tt);
		final double radiusVector = SunPosition.radiusVector(tt);
		final double meanEclipticObliquity = SunPosition.meanEclipticObliquity(tt);
		final double trueEclipticObliquity = SunPosition.trueEclipticObliquity(meanEclipticObliquity, nutation[1]);

		final double meanSiderealTime = TimeHelper.meanSiderealTime(ut);
		final double apparentSiderealTime = TimeHelper.apparentSiderealTime(meanSiderealTime, trueEclipticObliquity, nutation[0]);
		final double localMeanSiderealTime = TimeHelper.localMeanSiderealTime(apparentSiderealTime, location);
		final double localHourAngle = TimeHelper.localHourAngle(localMeanSiderealTime, coord.getRightAscension());
		if(Math.abs(localHourAngle - 11.105902) > 0.000001)
			throw new IllegalArgumentException("localHourAngle: " + (localHourAngle - 11.105902));

		//compute the sun position (rigth ascension and declination) with respect to the observer local position at the Earth surface:
		final double equatorialHorizontalParallax = calc.equatorialHorizontalParallax(radiusVector);
		final double latitude = StrictMath.toRadians(location.getLatitude());
		final double u = StrictMath.atan((1. - SunPosition.EARTH_FLATTENING) * StrictMath.tan(latitude));
		final double x = StrictMath.cos(u)
			+ (calc.getObserverElevation() / SunPosition.EARTH_EQUATORIAL_RADIUS) * StrictMath.cos(latitude);
		final double y = 0.99664719 * StrictMath.sin(u)
			+ (calc.getObserverElevation() / SunPosition.EARTH_EQUATORIAL_RADIUS) * StrictMath.sin(latitude);
		final double horizontalParallax = StrictMath.toRadians(equatorialHorizontalParallax);
		final double lha = StrictMath.toRadians(localHourAngle);
		final double decl = StrictMath.toRadians(coord.getDeclination());
		final double deltaRightAscension = StrictMath.toDegrees(StrictMath.atan2(
			-x * StrictMath.sin(horizontalParallax) * StrictMath.sin(lha),
			StrictMath.cos(decl) - x * StrictMath.sin(horizontalParallax) * StrictMath.cos(lha)
		));
		//calculate the topocentric Sun Right Ascension: α'
		final double rightAscensionTopocentric = coord.getRightAscension() + deltaRightAscension;
		if(Math.abs(rightAscensionTopocentric - 202.22704) > 0.00001)
			throw new IllegalArgumentException("rightAscensionTopocentric: " + (rightAscensionTopocentric - 202.22704));
		final double declinationTopocentric = StrictMath.toDegrees(StrictMath.atan2(
			((StrictMath.sin(decl) - y * StrictMath.sin(horizontalParallax)) * StrictMath.cos(StrictMath.toRadians(deltaRightAscension))),
			StrictMath.cos(decl) - y * StrictMath.sin(horizontalParallax) * StrictMath.cos(lha)
		));
		if(Math.abs(declinationTopocentric - -9.316179) > 0.000001)
			throw new IllegalArgumentException("declinationTopocentric: " + (declinationTopocentric - -9.316179));
		//calculate the topocentric local hour angle: H’
		final double localHourAngleTopocentric = localHourAngle - deltaRightAscension;
		if(Math.abs(localHourAngleTopocentric - 11.10629) > 0.00002)
			throw new IllegalArgumentException("localHourAngleTopocentric: " + (localHourAngleTopocentric - 11.10629));
		//calculate the topocentric elevation angle without atmospheric refraction correction: e0
		final double lhaTopocentric = StrictMath.toRadians(localHourAngleTopocentric);
		final double trueElevation = StrictMath.toDegrees(StrictMath.asin(
			StrictMath.sin(latitude) * StrictMath.sin(StrictMath.toRadians(declinationTopocentric))
				+ StrictMath.cos(latitude) * StrictMath.cos(StrictMath.toRadians(declinationTopocentric)) * StrictMath.cos(lhaTopocentric)
		));
		//calculate the atmospheric refraction correction: Δe
		final double deltaElevation = AtmosphereHelper.atmosphericRefractionCorrection(calc.getPressure(), calc.getTemperature(),
			trueElevation);
		//calculate the topocentric elevation angle: e
		final double elevationTopocentric = trueElevation + deltaElevation;
		//calculate the topocentric zenith angle: θ
		final double zenithTopocentric = 90. - elevationTopocentric;
		if(Math.abs(zenithTopocentric - 50.11162) > 0.00001)
			throw new IllegalArgumentException("zenithTopocentric: " + (zenithTopocentric - 50.11162));
		//calculate the topocentric astronomers azimuth angle (measured westward from south): Γ
		final double azimuthTopocentric = MathHelper.limitRangeDegree(StrictMath.toDegrees(StrictMath.atan2(
			StrictMath.sin(lhaTopocentric),
			StrictMath.cos(lhaTopocentric) * StrictMath.sin(latitude)
			- StrictMath.tan(StrictMath.toRadians(declinationTopocentric)) * StrictMath.cos(latitude)
		)));
		//calculate the topocentric azimuth angle (measured westward from north): M
		final double azimuthTopocentricNavigators = MathHelper.limitRangeDegree(azimuthTopocentric + 180.);
		if(Math.abs(azimuthTopocentricNavigators - 194.34024) > 0.00001)
			throw new IllegalArgumentException("azimuthTopocentricNavigators: " + (azimuthTopocentricNavigators - 194.34024));
		//calculate the incidence angle for a surface oriented in any direction, I
		final double zenith = StrictMath.toRadians(zenithTopocentric);
		final double slope = StrictMath.toRadians(calc.getSurfaceSlope());
		final double incidence = StrictMath.toDegrees(StrictMath.acos(
			StrictMath.cos(zenith) * StrictMath.cos(slope)
				+ StrictMath.sin(zenith) * StrictMath.sin(slope)
				* StrictMath.cos(StrictMath.toRadians(azimuthTopocentric - calc.getSurfaceAzimuthRotation()))
		));
		if(Math.abs(incidence - 25.18700) > 0.00001)
			throw new IllegalArgumentException("incidence: " + (incidence - 25.18700));
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
		final int expectedSeconds = getSeconds(expectedTime);
		final int actualSeconds = getSeconds(actualTime.toLocalTime().toString());
		Assertions.assertEquals(expectedSeconds, actualSeconds, 1,
			"Expected " + expectedTime + ", got " + actualTime.format(DateTimeFormatter.ofPattern("HH:MM:ss")));
	}

	private static int getSeconds(final String time){
		final String[] timeParts = time.split("\\:");
		return 60 * (60 * Integer.valueOf(timeParts[0]) + Integer.valueOf(timeParts[1])) + Integer.valueOf(timeParts[2]);
	}

}
