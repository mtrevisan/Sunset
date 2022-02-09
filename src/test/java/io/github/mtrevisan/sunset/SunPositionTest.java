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
import io.github.mtrevisan.sunset.test.AstroDay;
import io.github.mtrevisan.sunset.test.AstroLib;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.time.LocalDateTime;
import java.time.LocalTime;
import java.time.format.DateTimeFormatter;


@SuppressWarnings("ALL")
class SunPositionTest{

	@Test
	void sunPosition1(){
		final double jd = JulianDay.of(1957, 10, 4) + JulianDay.timeOf(LocalTime.of(19, 29));
		final double tt = JulianDay.centuryJ2000Of(jd);
		EquatorialCoordinate coord = SunPosition.sunPosition(tt);

		AstroDay astroDay = new AstroDay();
		AstroLib.computeAstroDay(jd, astroDay);

		Assertions.assertEquals("EquatorialCoordinate{α: 12h 42m 18.26s, δ: -4° 31' 37.18\"}", coord.toString());
	}

	@Test
	void sunPosition2(){
		final double jd = JulianDay.of(2020, 1, 1);
		final double tt = JulianDay.centuryJ2000Of(jd);
		EquatorialCoordinate coord = SunPosition.sunPosition(tt);

		Assertions.assertEquals("EquatorialCoordinate{α: 18h 42m 18.7s, δ: -23° 3' 3.69\"}", coord.toString());
	}


	@Test
	void sunPosition_stepByStep() throws SolarEventException{
		final double ut = JulianDay.of(2003, 10, 17)
			+ JulianDay.timeOf(LocalTime.of(19, 30, 30));
		final double jd = TimeHelper.universalTimeToTerrestrialTime(ut, 67.);
		final double tt = JulianDay.centuryJ2000Of(jd);

		final double geometricMeanLongitude = SunPosition.geometricMeanLongitude(tt);
		if(Math.abs(geometricMeanLongitude - 204.0182616917) > 0.0000000001)
			throw new IllegalArgumentException("geometricMeanLongitude: " + (geometricMeanLongitude - 204.0182616917));
		final double[] nutation = SunPosition.nutationCorrection(tt);
		if(Math.abs(nutation[0] - -0.25189948) > 0.00000001)
			throw new IllegalArgumentException("nutationInLongitude: " + (nutation[0] - -0.25189948));
		if(Math.abs(nutation[1] - 0.10499379) > 0.00000001)
			throw new IllegalArgumentException("nutationInObliquity: " + (nutation[1] - 0.10499379));
		final double radiusVector = SunPosition.radiusVector(tt);
		if(Math.abs(radiusVector - 0.9965422974) > 0.0000000001)
			throw new IllegalArgumentException("radiusVector: " + (radiusVector - 0.9965422974));
		final double aberration = SunPosition.aberrationCorrection(radiusVector);
		final double apparentGeometricLongitude = SunPosition.apparentGeometricLongitude(geometricMeanLongitude, nutation[0], aberration);
		if(Math.abs(apparentGeometricLongitude - 203.7606508546) > 0.0000000002)
			throw new IllegalArgumentException("apparentGeometricLongitude: " + (apparentGeometricLongitude - 203.7606508546));
		final double meanEclipticObliquity = SunPosition.meanEclipticObliquity(tt);
		final double trueEclipticObliquity = SunPosition.trueEclipticObliquity(meanEclipticObliquity, nutation[1]);
		if(Math.abs(trueEclipticObliquity - 23.543792) > 0.000001)
			throw new IllegalArgumentException("trueEclipticObliquity: " + (trueEclipticObliquity - 23.543792));
		final double geometricMeanLatitude = SunPosition.geometricMeanLatitude(tt);
		if(Math.abs(geometricMeanLatitude - 0.0001011219) > 0.0000000001)
			throw new IllegalArgumentException("geometricMeanLatitude: " + (geometricMeanLatitude - 0.0001011219));
		EquatorialCoordinate coord = EquatorialCoordinate.createFromEcliptical(geometricMeanLatitude, apparentGeometricLongitude,
			trueEclipticObliquity);
		if(Math.abs(coord.getRightAscension() - 201.97832) > 0.00001)
			throw new IllegalArgumentException("rightAscension: " + (coord.getRightAscension() - 201.97832));
		if(Math.abs(coord.getDeclination() - -9.26166) > 0.00001)
			throw new IllegalArgumentException("declination: " + (coord.getDeclination() - -9.26166));

		Assertions.assertEquals("EquatorialCoordinate{α: 13h 27m 54.8s, δ: -9° 15' 41.98\"}", coord.toString());
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
