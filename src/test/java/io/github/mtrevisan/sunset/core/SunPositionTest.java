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

import io.github.mtrevisan.sunset.coordinates.AtmosphericModel;
import io.github.mtrevisan.sunset.JulianDate;
import io.github.mtrevisan.sunset.coordinates.EclipticCoordinate;
import io.github.mtrevisan.sunset.coordinates.EquatorialCoordinate;
import io.github.mtrevisan.sunset.coordinates.GeographicLocation;
import io.github.mtrevisan.sunset.coordinates.HorizontalCoordinate;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.time.LocalTime;


@SuppressWarnings("ALL")
class SunPositionTest{

	@Test
	void test(){
		double jd = JulianDate.of(2005, 1, 1)
			+ JulianDate.timeOf(LocalTime.of(19, 29));
		EclipticCoordinate eclipticCoord = SunPosition.sunEclipticPosition(jd);

		Assertions.assertEquals("EclipticCoordinate{b: 0° 0' 0.55\", l: 191° 18' 9.93\", r: 1}", eclipticCoord.toString());
	}

	@Test
	void sunEclipticPositionSputnik(){
		double jd = JulianDate.of(1957, 10, 4)
			+ JulianDate.timeOf(LocalTime.of(19, 29));
		EclipticCoordinate eclipticCoord = SunPosition.sunEclipticPosition(jd);

		Assertions.assertEquals("EclipticCoordinate{b: 0° 0' 0.55\", l: 191° 18' 9.93\", r: 1}", eclipticCoord.toString());
	}

	@Test
	void sunEquatorialPositionSputnik(){
		double jd = JulianDate.of(1957, 10, 4)
			+ JulianDate.timeOf(LocalTime.of(19, 29));
		EclipticCoordinate eclipticCoord = SunPosition.sunEclipticPosition(jd);
		EquatorialCoordinate equatorialCoord = SunPosition.sunEquatorialPosition(eclipticCoord, jd);

		//https://www.eso.org/observing/etc/bin/gen/form?INS.MODE=swspectr+INS.NAME=SKYCALC
		//α☉ = 12h 41m 32s  δ☉ = -04° 28' 14"
		//http://cosinekitty.com/solar_system.html
		Assertions.assertEquals("EquatorialCoordinate{α: 12h 41m 33.28s, δ: -4° 28' 15.41\"}", equatorialCoord.toString());
	}

	@Test
	void sunTopocentricPosition(){
		double jd = JulianDate.of(1957, 10, 4)
			+ JulianDate.timeOf(LocalTime.of(19, 29));
		GeographicLocation location = GeographicLocation.create(45.714920, 12.194179, 23.);
		AtmosphericModel atmosphericModel = AtmosphericModel.create(1017., 18.);
		EclipticCoordinate eclipticCoord = SunPosition.sunEclipticPosition(jd);
		EquatorialCoordinate equatorialCoord = SunPosition.sunEquatorialPosition(eclipticCoord, jd);
		HorizontalCoordinate horizontalCoord = SunPosition.sunTopocentricPosition(location, atmosphericModel, eclipticCoord, jd);

		Assertions.assertEquals("HorizontalCoordinate{azi: 295° 18' 11.89\", alt: -28° 28' 53.12\", r: 1}", horizontalCoord.toString());
	}

}
