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

import io.github.mtrevisan.sunset.JulianDay;
import io.github.mtrevisan.sunset.coordinates.EclipticCoordinate;
import io.github.mtrevisan.sunset.coordinates.EquatorialCoordinate;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.time.LocalTime;


@SuppressWarnings("ALL")
class SunPositionTest{

	@Test
	void sunEclipticPositionSputnik(){
		double jd = JulianDay.of(1957, 10, 4)
			+ JulianDay.timeOf(LocalTime.of(19, 29));
		EclipticCoordinate coord = SunPosition.sunEclipticPosition(jd);

		Assertions.assertEquals("EclipticCoordinate{b: 0° 0' 0.55\", l: 191° 18' 9.93\", r: 1}", coord.toString());
	}

	@Test
	void sunEquatorialPositionSputnik(){
		double jd = JulianDay.of(1957, 10, 4)
			+ JulianDay.timeOf(LocalTime.of(19, 29));
		EclipticCoordinate eclipticCoord = SunPosition.sunEclipticPosition(jd);
		EquatorialCoordinate coord = SunPosition.sunEquatorialPosition(eclipticCoord, jd);

		//https://www.eso.org/observing/etc/bin/gen/form?INS.MODE=swspectr+INS.NAME=SKYCALC
		//α☉ = 12h 41m 32s  δ☉ = -04° 28' 14"
		//http://cosinekitty.com/solar_system.html
		Assertions.assertEquals("EquatorialCoordinate{α: 12h 41m 33.28s, δ: -4° 28' 15.41\"}", coord.toString());
	}

}
