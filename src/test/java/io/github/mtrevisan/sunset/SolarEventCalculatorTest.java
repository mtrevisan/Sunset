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
import java.time.format.DateTimeFormatter;


@SuppressWarnings("ALL")
class SolarEventCalculatorTest{

	@Test
	void sunPosition(){
		GNSSLocation location = GNSSLocation.create(45.65, 12.19);
		SolarEventCalculator calc = SolarEventCalculator.create(location);

		final double jd = JulianDay.of(1957, 10, 4) + (19. + 29. / 60.) / 24.;
		EquatorialCoordinate coord = SolarEventCalculator.sunPosition(jd);

		Assertions.assertEquals("EquatorialCoordinate{α: 12h 41m 33.57s, δ: -4° 28' 17.22\"}", coord.toString());
	}

	@Test
	void asd() throws SolarEventException{
		GNSSLocation location = GNSSLocation.create(
			SolarEventCalculator.toDegrees(-3, 4, 0.),
			SolarEventCalculator.toDegrees(37, 21, 33.));
		SolarEventCalculator calc = SolarEventCalculator.create(location);

		LocalDateTime datetime = calc.sunset(LocalDate.of(2017, 12, 22), Zenith.ASTRONOMICAL);

		assertTimeEquals("18:47:47", datetime);
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
			"Expected " + expectedTime + ", got " + actualTime.format(DateTimeFormatter.ofPattern("HH:MM")));
	}

	private static int getSeconds(final String time){
		final String[] timeParts = time.split("\\:");
		return 60 * (60 * Integer.valueOf(timeParts[0]) + Integer.valueOf(timeParts[1])) + Integer.valueOf(timeParts[2]);
	}

}
