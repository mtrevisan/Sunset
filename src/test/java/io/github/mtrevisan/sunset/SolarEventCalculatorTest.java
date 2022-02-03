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

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.time.LocalDate;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;


@SuppressWarnings("ALL")
class SolarEventCalculatorTest{

	@Test
	void sunPosition(){
		Location location = Location.create(45.65, 12.19);
		SolarEventCalculator calc = SolarEventCalculator.create(location);

		final double jd = SolarEventCalculator.julianDay(1957, 10, 4) + (19. + 29. / 60.) / 24.;
		EquatorialCoordinate coord = calc.sunPosition(jd);

		Assertions.assertEquals("EquatorialCoordinate{α: 12h 41m 33.57s, δ: -4° 28' 17.22\"}", coord.toString());
	}

	@Test
	void astronomicalSunset() throws SolarEventException{
		Location location = Location.create(45.65, 12.19);
		SolarEventCalculator calc = SolarEventCalculator.create(location);

		LocalDateTime datetime = calc.sunset(LocalDate.of(2022, 12, 25), Zenith.ASTRONOMICAL);

		assertTimeEquals("17:20", datetime);
	}

	@Test
	void civilSunset() throws SolarEventException{
		Location location = Location.create(45.65, 12.19);
		SolarEventCalculator calc = SolarEventCalculator.create(location);

		LocalDateTime datetime = calc.sunset(LocalDate.of(2022, 12, 25), Zenith.CIVIL);

		assertTimeEquals("16:06", datetime);
	}

	@Test
	void officialSunset() throws SolarEventException{
		Location location = Location.create(45.65, 12.19);
		SolarEventCalculator calc = SolarEventCalculator.create(location);

		LocalDateTime datetime = calc.sunset(LocalDate.of(2022, 12, 25), Zenith.OFFICIAL);

		assertTimeEquals("15:32", datetime);
	}


	private static void assertTimeEquals(final String expectedTime, final LocalDateTime actualTime){
		final int expectedMinutes = getMinutes(expectedTime);
		final int actualMinutes = getMinutes(actualTime.toLocalTime().toString());
		Assertions.assertEquals(expectedMinutes, actualMinutes, 1,
			"Expected " + expectedTime + ", got " + actualTime.format(DateTimeFormatter.ofPattern("HH:MM")));
	}

	private static int getMinutes(final String time){
		final String[] timeParts = time.split("\\:");
		if(timeParts[0].equals("00"))
			timeParts[0] = "24";
		return 60 * Integer.valueOf(timeParts[0]) + Integer.valueOf(timeParts[1]);
	}

}
