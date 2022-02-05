/**
 * Copyright (c) 2021 Mauro Trevisan
 * <p>
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 * <p>
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * <p>
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


public final class TimeHelper{

	private TimeHelper(){}


	/**
	 * Convert <a href="https://en.wikipedia.org/wiki/Universal_Time">Universal Time</a> to
	 * <a href="https://en.wikipedia.org/wiki/Terrestrial_Time">Terrestrial (Dynamical) Time</a>.
	 * <p>
	 * The Universal Time (UT), or Greenwich civil time, is based on the Earth’s rotation and counted from 0-hour at midnight; the unit is
	 * mean solar day.<br/>
	 * The Terrestrial (Dynamical) Time (TDT or TT) is the time scale of ephemerides for observations from the Earth surface. For practical
	 * purposes, TT must be realized by actual clocks in the Earth system.
	 * </p>
	 *
	 * @param ut	Julian Day of Universal Time. [day]
	 * @param deltaT   The difference between the Earth rotation time and the TT; it is derived from observation only and reported yearly
	 * 	in the Astronomical Almanac [s].
	 * @return	The Terrestrial Time [day].
	 */
	public static double universalTimeToTerrestrialTime(final double ut, final double deltaT){
		return ut + deltaT / JulianDay.SECONDS_IN_DAY;
	}

}
