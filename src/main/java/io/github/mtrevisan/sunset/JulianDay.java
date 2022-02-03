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

import java.time.LocalDate;
import java.time.LocalDateTime;


/**
 * <The Julian day is the Julian day number for the preceding noon plus the fraction of the day (counting 86400 seconds) since that instant.
 *
 * <p>
 * A Julian day number is defined as continuous number associated with the solar day and is zero at Greenwich mean noon on 1st of January
 * 4713 BC.
 * </p>
 */
public final class JulianDay{

	/**
	 * Calculate Julian Day at 0 UTC.
	 *
	 * @param date	The date.
	 * @return	The Julian Day [day].
	 */
	public static double of(final LocalDateTime date){
		final double jdNoon = of(date.toLocalDate());
		final double time = timeOf(date);
		return jdNoon + time;
	}

	public static double timeOf(final LocalDateTime date){
		return date.toLocalTime().toSecondOfDay() / (24. * 60. * 60.);
	}

	/**
	 * Calculate Julian Day at 0 UTC.
	 *
	 * @param date	The date.
	 * @return	The Julian Day [day].
	 */
	public static double of(final LocalDate date){
		return of(date.getYear(), date.getMonthValue(), date.getDayOfMonth());
	}

	/**
	 * Calculate Julian Day at 0 UTC.
	 *
	 * @param year	The year.
	 * @param month	The month (1 is January).
	 * @param day	The day.
	 * @return	The Julian Day [day].
	 */
	public static double of(int year, int month, final int day){
		if(month <= 2){
			year --;
			month += 12;
		}

		final double a = StrictMath.floor(year / 100.);
		final double b = 2. - a + StrictMath.floor(a / 4.);
		return StrictMath.floor(365.25 * (year + 4716)) + StrictMath.floor(30.6001 * (month + 1)) + day + b - 1524.5;
	}

	/**
	 * Calculated the Julian Century since JDE2451545, that is J2000.0.
	 *
	 * @param jd	The Julian Day [day].
	 * @return	The Julian Century.
	 */
	public static double centuryJ2000Of(final double jd){
		return (jd - 2451545.) / 36525.;
	}

	/**
	 * Calculated the Modified Julian Day starting at midnight since JDE2451545, that is J2000.0.
	 *
	 * @param jd	The Julian Day [day].
	 * @return	The Modified Julian Day starting at midnight.
	 */
	public static double mjdOf(final double jd){
		return jd - 2400000.5;
	}

}
