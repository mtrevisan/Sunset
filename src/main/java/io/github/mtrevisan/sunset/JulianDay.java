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

import java.time.LocalDate;
import java.time.LocalDateTime;
import java.time.LocalTime;


/**
 * <The Julian day is the Julian day number for the preceding noon plus the fraction of the day (counting 86400 seconds) since that instant.
 *
 * <p>
 * A Julian day number is defined as continuous number associated with the solar day and is zero at Greenwich mean noon on 1st of January
 * 4713 BC.
 * </p>
 */
public final class JulianDay{

	/** Last day of the Julian calendar, that is October 4, 1582, after which 10 day were added [day]. */
	private static final double JULIAN_CALENDAR_END = 2299160.;
	private static final int GREGOR_XIII_REFORM_YEAR = 1582;

	public static final double MJD = 2400000.5;
	/** 1.5 Jan 2000 UT - Julian epoch. */
	private static final double J2000 = 2451545.;

	public static final double CIVIL_SAECULUM = 36525.;
	public static final double CIVIL_MILLENNIUM = CIVIL_SAECULUM * 10.;
	//[°/h]
	public static final double DEGREES_PER_HOUR = 15.;
	public static final double HOURS_IN_DAY = 24.;
	public static final double SECONDS_IN_HOUR = 3600.;
	public static final double SECONDS_IN_DAY = SECONDS_IN_HOUR * HOURS_IN_DAY;


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
		return timeOf(date.toLocalTime());
	}

	public static double timeOf(final LocalTime time){
		return time.toSecondOfDay() / SECONDS_IN_DAY;
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
	public static double of(final int year, final int month, final int day){
		int yy = year;
		int mm = month;
		if(month < 3){
			yy --;
			mm += 12;
		}

		double jd =  (int)(365.25 * (yy + 4716)) + (int)(30.6001 * (mm + 1)) + day - 1524.5;

		if(jd > JULIAN_CALENDAR_END){
			//number of full centuries
			final int aa = (int)(yy / 100.);
			//days within the whole centuries (in the Julian Calendar) adding back days removed in the Gregorian Calendar
			jd += 2 - aa + (int)(aa / 4);
		}

		return jd;
	}

	/**
	 * Calculated the Julian Century since JDE2451545, that is J2000.0.
	 *
	 * @param jd	The Julian Day [day].
	 * @return	The Julian Century.
	 */
	public static double centuryJ2000Of(final double jd){
		return (jd - J2000) / CIVIL_SAECULUM;
	}

	/**
	 * Calculated the Julian Millennium since JDE2451545, that is J2000.0.
	 *
	 * @param jd	The Julian Day [day].
	 * @return	The Julian Century.
	 */
	public static double millenniumJ2000Of(final double jd){
		return (jd - J2000) / CIVIL_MILLENNIUM;
	}

	/**
	 * Calculated the Modified Julian Day starting at midnight since JDE2451545, that is J2000.0.
	 *
	 * @param jd	The Julian Day [day].
	 * @return	The Modified Julian Day starting at midnight.
	 */
	public static double mjdOf(final double jd){
		return jd - MJD;
	}


	public static boolean isLeapYear(final int year){
		return ((year & 0x03) == 0 && (year < GREGOR_XIII_REFORM_YEAR || (year % 100) != 0 || (year % 400) == 0));
	}

	public static double[] extractDateAndTime(final double jd){
		//calculate time of day [day]
		final double time = ((jd - 0.5) % 1);
		//calculate UT at 0h [day]
		final double date = (long)jd - 0.5;
		return new double[]{date, time};
	}

}
