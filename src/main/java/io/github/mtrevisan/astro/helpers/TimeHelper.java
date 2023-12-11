/*
 * Copyright (c) 2023 Mauro Trevisan
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
package io.github.mtrevisan.astro.helpers;

import io.github.mtrevisan.astro.coordinates.GeographicLocation;
import io.github.mtrevisan.astro.core.NutationCorrections;

import java.io.IOException;
import java.time.LocalDate;
import java.time.LocalDateTime;
import java.time.ZonedDateTime;
import java.util.List;


/**
 * @see <a href="https://www.stjarnhimlen.se/comp/time.html">Time scales</>
 * @see <a href="https://webspace.science.uu.nl/~gent0113/deltat/deltat.htm">Delta T: Terrestrial Time, Universal Time and algorithms for historical periods</>
 * @see <a href="http://www.internetsv.info/UniClock.html">Universal Time Calculator</a>
 * @see <a href="https://lweb.cfa.harvard.edu/~jzhao/times.html">Astronomical Times</a>
 */
public final class TimeHelper{

	/**
	 * NOTE: EVERY MONTH UPDATE THIS FILE!!!
	 *
	 * @see <a href="https://maia.usno.navy.mil/ser7/deltat.data">Delta T (ΔT) values</a>
	 */
	private static List<double[]> DELTA_T_DATA;
	static{
		try{
			DELTA_T_DATA = ResourceReader.readPlain("deltat.dat");
		}
		catch(final IOException ignored){}
	}


	/** Lunar acceleration parameter (as of 2002 of −25.858±0.003″/cy2) ["/cy^2]. */
	private static final double LUNAR_ACCELERATION = -25.858;

	private static final double[] GREENWICH_MEAN_SIDEREAL_TIME_COEFFS = {280.460_618_375, 360.985_647_366_29 * JulianDate.CIVIL_SAECULUM,
		0.000_387_933, -1. / 38_710_000.};


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
	 * @param ut	Julian Date of Universal Time. [day]
	 * @param deltaT   The difference between the Earth rotation time and the TT; it is derived from observation only and reported yearly
	 * 	in the Astronomical Almanac [s].
	 * @return	The Terrestrial Time [day].
	 * @see #deltaT(double)
	 */
	public static double universalTimeToTerrestrialTime(final double ut, final double deltaT){
		return ut + deltaT / JulianDate.SECONDS_PER_DAY;
	}

	/**
	 * Convert <a href="https://en.wikipedia.org/wiki/Terrestrial_Time">Terrestrial (Dynamical) Time</a> to
	 * <a href="https://en.wikipedia.org/wiki/Universal_Time">Universal Time</a>.
	 * <p>
	 * The Universal Time (UT), or Greenwich civil time, is based on the Earth’s rotation and counted from 0-hour at midnight; the unit is
	 * mean solar day.<br/>
	 * The Terrestrial (Dynamical) Time (TDT or TT) is the time scale of ephemerides for observations from the Earth surface. For practical
	 * purposes, TT must be realized by actual clocks in the Earth system.
	 * </p>
	 *
	 * @param tt	Julian Date of Terrestrial Time. [day]
	 * @param deltaT   The difference between the Earth rotation time and the TT; it is derived from observation only and reported yearly
	 * 	in the Astronomical Almanac [s].
	 * @return	The Universal Time [day].
	 * @see #deltaT(double)
	 */
	public static double terrestrialTimeToUniversalTime(final double tt, final double deltaT){
		return tt - deltaT / JulianDate.SECONDS_PER_DAY;
	}

	/**
	 * Convert <a href="https://en.wikipedia.org/wiki/Terrestrial_Time">Terrestrial (Dynamical) Time</a> to
	 * <a href="https://en.wikipedia.org/wiki/Universal_Time">Universal Time</a>.
	 * <p>
	 * The Universal Time (UT), or Greenwich civil time, is based on the Earth’s rotation and counted from 0-hour at midnight; the unit is
	 * mean solar day.<br/>
	 * The Terrestrial (Dynamical) Time (TDT or TT) is the time scale of ephemerides for observations from the Earth surface. For practical
	 * purposes, TT must be realized by actual clocks in the Earth system.
	 * </p>
	 *
	 * @param zonedDateTime	Julian Date of Terrestrial Time. [day]
	 * @return	The Universal Time.
	 * @see #deltaT(double)
	 */
	public static ZonedDateTime terrestrialTimeToUniversalTime(final ZonedDateTime zonedDateTime){
		return terrestrialTimeToUniversalTime(zonedDateTime, deltaT(zonedDateTime));
	}

	/**
	 * Convert <a href="https://en.wikipedia.org/wiki/Terrestrial_Time">Terrestrial (Dynamical) Time</a> to
	 * <a href="https://en.wikipedia.org/wiki/Universal_Time">Universal Time</a>.
	 * <p>
	 * The Universal Time (UT), or Greenwich civil time, is based on the Earth’s rotation and counted from 0-hour at midnight; the unit is
	 * mean solar day.<br/>
	 * The Terrestrial (Dynamical) Time (TDT or TT) is the time scale of ephemerides for observations from the Earth surface. For practical
	 * purposes, TT must be realized by actual clocks in the Earth system.
	 * </p>
	 *
	 * @param zonedDateTime	Julian Date of Terrestrial Time. [day]
	 * @param deltaT   The difference between the Earth rotation time and the TT; it is derived from observation only and reported yearly
	 * 	in the Astronomical Almanac [s].
	 * @return	The Universal Time.
	 * @see #deltaT(double)
	 */
	public static ZonedDateTime terrestrialTimeToUniversalTime(final ZonedDateTime zonedDateTime, final double deltaT){
		final double tt = JulianDate.of(zonedDateTime);
		final double utc = terrestrialTimeToUniversalTime(tt, deltaT);
		return JulianDate.dateTimeOf(utc)
			.plusSeconds(zonedDateTime.getOffset().getTotalSeconds())
			.atZone(zonedDateTime.getZone());
	}

	/**
	 * Convert <a href="https://en.wikipedia.org/wiki/Universal_Time">Universal Time</a> to
	 * <a href="https://en.wikipedia.org/wiki/Solar_time#Apparent_solar_time">Apparent Solar Time</a>.
	 *
	 * @param zonedDateTime	Julian Date of Universal Time. [day]
	 * @return	The Solar Time.
	 */
	public static LocalDateTime universalTimeToApparentSolarTime(final ZonedDateTime zonedDateTime, final GeographicLocation location){
		final double longitude = MathHelper.mod(location.getLongitude(), 360.) / 360. * JulianDate.SECONDS_PER_DAY;
		final int offset = zonedDateTime.getOffset().getTotalSeconds();
		return zonedDateTime.toLocalDateTime()
			.plusSeconds(Math.round(longitude - offset));
	}


	/**
	 * Estimate <code>∆T</code>, that is the difference between the Earth rotation time and the TT, for the given decimal year.
	 * <p>
	 *    This is based on Espenak and Meeus, "Five Millennium Canon of Solar Eclipses: -1999 to +3000" (NASA/TP-2006-214141) and updated by
	 *    Espenak in 2014 at <a href="https://www.eclipsewise.com/help/deltatpoly2014.html">Eclipsewise</a>.
	 * </p>
	 *
	 * @param zonedDateTime	Date and time.
	 * @return	Estimated <code>ΔT = TT – UT1</code> value [s].
	 *
	 * @see <a href="https://eclipse.gsfc.nasa.gov/SEhelp/deltaT.html">Delta T (ΔT) And Universal Time</a>
	 * @see <a href="https://eclipse.gsfc.nasa.gov/SEhelp/deltat2004.html">Historical Values Of Delta T (ΔT)</a>
	 */
	public static double deltaT(final ZonedDateTime zonedDateTime){
		final double year = decimalYear(zonedDateTime);
		return deltaT(year);
	}

	private static double decimalYear(final ZonedDateTime zonedDateTime){
		return zonedDateTime.getYear() + ((double)zonedDateTime.getDayOfYear() / (zonedDateTime.toLocalDate().isLeapYear()? 366: 365));
	}

	/**
	 * Estimate <code>∆T</code>, that is the difference between the Earth rotation time and the TT, for the given decimal year.
	 * <p>
	 *    This is based on Espenak and Meeus, "Five Millennium Canon of Solar Eclipses: -1999 to +3000" (NASA/TP-2006-214141) and updated by
	 *    Espenak in 2014 at <a href="https://www.eclipsewise.com/help/deltatpoly2014.html">Eclipsewise</a>.
	 * </p>
	 *
	 * @param localDate	Date and time.
	 * @return	Estimated <code>ΔT = TT – UT1</code> value [s].
	 *
	 * @see <a href="https://eclipse.gsfc.nasa.gov/SEhelp/deltaT.html">Delta T (ΔT) And Universal Time</a>
	 * @see <a href="https://eclipse.gsfc.nasa.gov/SEhelp/deltat2004.html">Historical Values Of Delta T (ΔT)</a>
	 */
	public static double deltaT(final LocalDate localDate){
		final double year = decimalYear(localDate);
		return deltaT(year);
	}

	private static double decimalYear(final LocalDate localDate){
		return localDate.getYear() + ((double)localDate.getDayOfYear() / (localDate.isLeapYear()? 366: 365));
	}

	/**
	 * Estimate <code>∆T</code>, that is the difference between the Earth rotation time and the TT, for the given decimal year.
	 * <p>
	 *    This is based on Espenak and Meeus, "Five Millennium Canon of Solar Eclipses: -1999 to +3000" (NASA/TP-2006-214141) and updated by
	 *    Espenak in 2014 at <a href="https://www.eclipsewise.com/help/deltatpoly2014.html">Eclipsewise</a>.
	 * </p>
	 *
	 * @param year	Decimal year.
	 * @return	Estimated <code>ΔT = TT – UT1</code> value [s].
	 *
	 * @see <a href="https://eclipse.gsfc.nasa.gov/SEhelp/deltaT.html">Delta T (ΔT) And Universal Time</a>
	 * @see <a href="https://eclipse.gsfc.nasa.gov/SEhelp/deltat2004.html">Historical Values Of Delta T (ΔT)</a>
	 */
	public static double deltaT(final double year){
		double deltaT = 0.;
		final double u;
		if(year <= -500. || year > 2150.){
			u = (year - 1820.) / 100.;
			deltaT = -20. + 32. * u * u;
		}
		else if(year <= 500.){
			u = year / 100.;
			deltaT = MathHelper.polynomial(u, new double[]{10583.6, -1014.41, 33.78311, -5.952053, -0.1798452, 0.022174192, 0.0090316521});
		}
		else if(year <= 1600.){
			u = (year - 1000.) / 100.;
			deltaT = MathHelper.polynomial(u, new double[]{1574.2, -556.01, 71.23472, 0.319781, -0.8503463, -0.005050998, 0.0083572073});
		}
		else if(year <= 1700.){
			u = (year - 1600.) / 100.;
			deltaT = MathHelper.polynomial(u, new double[]{120, -98.08, -153.2, 1. / 0.007129});
		}
		else if(year <= 1800.){
			u = (year - 1700.) / 100.;
			deltaT = MathHelper.polynomial(u, new double[]{8.83, 16.03, -59.285, 133.36, -1. / 0.01174});
		}
		else if(year <= 1860.){
			u = (year - 1800.) / 100.;
			deltaT = MathHelper.polynomial(u, new double[]{13.72, -33.2447, 68.612, 4111.6, -37436., 121272., -169900., 87500.});
		}
		else if(year <= 1900.){
			u = (year - 1860.) / 100.;
			deltaT = MathHelper.polynomial(u, new double[]{7.62, 57.37, -2517.54, 16806.68, - 44736.24, 1. / 0.0000233174});
		}
		else if(year <= 1920.){
			u = (year - 1900.) / 100.;
			deltaT = MathHelper.polynomial(u, new double[]{-2.79, 149.4119, -598.939, 6196.6, -19700.});
		}
		else if(year <= 1941.){
			u = (year - 1920.) / 100.;
			deltaT = MathHelper.polynomial(u, new double[]{21.20, 84.493, -761.00, 2093.6});
		}
		else if(year <= 1961.){
			u = (year - 1950.) / 100.;
			deltaT = MathHelper.polynomial(u, new double[]{29.07, 40.7, -1. / 0.0233, 1. / 0.002547});
		}
		else if(year < 1973. + 2./12.){
			u = (year - 1975.) / 100.;
			deltaT = MathHelper.polynomial(u, new double[]{45.45, 106.7, -1. / 0.026, -1. / 0.000718});
		}
		else if(year <= 2023. + 10./12.)
			//interpolate the table between 19730201 and 20231001
			deltaT = interpolateModern(year);
//		else if(year <= 1986.){
//			u = (year - 1975.) / 100.;
//			deltaT = MathHelper.polynomial(u, new double[]{45.45, 106.7, -1. / 0.026, -1. / 0.000718});
//		}
//		else if(year <= 2005.){
//			u = (year - 2000.) / 100.;
//			deltaT = MathHelper.polynomial(u, new double[]{63.86, 33.45, -603.74, 1727.5, 65181.4, 237359.9});
//		}
//		else if(year < 2015.){
//			u = (year - 2005.) / 100.;
//			deltaT = MathHelper.polynomial(u, new double[]{64.69, 29.30});
//		}
		else if(year <= 3000.){
			u = (year - 2015.) / 100.;
			deltaT = MathHelper.polynomial(u, new double[]{67.62, 36.45, 39.755});
		}

		if(year < 1955.5){
			final double yat = (year - 1955.5) / 100.;
			deltaT -= 0.91072 * yat * yat * (LUNAR_ACCELERATION + 26.);
		}
		return deltaT;
	}

	private static double interpolateModern(final double year){
		int low = 0;
		int high = DELTA_T_DATA.size() - 1;
		while(low <= high){
			final int mid = (low + high) / 2;

			final double[] dataCurrent = DELTA_T_DATA.get(mid);
			final double[] dataNext = DELTA_T_DATA.get(mid + 1);
			final LocalDate currentDate = LocalDate.of((int)dataCurrent[0], (int)dataCurrent[1], (int)dataCurrent[2]);
			final LocalDate nextDate = LocalDate.of((int)dataNext[0], (int)dataNext[1], (int)dataNext[2]);
			final double yearCurrent = currentDate.getYear() + currentDate.getDayOfYear() / (currentDate.isLeapYear()? 366.: 365.);
			final double yearNext = nextDate.getYear() + nextDate.getDayOfYear() / (nextDate.isLeapYear()? 366.: 365.);

			if(yearCurrent <= year && year <= yearNext){
				final LocalDate targetDate = LocalDate.of((int)year, 1, 1)
					.plusDays((long)(MathHelper.frac(year) * (JulianDate.isLeapYear((int)year)? 366.: 365.)));

				final double startValue = dataCurrent[3];
				final double endValue = dataNext[3];
				final double daysBetween = (double)java.time.temporal.ChronoUnit.DAYS.between(currentDate, nextDate);
				final double daysToTarget = (double)java.time.temporal.ChronoUnit.DAYS.between(currentDate, targetDate);
				return startValue + (endValue - startValue) * (daysToTarget / daysBetween);
			}
			else if(year < yearCurrent)
				high = mid - 1;
			else
				low = mid + 1;
		}
		throw new IllegalArgumentException("Date not in ∆T table: cannot be!");
	}


	/**
	 * Calculate Greenwich Mean Sidereal Time, <code>Θ_{GMST}</code>.
	 *
	 * @param jce	Julian Century of Terrestrial Time from J2000.0.
	 * @return	mean Sidereal time at Greenwich [deg].
	 */
	public static double greenwichMeanSiderealTime(double jce){
//		//FIXME
//		final double[] dateAndTime = JulianDate.extractDateAndTime(jce * JulianDate.CIVIL_SAECULUM + JulianDate.J2000);
//		final double ut0 = JulianDate.centuryJ2000Of(dateAndTime[0]);
//		//[s]
//		final double t = dateAndTime[1];
//
//		//Greenwich Sidereal Time at midnight [day]
//		final double h0 = MathHelper.polynomial(ut0, new double[]{24_110.548_41, 8_640_184.812_866, 0.093_104, -6.2e-6}) / JulianDate.SECONDS_PER_DAY;
//
//		final double earthSiderealRotationRate = earthSiderealRotationRate(ut0);
//		/*
//		This is the difference between UT1 (time using the mean rotating Earth as a clock) and UTC (time that runs at the same rate as
//		an atomic clock, but with leap seconds occasionally inserted to keep UT1 and UTC in sync).
//		Since dUT1 = |UT1 − UTC| < 0.9 s, it is ignored.
//		See <a href="https://en.wikipedia.org/wiki/DUT1">DUT1</a>.
//		*/
//		//[day]
//		final double dUT1 = 0.;
//		//[day]
//		final double h = MathHelper.frac(MathHelper.limitRangeDay(h0 + earthSiderealRotationRate * (t - dUT1)));
//		return h * JulianDate.HOURS_PER_DAY * JulianDate.DEGREES_PER_HOUR;

		//alternative:
		return MathHelper.mod(MathHelper.polynomial(jce, GREENWICH_MEAN_SIDEREAL_TIME_COEFFS), 360.);
	}

	/**
	 * Calculate the Earth's sidereal rotation rate.
	 *
	 * @param ut	Julian Century of Universal Time from J2000.0.
	 * @return	Earth's sidereal rotation rate [sidereal day/UT s].
	 */
	private static double earthSiderealRotationRate(final double ut){
		//FIXME
		//[rad/s]
		final double rate = 7.292_115_855_3e-5 + 4.3e-15 * ut;
		return rate * JulianDate.SECONDS_PER_DAY / MathHelper.TWO_PI;
	}

	/**
	 * Calculate Greenwich Apparent Sidereal Time, <code>Θ_{GAST}</code>.
	 *
	 * @param jce	Julian Century of Terrestrial Time from J2000.0.
	 * @param nutationCorrections	Nutation corrections.
	 * @param trueEclipticObliquity	Obliquity of the ecliptic, corrected for nutation [rad].
	 * @return apparent Sidereal time at Greenwich [deg].
	 *
	 * @see #greenwichApparentSiderealTime(double, double, double, double)
	 */
	public static double greenwichApparentSiderealTime(final double jce, final NutationCorrections nutationCorrections,
			final double trueEclipticObliquity){
		final double moonLongitudeAscendingNode = NutationCorrections.moonLongitudeAscendingNode(jce);
		final double greenwichMeanSiderealTime = greenwichMeanSiderealTime(jce);
		return greenwichApparentSiderealTime(greenwichMeanSiderealTime, nutationCorrections.getDeltaPsi(), trueEclipticObliquity,
			moonLongitudeAscendingNode);
	}

	/**
	 * Calculate Greenwich Apparent Sidereal Time, <code>Θ_{GAST}</code>.
	 *
	 * @param greenwichMeanSiderealTime	Greenwich Mean Sidereal Time [deg].
	 * @param deltaPsi	Nutation in longitude [deg].
	 * @param trueEclipticObliquity	Obliquity of the ecliptic, corrected for nutation [rad].
	 * @param moonLongitudeAscendingNode	Longitude of the ascending node of the Moon [rad].
	 * @return apparent Sidereal time at Greenwich [deg].
	 *
	 * @see <a href="https://www.research.unipd.it/retrieve/e14fb26f-e069-3de1-e053-1705fe0ac030/PhD_Thesis_Zoccarato.pdf">Determinazione Orbitale Precisa di satelliti in orbita LEO per la Radio Occultazione attraverso sistemi GNSS</a>
	 */
	public static double greenwichApparentSiderealTime(final double greenwichMeanSiderealTime, final double deltaPsi,
			final double trueEclipticObliquity, final double moonLongitudeAscendingNode){
		final double equationOfTheEquinoxes = equationOfTheEquinoxes(deltaPsi, trueEclipticObliquity, moonLongitudeAscendingNode);
		return MathHelper.mod(greenwichMeanSiderealTime + equationOfTheEquinoxes, 360.);
	}

	private static double equationOfTheEquinoxes(final double deltaPsi, final double trueEclipticObliquity,
			final double moonLongitudeAscendingNode){
		return StrictMath.toDegrees(deltaPsi) * StrictMath.cos(trueEclipticObliquity)
			+ (0.002_64 * StrictMath.sin(moonLongitudeAscendingNode)
			+ 0.000_063 * StrictMath.sin(moonLongitudeAscendingNode * 2.)) / JulianDate.SECONDS_PER_HOUR;
	}

	/**
	 * Calculate Local Mean Sidereal Time, <code>Θ_{LMST}</code>.
	 *
	 * @param meanSiderealTime	Greenwich Mean Sidereal Time [deg].
	 * @return	The apparent local Sidereal time at Greenwich [deg].
	 */
	public static double localSiderealTime(final double meanSiderealTime, final GeographicLocation location){
		return meanSiderealTime + location.getLongitude();
	}

	/**
	 * Calculate the hour angle of a body, <code>H</code>.
	 *
	 * @param localSiderealTime	Local Sidereal Time [deg].
	 * @param rightAscension	Right ascension [deg].
	 * @return	The hour angle [deg].
	 */
	public static double localHourAngle(final double localSiderealTime, final double rightAscension){
		return MathHelper.mod(localSiderealTime - rightAscension, -180., 180.);
	}

}
