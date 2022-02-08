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

import io.github.mtrevisan.sunset.coordinates.GNSSLocation;

import java.time.LocalTime;


public final class TimeHelper{

	/** Lunar acceleration parameter ["/cy^2]. */
	private static final double LUNAR_ACCELERATION = -25.7376;

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
	 * @see #deltaT(int)
	 */
	public static double universalTimeToTerrestrialTime(final double ut, final double deltaT){
		return ut + deltaT / JulianDay.SECONDS_IN_DAY;
	}

	/**
	 * Calculate the predicted difference between the Earth rotation time and the TT.
	 *
	 * @param year	The year.
	 * @return	The predicted {@code TT – UT1}.
	 */
	public static double deltaT(final int year){
		double yat2 = 0.;
		if(year < 1955.5){
			final double yat = (year - 1955.5) / 100.;
			yat2 = 0.91072 * yat * yat;
		}

		double deltaT = 0;
		final double u;
		if(year <= -500 || year > 2150){
			u = (year - 1820) / 100.;
			deltaT = -20. + 32. * u * u;
		}
		else if(year > -500 && year <= 500){
			u = year / 100.;
			deltaT = MathHelper.eval(u, new double[]{10583.6, -1014.41, 33.78311, -5.952053, -0.1798452, 0.022174192, 0.0090316521});
		}
		else if(year > 500 && year <= 1600){
			u = (year - 1000) / 100.;
			deltaT = MathHelper.eval(u, new double[]{1574.2, -556.01, 71.23472, 0.319781, -0.8503463, -0.005050998, 0.0083572073});
		}
		else if(year > 1600 && year <= 1700){
			u = (year - 1600) / 100.;
			deltaT = MathHelper.eval(u, new double[]{120, -98.08, -153.2, 1. / 0.007129});
		}
		else if(year > 1700 && year <= 1800){
			u = (year - 1700) / 100.;
			deltaT = MathHelper.eval(u, new double[]{8.83, 16.03, -59.285, 133.36, -1. / 0.01174});
		}
		else if(year > 1800 && year <= 1860){
			u = (year - 1800) / 100.;
			deltaT = MathHelper.eval(u, new double[]{13.72, -33.2447, 68.612, 4111.6, -37436, 121272, -169900, 87500});
		}
		else if(year > 1860 && year <= 1900){
			u = (year - 1860) / 100.;
			deltaT = MathHelper.eval(u, new double[]{7.62, 57.37, -2517.54, 16806.68, - 44736.24, 1. / 0.0000233174});
		}
		else if(year > 1900 && year <= 1920){
			u = (year - 1900) / 100.;
			deltaT = MathHelper.eval(u, new double[]{-2.79, 149.4119, -598.939, 6196.6, -19700.});
		}
		else if(year > 1920 && year <= 1941){
			u = (year - 1920) / 100.;
			deltaT = MathHelper.eval(u, new double[]{21.20, 84.493, -761.00, 2093.6});
		}
		else if(year > 1941 && year <= 1961){
			u = (year - 1950) / 100.;
			deltaT = MathHelper.eval(u, new double[]{29.07, 40.7, -1. / 0.0233, 1. / 0.002547});
		}
		else if(year > 1961 && year <= 1986){
			u = (year - 1975) / 100.;
			deltaT = MathHelper.eval(u, new double[]{45.45, 106.7, -1. / 0.026, -1. / 0.000718});
		}
		else if(year > 1986 && year <= 2005){
			u = (year - 2000) / 100.;
			deltaT = MathHelper.eval(u, new double[]{63.86, 33.45, -603.74, 1727.5, 65181.4, 237359.9});
		}
		else if(year > 2005 && year <= 2050){
			u = (year - 2000) / 100.;
			deltaT = MathHelper.eval(u, new double[]{62.92, 32.217, 55.89});
		}
		else if(year > 2050 && year <= 2150){
			u = (year - 1820) / 100.;
			deltaT = MathHelper.eval(u, new double[]{-205.72, 56.28, 32});
		}
		return deltaT - yat2 * (LUNAR_ACCELERATION + 26.);
	}


	/**
	 * Calculate Greenwich Mean Sidereal Time, ΘGMST.
	 *
	 * @param ut	Julian Day of Universal Time from J2000.0.
	 * @return	mean Sidereal time at Greenwich [°].
	 */
	static double meanSiderealTime(final double ut){
		final double[] dateAndTime = JulianDay.extractDateAndTime(ut);
		final double ut0 = JulianDay.centuryJ2000Of(dateAndTime[0]);
		//[s]
		final double t = dateAndTime[1] * JulianDay.SECONDS_IN_DAY;

		//Greenwich Sidereal Time at midnight [day]
		final double h0 = MathHelper.eval(ut0, new double[]{24110.54841, 8640184.812866, 0.093104, -6.2e-6}) / JulianDay.SECONDS_IN_DAY;

		final double earthSiderealRotationRate = earthSiderealRotationRate(ut0);
		/*
		This is the difference between UT1 (time using the mean rotating Earth as a clock) and UTC (time that runs at the same rate as
		an atomic clock, but with leap seconds occasionally inserted to keep UT1 and UTC in sync).
		Since dUT1 = |UT1 − UTC| < 0.9 s, it is ignored.
		See <a href="https://en.wikipedia.org/wiki/DUT1">DUT1</a>.
		*/
		//[s]
		final double dUT1 = 0.;
		//[day]
		final double h = MathHelper.frac(MathHelper.limitRangeDay(h0 + earthSiderealRotationRate * (t - dUT1)));
		return h * JulianDay.HOURS_IN_DAY * JulianDay.DEGREES_PER_HOUR;
		//alternative:
		//return MathHelper.limitRangeDegree(MathHelper.eval(JulianDay.centuryJ2000Of(ut), new double[]{280.46061837, 360.98564736629 * JulianDay.CIVIL_SAECULUM, 0.000387933, -1. / 38710000.}));
	}

	/**
	 * Calculate the Earth's sidereal rotation rate.
	 *
	 * @param ut	Julian Century of Universal Time from J2000.0.
	 * @return	Earth's sidereal rotation rate [sidereal day/UT s].
	 */
	private static double earthSiderealRotationRate(final double ut){
		//[rad/s]
		final double rate = 7.2921158553e-5 + 4.3e-15 * ut;
		return rate / (2. * StrictMath.PI);
	}

	/**
	 * Calculate Greenwich Apparent Sidereal Time, ΘGAST.
	 *
	 * @param meanSiderealTime	Greenwich Mean Sidereal Time [°].
	 * @param trueEclipticObliquity	Obliquity of the ecliptic, corrected for nutation [°].
	 * @param deltaPsi	Nutation in longitude [°].
	 * @return	apparent Sidereal time at Greenwich [°].
	 */
	static double apparentSiderealTime(final double meanSiderealTime, final double trueEclipticObliquity, final double deltaPsi){
		final double equationOfTheEquinoxes = deltaPsi * StrictMath.cos(StrictMath.toRadians(trueEclipticObliquity));
		return MathHelper.limitRangeDegree(meanSiderealTime + equationOfTheEquinoxes);
	}

	/**
	 * Calculate Local Mean Sidereal Time, ΘLMST.
	 *
	 * @param meanSiderealTime	Greenwich Mean Sidereal Time [°].
	 * @return	The apparent local Sidereal time at Greenwich [°].
	 */
	static double localMeanSiderealTime(final double meanSiderealTime, final GNSSLocation location){
		return meanSiderealTime + location.getLongitude();
	}

	/**
	 * Calculate the hour angle of a body, H.
	 *
	 * @param localSiderealTime	Local Sidereal Time [°].
	 * @param rightAscension	Right ascension [°].
	 * @return	The hour angle [°].
	 */
	static double localHourAngle(final double localSiderealTime, final double rightAscension){
		return localSiderealTime - rightAscension;
	}

}
