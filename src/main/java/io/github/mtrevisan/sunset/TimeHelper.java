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

}
