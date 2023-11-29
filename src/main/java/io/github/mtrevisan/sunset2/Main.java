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
package io.github.mtrevisan.sunset2;

import io.github.mtrevisan.sunset2.coordinates.GeographicalLocation;

import java.time.Instant;
import java.time.ZoneId;
import java.time.ZonedDateTime;


public class Main{

	public static void main(final String[] args){
		final double latitude = 45.651727;
		final double longitude = 12.212561;
		final double elevation = 16.;

		SolarEventCalculator solarEventCalculator = new SolarEventCalculator(GeographicalLocation.of(latitude, longitude, elevation),
			ZoneId.of("Europe/Rome"));
		ZonedDateTime officialSunset = solarEventCalculator.sunset(Zenith.CIVIL, ZonedDateTime.now());
		System.out.println(officialSunset.toLocalTime());


//		final ZoneId zoneId = ZoneId.of("Europe/Rome");
//		final ZonedDateTime currentTime = ZonedDateTime.now(zoneId);
//
//		calc(currentTime.toEpochSecond(), latitude, longitude, elevation, zoneId);
	}

	private static void calc(long currentTimestamp, double f, double l_w, double elevation, ZoneId zoneId){
		System.out.println("Latitude f = " + f);
		System.out.println("Longitude l_w = " + l_w);
		System.out.println("Now = " + ZonedDateTime.ofInstant(Instant.ofEpochSecond(currentTimestamp), zoneId));

		double jDate = ts2j(currentTimestamp);
		System.out.println("Julian date j_date = " + String.format("%.3f days", jDate));

		double n = Math.ceil(jDate - l_w / 360.);
		System.out.println("Mean solar time J_ = " + String.format("%.9f days", n));

		double mDegrees = Math.toDegrees(Math.toRadians(357.5291) + 0.98560028 * n) % 360;
		System.out.println("Solar mean anomaly M = " + mDegrees);

		double cDegrees = 1.9148 * Math.sin(Math.toRadians(mDegrees)) + 0.02 * Math.sin(Math.toRadians(2 * mDegrees)) + 0.0003 * Math.sin(Math.toRadians(3 * mDegrees));
		System.out.println("Equation of the center C = " + cDegrees);

		double lDegrees = (mDegrees + cDegrees + 180. + 102.9372) % 360;
		System.out.println("Ecliptic longitude L = " + lDegrees);

		double jTransit = 2451545. + n + 0.0053 * Math.sin(Math.toRadians(mDegrees)) - 0.0069 * Math.sin(Math.toRadians(2 * lDegrees));
		System.out.println("Solar transit time J_trans = " + ZonedDateTime.ofInstant(Instant.ofEpochSecond((long)jTransit), zoneId));

		double sinD = Math.sin(Math.toRadians(lDegrees)) * Math.sin(Math.toRadians(23.4397));
		double cosD = Math.cos(Math.asin(sinD));

		double someCos = (Math.sin(Math.toRadians(- 0.833 - 2.076 * Math.sqrt(elevation) / 60.)) - Math.sin(Math.toRadians(f)) * sinD) / (Math.cos(Math.toRadians(f)) * cosD);
		double w0Degrees = Math.toDegrees(Math.acos(someCos));

		double jRise = jTransit - w0Degrees / 360;
		double jSet = jTransit + w0Degrees / 360;

		System.out.println("Sunrise j_rise = " + ZonedDateTime.ofInstant(Instant.ofEpochSecond((long)jRise), zoneId));
		System.out.println("Sunset j_set = " + ZonedDateTime.ofInstant(Instant.ofEpochSecond((long)jSet), zoneId));
	}

	private static double ts2j(long ts){
		return (ts / 86400.) + 2440587.5;
	}

}
