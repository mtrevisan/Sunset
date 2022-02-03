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


/** Defines the solar declination used in computing the sunrise/sunset. */
public enum Zenith{
	/** Official sunrise/set is when the sun is 50' below the horizon (to account for refraction). */
	OFFICIAL(0.5),
	/** Civil sunrise/set (dawn/dusk) is when the sun is 6 degrees below the horizon. */
	CIVIL(6.),
	/** Nautical sunrise/set is when the sun is 12 degrees below the horizon. */
	NAUTICAL(12.),
	/** Astronomical sunrise/set is when the sun is 18 degrees below the horizon. */
	ASTRONOMICAL(18.);


	/** Solar declination [°]. */
	private final double radians;


	Zenith(final double degrees){
		radians = convertDegreesToRadians(90. + degrees);
	}

	public double getRadians(){
		return radians;
	}

	private double convertDegreesToRadians(final double degrees){
		return degrees * (Math.PI / 180.);
	}

}
