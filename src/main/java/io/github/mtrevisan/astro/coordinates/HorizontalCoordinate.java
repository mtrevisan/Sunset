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
package io.github.mtrevisan.astro.coordinates;

import io.github.mtrevisan.astro.helpers.StringHelper;

import java.text.DecimalFormat;


/**
 * @see <a href="https://en.wikipedia.org/wiki/Horizontal_coordinate_system">Horizontal coordinate system</>
 */
public final class HorizontalCoordinate{

	//[rad]
	private final double azimuth;
	//[rad]
	private final double altitude;
	//[AU]
	private final double radius;


	/**
	 * Creates a new instance with the given parameters.
	 *
	 * @param azimuth	The azimuth [rad].
	 * @param altitude	The altitude [rad].
	 * @param radius	The radius [AU].
	 * @return	An instance.
	 */
	public static HorizontalCoordinate create(final double azimuth, final double altitude, final double radius){
		return new HorizontalCoordinate(azimuth, altitude, radius);
	}


	private HorizontalCoordinate(final double azimuth, final double altitude, final double radius){
		if(Double.isNaN(azimuth) || Double.isInfinite(azimuth))
			throw new IllegalArgumentException("Not finite azimuth: " + azimuth);
		if(Double.isNaN(altitude) || Double.isInfinite(altitude))
			throw new IllegalArgumentException("Not finite altitude: " + altitude);
		if(Double.isNaN(radius) || Double.isInfinite(radius))
			throw new IllegalArgumentException("Not finite radius: " + radius);

		this.azimuth = azimuth;
		this.altitude = altitude;
		this.radius = radius;
	}


	/**
	 * The azimuth.
	 *
	 * @return	The azimuth [rad].
	 */
	public double getAzimuth(){
		return azimuth;
	}

	/**
	 * The altitude.
	 *
	 * @return	The altitude [rad].
	 */
	public double getAltitude(){
		return altitude;
	}

	/**
	 * The radius.
	 *
	 * @return	The radius [AU].
	 */
	public double getRadius(){
		return radius;
	}

	@Override
	public String toString(){
		final DecimalFormat formatter = StringHelper.decimalFormat(2);
		return "HorizontalCoordinate{"
			+ "azi: " + StringHelper.degreeToDegMinSecString(StrictMath.toDegrees(azimuth), 2)
			+ ", alt: " + StringHelper.degreeToDegMinSecString(StrictMath.toDegrees(altitude), 2)
			+ ", r: " + formatter.format(radius)
			+ '}';
	}

}
