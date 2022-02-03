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


/**
 * @see <a href="https://en.wikipedia.org/wiki/Equatorial_coordinate_system">Equatorial Coordinate System</>
 */
public final class EquatorialCoordinate{

	//[°]
	private final double rightAscension;
	//[°]
	private final double declination;


	/**
	 * Creates a new instance with the given parameters.
	 *
	 * @param rightAscension	The right ascension [°].
	 * @param declination	The declination [°].
	 * @return	An instance.
	 */
	public static EquatorialCoordinate create(final double rightAscension, final double declination){
		return new EquatorialCoordinate(rightAscension, declination);
	}


	private EquatorialCoordinate(final double rightAscension, final double declination){
		this.rightAscension = rightAscension;
		this.declination = declination;
	}


	/**
	 * The right ascension.
	 *
	 * @return	The right ascension.
	 */
	public double getRightAscension(){
		return rightAscension;
	}

	/**
	 * The declination.
	 *
	 * @return	The declination.
	 */
	public double getLongitude(){
		return declination;
	}

	@Override
	public String toString(){
		return "EquatorialCoordinate{"
			+ "α: " + StringHelper.degreeToHMSString(rightAscension, 2)
			+ ", δ: " + StringHelper.degreeToDegMinSecString(declination, 2)
			+ '}';
	}

}
