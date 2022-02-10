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
package io.github.mtrevisan.sunset.coordinates;

import io.github.mtrevisan.sunset.MathHelper;
import io.github.mtrevisan.sunset.StringHelper;


/**
 * Describes a (polar) celestial coordinate system which projects the Earth equator and poles onto the celestial sphere
 * using right ascension and declination in the reference frame J2000 as epoch.
 *
 * <p>
 * See also: <a href="https://en.wikipedia.org/wiki/Celestial_coordinate_system">Celestial coordinate system</a>.
 * The effect of precession will have an impact on right ascension and declination if times far away
 * from year 2000 are considered.
 * </p>
 *
 * @see <a href="https://en.wikipedia.org/wiki/Equatorial_coordinate_system">Equatorial coordinate system</>
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

	/**
	 * Creates a new instance from ecliptical coordinates.
	 *
	 * @param eclipticLatitude   Geometric ecliptic latitude [°].
	 * @param eclipticLongitude   Geometric ecliptic longitude [°].
	 * @param eclipticObliquity   Obliquity of the ecliptic [°].
	 * @return	An instance.
	 */
	public static EquatorialCoordinate createFromEcliptical(double eclipticLatitude, double eclipticLongitude, double eclipticObliquity){
		eclipticLatitude = StrictMath.toRadians(eclipticLatitude);
		eclipticLongitude = StrictMath.toRadians(eclipticLongitude);
		eclipticObliquity = StrictMath.toRadians(eclipticObliquity);

		final double sinLat = StrictMath.sin(eclipticLatitude);
		final double cosLat = StrictMath.cos(eclipticLatitude);
		final double sinLon = StrictMath.sin(eclipticLongitude);
		final double cosLon = StrictMath.cos(eclipticLongitude);
		final double cosObl = StrictMath.cos(eclipticObliquity);
		final double sinObl = StrictMath.sin(eclipticObliquity);

		final double rightAscension = MathHelper.limitRangeDegree(StrictMath.toDegrees(
			StrictMath.atan2(sinLon * cosObl - sinLat * sinObl / cosLat, cosLon)
		));
		final double declination = StrictMath.toDegrees(
			StrictMath.asin(sinLat * cosObl + cosLat * sinObl * sinLon)
		);

		return create(rightAscension, declination);
	}


	private EquatorialCoordinate(final double rightAscension, final double declination){
		if(Double.isNaN(rightAscension) || Double.isInfinite(rightAscension))
			throw new IllegalArgumentException("Not finite right ascension: " + rightAscension);
		if(Double.isNaN(declination) || Double.isInfinite(declination))
			throw new IllegalArgumentException("Not finite declination: " + declination);

		this.rightAscension = rightAscension;
		this.declination = declination;
	}


	/**
	 * The right ascension.
	 *
	 * @return	The right ascension [°].
	 */
	public double getRightAscension(){
		return rightAscension;
	}

	/**
	 * The declination.
	 *
	 * @return	The declination [°].
	 */
	public double getDeclination(){
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
