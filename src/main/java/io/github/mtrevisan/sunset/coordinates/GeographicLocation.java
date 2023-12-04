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

import io.github.mtrevisan.sunset.StringHelper;


/**
 * Describes a geographical position.
 */
public final class GeographicLocation{

	//[deg]
	private final double latitude;
	//[deg]
	private final double longitude;
	//[m]
	private final double altitude;

	private AtmosphericModel atmosphere;


	/**
	 * Creates a new instance with the given parameters.
	 *
	 * @param latitude	The latitude of this location [deg]. North latitude is positive, south negative.
	 * @param longitude	The longitude of this location [deg]. East longitude is positive, west negative.
	 * @return	An instance.
	 */
	public static GeographicLocation create(final double latitude, final double longitude){
		return create(latitude, longitude, 0.);
	}

	/**
	 * Creates a new instance with the given parameters.
	 *
	 * @param latitude	The latitude of this location [°]. North latitude is positive, south negative.
	 * @param longitude	The longitude of this location [°]. East longitude is positive, west negative.
	 * @param altitude	The altitude of this location [m].
	 * @return	An instance.
	 */
	public static GeographicLocation create(final double latitude, final double longitude, final double altitude){
		return new GeographicLocation(latitude, longitude, altitude);
	}


	private GeographicLocation(final double latitude, final double longitude, final double altitude){
		if(!Double.isFinite(latitude))
			throw new IllegalArgumentException("Latitude must be a finite value: " + latitude);
		if(!Double.isFinite(longitude))
			throw new IllegalArgumentException("Longitude must be a finite value: " + longitude);
		if(latitude > 90. || latitude < -90.)
			throw new IllegalArgumentException("Degrees out of range -90 <= latitude <= +90: " + latitude);
		if(longitude >= 180. || longitude < -180.)
			throw new IllegalArgumentException("Degrees out of range -180 <= longitude < +180: " + longitude);
		if(altitude < -11_000. || altitude >= 9_000.)
			throw new IllegalArgumentException("Meters out of range -11000 <= altitude < +9000: " + altitude);

		this.latitude = latitude;
		this.longitude = longitude;
		this.altitude = altitude;
	}

	/**
	 * @param atmosphere	The atmospheric model of the location.
	 * @return	This instance.
	 */
	public GeographicLocation withAtmosphere(final AtmosphericModel atmosphere){
		this.atmosphere = atmosphere;

		return this;
	}


	/**
	 * The latitude of the location.
	 *
	 * @return	The latitude [deg].
	 */
	public double getLatitude(){
		return latitude;
	}

	/**
	 * The longitude of the location.
	 *
	 * @return	The longitude [deg].
	 */
	public double getLongitude(){
		return longitude;
	}

	/**
	 * The altitude of the location.
	 *
	 * @return	The altitude [m].
	 */
	public double getAltitude(){
		return altitude;
	}

	/**
	 * The atmospheric model of the location.
	 *
	 * @return	The atmospheric model.
	 */
	public AtmosphericModel getAtmosphere(){
		return atmosphere;
	}


	@Override
	public String toString(){
		return "Location{"
			+ "lat: " + StringHelper.degreeToDegMinSecString(latitude, 2)
			+ ", lon: "  + StringHelper.degreeToDegMinSecString(longitude, 2)
			+ ", alt: "  + altitude + " m"
			+ '}';
	}

}
