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

import java.math.BigDecimal;


public final class Location{

	private final BigDecimal latitude;
	private final BigDecimal longitude;


	/**
	 * Creates a new instance of <code>Location</code> with the given parameters.
	 *
	 * @param latitude	The latitude of this location [°]. North latitude is positive, south negative.
	 * @param longitude	The longitude of this location [°]. East longitude is positive, west negative.
	 * @return	An instance of a location.
	 */
	public static Location create(final String latitude, final String longitude){
		return new Location(latitude, longitude);
	}

	/**
	 * Creates a new instance of <code>Location</code> with the given parameters.
	 *
	 * @param latitude	The latitude of this location [°]. North latitude is positive, south negative.
	 * @param longitude	The longitude of this location [°]. East longitude is positive, west negative.
	 * @return	An instance of a location.
	 */
	public static Location create(final double latitude, final double longitude){
		return new Location(latitude, longitude);
	}


	private Location(final String latitude, final String longitude){
		this.latitude = new BigDecimal(latitude);
		this.longitude = new BigDecimal(longitude);
	}

	private Location(final double latitude, final double longitude){
		this.latitude = new BigDecimal(latitude);
		this.longitude = new BigDecimal(longitude);
	}


	/**
	 * The latitude of the location.
	 *
	 * @return	The latitude.
	 */
	public BigDecimal getLatitude(){
		return latitude;
	}

	/**
	 * The longitude of the location.
	 *
	 * @return	The longitude.
	 */
	public BigDecimal getLongitude(){
		return longitude;
	}

}
