package io.github.mtrevisan.sunset2.coordinates;


/**
 * Describes a geographical position.
 */
public class GeographicalLocation{

	//[deg]
	private final double latitude;
	//[deg]
	private final double longitude;
	//[m]
	private final double altitude;


	/**
	 * Creates a new instance of <code>Location</code> with the given parameters.
	 *
	 * @param latitude	The latitude, in degrees, of this location. North latitude is positive, south negative.
	 * @param longitude	The longitude, in degrees of this location. East longitude is positive, west negative.
	 */
	public static GeographicalLocation of(final double latitude, final double longitude){
		return new GeographicalLocation(latitude, longitude, 0.);
	}

	/**
	 * Creates a new instance of <code>Location</code> with the given parameters.
	 *
	 * @param latitude	The latitude, in degrees, of this location. North latitude is positive, south negative.
	 * @param longitude	The longitude, in degrees of this location. East longitude is positive, west negative.
	 * @param altitude	The altitude, in meters of this location.
	 */
	public static GeographicalLocation of(final double latitude, final double longitude, final double altitude){
		return new GeographicalLocation(latitude, longitude, altitude);
	}


	private GeographicalLocation(final double latitude, final double longitude, final double altitude){
		this.latitude = latitude;
		this.longitude = longitude;
		this.altitude = altitude;
	}

	public double getLatitude(){
		return latitude;
	}

	public double getLongitude(){
		return longitude;
	}

}
