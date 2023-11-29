package io.github.mtrevisan.sunset2;


/**
 * Defines the solar declination.
 */
public enum Zenith{
	/** Official sunrise/set is when the sun is 50' below the horizon (to account for refraction). */
	OFFICIAL(-50. / 60.),
	/** Civil sunrise/set (dawn/dusk) is when the sun is 6 degrees below the horizon. */
	CIVIL(-6.),
	/** Nautical sunrise/set is when the sun is 12 degrees below the horizon. */
	NAUTICAL(-12.),
	/** Astronomical sunrise/set is when the sun is 18 degrees below the horizon. */
	ASTRONOMICAL(-18.);


	/** Solar declination [rad]. */
	private final double angle;


	Zenith(final double degrees){
		angle = degrees;
	}

	public double getAngle(){
		return angle;
	}

}
