package io.github.mtrevisan.astro.core;


/**
 * Interface defining methods for calculating or retrieving the solar elevation angle.
 * <p>
 * The solar elevation is often used in astronomical or geographic computations to determine the height of the Sun
 * above the horizon.
 * </p>
 */
public interface ZenithInterface{

	/**
	 * @return	The elevation [rad].
	 */
	double getElevation();

}
