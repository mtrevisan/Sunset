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


/**
 * @see <a href="https://en.wikipedia.org/wiki/Orbital_elements">Orbital elements</>
 * @see <a href="https://en.wikipedia.org/wiki/Mean_motion">Mean motion</>
 * @see <a href="http://www.astro.utoronto.ca/~astrolab/files/AST326_Lab4_2017.pdf">Astrometric orbit determination</a>
 * @see <a href="https://phas.ubc.ca/~newhouse/p210/orbits/cometreport.pdf">Calculating celestial coordinates from orbital parameters using Javascript</a>
 */
public class OrbitalElements{

	//epoch of orbital elements as Julian Ephemeris Date, t [JD].
	public double t;

	//shape and size of the ellipse:
	/** Shape of the ellipse (0 = circular, 1 = parabolic, > 1 = hyperbolic), e. */
	public double eccentricity;
	/** Semimajor axis, a [AU]. */
	public double semimajorAxis;

	//orientation of the orbital plane in which the ellipse is embedded:
	/**
	 * Vertical tilt of the ellipse with respect to the reference plane, measured at the ascending node (where the orbit passes upward
	 * through the reference plane, the green angle i in the diagram), i [rad].
	 * <p>Tilt angle is measured perpendicular to line of intersection between orbital plane and reference plane.</p>
	 */
	public double inclination;
	/**
	 * Horizontally orients the ascending node of the ellipse (where the orbit passes upward through the reference plane) with respect to
	 * the reference frame's vernal point, Ω [rad].
	 */
	public double longitudeAscendingNode;

	/** Orientation of the ellipse in the orbital plane, as an angle measured from the ascending node to the periapsis, ω [rad]. */
	public double argumentOfPerihelion;
	/**
	 * Fraction of an elliptical orbit's period that has elapsed since the orbiting body passed periapsis, expressed as an angle which can
	 * be used in calculating the position of that body in the classical two-body problem. It is the angular distance from the pericenter
	 * which a fictitious body would have if it moved in a circular orbit, with constant speed, in the same orbital period as the actual
	 * body in its elliptical orbit, M [rad].
	 * <p>
	 * True anomaly, v, is meanAnomaly + equationOfCenter.
	 * </p>
	 */
	public double meanAnomaly;

	/**
	 *  Angular speed required for a body to complete one orbit, assuming constant speed in a circular orbit which completes in the same
	 *  time as the variable speed, elliptical orbit of the actual body, n [rad/day].
	 */
	public double meanMotion;

}
