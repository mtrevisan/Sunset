/**
 * Copyright (c) 2021 Mauro Trevisan
 * <p>
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 * <p>
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * <p>
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


public final class AtmosphereHelper{

	/** [°C] */
	static final double ABSOLUTE_ZERO = 273.15;


	private AtmosphereHelper(){}


	/**
	 * Calculate the atmospheric refraction correction, Δe.
	 *
	 * @param pressure	The pressure [hPa].
	 * @param temperature	The temperature [°C].
	 * @param e0	The topocentric elevation angle without atmospheric refraction correction [°].
	 * @return	The correction [°].
	 */
	static double atmosphericRefractionCorrection(final double pressure, final double temperature, final double e0){
		return (pressure / 1010.)
			* ((ABSOLUTE_ZERO + 10.) / (ABSOLUTE_ZERO + temperature))
			* (1.02 / (60. * StrictMath.tan(StrictMath.toRadians(e0 + 10.3 / (e0 + 5.11)))));
	}

}
