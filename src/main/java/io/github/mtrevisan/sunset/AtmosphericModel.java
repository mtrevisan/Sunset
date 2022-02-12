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
package io.github.mtrevisan.sunset;


/**
 * @see <a href="https://en.wikipedia.org/wiki/International_Standard_Atmosphere">International Standard Atmosphere</>
 */
public final class AtmosphericModel{

	/** [°C] */
	private static final double ABSOLUTE_ZERO = 273.15;


	//[hPa]
	private final double pressure;
	//[°C]
	private final double temperature;


	/**
	 * Creates a new instance with the given parameters.
	 *
	 * @param pressure	The pressure [hPa].
	 * @param temperature	The temperature [°C].
	 * @return	An instance.
	 */
	public static AtmosphericModel create(final double pressure, final double temperature){
		return new AtmosphericModel(pressure, temperature);
	}


	private AtmosphericModel(final double pressure, final double temperature){
		if(Double.isNaN(pressure) || Double.isInfinite(pressure))
			throw new IllegalArgumentException("Not finite pressure: " + pressure);
		if(Double.isNaN(temperature) || Double.isInfinite(temperature) || temperature < ABSOLUTE_ZERO)
			throw new IllegalArgumentException("Not finite temperature: " + temperature);

		this.pressure = pressure;
		this.temperature = temperature;
	}


	/**
	 * The pressure.
	 *
	 * @return	The pressure [hPa].
	 */
	public double getPRessure(){
		return pressure;
	}

	/**
	 * The temperature.
	 *
	 * @return	The temperature [°C].
	 */
	public double getTemperature(){
		return temperature;
	}


	/**
	 * Calculate the atmospheric refraction correction, Δe.
	 *
	 * @param elevation	The topocentric elevation angle without atmospheric refraction correction [rad].
	 * @return	The correction [rad].
	 */
	public double atmosphericRefractionCorrection(final double elevation){
		final double elev = StrictMath.toDegrees(elevation);
		return StrictMath.toRadians((pressure / 1010.)
			* ((ABSOLUTE_ZERO + 10.) / (ABSOLUTE_ZERO + temperature))
			* (1.02 / (60. * StrictMath.tan(StrictMath.toRadians(elev + 10.3 / (elev + 5.11)))))
		);
	}


	@Override
	public String toString(){
		return "AtmosphericModel{"
			+ "pressure: " + StringHelper.degreeToHMSString(StrictMath.toDegrees(pressure), 2) + " hPa"
			+ ", temperature: " + StringHelper.degreeToDegMinSecString(temperature, 1) + " °C"
			+ '}';
	}

}
