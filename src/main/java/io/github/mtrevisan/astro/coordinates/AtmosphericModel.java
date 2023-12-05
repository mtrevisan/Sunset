/*
 * Copyright (c) 2023 Mauro Trevisan
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
package io.github.mtrevisan.astro.coordinates;

import io.github.mtrevisan.astro.helpers.StringHelper;


/**
 * @see <a href="https://en.wikipedia.org/wiki/International_Standard_Atmosphere">International Standard Atmosphere</>
 */
public final class AtmosphericModel{

	/** [hPa] */
	public static final double STANDARD_PRESSURE = 1013.25;

	/** [°C] */
	private static final double ABSOLUTE_ZERO = -273.15;


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
	public double getPressure(){
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
	 * Calculate the atmospheric refraction correction, <code>Δe</code>.
	 *
	 * @param elevation	The observed topocentric elevation angle without atmospheric refraction correction [rad].
	 * @return	The correction [rad].
	 *
	 * @see <a href="https://digitalcommons.mtu.edu/etdr/697/">Evaluating the Effectiveness of Current Atmospheric Refraction Models in Predicting Sunrise and Sunset Times</a>
	 * @see <a href="https://github.com/Starlink/starjava/blob/master/pal/src/main/uk/ac/starlink/pal/Pal.java">Pal.java</a>
	 */
	public double atmosphericRefractionCorrection(final double elevation){
		final double elev = StrictMath.toDegrees(elevation);
		//Bennett-NA's equation [deg]
		final double deltaElevation = 0.0167 / StrictMath.tan(StrictMath.toRadians(elev + 7.32 / (elev + 4.32)))
			* (pressure / 1010.) * ((10. - ABSOLUTE_ZERO) / (temperature - ABSOLUTE_ZERO));
		return StrictMath.toRadians(deltaElevation);
	}


	@Override
	public String toString(){
		return "AtmosphericModel{"
			+ "pressure: " + StringHelper.decimalFormat(2).format(pressure) + " hPa"
			+ ", temperature: " + StringHelper.decimalFormat(1).format(temperature) + " °C"
			+ '}';
	}

}
