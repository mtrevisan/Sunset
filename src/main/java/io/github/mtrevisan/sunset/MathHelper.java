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


public final class MathHelper{

	private static final double TWO_PI = StrictMath.PI * 2.;


	private MathHelper(){}


	public static double frac(final double x){
		return x - (int)x;
	}

	// Reduces an angle in degrees to the range 0 to 2.0 * M_PI.
	public static double mod2pi(final double rad){
		return rad - TWO_PI * StrictMath.floor(rad / TWO_PI);
	}

	public static double toDegrees(final int degree, final int minute, final double second){
		return degree + (minute + second / 60.) / 60.;
	}

	public static double limitRangeDegree(double degree){
		degree %= 360.;
		return (degree < 0.? degree + 360.: degree);
	}

	public static double limitRangeHour(double degree){
		degree %= 24;
		return (degree < 0.? degree + 24.: degree);
	}

	public static double limitRangeDay(double value){
		value %= 1;
		return (value < 0.? value + 1.: value);
	}

	public static double degToHrs(final double degrees){
		return degrees / JulianDay.DEGREES_PER_HOUR;
	}

	/**
	 * Use Horner's method to compute and return the polynomial evaluated at {@code x}:<br/>
	 * {@code p[0] + p[1] * x^1 + p[2] * x^2 + ... + p[n-1] * x^n-1}
	 *
	 * @param x	The value at which to calculate the polynomial.
	 * @param p	The polynomial coefficients.
	 * @return	The value of the polynomial.
	 */
	public static double eval(final double x, final double[] p){
		double result = 0.;
		for(int i = p.length - 1; i >= 0; i --)
			result = p[i] + result * x;
		return result;
	}

}
