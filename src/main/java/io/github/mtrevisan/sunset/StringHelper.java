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

import org.slf4j.helpers.MessageFormatter;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;


/**
 * A collection of convenience methods for working with {@link String} objects.
 */
public final class StringHelper{

	private StringHelper(){}


	public static String degreeToHMSString(double degree, final int precision){
		degree /= 15.;
		final int hour = (int)degree;
		degree -= hour;
		degree *= 60.;
		final int minute = (int)degree;
		degree -= minute;
		degree *= 60.;
		final double second = degree;
		final DecimalFormat df = decimalFormat(precision);
		return hour + "h " + minute + "m " + df.format(second) + "s";
	}

	@SuppressWarnings("NumericCastThatLosesPrecision")
	public static String degreeToDegMinSecString(final double degree, final int precision){
		final int hour = (int)degree;
		final double deg = Math.abs((degree - hour) * 60.);
		final int minute = (int)deg;
		final double second = (deg - minute) * 60.;
		final DecimalFormat df = decimalFormat(precision);
		return hour + "° " + minute + "' " + df.format(second) + "\"";
	}

	private static DecimalFormat decimalFormat(final int precision){
		final DecimalFormat df = (DecimalFormat)NumberFormat.getNumberInstance(Locale.US);
		df.setMaximumFractionDigits(precision);
		return df;
	}

}
