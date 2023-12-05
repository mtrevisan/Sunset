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
package io.github.mtrevisan.astro;

import java.io.Serial;


/**
 * Thrown if the Sun never rises or never sets.
 */
public class SunlightPhaseException extends Exception{

	@Serial
	private static final long serialVersionUID = 3660516823104741068L;


	/** Type of error. */
	private final SunlightPhaseError error;


	/**
	 * Constructs a new exception with the specified error.
	 *
	 * @param error	The error.
	 * @return	An instance of this exception.
	 */
	public static SunlightPhaseException create(final SunlightPhaseError error){
		return new SunlightPhaseException(error);
	}


	/**
	 * Constructs a new exception with the specified error.
	 *
	 * @param error	The error.
	 */
	private SunlightPhaseException(final SunlightPhaseError error){
		this.error = error;
	}

}
