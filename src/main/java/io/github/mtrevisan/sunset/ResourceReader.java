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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;


public final class ResourceReader{

	private ResourceReader(){}


	public static Map<String, Collection<Double[]>> read(final String filename) throws IOException{
		final ClassLoader classLoader = ResourceReader.class.getClassLoader();
		final InputStream is = classLoader.getResourceAsStream(filename);
		if(is == null)
			throw new IllegalArgumentException("file not found! " + filename);

		try(
				final InputStreamReader sr = new InputStreamReader(is, StandardCharsets.UTF_8);
				final BufferedReader reader = new BufferedReader(sr)){
			final Map<String, Collection<Double[]>> result = new HashMap<>(0);
			String key = null;
			String line;
			Collection<Double[]> values = new LinkedList<>();
			while((line = reader.readLine()) != null){
				if(line.isEmpty()){
					result.put(key, values);
					key = null;

					continue;
				}
				if(key == null){
					key = line;
					values = new LinkedList<>();

					continue;
				}

				final String[] parameters = line.split(" ");
				values.add(Arrays.stream(parameters)
					.map(Double::valueOf)
					.toArray(Double[]::new));
			}
			if(key != null)
				result.put(key, values);
			return result;
		}
	}

}
