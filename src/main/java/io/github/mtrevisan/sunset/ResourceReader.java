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

import org.apache.commons.lang3.StringUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;


public final class ResourceReader{

	private ResourceReader(){}


	public enum VariableIndex{
		//semi-major axis [UA]
		A(1),
		//mean longitude [rad]
		L(2),
		//e * cos(v)
		K(3),
		//e * sin(v)
		H(4),
		//sin(i/2) * cos(W)
		Q(5),
		//sin(i/2) * sin(W)
		P(6);

		private int index;
		VariableIndex(final int index){
			this.index = index;
		}

		public static VariableIndex of(final int idx){
			for(final VariableIndex vi : values())
				if(vi.index == idx)
					return vi;
			return null;
		}
	}

	public static final class VSOP2013Data{
		int timePower;
		VSOP2013Coeffs[] coeffs;
	}

	public static final class VSOP2013Coeffs{
		final short[] iphi = new short[17];
		double sine;
		double cosine;
	}

	public static final class Orbit{
		//epoch of orbital elements as Julian Ephemeris Date [JD]
		double t;
		//[AU]
		double periapseDistance;
		//0 = circular, 1 = parabolic, > 1 = hyperbolic
		double eccentricity;
		//[rad]
		double inclinationReferencePlane;
		//[rad]
		double argumentPeriapse;
		//[rad]
		double longitudeAscendingNode;
		//[rad]
		double meanAnomalyAtEpoch;
		//[rad/day]
		double meanMotion;
	}


	public static Map<VariableIndex, List<VSOP2013Data>> readData(final String filename) throws IOException{
		final ClassLoader classLoader = ResourceReader.class.getClassLoader();
		final InputStream is = classLoader.getResourceAsStream(filename);
		if(is == null)
			throw new IllegalArgumentException("file not found! " + filename);

		final Map<VariableIndex, List<VSOP2013Data>> datas = new HashMap<>();
		try(
				final InputStreamReader sr = new InputStreamReader(is, StandardCharsets.UTF_8);
				final BufferedReader reader = new BufferedReader(sr)){
			String line;
			while((line = reader.readLine()) != null){
				final String[] header = StringUtils.split(line, ' ');

				//3 = Earth
				//final int planet = Integer.parseInt(header[1]);

				final VSOP2013Data data = new VSOP2013Data();
				final VariableIndex variable = VariableIndex.of(Integer.parseInt(header[2]));
				data.timePower = Integer.parseInt(header[3]);
				final int count = Integer.parseInt(header[4]);
				data.coeffs = new VSOP2013Coeffs[count];
				for(int i = 0; i < count; i ++){
					line = reader.readLine();

					final String[] terms = StringUtils.split(line, ' ');

					//read coefficients
					final VSOP2013Coeffs coeffs = new VSOP2013Coeffs();
					coeffs.iphi[0] = Short.parseShort(terms[1]);
					coeffs.iphi[1] = Short.parseShort(terms[2]);
					coeffs.iphi[2] = Short.parseShort(terms[3]);
					coeffs.iphi[3] = Short.parseShort(terms[4]);

					coeffs.iphi[4] = Short.parseShort(terms[5]);
					coeffs.iphi[5] = Short.parseShort(terms[6]);
					coeffs.iphi[6] = Short.parseShort(terms[7]);
					coeffs.iphi[7] = Short.parseShort(terms[8]);
					coeffs.iphi[8] = Short.parseShort(terms[9]);

					coeffs.iphi[9] = Short.parseShort(terms[10]);
					coeffs.iphi[10] = Short.parseShort(terms[11]);
					coeffs.iphi[11] = Short.parseShort(terms[12]);
					coeffs.iphi[12] = Short.parseShort(terms[13]);

					coeffs.iphi[13] = Short.parseShort(terms[14]);

					coeffs.iphi[14] = Short.parseShort(terms[15]);
					coeffs.iphi[15] = Short.parseShort(terms[16]);
					coeffs.iphi[16] = Short.parseShort(terms[17]);

					coeffs.sine = Double.parseDouble(terms[18] + "e" + terms[19]);
					coeffs.cosine = Double.parseDouble(terms[20] + "e" + terms[21]);

					data.coeffs[i] = coeffs;
				}

				datas.putIfAbsent(variable, new ArrayList<>(1))
					.add(data);
			}

			return datas;
		}
	}

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
