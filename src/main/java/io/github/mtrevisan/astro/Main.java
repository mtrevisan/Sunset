package io.github.mtrevisan.astro;

import io.github.mtrevisan.astro.core.Season;
import io.github.mtrevisan.astro.core.Zenith;
import io.github.mtrevisan.astro.coordinates.AtmosphericModel;
import io.github.mtrevisan.astro.coordinates.GeographicLocation;
import io.github.mtrevisan.astro.core.SunlightPhase;
import io.github.mtrevisan.astro.core.EarthCalculator;
import io.github.mtrevisan.astro.helpers.TimeHelper;

import java.time.ZoneId;
import java.time.ZonedDateTime;
import java.time.format.DateTimeFormatter;


public class Main{

	//https://www.agopax.it/Libri_astronomia/pdf/Astronomical%20Algorithms.pdf
	//https://gml.noaa.gov/grad/solcalc/calcdetails.html
	//https://www.nrel.gov/docs/fy08osti/34302.pdf
	//https://www.sunearthtools.com/dp/tools/pos_sun.php
	public static void main(final String[] args){
		final int year = 2023;
		final ZoneId zoneId = ZoneId.of("Europe/Rome");
		//Padova 3 °C, Treviso 4 °C, Venezia 1015.5 hPa 4 °C, Verona 3 °C, Vicenza 3 °C
		AtmosphericModel atmosphere = AtmosphericModel.create(1021.1, 4.);
		//16:53
//		final GeographicLocation location = GeographicLocation.create(45.4421304114541, 10.622970782783327, 0.)
		//16:52
//		final GeographicLocation location = GeographicLocation.create(45.63956478260069, 13.100940416786331, 0.)
		//16:55
//		final GeographicLocation location = GeographicLocation.create(44.79124286512801, 12.398222626853773, 0.)
		//16:49
//		final GeographicLocation location = GeographicLocation.create(46.68058098966078, 12.477496139552612, 0.)
		//18:42
//		final GeographicLocation location = GeographicLocation.create(45.827049477073956, 10.857018789615736, 16.)
		//16:53
		final GeographicLocation location = GeographicLocation.create(45.714920, 12.194179, 16.)
			.withAtmosphere(atmosphere);

//		final DecimalFormat decimalFormatter = (DecimalFormat)NumberFormat.getNumberInstance(Locale.US);


		ZonedDateTime winterSolstice = EarthCalculator.season(year, zoneId, Season.WINTER);
System.out.println("ws " + DateTimeFormatter.ISO_LOCAL_DATE_TIME.format(TimeHelper.terrestrialTimeToUniversalTime(winterSolstice)) + " UTC");

//ZonedDateTime date = winterSolstice
//	.atZone(ZoneId.of("Etc/UTC"));
//SunriseResult result = SPA.calculateSunriseTransitSet(date, 45.714920, 12.194179, deltaT, SPA.Horizon.CIVIL_TWILIGHT);
//if(result instanceof SunriseResult.RegularDay regular)
//	System.out.println("spa " + regular.sunset());
		EarthCalculator calculator = EarthCalculator.create(location);
		SunlightPhase sunlightPhase = calculator.sunlightPhase(winterSolstice.plusMonths(3), Zenith.CIVIL);
		if(sunlightPhase instanceof SunlightPhase.RegularDay event){
			final ZonedDateTime utcSunset = TimeHelper.terrestrialTimeToUniversalTime(event.sunset());
			System.out.println(DateTimeFormatter.ISO_LOCAL_TIME.format(utcSunset) + " UTC");
System.out.println(DateTimeFormatter.ISO_LOCAL_TIME.format(TimeHelper.universalTimeToApparentSolarTime(utcSunset, location)) + " LMT");
		}
	}

}
