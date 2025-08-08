package io.github.mtrevisan.astro;

import io.github.mtrevisan.astro.core.*;
import io.github.mtrevisan.astro.coordinates.AtmosphericModel;
import io.github.mtrevisan.astro.coordinates.GeographicLocation;
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
		final int year = 2025;
		final ZoneId zoneId = ZoneId.of("Europe/Rome");
		final AtmosphericModel atmosphere = AtmosphericModel.create(1021.1, 4.);
		//cortina: 16:56 (2000), 16:58 (1900-1800), 16:57 (1700), 16:56 (1600), 16:55 (1500), 16:53 (1400), 16:50 (1300), 16:47 (1200), 16:41 (1100), 16:33 (1000)
//		final GeographicLocation location = GeographicLocation.create(46.54056257185929, 12.135682113077916, 1224.)
		//portogruaro: 16:52 (2000), 16:53 (1900), 16:54 (1800-1700), 16:52 (1600), 16:51 (1500), 16:49 (1400), 16:46 (1300), 16:42 (1200), 16:37 (1100), 16:29 (1000)
//		final GeographicLocation location = GeographicLocation.create(45.776955464421896, 12.840035227634429, 7.)
		//16:53
		final GeographicLocation location = GeographicLocation.create(45.714920, 12.194179, 16.)
			.withAtmosphere(atmosphere);

//		final DecimalFormat decimalFormatter = (DecimalFormat)NumberFormat.getNumberInstance(Locale.US);

		final Season season = Season.WINTER;
		final ZonedDateTime winterSolstice = EarthCalculator.season(year, zoneId, season);
System.out.println(season + " " + DateTimeFormatter.ISO_LOCAL_DATE_TIME.format(TimeHelper.terrestrialTimeToUniversalTime(winterSolstice)) + " UTC");

		final EarthCalculator calculator = EarthCalculator.create(location);
		//"un'ora dopo il tramonto" corrisponderebbe a circa -15.83°
//		final ZenithInterface zenith = Zenith.NAUTICAL;
		final ZenithInterface zenith = () -> StrictMath.toRadians(-15.83);
		final SunlightPhase sunlightPhase = calculator.sunlightPhase(winterSolstice, zenith);
		if(sunlightPhase instanceof SunlightPhase.RegularDay event){
			final ZonedDateTime utcSunset = TimeHelper.terrestrialTimeToUniversalTime(event.sunset());
			System.out.println(DateTimeFormatter.ISO_LOCAL_TIME.format(utcSunset) + " UTC");
System.out.println(DateTimeFormatter.ISO_LOCAL_TIME.format(TimeHelper.universalTimeToApparentSolarTime(utcSunset, location.getLongitude())) + " LMT");
		}
	}

}
