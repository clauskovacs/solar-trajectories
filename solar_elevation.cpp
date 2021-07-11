// Calculate the solar trajectories according to http://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html
#include <math.h>

#include "solar_elevation.hpp"

double julian_date(float hour, float minute, float seconds, float day, float month, float year)
{
	// for TT.MM.YYYY >= 15.10.1582 !
	double JD, y, m, d, h, b, a;

	if (month > 2)
	{
		y = year;
		m = month;
	}
	else
	{
		y = year - 1;
		m = month + 12;
	}

	d = day;
	h = hour / 24 + minute / 1440 + seconds / 86400;

	a = floor(y/100);

	b = 2 - a + floor(a/4);

	JD = floor(365.25 * (y + 4716)) + floor(30.6001 * (m + 1)) + d + h + b - 1524.5;

	return JD;
}


sun_position calculate_sun_position(float latitude, float longitude, float hours, float minutes, float seconds, float day, float month, float year, float timezone)
{
	sun_position computed_position;

	float time_past_local_midnight = (minutes / 3600) / 24 + (minutes / 60) / 24 + hours / 24;

	float JD = julian_date(hours, minutes, seconds, day, month, year);

	float julian_century = (JD - 2451545)/36525;

	float geom_mean_long_sun_deg = 280.46646 + julian_century * (36000.76983 + julian_century * 0.0003032);

	while (geom_mean_long_sun_deg > 360)
		geom_mean_long_sun_deg -= 360;

	float geom_mean_anom_sun_deg = 357.52911 + julian_century * (35999.05029 - 0.0001537 * julian_century);

	float eccent_earth_orbit = 0.016708634 - julian_century * (0.000042037 + 0.0000001267 * julian_century);

	float fr = M_PI/180;

	float sun_eq_of_ctr = sin(fr * (geom_mean_anom_sun_deg)) * (1.914602 - julian_century * (0.004817 + 0.000014 * julian_century)) + sin(fr * (2 * geom_mean_anom_sun_deg)) * (0.019993 - 0.000101 * julian_century) + sin(fr * (3 * geom_mean_anom_sun_deg)) * 0.000289;

	float sun_true_long_deg = geom_mean_long_sun_deg + sun_eq_of_ctr;

//	float sun_true_anom_deg = geom_mean_anom_sun_deg + sun_eq_of_ctr;

//	float sun_rad_vec_AUs =(1.000001018 * (1 - eccent_earth_orbit * eccent_earth_orbit))/(1 +  eccent_earth_orbit * cos(fr * (sun_true_anom_deg)));

	float sun_app_long_deg = sun_true_long_deg - 0.00569 - 0.00478 * sin(fr * (125.04 - 1934.136 * julian_century));

	float mean_obliq_ecliptic_deg = 23 + (26 + ((21.448 - julian_century * (46.815 + julian_century * (0.00059 - julian_century * 0.001813))))/60)/60;

	float obliq_corr_deg = mean_obliq_ecliptic_deg + 0.00256 * cos(fr * (125.04 - 1934.136 * julian_century));

//	float sun_rt_ascen_deg = atan2(cos(fr * (obliq_corr_deg))*sin(fr * (sun_app_long_deg)), cos(fr * (sun_app_long_deg))) / fr;

	float sun_declin_deg = asin(sin(fr * (obliq_corr_deg)) * sin(fr * (sun_app_long_deg))) / fr;

	float var_y = tan(fr * (obliq_corr_deg/2)) * tan(fr * (obliq_corr_deg/2));

	float eq_of_time_min = 4 * (var_y * sin(2 * fr * (geom_mean_long_sun_deg)) - 2 * eccent_earth_orbit * sin(fr * (geom_mean_anom_sun_deg)) + 4 * eccent_earth_orbit * var_y * sin(fr * (geom_mean_anom_sun_deg)) * cos(2 * fr * (geom_mean_long_sun_deg)) - 0.5 * var_y * var_y * sin(4 * fr * (geom_mean_long_sun_deg)) - 1.25 * eccent_earth_orbit * eccent_earth_orbit * sin(2 * fr * (geom_mean_anom_sun_deg))) / fr;

//	float HA_sunrise_deg = (acos(cos(fr * (90.833)) / (cos(fr * (latitude))*cos(fr * (sun_declin_deg))) - tan(fr * (latitude)) * tan(fr * (sun_declin_deg)))) / fr;

	float true_solar_time_min = (time_past_local_midnight * 1440) + eq_of_time_min + (4 * longitude) - (60 * timezone);

	while (true_solar_time_min > 1440)
		true_solar_time_min -= 1440;

	float factor = 0;

	if ((true_solar_time_min / 4 ) < 0)
		 factor = 180;
	else
		 factor = -180;

	float hour_angle_deg = (true_solar_time_min/4) + factor;

	float solar_zenith_angle_deg = (acos(sin(fr * (latitude)) * sin(fr * (sun_declin_deg)) + cos(fr * (latitude)) * cos(fr * (sun_declin_deg)) * cos(fr * (hour_angle_deg)))) / fr;

	float solar_elevation_angle_deg = 90 - solar_zenith_angle_deg;

	float approx_atmospheric_ref_deg = 0;

	if(solar_elevation_angle_deg > 85)
	{
		approx_atmospheric_ref_deg = 0;
	}
	else
	{

		if (solar_elevation_angle_deg > 5)
		{
			approx_atmospheric_ref_deg = (58.1 / tan(fr * (solar_elevation_angle_deg)) - 0.07 / pow(tan(fr * (solar_elevation_angle_deg)),3) + 0.000086 / pow(tan(fr * (solar_elevation_angle_deg)),5)) / 3600;
		}
		else
		{

			if (solar_elevation_angle_deg > - 0.757)
			{
				approx_atmospheric_ref_deg = (1735 + solar_elevation_angle_deg * (-518.2 + solar_elevation_angle_deg * (103.4 + solar_elevation_angle_deg * (-12.79 + solar_elevation_angle_deg * 0.711)))) / 3600;
			}
			else
			{
				approx_atmospheric_ref_deg = (-20.772 / tan(fr * (solar_elevation_angle_deg))) / 3600;
			}
		}
	}

	float solar_elevation_corrected_for_atm_refraction_deg = solar_elevation_angle_deg + approx_atmospheric_ref_deg;

	float solar_azimuth_angle_deg_cw_from_n = 0;

	if (hour_angle_deg > 0)
	{
		solar_azimuth_angle_deg_cw_from_n = (acos(((sin(fr * (latitude)) * cos(fr * (solar_zenith_angle_deg))) - sin(fr * (sun_declin_deg))) / (cos(fr * (latitude)) * sin(fr * (solar_zenith_angle_deg))))) / fr + 180;
	}
	else
	{
		solar_azimuth_angle_deg_cw_from_n = 540 - (acos(((sin(fr * (latitude)) * cos(fr * (solar_zenith_angle_deg))) - sin(fr * (sun_declin_deg))) / (cos(fr * (latitude)) * sin(fr * (solar_zenith_angle_deg))))) / fr;
	}

	while(solar_azimuth_angle_deg_cw_from_n > 360)
		solar_azimuth_angle_deg_cw_from_n -= 360;

	computed_position.elevation_angle = solar_elevation_corrected_for_atm_refraction_deg;
	computed_position.azimuth_angle = solar_azimuth_angle_deg_cw_from_n;

	return computed_position;
}
