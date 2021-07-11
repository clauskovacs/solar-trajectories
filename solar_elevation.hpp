#ifndef __SOLAR_ELEVATION_H_INCLUDED__
#define __SOLAR_ELEVATION_H_INCLUDED__

double julian_date(float hour, float minute, float seconds, float day, float month, float year);

struct sun_position
{
	float elevation_angle;
	float azimuth_angle;
};

sun_position calculate_sun_position(float latitude, float longitude, float hours, float minutes, float seconds, float day, float month, float year, float timezone);

#endif
