# solar-trajectories
Calculate the position (azimuth and elevation) of the sun depending on the position on earth and the date.

## Description
This program uses the astronomical algorithms by Jean Meeus to calculate the position of the sun [[1]](#1). This is done via the function `calculate_sun_position` which returns the **elevation** and **azimuth** of the sun depending on the **latitude**, **longitude**, **hours**, **minutes**, **seconds**, **day**, **month**, **year** and **timezone**. The nine latter mentioned variables represent the input variables of the function, i.e., they define the position of the observer in space and time on earth for which the sun position should be determined.

## Additional Information
Additional information can be found at <http://phys.ik.cx/programming/cpp/solar_elevation/index.php>.

The algorithm in action is visualised through a video at <https://www.youtube.com/watch?v=9oRLf2-qdus>.

## References
<a id="1">[1]</a> 
<https://gml.noaa.gov/grad/solcalc/calcdetails.html>
