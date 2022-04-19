First, install the System.ValueTuple package from NuGet.

	Install-Package "System.ValueTuple"

All spherical calculations are based on vectors, their dot, cross, triple products, and positions relative to planes.

A vector must implement **ICartesian** interface:

	public interface ICartesian
	{
		double X { get; }
		double Y { get; }
		double Z { get; }
		
		void SetCartesian(double x, double y, double z);
	}

Some extension methods for the **ICartesian** interface are generic and can return a result of any type implementing **ICartesian**. Interface method **SetCartesian** is used there for initializing the result with X,Y,Z values. Classes implementing **ICartesian** may have some other properties depending on X,Y,Z. **SetCartesian** is a place where to maintain consistency of these properties with X,Y,Z. aleprojects.spherical namespace provides **Cartesian** class representing a basic vector. 

A location on a sphere must implement the **IGeoCoordinate** interface:

	public interface IGeoCoordinate
	{
		double Lat { get; } // Latitude in radians.
		double Lon { get; } // Longitude in radians.
		
		void SetSpherical(double lat, double lon);
	}

Purpose of **SetSpherical** is the same as of **SetCartesian**.
aleprojects.spherical namespace provides **GeoCoordinate** class that implements both **ICartesian** and **IGeoCoordinate** interfaces.

# Calculations

The SphericalExtension static class provides extension methods for the ICartesian and IGeoCoordinate interfaces. Examples:

	Cartesian vectorI = new Cartesian(1.0, 0.0, 0.0, false);
	Cartesian vectorJ = new Cartesian(0.0, 1.0, 0.0, false);
	Cartesian N = vectorI.CrossProduct<Cartesian>(vectorJ);
	GeoCoordinate L = vectorJ.CrossProduct<GeoCoordinate>(vectorI);
	double angle = vectorI.Angle(vectorJ);
	
	GeoCoordinate location1 = new GeoCoordinate(33.742763, -84.393547);
	GeoCoordinate location2 = new GeoCoordinate(32.392110, -86.304144);
	double distance = location1.DistanceTo(location2);

I'll not focus on such methods like DotProduct, CrossProduct, Angle, DistanceTo, etc, because their use and purpose are obvious. But some methods for typical geospatial tasks are worth mentioning and explaining.

### InsideTriangle extension method
	public static bool InsideTriangle(this ICartesian cartesian, ICartesian vertex1, ICartesian vertex2, ICartesian vertex3)
	public static bool InsideTriangle(this ICartesian cartesian, double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3)

This method checks if the vector **cartesian** is inside the spherical triangle with the **vertex1**, **vertex2**, **vertex3** vertices. Or with the vertices given by (x1,y1,z1), (x2,y2,z2), (x3,y3,z3). The method calculates triple products of the **cartesian** with each pair of vectors. The first criterion is that all results are of one sign. If the vector is inside, an opposite vector will produce triple products of one sign too. So, the second criterion is a plane formed by the triangle vertices ends and moved to the coordinate origin. End of the  **cartesian** vector must be on the one side with any of the vertices relative to this plane. This method is a cornerstone of the sphere triangular grid.

### PolygonInfo method
    public static double[] PolygonInfo(IEnumerable<ICartesian> polygon)

Returns an array of angles between planes formed by vertices of a polygon. Any negative angle in the returned array indicates that the polygon is not a convex.
Example:

    string line = "svx}FhefvOa_HfsDykA~zKleIv_C`_Bo`Ggr@}wFhrCglA|kBlyNalJhvIkeJquEwkAqcLbrA{bHbrDy{D";
    var polygon = SphericalExtension.DecodeGooglePolyline<GeoCoordinate>(line);
    var info = SphericalExtension.PolygonInfo(polygon);

![Polygon vertices and angles](https://aleprojects.com/upload/images/polygon-info.png)

For this polygon **PolygonInfo** will return the following array:

![Test section](https://aleprojects.com/upload/images/polygon-info-values.png)

Vertices from 1 to 4 make this polygon not a convex. *This does not mean that if to remove these vertices, the polygon will become a convex*. If some value of this array is zero, this means that 3 points are on one geodesical line and the middle one can be safely removed. 

### InflateConvex method

This method inflates or deflates a convex polygon. 

    public static List<T> InflateConvex<T>(IEnumerable<ICartesian> polygon, IList<double> angles)

The idea behind calculations is to rotate planes formed by vertices vectors to an angle which is a distance on a sphere divided by sphere radius. **InflateConvex** rotates normal vectors N1 = V1 x V2, N2 = V2 x V3, etc. that represent planes. Next step is to find new vertices of the inflated convex. They are cross products of the rotated normal vectors N1 x N2, etc.

![Inflate convex polygon](https://aleprojects.com/upload/images/inflate-convex.jpg)

**angles** parameter must contain at least 1 element. **InflateConvex** takes the next value in **angles** for the next plane and uses the last value for the rest planes if the length of **angles** is less than the number of the convex vertices. 


### SectionsIntersect method

	public static int SectionsIntersect(ICartesian vertex1S1, ICartesian vertex2S1, ICartesian vertex1S2, ICartesian vertex2S2)
	public static int SectionsIntersect(ICartesian vertex1S1, ICartesian vertex2S1, ICartesian vertex1S2, ICartesian vertex2S2, double tolerance)

This method checks that two sections (geodesical lines) on a sphere intersect. First, this method takes a plane formed by **vertex1S1** and **vertex2S1** and checks positions of **vertex1S2** and **vertex2S2** relative to the plane. Then the same with a plane by **vertex1S2** and **vertex2S2** and the positions of **vertex1S1** and **vertex2S1**. The method returns a negative number if the sections don't intersect, positive if intersect, zero if some section end belongs to another section. Overloaded method with **tolerance** parameter tests if the sections are close enough to each other to consider them intersecting. Tolerance is an angle in radians at the coordinates origin, real distance on the surface is tolerance\*sphere_radius.

### InsidePolygon extension method

	public static bool InsidePolygon(this ICartesian cartesian, IEnumerable<ICartesian> polygon)

This method checks if a vector **cartesian** is inside a polygon including its borders using winding number algorithm. The polygon can be of any shape. If the polygon is a convex, it is better to use **InsideConvex** method. 

### TestSection extension method

	public static double TestSection(this ICartesian cartesian, ICartesian vertex1, ICartesian vertex2)

This method checks if a vector **cartesian** is between two vectors **vertex1** and **vertex2**. First, it builds a plane using **vertex1** and **vertex2**. Then it takes an angle between the vector **cartesian** and this plane and corrects it with PI/2 angle. The method returns this angle, its positive value indicates that **cartesian** is between the vectors, negative if not. The next step is to build another plane containing **cartesian** and perpendicular to the plane obtained in the first step. Finally, it checks positions of **vertex1** and **vertex2** relative to this plane. If they are on the opposite sides, the result is positive.

![Test section](https://aleprojects.com/upload/images/test_section.jpg)

This drawing demonstrates **TestSection** method. If to call C.TestSection(V1,V2) with C,V1,V2 as they drawn, it will return positive alpha angle indicating success. Alpha\*sphere_radius is a distance from the point C to the arc on the sphere. If C is out of the arc, method will return negative value of the angle.

### TestPolyline extension method

	public static int TestPolyline(this ICartesian cartesian, IEnumerable<ICartesian> polyline, double tolerance, bool reverse, Func<ICartesian, ICartesian, ICartesian, double, double, bool> userTest, out double angleFromPolyline)
	public static PolylineTestResult TestPolyline(this ICartesian cartesian, IEnumerable<ICartesian> polyline, double tolerance, bool reverse, Func<ICartesian, ICartesian, ICartesian, double, double, bool> userTest, PolylineTestResult result)

This method works in a similar way as **TestSection**. It applies **TestSection** to all sections of the polyline with some additional checks. First additional check is how far from the polyline the point can stand. The **tolerance** parameter is an angle in radians at the coordinates origin, it represents maximum distance as tolerance\*sphere_radius. Second optional check is the **userTest** delegate, if it is not null.  If the task is only to check the position of a point relative to a polyline, no need to use delegate. If the task is to guide moving point along directional route polyline, then the delegate is required.

Example:

![Test polyline](https://aleprojects.com/upload/images/test_polyline_1.png) 

	List<GeoCoordinate> polyline = SphericalExtension.DecodeGooglePolyline<GeoCoordinate>("m{k~F`c~uOOuCCqCEiCSkD_@gC");
	GeoCoordinate location = new GeoCoordinate(41.84524, -87.64902);
	double heading = 85.0 * Math.PI / 180.0; // 85 degrees, almost to the East
	double speed = 40; // units don't matter in this example, let it be mph.
	double tolerance = 30; // 30 meters
	
	Func<ICartesian, ICartesian, ICartesian, double, double, bool> geometryTest = (v0, v1, v2, angle, azimuth) =>
	{
		// the higher speed, the less fluctuations of a heading calculated from GPS data. 
		double diff = speed < 10.0 ? Math.PI / 4.0 : Math.PI / 10.0;
		double azimuth2;

		return Math.Abs(heading.HeadingsDiff(azimuth)) < diff ||
			(v0 != null &&
			 (Math.Abs(heading.HeadingsDiff(azimuth2 = v0.AzimuthOnSphere(v1))) < diff ||
			 heading.HeadingInRange(azimuth, azimuth2)));
	};
	
	PolylineTestResult testResult = location.TestPolyline(
		polyline, 
		tolerance / SphericalExtension.EARTH_MEAN_RADIUS, 
		false, 
		geometryTest, 
		null);


In most cases the **geometryTest** delegate is called with v0 = null and v1, v2 corresponding to the section of the polyline most closest to the point. This is the case of point C on the drawing below. In this case it is enough to compare current heading with the azimuth from v1 to v2. The **azumuth** parameter of the delegate holds this azimuth. **HeadingsDiff** is the extension method of the SphericalExtesion static class. The **angle** parameter if multiplied by SphericalExtension.EARTH_MEAN_RADIUS represents a distance to the section.

![Test polyline](https://aleprojects.com/upload/images/test_polyline2.jpg)

Sometimes it is possible that the point is close enough to some vertex of polyline, but failed the **TestSection** for both sections. Point C' on the drawing demonstrates this situation. When such a thing happens, v0 parameter is not null. In this case we need to make additional test for the heading. First, we compare the heading with the azimuth of v0,v1 (points marked 0, 1 on the drawing). If this test fails, we check if the heading is within range of the azimuths of neighbour sections. **HeadingInRange** is the standard extension method. 

Finally, if the polyline test succeeds, it return non-null object of **PolylineTestResult** type. **SectionIndex** property contains index of the section found, the same as the index of the first vertex of the section. There is useful method **DistanceToEnd** returning a distance along the polyline from the point to the end.

If we apply **TestPolyline** to the next received point with the same polyline and with the result for the previous point, new test will begin from **SectionIndex** vertex of the polyline.

	testResult = nextLocation.TestPolyline(
		polyline, 
		tolerance / SphericalExtension.EARTH_MEAN_RADIUS, 
		false, 
		geometryTest, 
		testResult); // <----


If the new test succeeds, testResult will be updated with new values. If the test fails, **Fails** property will be incremented by 1. How to handle non-zero **Fails** is up to developer. This may happen if bad GPS data is received, but desired behavior is not to leave a route immediately.


# Sunrise and sunset

### Sunrise extension method

	public static DateTime Sunrise(this IGeoCoordinate location, DateTimeOffset date)

This method calculates a time of sunrise for a given **location** and **date**. The method uses only date.Year, date.Month, date.Day, and date.Offset components of the **date**. date.Offset is important for the calculation and considered as the time zone. The method returns DateTime.MaxValue if the Sun never rises at this location, DateTime.MinValue if never sets. In other cases, it returns local  sunrise time.

### Sunset extension method

	public static DateTime Sunset(this IGeoCoordinate location, DateTimeOffset date)

This method calculates a time of sunset for a given **location** and **date**. The method uses only date.Year, date.Month, date.Day, and date.Offset components of the **date**. date.Offset is important for the calculation and considered as a time zone. The method returns DateTime.MaxValue if the Sun never rises at this location, DateTime.MinValue if never sets. In other cases, it returns local sunset time.

### NearestSunriseSunset

	public static (DateTime sunrise, DateTime sunset) NearestSunriseSunset(IGeoCoordinate location, DateTimeOffset date)

Returns times of nearest past and future events. If specified time is after sunset for this date, the method returns time of sunrise for the next day. If specified time is before sunrise for this date, the method returns time of sunset for the previous day. Otherwise, it returns time of sunrise and sunset for the specified date. This makes easy to check if the current time is dark or not, and to switch your app to the dark theme. Example:

	var time = DateTimeOffset.Now;
	var (sunrise, sunset) = location.NearestSunriseSunset(time);
	if (time.DateTime > sunrise && time.DateTime < sunset) 
	{
		// Day time
	}
	else
	{
		// Dark time
	}

if the Sun never rises, the method returns (DateTime.MaxValue, DateTime.MinValue). If the Sun never sets, it returns (DateTime.MinValue, DateTime.MaxValue).

Extension method **DayTime** does the same as the example above. 

	public static bool DayTime(IGeoCoordinate location, DateTimeOffset date)


