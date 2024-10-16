# Sphere triangular grid 

The sphere triangular grid is a recursive decomposition of a sphere surface to spherical triangles (tiles) almost uniformly covering it everywhere including poles, unlike the coordinate grid. At zero level, the sphere is decomposed to 8 tiles (faces of a octahedron inscribed into the sphere), 3 bits required to index them. 

![Octahedron inside sphere](https://aleprojects.com/upload/images/octahedron-inside-sphere.png)

On every next step (level), each tile is splitted to 4 smaller tiles, 2 bits required to index them. 

![Sphere triangle tile](https://aleprojects.com/upload/images/sphere-triangle.png)

64-bit integer can address a sphere tile of up to 30 level. Such key can be assigned to sphere objects and will allow to retrieve them from a database or hash-table very quickly. 

Example

```
long quadKey = SphereGridHelper.BuildQuadKey(48.4432, 1.9574, 20);
```

This code calculates a key for the location with latitude = 48.4432 and longitude = 1.9574 at the 20th grid level. The result is `291816803624026112`. In binary presentation it is `0 000 01 00 00 00 11 00 10 11 11 01 11 01 01 10 00 11 01 11 10 00 00000000000000000000`.
The highest bit is the sign bit, it is always 0. Then 3 zero bits follow which are the index of the first topmost triangle. Next two bits 01 are the index of the triangle on the next grid level. It is 1 since latitude is greater than 45 degrees. Next pairs of bits are indexes of triangles at the next levels.

The **SphereGridHelper** static class provides coordinate system, topmost triangular tiles at zero level, calculates keys at given level for a location.

The **SphereGridTile** class represents triangular tile and allows to cover regions like  circles, polygons, and polylines with tiles.

Example

```
Cartesian location = new Cartesian(43.60943, 1.440064);
double angle = 2000.0 / SphericalExtension.EARTH_MEAN_RADIUS;
int level = 14;
bool joinTiles = true;

List<SphereGridTiles> tiles = SphereGridTile.CoverCircleByTiles(location, angle, level, joinTiles);
```

This is the visualization of tiles at the 14th level covering the circle with the center in (43.60943, 1.440064) and radius of 2000 m. 

![Tiles covering a circle](https://aleprojects.com/upload/images/tiles-circle.jpg)

Suppose we have some database table:

```
public class MyDataContext : DataContext
{
    public Table<Node> Nodes;
    
    public MyDataContext() : base("connection string goes here") { }
}

[Table(Name = "Nodes")]
public class Node : IGeoCoordinate
{
    [Column]
    public long QuadKey { get; set; }
    
    [Column]
    public double Latitude { get; set; }
    
    [Column]
    public double Longitude { get; set; }
   
}
```

This table is indexed by the QuadKey field and has millions of rows. With the list of tiles it is possible to query database objects for the circle quickly.

```
List<Node> result = new List<Node>();

using (MyDataContext context = new MyDataContext())
{
    long lower = 0;
    long upper = 0;

    // this query is translated to
    // SELECT * FROM Nodes WHERE QuadKey >= @lower AND QuadKey <= @upper
    var query = context.Nodes.Where(n => n.QuadKey >= lower && n.QuadKey <= upper);

    foreach (SphereGridTile tile in tiles)
    {
        lower = tile.QuadKey;
        upper = tile.QuadKeyUpperValue;

        result.AddRange(query);
    }
}

// optionally, it is possible to check precise distance

result.RemoveAll(n => n.DistanceTo(location.Latitude, location.Longitude) > angle * SphericalExtension.EARTH_MEAN_RADIUS);
```

Joining tiles by the **SphereGridTile.CoverCircleByTiles** helps to reduce number of database requests. 

How to find a balance between desired precision and number of covering tiles for a circle of given radius. From one hand, it is possible to select high grid level (small tiles), but for a large circle number of tiles will be very large. From other hand, it is possible to select low grid level (ultimately 0), but precision will be very low and the query will return too many rows. 

**SphereGridHelper** class has the function **LevelForCircleToTriangleRatio**. 

```
public static int LevelForCircleToTriangleRatio(double angle, double circleToTriangleRatio)
```

This function takes an angle representing a circle radius and **circleToTriangleRatio** which value is a circle to triangle surfaces ratio. First, this function takes 1/8 of the unit sphere surface and compares it with 2π(1-cos(angle))/circleToTriangleRatio value, where 2π(1-cos(angle)) is a surface of a circle on the unit sphere. While this value is larger than triangle surface (initially 1/8 of the sphere), it continues to divide triangle surface by 4 and compare. Finally, it returns a level when triangle surface is less than circle surface at least **circleToTriangleRatio** times.

Pictures below demonstrates coverage of circles with 4000m and 1000m radius and 4 and 64 **circleToTriangleRatio**.

```
Cartesian location = new Cartesin(48.854536, 2.29818);
double angle = 4000.0 / SphericalExtension.EARTH_MEAN_RADIUS;
int level = SphericalGridHelper.LevelForCircleToTriangleRatio(angle, 4.0); // or 64.0
bool joinTiles = true;

List<SphereGridTiles> tiles = SphereGridTile.CoverCircleByTiles(location, angle, level, joinTiles);
```

![Tiles covering a circle](https://aleprojects.com/upload/images/circle-to-triangle-4.jpg)

4000 m and 1000 m circles with 4 circle-to-triangle ratio.

![Tiles covering a circle](https://aleprojects.com/upload/images/circle-to-triangle-64.jpg)

4000 m and 1000 m circles with 64 circle-to-triangle ratio.

With **circleToTriangleRatio** = 64 coverage is more precise, but produces more tiles.  

Everything the same is possible for polygons and polylines.

```
public static List<SphereGridTile> CoverPolylineByTiles<T>(IEnumerable<T> polyline, int level, double tolerance, bool join)
	where T : ICartesian

public static List<SphereGridTile> CoverPolygonByTiles<T>(IEnumerable<T> polygon, int level, bool join)
	where T : ICartesian    
```

There is additional parameter **tolerance** for polylines. This is an angle, equal to a distance if multiplied by a sphere radius. Polyline is 1-dimensional, "very thin" object and it can be located very close to the triangle side. Typical task is to cover a polyline by tiles and to find objects along the route. Without tolerance, objects from one side will be missed in this case.