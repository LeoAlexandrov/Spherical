# Spherical triangular mesh 

The spherical triangular mesh is a recursive decomposition of a spherical surface to spherical triangles (cells) uniformly covering it everywhere including poles, unlike the coordinate grid.

At zero level, the sphere is decomposed to 20 cells (the faces of an inscribed regular icosahedron projected to the sphere surface), 5 bits required to index them. For the Earth, each cell edge has the length equal to 7053644 meters.

On every next step (level), each cell is split into 4 smaller tiles by connecting the midpoints of its edges, 2 bits required to index them. 

![Octahedron inside sphere](https://media.aleprojects.com/github/Icosahedron.png)


A 64-bit integer can address cells of up to 28 level. Such a key can be assigned to sphere objects and will allow to retrieve them from a database or hash-table very quickly. 

The `Mesh` static helper class provides a coordinate system, topmost triangular cells at zero level, calculates keys at given level for a location.

The `MeshCell` class represents a triangular cell and allows to cover regions like circles, polygons, and polylines by cells.

```
public struct MeshCell : IComparable<MeshCell>, IEquatable<MeshCell>
{
	public long Key { get; private set; }
	public int Level { get; private set; }
	public CartesianValue Vertex1 { get; private set; }
	public CartesianValue Vertex2 { get; private set; }
	public CartesianValue Vertex3 { get; private set; }

	// other computed properties and methods
}
```

### Example 1

```
long key = Mesh.BuildKey(48.4432, 1.9574, 20);
```

This code calculates a key for the location with latitude = 48.4432 and longitude = 1.9574 at the 20th subdivision level. The result is `115945573754077184`. In binary presentation `0 00000 01 10 01 10 11 11 10 10 11 11 10 00 10 10 01 00 00 11 10 01 00 0000 0000 0000 0000`.  
The highest bit is the sign bit, always 0. The next 5 zero bits are the index of the first topmost triangle (icosahedron face). Next two bits 01 are the index of the triangle on the next grid level. Next pairs of bits are indexes of triangles at the next levels.

### Example 2

Let's cover a circle of some radius with the center in some location by cells at some level of subdivision.

```csharp
using AleProjects.Spherical;
using AleProjects.Spherical.Mesh;

double latitude = 47.987783;
double longitude = 0.240283;
double radius = 1000.0;
double angle = radius / SphericalExtension.EARTH_MEAN_RADIUS; // radians
int level = 13;                                               // level of subdivision, ~ 861 meters edge length
bool joinCells = true;                                        // join cells to reduce their number
var location = new Cartesian(latitude, longitude);

List<MeshCell> cells = MeshCell.CoverCircle(location, angle, level, joinCells);
```

This is how the cells look at the 13th level of subdivision covering the circle with the center in (47.987783, 0.240283) and radius of 1000 m. 

![Cells covering a circle](https://media.aleprojects.com/github/CoverCircle-1000m-13lev.jpg)

Suppose we have a database:

```csharp
using Microsoft.EntityFrameworkCore;
using AleProjects.Spherical;

public class Pokemon : IGeoCoordinate
{
	public int Id { get; set; }
	public double Longitude { get; set; }
	public double Latitude { get; set; }
	public long SpatialKey { get; set; }
	public bool IsActive { get; set; }
}

...

public class PokemonsDbContext : DbContext
{
	public DbSet<Pokemon> Pokemons { get; set; }
}
```

and we need to find all active pokemons around a user. So

```csharp
using AleProjects.Spherical;

...
// dbContext here is an instance of PokemonsDbContext,
// cells, location, and angle are defined as in the code above

var ranges = cells
	.Select(c => new { Lower = c.Key, Upper = c.KeyUpperValue })
	.ToArray();

var query = dbContext.Pokemons
	.AsNoTracking()
	.Where(p => p.SpatialKey >= ranges[0].Lower && p.SpatialKey <= ranges[0].Upper)
	.Where(p => p.IsActive);

for (int i = 1; i < ranges.Length; i++)
{
	var range = ranges[i];

	query = query.Union(
		dbContext.Pokemons
			.AsNoTracking()
			.Where(p => p.SpatialKey >= range.Lower && p.SpatialKey <= range.Upper)
			.Where(p => p.IsActive)
	);
}

var pokemons = await query.ToArrayAsync();

// optionally, it is possible to filter out more precisely by the distance  
// pokemons.RemoveAll(p => p.DistanceTo(location.Latitude, location.Longitude) > angle * SphericalExtension.EARTH_MEAN_RADIUS);
```

What is the `KeyUpperValue` and what is it for? It is a value of the key at a given subdivision level padded with 1 bits to the right. And the condition `p.SpatialKey >= range.Lower && p.SpatialKey <= range.Upper` includes all pokemons within the given cell even if they have the key for higher subdivision level.

Of course, table column `SpatialKey` must be indexed.

**Important**: `pokemon.SpatialKey` must be built for the subdivision level equal or higher than the level of covering cells used in the query.

### Example 3

Let's make server side clustering of pokemons. We will cluster pokemons at the current map zoom level and then calculate the number of pokemons in each cluster. Just general approach, not ready-to-use code.


```csharp
using AleProjects.Spherical;
using AleProjects.Spherical.Mesh;


// define cluster struct
public struct MapCluster
{
	public string Id { get; set; }
	public long SpatialKey { get; set; }
	public double Longitude { get; set; }
	public double Latitude { get; set; }
	public int Count { get; set; }
}


var SW = new Cartesian( .... ); // South-West corner of the map
var NE = new Cartesian( .... ); // North-East corner of the map

var polygon = [
	SW,
	new Cartesian(NE.Latitude, SW.Longitude),
	NE,
	new Cartesian(SW.Latitude, NE.Longitude)
];

// here is Gotcha.
// If the map is zoomed out so much (zoom < 6), map view area significantly differs
// from the above polygon on the sphere surface due to mercator projection.
// workaround is to add some more points to the polygon, for example mid points of top and bottom sides of the map view.

int level = MapZoom - 3; // adjust subdivision level to the map zoom, 3 is just a heuristic value.
bool joinCells = false;  // no joining, otherwise cluster areas will be not uniform.

var cells = Mesh.CoverPolygon(polygon, level, joinCells);

long mask = Mesh.KeyMask(level - 1);

var ranges = cells
	.Select(c => new { Lower = c.Key, Upper = c.KeyUpperValue })
	.ToArray();

var query = dbContext.Pokemons
	.AsNoTracking()
	.Where(p => p.SpatialKey >= ranges[0].Lower && p.SpatialKey <= ranges[0].Upper)
	.Where(p => p.IsActive)
	.GroupBy(t => t.SpatialKey & mask)
	.Select(g => new MapCluster()
		{
			Id = g.Max(t => t.Id),
			SpatialKey = g.Key,
			Count = g.Count(),
			Longitude = g.Average(t => t.Longitude),
			Latitude = g.Average(t => t.Latitude)
		});

for (int i = 1; i < ranges.Length; i++)
{
	var range = ranges[i];

	query = query.Union(
		dbContext.Pokemons
			.AsNoTracking()
			.Where(p => p.SpatialKey >= range.Lower && p.SpatialKey <= range.Upper)
			.Where(p => p.IsActive)
			.GroupBy(t => t.SpatialKey & mask)
			.Select(g => new MapCluster()
				{
					Id = g.Max(t => t.Id),
					SpatialKey = g.Key,
					Count = g.Count(),
					Longitude = g.Average(t => t.Longitude),
					Latitude = g.Average(t => t.Latitude)
				}));
}

var clusters = await query.ToArrayAsync();
```

Now we have clusters with the number of pokemons and average location for the map marker placement. We can serialize them to json, send back to the client, and put markers on the map.

If the `Count` property is 1, we can assume that the cluster contains only one pokemon and use its coordinates as the real location, use different marker icon, and attach **onclick** event opening something like `/api/pokemons/${id}`.

How the clusters look at zoom level 9 and subdivision level 6.

![Clusters at zoom level 9](https://media.aleprojects.com/github/Clusters-9zoom.jpg)

**Important**: `pokemon.SpatialKey` must be built for the subdivision level equal or higher than the level of covering cells used in the query.

### Example 4

About indexing continuous objects like polylines with the mesh. Below is the picture demonstrating selected OpenStreetMap ways along some route.

![Cover polylines](https://media.aleprojects.com/github/CoverPolyline.jpg)

How it works:

- First, we extract nodes (a point with Id, Latitude, and Longitude) and ways (a collection of node references + Id, name, and other properties) from the OpenStreetMap `pbf` files.
- Then, the way is a polyline, we cover it by cells at some subdivision level, I took level 14 (`Mesh.CoverPolyline(wayNodes, 14, 0.0, false)`). Save cell keys in the database table `WaysCells` with columns `Id`, `WayId`, and `Key`.
- Then, save way properties (including Id as primary key) and its nodes as Google-encoded polyline (for simplicity) in the database table `Ways`.

Tables `Ways` and `WaysCells` are in one-to-many relationship by the `Ways.Id` and `WaysCells.WayId` columns. Each way has at least one covering cell. At this step we can select all ways at given location:

```csharp
using AleProjects.Spherical;
using AleProjects.Spherical.Mesh;
using Dapper;

int level = 14;
long key = Mesh.BuildKey(location, level);
string sql = $"SELECT w.* FROM Ways as w INNER JOIN WaysCells as wc ON w.Id=wc.WayId WHERE wc.Key={key};"

var ways = await conn.QueryAsync<Way>(sql);
```

If we have some route and we'd like to get all ways along it as shown in the picture:

```csharp
using AleProjects.Spherical;
using AleProjects.Spherical.Mesh;
using Dapper;


List<Cartesian> route = [ ... ];
int level = 14;
double tolerance = 100.0;
bool joinCells = false;

var cells = Mesh.CoverPolyline(route, level, tolerance, joinCells);
var keys = cells.Select(c => c.Key.ToString());

string sql = $"SELECT DISTINCT w.* FROM Ways as w INNER JOIN WaysCells as wc ON w.Id=wc.WayId WHERE wc.Key IN ({string.Join(',', keys)};";

var ways = await conn.QueryAsync<Way>(sql);
```

What is `double tolerance = 100.0;`?  
Think of it as polyline 'thickness'. More exactly, half-thickness because it works for both polyline sides.  
What happens if to use tolerance = 0.0 :

![Cover polyline tolerance](https://media.aleprojects.com/github/CoverPolyline_Tolerance_Example.jpg)

Cells marked red will not be included.