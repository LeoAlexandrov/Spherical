This repository contains two C# code files with a set of interfaces and classes for calculations on a sphere and for work with sphere triangular grid.

_This is a new version of the library, it is incompatible with previous one which stays in **v1** branch._

The first file AleSpherical.cs defines namespace aleprojects.spherical that provides interfaces ICartesian (vector), IGeoCoordinate (location on a sphere), classes implementing them, and static class SphericalExtension with extension methods. Extension methods allow to solve such geospatial tasks like to find a distance or azimuth from one point to another, to check if a point is inside a spherical polygon or belongs to a polyline, to inflate a spherical polygon, to calculate sunrise and sunset times (for the Earth), polyline encoding and decoding, etc. [Learn more](/doc/AleSpherical.md)

The second file AleSphericalGrid.cs defines namespace aleprojects.spherical.grid that provides static helper class SphereGridHelper for creating and handling of a sphere triangular grid, and  SphereGridTile class representing a grid tile. The SphereGridTile class has methods that allow to find tiles at the specified grid level covering continuous objects like polygons, polylines, and circles. [Learn more](/doc/AleSphericalGrid.md)

![Polygon covered by grid tiles](https://aleprojects.com/upload/images/tiles-polygon.jpg) ![Polyline covered by grid tiles](https://aleprojects.com/upload/images/tiles-polyline.jpg) ![Circle covered by grid tiles](https://aleprojects.com/upload/images/tiles-circle0.jpg)
