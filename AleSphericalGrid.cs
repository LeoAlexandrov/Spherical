/* (c) AleProjects.com, 2018
 * v.1.0
 * 
 * MIT License
 * 
 * Convention in variable names: 
 *      Lon, Lat (short) - radians; 
 *      Longitude, Latitude (full) - degrees;
 */

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace AleProjects.Spherical.Grid
{

	public interface IQuadKeyIdentified : IGeoCoordinate
	{
		long QuadKey { get; set; }
	}



	/// <summary>
	/// Helper class for building and handling of the sphere grid.
	/// </summary>
	public static class SphereGridHelper
	{
		// ABSOLUTE_MAX_LEVEL = 30 = ((64 - 1) - 3) / 2
		// "64" bits in long type,  
		// "1" is a sign bit, 
		// "3" three bits for 0..7 value indexing primary tile which is 1/8 of sphere, 
		// "2" two bits to represent 0..3 value which is index of sub-tile inside parent tile.

		/// <summary>
		/// Maximum level of the sphere grid.
		/// </summary>
		public const int ABSOLUTE_MAX_LEVEL = ((64 - 1) - 3) / 2;

		#region unit-vector
		/// <summary>
		/// Implements unit vectors for coordinate system.
		/// </summary>
		public struct UnitVector : ICartesian
		{
			public double X { get; private set; }
			public double Y { get; private set; }
			public double Z { get; private set; }

			public UnitVector(double x, double y, double z)
			{
				X = 0;
				Y = 0;
				Z = 0;

				if (x > 0.0) X = 1.0;
				else if (x < 0.0) X = -1.0;
				else if (y > 0.0) Y = 1.0;
				else if (y < 0.0) Y = -1.0;
				else if (z > 0.0) Z = 1.0;
				else if (z < 0.0) Z = -1.0;
				else throw new ArgumentException();
			}

			public void SetCartesian(double x, double y, double z)
			{
				X = 0;
				Y = 0;
				Z = 0;

				if (x > 0.0) X = 1.0;
				else if (x < 0.0) X = -1.0;
				else if (y > 0.0) Y = 1.0;
				else if (y < 0.0) Y = -1.0;
				else if (z > 0.0) Z = 1.0;
				else if (z < 0.0) Z = -1.0;
				else throw new ArgumentException();
			}
		}
		#endregion

		#region unit-vectors
		// Unit vectors for the coordinate system.
		public static readonly UnitVector PositiveI = new UnitVector(1.0, 0.0, 0.0);
		public static readonly UnitVector PositiveJ = new UnitVector(0.0, 1.0, 0.0);
		public static readonly UnitVector PositiveK = new UnitVector(0.0, 0.0, 1.0);
		public static readonly UnitVector NegativeI = new UnitVector(-1.0, 0.0, 0.0);
		public static readonly UnitVector NegativeJ = new UnitVector(0.0, -1.0, 0.0);
		public static readonly UnitVector NegativeK = new UnitVector(0.0, 0.0, -1.0);
		#endregion

		#region primary-tiles
		/// <summary>
		/// Represents sphere primary tiles (1/8 of sphere, face of octahedron projected in sphere). 
		/// </summary>
		public static readonly UnitVector[][] PrimaryTiles =
		{
			// North hemisphere
			new UnitVector[] { PositiveI, PositiveJ, PositiveK },
			new UnitVector[] { PositiveJ, NegativeI, PositiveK },
			new UnitVector[] { NegativeI, NegativeJ, PositiveK },
			new UnitVector[] { NegativeJ, PositiveI, PositiveK },
			// South hemisphere
			new UnitVector[] { PositiveJ, PositiveI, NegativeK },
			new UnitVector[] { NegativeI, PositiveJ, NegativeK },
			new UnitVector[] { NegativeJ, NegativeI, NegativeK },
			new UnitVector[] { PositiveI, NegativeJ, NegativeK }
		};
		#endregion

		#region vector3D for internal use
		/// <summary>
		/// Represents vector and used for QuadKey building
		/// </summary>
		private struct Vector3D : ICartesian
		{
			public double X { get; private set; }
			public double Y { get; private set; }
			public double Z { get; private set; }

			public Vector3D(double x, double y, double z)
			{
				double len = Math.Sqrt(x * x + y * y + z * z);

				X = x / len;
				Y = y / len;
				Z = z / len;
			}

			public void SetCartesian(double x, double y, double z)
			{
				X = x;
				Y = y;
				Z = z;
			}
		}
		#endregion


		/// <summary>
		/// Finds a primary tile (1/8 part of a sphere, face of octahedron projected on a sphere surface) containing specified location on a sphere.
		/// </summary>
		/// <param name="location">Location on the sphere.</param>
		/// <returns>Index of primary sphere tile.</returns>
		public static int FindPrimaryTile(ICartesian location)
		{
			for (int i = 0; i < PrimaryTiles.Length; i++)
				if (location.InsideTriangle(PrimaryTiles[i][0], PrimaryTiles[i][1], PrimaryTiles[i][2])) return i;

			throw new InvalidOperationException();
		}

		/// <summary>
		/// Finds a primary tile (1/8 part of a sphere, face of octahedron projected on a sphere surface) containing specified location on a sphere.
		/// </summary>
		/// <param name="lat">Latitude of the location.</param>
		/// <param name="lon">Longitude of the location.</param>
		/// <param name="degrees">When true, latitude and longitude assumed in degrees.</param>
		/// <returns></returns>
		public static int FindPrimaryTile(double lat, double lon, bool degrees = true)
		{
			return FindPrimaryTile(new Cartesian(lat, lon, degrees));
		}

		/// <summary>
		/// Calculates QuadKey of a tile at specified grid level containing a location with given latitude and longitude.
		/// </summary>
		/// <param name="lat">Latitude of location.</param>
		/// <param name="lon">Longitude of location.</param>
		/// <param name="level">Grid level.</param>
		/// <returns>QuadKey which uniquely identifies tile containing given location.</returns>
		public static long BuildQuadKey(double lat, double lon, int level)
		{
			if (level < 0 ||
				level > ABSOLUTE_MAX_LEVEL) throw new ArgumentException();

			Vector3D p = new Vector3D(Math.Cos(lat) * Math.Cos(lon), Math.Cos(lat) * Math.Sin(lon), Math.Sin(lat));

			Vector3D vertex1 = new Vector3D();
			Vector3D vertex2 = new Vector3D();
			Vector3D vertex3 = new Vector3D();
			Vector3D v1, v2, v3;
			long result = 0;

			for (int i = 0; i < PrimaryTiles.Length; i++)
				if (p.InsideTriangle(PrimaryTiles[i][0], PrimaryTiles[i][1], PrimaryTiles[i][2]))
				{
					vertex1.SetCartesian(PrimaryTiles[i][0].X, PrimaryTiles[i][0].Y, PrimaryTiles[i][0].Z);
					vertex2.SetCartesian(PrimaryTiles[i][1].X, PrimaryTiles[i][1].Y, PrimaryTiles[i][1].Z);
					vertex3.SetCartesian(PrimaryTiles[i][2].X, PrimaryTiles[i][2].Y, PrimaryTiles[i][2].Z);
					result = (long)i << (64 - 1 - 3);
					break;
				}

			for (int i = 0; i < level; i++)
			{
				v1 = new Vector3D((vertex1.X + vertex3.X) / 2.0, (vertex1.Y + vertex3.Y) / 2.0, (vertex1.Z + vertex3.Z) / 2.0);
				v2 = new Vector3D((vertex3.X + vertex2.X) / 2.0, (vertex3.Y + vertex2.Y) / 2.0, (vertex3.Z + vertex2.Z) / 2.0);
				v3 = new Vector3D((vertex2.X + vertex1.X) / 2.0, (vertex2.Y + vertex1.Y) / 2.0, (vertex2.Z + vertex1.Z) / 2.0);

				if (p.InsideTriangle(vertex1, v3, v1))
				{
					vertex2 = v3;
					vertex3 = v1;
					result |= (long)0x0000000000000000 >> (i * 2);
					continue;
				}

				if (p.InsideTriangle(v1, v2, vertex3))
				{
					vertex1 = v1;
					vertex2 = v2;
					result |= 0x0400000000000000 >> (i * 2);
					continue;
				}

				if (p.InsideTriangle(v3, vertex2, v2))
				{
					vertex1 = v3;
					vertex3 = v2;
					result |= 0x0800000000000000 >> (i * 2);
					continue;
				}

				if (p.InsideTriangle(v1, v2, v3))
				{
					vertex1 = v1;
					vertex2 = v2;
					vertex3 = v3;
					result |= 0x0c00000000000000 >> (i * 2);
					continue;
				}

				throw new InvalidOperationException();
			}

			return result;
		}

		/// <summary>
		/// Calculates QuadKey of a tile at specified grid level containing a location with given latitude and longitude.
		/// </summary>
		/// <param name="location">Location on sphere.</param>
		/// <param name="level">Grid level.</param>
		/// <returns>QuadKey which uniquely identifies tile containing given location.</returns>
		public static long BuildQuadKey(IGeoCoordinate location, int level)
		{
			return BuildQuadKey(location.Lat, location.Lon, level);
		}

		/// <summary>
		/// Returns bit mask for a Quadkey of a tile at specified grid level.
		/// </summary>
		/// <param name="level">Grid level.</param>
		/// <returns>Bit mask.</returns>
		public static long QuadKeyMask(int level)
		{
			long result = 0x70_00_00_00_00_00_00_00;
			long mask = 0x0C_00_00_00_00_00_00_00;

			while (level-- > 0)
			{
				result |= mask;
				mask >>= 2;
			}

			return result;
		}

		/// <summary>
		/// Returns a quadkey for a primary (most top) tile by its index.
		/// </summary>
		/// <param name="tileNumber">Index of primary tile.</param>
		/// <returns>quadkey for primary tile.</returns>
		public static long QuadKeyForPrimaryTile(int tileNumber)
		{
			return (long)(tileNumber & 7) << (64 - 3 - 1);
		}

		/// <summary>
		/// Decomposes QuadKey to triangles indexes at each grid level. 
		/// </summary>
		/// <param name="quadKey">Quadkey to decompose.</param>
		/// <returns>Decomposed quadkey.</returns>
		public static int[] QuadKeyParts(long quadKey)
		{
			int[] keys = new int[ABSOLUTE_MAX_LEVEL + 1];

			for (int i = ABSOLUTE_MAX_LEVEL; i > 0; i--)
			{
				keys[i] = (int)(quadKey & 3);
				quadKey >>= 2;
			}

			keys[0] = (int)(quadKey & 7);

			return keys;
		}

		/// <summary>
		/// Returns grid level where square of a tile relates to square of a circle in given ratio (approximately).  
		/// </summary>
		/// <param name="angle">Represents the circle radius.</param>
		/// <param name="triangleToCircleRatio">Tile to circle ratio.</param>
		/// <returns>Grid level.</returns>
		public static int LevelForCircleToTriangleRatio(double angle, double circleToTriangleRatio)
		{
			double circleArea = 2.0 * Math.PI * (1.0 - Math.Cos(angle)); // A = 2pi* r^2 *(1-cos(angle))
			double triangleArea = 4.0 * Math.PI / 8.0; // A = 4pi*r^2
			int level = 0;

			while (circleArea < triangleArea * circleToTriangleRatio && level < ABSOLUTE_MAX_LEVEL)
			{
				triangleArea /= 4.0;
				level++;
			}

			return level;
		}
	}


	/// <summary>
	/// Sphere grid tile. Tiles are triangles uniformly (almost) covering the whole sphere. 
	/// At first step tile is 1/8 of sphere, on every next level tile is splitted to 4 sub-tiles.
	/// </summary>
	public class SphereGridTile : IComparable<SphereGridTile>, IEquatable<SphereGridTile>
	{
		protected SphereGridTile _Parent = null;

		/// <summary>
		/// Unique identifier of this tile.
		/// </summary>
		public long QuadKey { get; protected set; }

		/// <summary>
		/// Sphere grid level of this tile.
		/// </summary>
		public int Level { get; protected set; }

		/// <summary>
		/// First vertex of the tile triangle.
		/// </summary>
		public Cartesian Vertex1 { get; protected set; }

		/// <summary>
		/// Second vertex of the tile triangle.
		/// </summary>
		public Cartesian Vertex2 { get; protected set; }

		/// <summary>
		/// Third vertex of the tile triangle.
		/// </summary>
		public Cartesian Vertex3 { get; protected set; }

		/// <summary>
		/// Parent tile at the upper level.
		/// </summary>
		public SphereGridTile Parent
		{
			get
			{
				if (_Parent == null && Level > 0)
				{
					_Parent = new SphereGridTile(QuadKey, Level - 1);
				}

				return _Parent;
			}

			protected set => _Parent = value;
		}

		/// <summary>
		/// Most top (primary) tile containing this tile.
		/// </summary>
		public SphereGridTile Top
		{
			get
			{
				SphereGridTile result = this;

				while (result.Parent != null) result = result.Parent;

				return result;
			}
		}

		/// <summary>
		/// Index of this tile inside upper parent tile.
		/// </summary>
		public int LevelKey
		{
			get
			{
				long result = QuadKey >> ((SphereGridHelper.ABSOLUTE_MAX_LEVEL - Level) * 2);

				return Level == 0 ? (int)(result & 7) : (int)(result & 3);
			}
		}

		/// <summary>
		/// QuadKey of the parent tile. 
		/// </summary>
		public long ParentQuadKey
		{
			get => Level > 0 ? QuadKey & (long)(0x7f_ff_ff_ff_ff_ff_ff_ff & (0xff_ff_ff_ff_ff_ff_ff_ff << ((SphereGridHelper.ABSOLUTE_MAX_LEVEL - Level + 1) * 2))) : 0;
		}

		/// <summary>
		/// Upper value of QuadKey. Allows to select objects inside this tile using filter &gt;= QuadKey and &lt;= QuadKeyUpperValue. 
		/// </summary>
		public long QuadKeyUpperValue
		{
			get => QuadKey | (0x4f_ff_ff_ff_ff_ff_ff_ff >> (Level * 2 + 2));
		}

		/// <summary>
		/// Hidden basic constructor. 
		/// </summary>
		protected SphereGridTile()
		{
		}

		/// <summary>
		/// Creates a tile at the specified grid level for a given location.
		/// </summary>
		/// <param name="location">Location to buid tile for.</param>
		/// <param name="level">Grid level.</param>
		public SphereGridTile(ICartesian location, int level)
		{
			if (location == null ||
				level < 0 ||
				level > SphereGridHelper.ABSOLUTE_MAX_LEVEL) throw new ArgumentException();

			int k = SphereGridHelper.FindPrimaryTile(location);
			long quadkey = (long)k << (64 - 1 - 3); // 1 - sign, 3 - bits for 0..7 values

			if (level > 0)
			{
				SphereGridTile tile = new SphereGridTile()
				{
					QuadKey = quadkey,
					Vertex1 = new Cartesian(SphereGridHelper.PrimaryTiles[k][0]),
					Vertex2 = new Cartesian(SphereGridHelper.PrimaryTiles[k][1]),
					Vertex3 = new Cartesian(SphereGridHelper.PrimaryTiles[k][2]),
				};

				SphereGridTile nextTile;

				for (int i = 0; i < level; i++)
					if ((nextTile = tile.SplitAndFind(location)) != null)
					{
						nextTile.Parent = tile;
						tile = nextTile;
					}
					else throw new InvalidOperationException();

				QuadKey = tile.QuadKey;
				Level = level;
				Vertex1 = tile.Vertex1;
				Vertex2 = tile.Vertex2;
				Vertex3 = tile.Vertex3;
				Parent = tile.Parent;
			}
			else
			{
				QuadKey = quadkey;
				Vertex1 = new Cartesian(SphereGridHelper.PrimaryTiles[k][0]);
				Vertex2 = new Cartesian(SphereGridHelper.PrimaryTiles[k][1]);
				Vertex3 = new Cartesian(SphereGridHelper.PrimaryTiles[k][2]);
			}
		}

		/// <summary>
		/// Creates a tile at the specified grid level using existing QuadKey.
		/// </summary>
		/// <param name="quadKey">Existing QuadKey.</param>
		/// <param name="level">Grid level.</param>
		public SphereGridTile(long quadKey, int level)
		{
			if (quadKey < 0 ||
				level < 0 ||
				level > SphereGridHelper.ABSOLUTE_MAX_LEVEL) throw new ArgumentException();


			if (level > 0)
			{
				int[] keys = SphereGridHelper.QuadKeyParts(quadKey);
				int k = keys[0];

				SphereGridTile tile = new SphereGridTile()
				{
					Vertex1 = new Cartesian(SphereGridHelper.PrimaryTiles[k][0]),
					Vertex2 = new Cartesian(SphereGridHelper.PrimaryTiles[k][1]),
					Vertex3 = new Cartesian(SphereGridHelper.PrimaryTiles[k][2]),
					QuadKey = quadKey & 0x70_00_00_00_00_00_00_00
				};

				SphereGridTile nextTile;
				SphereGridTile[] tiles;

				for (int i = 1; i <= level; i++)
				{
					tiles = tile.Split();
					nextTile = tiles[keys[i]];
					nextTile.Parent = tile;
					tile = nextTile;
				}

				QuadKey = tile.QuadKey;
				Parent = tile.Parent;
				Level = level;
				Vertex1 = tile.Vertex1;
				Vertex2 = tile.Vertex2;
				Vertex3 = tile.Vertex3;
			}
			else
			{
				int k = (int)(quadKey >> (64 - 1 - 3)); // 1 - sign, 3 - bits for 0..7 values

				QuadKey = quadKey;
				Vertex1 = new Cartesian(SphereGridHelper.PrimaryTiles[k][0]);
				Vertex2 = new Cartesian(SphereGridHelper.PrimaryTiles[k][1]);
				Vertex3 = new Cartesian(SphereGridHelper.PrimaryTiles[k][2]);
			}
		}


		/// <summary>
		/// Splits this tile and finds a sub-tile containing given location.
		/// </summary>
		/// <param name="location">Location.</param>
		/// <returns>Sub-tile containing location.</returns>
		public SphereGridTile SplitAndFind(ICartesian location)
		{
			Cartesian v1 = Vertex1.DivideLineInRatio<Cartesian>(Vertex3, 0.5, true);
			Cartesian v2 = Vertex3.DivideLineInRatio<Cartesian>(Vertex2, 0.5, true);
			Cartesian v3 = Vertex2.DivideLineInRatio<Cartesian>(Vertex1, 0.5, true);
			int level = this.Level + 1;

			if (location.InsideTriangle(this.Vertex1, v3, v1))
				return new SphereGridTile()
				{
					QuadKey = this.QuadKey | (0x0000000000000000L >> (level * 2)),
					Level = level,
					Vertex1 = this.Vertex1,
					Vertex2 = v3,
					Vertex3 = v1,
				};

			if (location.InsideTriangle(v1, v2, this.Vertex3))
				return new SphereGridTile()
				{
					QuadKey = this.QuadKey | (0x1000000000000000L >> (level * 2)),
					Level = level,
					Vertex1 = v1,
					Vertex2 = v2,
					Vertex3 = this.Vertex3,
				};

			if (location.InsideTriangle(v3, this.Vertex2, v2))
				return new SphereGridTile()
				{
					QuadKey = this.QuadKey | (0x2000000000000000L >> (level * 2)),
					Level = level,
					Vertex1 = v3,
					Vertex2 = this.Vertex2,
					Vertex3 = v2,
				};

			if (location.InsideTriangle(v1, v2, v3))
				return new SphereGridTile()
				{
					QuadKey = this.QuadKey | (0x3000000000000000L >> (level * 2)),
					Level = level,
					Vertex1 = v1,
					Vertex2 = v2,
					Vertex3 = v3,
				};

			return null;
		}

		/// <summary>
		/// Splits this tile and returns sub-tiles.
		/// </summary>
		/// <returns>Array of sub-tiles at the next grid level.</returns>
		public SphereGridTile[] Split()
		{
			Cartesian v1 = Vertex1.DivideLineInRatio<Cartesian>(Vertex3, 0.5, true);
			Cartesian v2 = Vertex3.DivideLineInRatio<Cartesian>(Vertex2, 0.5, true);
			Cartesian v3 = Vertex2.DivideLineInRatio<Cartesian>(Vertex1, 0.5, true);
			int level = this.Level + 1;

			return new SphereGridTile[]
			{
				new SphereGridTile()
				{
					QuadKey = this.QuadKey | (0x0000000000000000L >> (level * 2)),
					Level = level,
					Vertex1 = this.Vertex1,
					Vertex2 = v3,
					Vertex3 = v1
				},
				new SphereGridTile()
				{
					QuadKey = this.QuadKey | (0x1000000000000000L >> (level * 2)),
					Level = level,
					Vertex1 = v1,
					Vertex2 = v2,
					Vertex3 = this.Vertex3
				},
				new SphereGridTile()
				{
					QuadKey = this.QuadKey | (0x2000000000000000L >> (level * 2)),
					Level = level,
					Vertex1 = v3,
					Vertex2 = this.Vertex2,
					Vertex3 = v2
				},
				new SphereGridTile()
				{
					QuadKey = this.QuadKey | (0x3000000000000000L >> (level * 2)),
					Level = level,
					Vertex1 = v1,
					Vertex2 = v2,
					Vertex3 = v3
				}
			};
		}

		/// <summary>
		/// Checks if this tile is covered (fully or partially) by a circle with given center and radius represented as an angle.  
		/// </summary>
		/// <param name="center">Center of the circle.</param>
		/// <param name="angle">Angle representing circle radius.</param>
		/// <returns>True if the circle covers this tile.</returns>
		public bool CoveredByCircle(ICartesian center, double angle)
		{
			double circleCosine = Math.Cos(angle);
			double phi;

			// any vertex inside circle or circle center inside triangle
			if (center.Cosine(Vertex1) > circleCosine ||
				center.Cosine(Vertex2) > circleCosine ||
				center.Cosine(Vertex3) > circleCosine ||
				center.InsideTriangle(Vertex1, Vertex2, Vertex3)) return true;

			// check side Vertex1-Vertex2
			Cartesian vector = Vertex1.CrossProduct<Cartesian>(Vertex2, true);

			double circlePlaneA = center.Y * vector.Z - center.Z * vector.Y;
			double circlePlaneB = -(center.X * vector.Z - center.Z * vector.X);
			double circlePlaneC = center.X * vector.Y - center.Y * vector.X;

			if ((circlePlaneA * Vertex1.X + circlePlaneB * Vertex1.Y + circlePlaneC * Vertex1.Z) *
				(circlePlaneA * Vertex2.X + circlePlaneB * Vertex2.Y + circlePlaneC * Vertex2.Z) < 0.0) // from opposite plane sides
			{
				phi = Math.Acos(vector.Cosine(center));
				if (phi > Math.PI / 2) phi = Math.PI - phi;
				if (angle + phi > Math.PI / 2) return true;
			}

			// check side Vertex2-Vertex3
			vector = Vertex2.CrossProduct<Cartesian>(Vertex3, true);

			circlePlaneA = center.Y * vector.Z - center.Z * vector.Y;
			circlePlaneB = -(center.X * vector.Z - center.Z * vector.X);
			circlePlaneC = center.X * vector.Y - center.Y * vector.X;

			if ((circlePlaneA * Vertex2.X + circlePlaneB * Vertex2.Y + circlePlaneC * Vertex2.Z) *
				(circlePlaneA * Vertex3.X + circlePlaneB * Vertex3.Y + circlePlaneC * Vertex3.Z) < 0.0)
			{
				phi = Math.Acos(vector.Cosine(center));
				if (phi > Math.PI / 2) phi = Math.PI - phi;
				if (angle + phi > Math.PI / 2) return true;
			}

			// check side Vertex3-Vertex1
			vector = Vertex3.CrossProduct<Cartesian>(Vertex1, true);

			circlePlaneA = center.Y * vector.Z - center.Z * vector.Y;
			circlePlaneB = -(center.X * vector.Z - center.Z * vector.X);
			circlePlaneC = center.X * vector.Y - center.Y * vector.X;

			if ((circlePlaneA * Vertex3.X + circlePlaneB * Vertex3.Y + circlePlaneC * Vertex3.Z) *
				(circlePlaneA * Vertex1.X + circlePlaneB * Vertex1.Y + circlePlaneC * Vertex1.Z) < 0.0)
			{
				phi = Math.Acos(vector.Cosine(center));
				if (phi > Math.PI / 2) phi = Math.PI - phi;
				return angle + phi > Math.PI / 2;
			}

			return false;
		}

		/// <summary>
		/// Checks if this tile is covered (fully or partially) by a polyline.
		/// </summary>
		/// <param name="polyline">Vertices of the polyline.</param>
		/// <param name="tolerance">Test tolerance in radians.</param>
		/// <returns>True if the polyline covers this tile.</returns>
		public bool CoveredByPolyline(IEnumerable<ICartesian> polyline, double tolerance)
		{
			ICartesian vertex = polyline.First();

			foreach (ICartesian v in polyline.Skip(1))
			{
				if (SphericalExtension.SectionsIntersect(Vertex1, Vertex2, vertex, v, tolerance) >= 0 ||
					SphericalExtension.SectionsIntersect(Vertex2, Vertex3, vertex, v, tolerance) >= 0 ||
					SphericalExtension.SectionsIntersect(Vertex3, Vertex1, vertex, v, tolerance) >= 0)
					return true;

				vertex = v;
			}

			return EnclosesPolygon(polyline);
		}

		/// <summary>
		/// Checks if this tile is covered (fully or partially) by a polygon.
		/// </summary>
		/// <param name="polygon">Vertices of the polygon.</param>
		/// <returns>True if the polygon covers this tile.</returns>
		public bool CoveredByPolygon(IEnumerable<ICartesian> polygon)
		{
			ICartesian vertex = polygon.First();
			ICartesian last = polygon.Last();

			if (SphericalExtension.SectionsIntersect(Vertex1, Vertex2, vertex, last) >= 0 ||
				SphericalExtension.SectionsIntersect(Vertex2, Vertex3, vertex, last) >= 0 ||
				SphericalExtension.SectionsIntersect(Vertex3, Vertex1, vertex, last) >= 0) return true;

			foreach (ICartesian v in polygon.Skip(1))
			{
				if (SphericalExtension.SectionsIntersect(Vertex1, Vertex2, vertex, v) >= 0 ||
					SphericalExtension.SectionsIntersect(Vertex2, Vertex3, vertex, v) >= 0 ||
					SphericalExtension.SectionsIntersect(Vertex3, Vertex1, vertex, v) >= 0) return true;

				vertex = v;
			}

			return EnclosesPolygon(polygon) ||
				Vertex1.InsidePolygon(polygon) ||
				Vertex2.InsidePolygon(polygon) ||
				Vertex3.InsidePolygon(polygon);
		}

		/// <summary>
		/// Checks if this tile fully encloses a circle with given center and radius represented as an angle.
		/// </summary>
		/// <param name="center">Center of the circle.</param>
		/// <param name="angle">Angle representing circle radius.</param>
		/// <returns>True if this tile fully encloses the circle.</returns>
		public bool EnclosesCircle(ICartesian center, double angle)
		{
			double circleCosine = Math.Cos(angle);

			// any vertex outside circle or circle center outside triangle
			if (center.Cosine(Vertex1) > circleCosine ||
				center.Cosine(Vertex2) > circleCosine ||
				center.Cosine(Vertex3) > circleCosine ||
				!center.InsideTriangle(Vertex1, Vertex2, Vertex3)) return false;

			// check side Vertex1-Vertex2
			Cartesian vector = Vertex1.CrossProduct<Cartesian>(Vertex2, true);
			double phi = Math.Acos(vector.Cosine(center));
			if (phi > Math.PI / 2) phi = Math.PI - phi;
			if (angle + phi > Math.PI / 2) return false;

			// check side Vertex2-Vertex3
			vector = Vertex2.CrossProduct<Cartesian>(Vertex3, true);
			phi = Math.Acos(vector.Cosine(center));
			if (phi > Math.PI / 2) phi = Math.PI - phi;
			if (angle + phi > Math.PI / 2) return false;

			// check side Vertex3-Vertex1
			vector = Vertex3.CrossProduct<Cartesian>(Vertex1, true);
			phi = Math.Acos(vector.Cosine(center));
			if (phi > Math.PI / 2) phi = Math.PI - phi;
			return angle + phi <= Math.PI / 2; // <= because 'false' above
		}

		/// <summary>
		/// Checks if this tile fully encloses a polygon. Can be used for polylines too.
		/// </summary>
		/// <param name="polygon">Vertices of the polygon or polyline.</param>
		/// <returns>True if encloses.</returns>
		public bool EnclosesPolygon(IEnumerable<ICartesian> polygon)
		{
			foreach (var vertex in polygon)
				if (!vertex.InsideTriangle(Vertex1, Vertex2, Vertex3)) return false;

			return true;
		}

		/// <summary>
		/// Joins tiles to upper possible level.
		/// </summary>
		/// <param name="tiles">List of the tiles to join.</param>
		/// <param name="level">Grid level.</param>
		/// <param name="sort">Sort if necessary.</param>
		protected static void JoinTiles(List<SphereGridTile> tiles, int level, bool sort = true)
		{
			if (sort) tiles.Sort();

			bool hasJoined = true;
			int k = 0;

			for (int i = level; i > 0 && hasJoined; i--)
			{
				hasJoined = false;

				for (int j = 0; j < tiles.Count; j++)
					if (tiles[j].Level == i &&
						j < tiles.Count - 1 &&
						tiles[j + 1].Level == i &&
						tiles[j + 1].ParentQuadKey == tiles[j].ParentQuadKey)
					{
						k++;

						if (k == 3)
						{
							j -= 2;
							tiles[j] = tiles[j].Parent;
							tiles.RemoveRange(j + 1, 3);
							k = 0;
							hasJoined = true;
						}
					}
					else k = 0;
			}
		}

		/// <summary>
		/// Returns a list of tiles covering (fully or partially) a circle with given center and radius represented as an angle. 
		/// </summary>
		/// <param name="center">Center of circle.</param>
		/// <param name="angle">Angle representing circle radius.</param>
		/// <param name="level">Maximum grid level of the tiles.</param>
		/// <param name="join">Join tiles if possible.</param>
		/// <returns>List of tiles covering circle.</returns>
		public static List<SphereGridTile> CoverCircleByTiles(ICartesian center, double angle, int level, bool join)
		{
			if (level < 0 ||
				level > SphereGridHelper.ABSOLUTE_MAX_LEVEL) throw new ArgumentException();

			SphereGridTile tile = new SphereGridTile(center, 0);
			SphereGridTile t;
			List<SphereGridTile> result = new List<SphereGridTile>() { tile };
			int n;
			int k = tile.LevelKey;

			if (!tile.EnclosesCircle(center, angle))
				for (int i = 0; i < SphereGridHelper.PrimaryTiles.Length; i++)
					if (i != k && (t = new SphereGridTile(SphereGridHelper.QuadKeyForPrimaryTile(i), 0)).CoveredByCircle(center, angle))
						result.Add(t);

			for (int i = 0; i < level; i++)
			{
				n = result.Count;

				for (int j = 0; j < n; j++)
					result.AddRange(result[j].Split());

				result.RemoveRange(0, n);
				result.RemoveAll(tl => !tl.CoveredByCircle(center, angle));
			}

			if (join && level > 0 && result.Count > 1) JoinTiles(result, level);

			return result;
		}

		/// <summary>
		/// Returns a list of tiles covering (fully or partially) a polyline. 
		/// </summary>
		/// <param name="polyline">Vertices of the polyline.</param>
		/// <param name="level">Maximum grid level of tiles.</param>
		/// <param name="join">Join tiles if possible.</param>
		/// <returns>List of tiles covering the polyline.</returns>
		public static List<SphereGridTile> CoverPolylineByTiles(IEnumerable<ICartesian> polyline, int level, double tolerance, bool join)
		{
			if (level < 0 ||
				level > SphereGridHelper.ABSOLUTE_MAX_LEVEL) throw new ArgumentException();

			SphereGridTile tile = new SphereGridTile(polyline.First(), 0);
			SphereGridTile t;
			List<SphereGridTile> result = new List<SphereGridTile>() { tile };
			int n;
			int k = tile.LevelKey;

			for (int i = 0; i < SphereGridHelper.PrimaryTiles.Length; i++)
				if (i != k && (t = new SphereGridTile(SphereGridHelper.QuadKeyForPrimaryTile(i), 0)).CoveredByPolyline(polyline, tolerance))
					result.Add(t);

			for (int i = 0; i < level; i++)
			{
				n = result.Count;

				for (int j = 0; j < n; j++)
					result.AddRange(result[j].Split());

				result.RemoveRange(0, n);
				result.RemoveAll(tl => !tl.CoveredByPolyline(polyline, tolerance));
			}

			if (join && level > 0 && result.Count > 1) JoinTiles(result, level);

			return result;
		}

		/// <summary>
		/// Returns a list of tiles covering (fully or partially) a polygon. 
		/// </summary>
		/// <param name="polygon">Vertices of the polygon.</param>
		/// <param name="level">Maximum grid level of tiles.</param>
		/// <param name="join">Join tiles if possible.</param>
		/// <returns>List of tiles covering the polygon.</returns>
		public static List<SphereGridTile> CoverPolygonByTiles(IEnumerable<ICartesian> polygon, int level, bool join)
		{
			if (level < 0 ||
				level > SphereGridHelper.ABSOLUTE_MAX_LEVEL) throw new ArgumentException();

			SphereGridTile tile = new SphereGridTile(polygon.First(), 0);
			SphereGridTile t;
			List<SphereGridTile> result = new List<SphereGridTile>() { tile };
			int n;
			int k = tile.LevelKey;

			if (!tile.EnclosesPolygon(polygon))
				for (int i = 0; i < SphereGridHelper.PrimaryTiles.Length; i++)
					if (i != k && (t = new SphereGridTile(SphereGridHelper.QuadKeyForPrimaryTile(i), 0)).CoveredByPolygon(polygon))
						result.Add(t);

			for (int i = 0; i < level; i++)
			{
				n = result.Count;

				for (int j = 0; j < n; j++)
					result.AddRange(result[j].Split());

				result.RemoveRange(0, n);
				result.RemoveAll(tl => !tl.CoveredByPolygon(polygon));
			}

			if (join && level > 0 && result.Count > 1) JoinTiles(result, level);

			return result;
		}

		public int CompareTo(SphereGridTile other)
		{
			return QuadKey.CompareTo(other.QuadKey);
		}

		public bool Equals(SphereGridTile other)
		{
			return other != null && QuadKey.Equals(other.QuadKey);
		}

		public override int GetHashCode()
		{
			return QuadKey.GetHashCode();
		}

		public override string ToString()
		{
			int[] keys = SphereGridHelper.QuadKeyParts(QuadKey);

			StringBuilder result = new StringBuilder(keys[0].ToString(), 128);

			for (int i = 1; i <= SphereGridHelper.ABSOLUTE_MAX_LEVEL; i++)
				result.Append(",").Append(keys[i]);

			return result.ToString();
		}
	}

}
