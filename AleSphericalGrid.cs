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

	public interface IQuadKeyIdentified : ICartesian
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

		#region unit-vectors
		// Unit vectors for the coordinate system.
		public static readonly CartesianValue PositiveI = new CartesianValue(1.0, 0.0, 0.0);
		public static readonly CartesianValue PositiveJ = new CartesianValue(0.0, 1.0, 0.0);
		public static readonly CartesianValue PositiveK = new CartesianValue(0.0, 0.0, 1.0);
		public static readonly CartesianValue NegativeI = new CartesianValue(-1.0, 0.0, 0.0);
		public static readonly CartesianValue NegativeJ = new CartesianValue(0.0, -1.0, 0.0);
		public static readonly CartesianValue NegativeK = new CartesianValue(0.0, 0.0, -1.0);
		#endregion

		#region primary-tiles
		/// <summary>
		/// Represents sphere primary tiles (1/8 of sphere, face of octahedron projected in sphere). 
		/// </summary>
		public static readonly CartesianValue[][] PrimaryTiles =
		{
			// North hemisphere
			new CartesianValue[] { PositiveI, PositiveJ, PositiveK },
			new CartesianValue[] { PositiveJ, NegativeI, PositiveK },
			new CartesianValue[] { NegativeI, NegativeJ, PositiveK },
			new CartesianValue[] { NegativeJ, PositiveI, PositiveK },
			// South hemisphere
			new CartesianValue[] { PositiveJ, PositiveI, NegativeK },
			new CartesianValue[] { NegativeI, PositiveJ, NegativeK },
			new CartesianValue[] { NegativeJ, NegativeI, NegativeK },
			new CartesianValue[] { PositiveI, NegativeJ, NegativeK }
		};
		#endregion



		/// <summary>
		/// Finds a primary tile (1/8 part of a sphere, face of octahedron projected on a sphere surface) containing specified location on a sphere.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="location">Location on the sphere.</param>
		/// <returns>Index of primary sphere tile.</returns>
		public static int FindPrimaryTile<T>(T location)
			where T : ICartesian
		{
			for (int i = 0; i < PrimaryTiles.Length; i++)
				if (SphericalExtension._InsideTriangle(location.X, location.Y, location.Z, PrimaryTiles[i][0], PrimaryTiles[i][1], PrimaryTiles[i][2]))
					return i;

			throw new InvalidOperationException();
		}

		/// <summary>
		/// Finds a primary tile (1/8 part of a sphere, face of octahedron projected on a sphere surface) containing specified location on a sphere.
		/// </summary>
		/// <param name="lat">Latitude of the location in degrees.</param>
		/// <param name="lon">Longitude of the locationin degrees.</param>
		/// <returns></returns>
		public static int FindPrimaryTile(double lat, double lon)
		{
			CartesianValue p = new CartesianValue(lat, lon);

			return FindPrimaryTile(p);
		}

		/// <summary>
		/// Calculates QuadKey of a tile at specified grid level containing a location with given latitude and longitude.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="lat">Latitude of location.</param>
		/// <param name="lon">Longitude of location.</param>
		/// <param name="level">Grid level.</param>
		/// <returns>QuadKey which uniquely identifies tile containing given location.</returns>
		public static long BuildQuadKey<T>(T location, int level)
			where T : ICartesian
		{
			if (level < 0 || level > ABSOLUTE_MAX_LEVEL)
				throw new ArgumentOutOfRangeException(nameof(level));

			double x = location.X;
			double y = location.Y;
			double z = location.Z;

			CartesianValue vertex1 = new CartesianValue();
			CartesianValue vertex2 = new CartesianValue();
			CartesianValue vertex3 = new CartesianValue();
			CartesianValue v1 = new CartesianValue();
			CartesianValue v2 = new CartesianValue();
			CartesianValue v3 = new CartesianValue();

			long result = 0;

			for (int i = 0; i < PrimaryTiles.Length; i++)
				if (SphericalExtension._InsideTriangle(x, y, z, PrimaryTiles[i][0], PrimaryTiles[i][1], PrimaryTiles[i][2]))
				{
					vertex1.SetCartesian(PrimaryTiles[i][0].X, PrimaryTiles[i][0].Y, PrimaryTiles[i][0].Z);
					vertex2.SetCartesian(PrimaryTiles[i][1].X, PrimaryTiles[i][1].Y, PrimaryTiles[i][1].Z);
					vertex3.SetCartesian(PrimaryTiles[i][2].X, PrimaryTiles[i][2].Y, PrimaryTiles[i][2].Z);
					result = (long)i << (64 - 1 - 3);
					break;
				}

			for (int i = 0; i < level; i++)
			{
				v1.SetNormalized(vertex1.X + vertex3.X, vertex1.Y + vertex3.Y, vertex1.Z + vertex3.Z);
				v2.SetNormalized(vertex3.X + vertex2.X, vertex3.Y + vertex2.Y, vertex3.Z + vertex2.Z);
				v3.SetNormalized(vertex2.X + vertex1.X, vertex2.Y + vertex1.Y, vertex2.Z + vertex1.Z);

				if (SphericalExtension._InsideTriangle(x, y, z, vertex1, v3, v1))
				{
					vertex2 = v3;
					vertex3 = v1;
					result |= (long)0x0000000000000000 >> (i * 2);
					continue;
				}

				if (SphericalExtension._InsideTriangle(x, y, z, v1, v2, vertex3))
				{
					vertex1 = v1;
					vertex2 = v2;
					result |= 0x0400000000000000 >> (i * 2);
					continue;
				}

				if (SphericalExtension._InsideTriangle(x, y, z, v3, vertex2, v2))
				{
					vertex1 = v3;
					vertex3 = v2;
					result |= 0x0800000000000000 >> (i * 2);
					continue;
				}

				vertex1 = v1;
				vertex2 = v2;
				vertex3 = v3;
				result |= 0x0c00000000000000 >> (i * 2);
			}

			return result;
		}

		/// <summary>
		/// Calculates QuadKey of a tile at specified grid level containing a location with given latitude and longitude.
		/// </summary>
		/// <param name="location">Location on sphere.</param>
		/// <param name="level">Grid level.</param>
		/// <returns>QuadKey which uniquely identifies tile containing given location.</returns>
		public static long BuildQuadKey(double latitude, double longitude, int level)
		{
			CartesianValue location = new CartesianValue(latitude, longitude, true);
			return BuildQuadKey(location, level);
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
		/// <param name="tileNumber">Index of the primary tile.</param>
		/// <returns>Quadkey for the primary tile.</returns>
		public static long QuadKeyForPrimaryTile(int tileNumber)
		{
			return (long)(tileNumber & 7) << (64 - 3 - 1);
		}

		/// <summary>
		/// Decomposes QuadKey to triangle indexes at each grid level. 
		/// </summary>
		/// <param name="quadKey">Quadkey to decompose.</param>
		/// <returns>Decomposed quadkey.</returns>

#if NETCOREAPP2_1_OR_GREATER
		public static void QuadKeyParts(long quadKey, Span<int> keys)
#else
		public static void QuadKeyParts(long quadKey, int[] keys)
#endif
		{
			//= new int[ABSOLUTE_MAX_LEVEL + 1];

			for (int i = ABSOLUTE_MAX_LEVEL; i > 0; i--)
			{
				keys[i] = (int)(quadKey & 3);
				quadKey >>= 2;
			}

			keys[0] = (int)(quadKey & 7);
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
		public CartesianValue Vertex1 { get; protected set; }

		/// <summary>
		/// Second vertex of the tile triangle.
		/// </summary>
		public CartesianValue Vertex2 { get; protected set; }

		/// <summary>
		/// Third vertex of the tile triangle.
		/// </summary>
		public CartesianValue Vertex3 { get; protected set; }

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

				while (result.Parent != null)
					result = result.Parent;

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
			//get => QuadKey | (0x4f_ff_ff_ff_ff_ff_ff_ff >> (Level * 2 + 2));
			get => QuadKey | (0x7f_ff_ff_ff_ff_ff_ff_ff >> (Level * 2 + 3));
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
			if (location == null)
				throw new ArgumentNullException(nameof(location));

			if (level < 0 || level > SphereGridHelper.ABSOLUTE_MAX_LEVEL)
				throw new ArgumentOutOfRangeException(nameof(level));

			CartesianValue locationValue = new CartesianValue(location);
			int k = SphereGridHelper.FindPrimaryTile(locationValue);
			long quadkey = (long)k << (64 - 1 - 3); // 1 - sign bit, 3 - bits for 0..7 values

			if (level > 0)
			{
				SphereGridTile tile = new SphereGridTile()
				{
					QuadKey = quadkey,
					Vertex1 = new CartesianValue(SphereGridHelper.PrimaryTiles[k][0]),
					Vertex2 = new CartesianValue(SphereGridHelper.PrimaryTiles[k][1]),
					Vertex3 = new CartesianValue(SphereGridHelper.PrimaryTiles[k][2]),
				};

				SphereGridTile nextTile;

				for (int i = 0; i < level; i++)
					if ((nextTile = tile.SplitAndFind(locationValue)) != null)
					{
						nextTile.Parent = tile;
						tile = nextTile;
					}

				QuadKey = tile.QuadKey;
				Level = level;
				Parent = tile.Parent;
				Vertex1 = tile.Vertex1;
				Vertex2 = tile.Vertex2;
				Vertex3 = tile.Vertex3;
			}
			else
			{
				QuadKey = quadkey;
				Vertex1 = new CartesianValue(SphereGridHelper.PrimaryTiles[k][0]);
				Vertex2 = new CartesianValue(SphereGridHelper.PrimaryTiles[k][1]);
				Vertex3 = new CartesianValue(SphereGridHelper.PrimaryTiles[k][2]);
			}
		}

		/// <summary>
		/// Creates a tile at the specified grid level using existing QuadKey.
		/// </summary>
		/// <param name="quadKey">Existing QuadKey.</param>
		/// <param name="level">Grid level.</param>
		public SphereGridTile(long quadKey, int level)
		{
			if (quadKey < 0)
				throw new ArgumentOutOfRangeException(nameof(quadKey));

			if (level < 0 || level > SphereGridHelper.ABSOLUTE_MAX_LEVEL)
				throw new ArgumentOutOfRangeException(nameof(level));


			if (level > 0)
			{
#if NETCOREAPP2_1_OR_GREATER
				Span<int> keys = stackalloc int[SphereGridHelper.ABSOLUTE_MAX_LEVEL + 1];
#else
				int[] keys = new int[SphereGridHelper.ABSOLUTE_MAX_LEVEL + 1];
#endif
				SphereGridHelper.QuadKeyParts(quadKey, keys);
				int k = keys[0];

				SphereGridTile tile = new SphereGridTile()
				{
					Vertex1 = new CartesianValue(SphereGridHelper.PrimaryTiles[k][0]),
					Vertex2 = new CartesianValue(SphereGridHelper.PrimaryTiles[k][1]),
					Vertex3 = new CartesianValue(SphereGridHelper.PrimaryTiles[k][2]),
					QuadKey = quadKey & 0x70_00_00_00_00_00_00_00
				};

				SphereGridTile nextTile;
				SphereGridTile[] tiles = new SphereGridTile[4];


				for (int i = 1; i <= level; i++)
				{
					tile.Split(tiles);
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
				int k = (int)(quadKey >> (64 - 1 - 3)); // 1 - sign bit, 3 - bits for 0..7 values

				QuadKey = quadKey;
				Vertex1 = new CartesianValue(SphereGridHelper.PrimaryTiles[k][0]);
				Vertex2 = new CartesianValue(SphereGridHelper.PrimaryTiles[k][1]);
				Vertex3 = new CartesianValue(SphereGridHelper.PrimaryTiles[k][2]);
			}
		}


		/// <summary>
		/// Splits this tile and finds a sub-tile containing given location.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="location">Location.</param>
		/// <returns>Sub-tile containing location.</returns>
		public SphereGridTile SplitAndFind<T>(T location)
			where T : ICartesian
		{
			var (x, y, z) = SphericalExtension._Normalized(Vertex1.X + Vertex3.X, Vertex1.Y + Vertex3.Y, Vertex1.Z + Vertex3.Z);
			CartesianValue v1 = new CartesianValue(x, y, z);
			(x, y, z) = SphericalExtension._Normalized(Vertex3.X + Vertex2.X, Vertex3.Y + Vertex2.Y, Vertex3.Z + Vertex2.Z);
			CartesianValue v2 = new CartesianValue(x, y, z);
			(x, y, z) = SphericalExtension._Normalized(Vertex2.X + Vertex1.X, Vertex2.Y + Vertex1.Y, Vertex2.Z + Vertex1.Z);
			CartesianValue v3 = new CartesianValue(x, y, z);

			int level = this.Level + 1;

			x = location.X;
			y = location.Y;
			z = location.Z;

			if (SphericalExtension._InsideTriangle(x, y, z, this.Vertex1, v3, v1))
				return new SphereGridTile()
				{
					QuadKey = this.QuadKey | (0x0000000000000000L >> (level * 2)),
					Level = level,
					Vertex1 = this.Vertex1,
					Vertex2 = v3,
					Vertex3 = v1,
				};

			if (SphericalExtension._InsideTriangle(x, y, z, v1, v2, this.Vertex3))
				return new SphereGridTile()
				{
					QuadKey = this.QuadKey | (0x1000000000000000L >> (level * 2)),
					Level = level,
					Vertex1 = v1,
					Vertex2 = v2,
					Vertex3 = this.Vertex3,
				};

			if (SphericalExtension._InsideTriangle(x, y, z, v3, this.Vertex2, v2))
				return new SphereGridTile()
				{
					QuadKey = this.QuadKey | (0x2000000000000000L >> (level * 2)),
					Level = level,
					Vertex1 = v3,
					Vertex2 = this.Vertex2,
					Vertex3 = v2,
				};

			return new SphereGridTile()
			{
				QuadKey = this.QuadKey | (0x3000000000000000L >> (level * 2)),
				Level = level,
				Vertex1 = v1,
				Vertex2 = v2,
				Vertex3 = v3,
			};
		}

		/// <summary>
		/// Splits this tile and returns sub-tiles.
		/// </summary>
		/// <returns>Array of sub-tiles at the next grid level.</returns>
		public SphereGridTile[] Split()
		{
			var (x, y, z) = SphericalExtension._Normalized(Vertex1.X + Vertex3.X, Vertex1.Y + Vertex3.Y, Vertex1.Z + Vertex3.Z);
			CartesianValue v1 = new CartesianValue(x, y, z);
			(x, y, z) = SphericalExtension._Normalized(Vertex3.X + Vertex2.X, Vertex3.Y + Vertex2.Y, Vertex3.Z + Vertex2.Z);
			CartesianValue v2 = new CartesianValue(x, y, z);
			(x, y, z) = SphericalExtension._Normalized(Vertex2.X + Vertex1.X, Vertex2.Y + Vertex1.Y, Vertex2.Z + Vertex1.Z);
			CartesianValue v3 = new CartesianValue(x, y, z);

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

		public void Split(SphereGridTile[] list)
		{
			var (x, y, z) = SphericalExtension._Normalized(Vertex1.X + Vertex3.X, Vertex1.Y + Vertex3.Y, Vertex1.Z + Vertex3.Z);
			CartesianValue v1 = new CartesianValue(x, y, z);
			(x, y, z) = SphericalExtension._Normalized(Vertex3.X + Vertex2.X, Vertex3.Y + Vertex2.Y, Vertex3.Z + Vertex2.Z);
			CartesianValue v2 = new CartesianValue(x, y, z);
			(x, y, z) = SphericalExtension._Normalized(Vertex2.X + Vertex1.X, Vertex2.Y + Vertex1.Y, Vertex2.Z + Vertex1.Z);
			CartesianValue v3 = new CartesianValue(x, y, z);

			int level = this.Level + 1;

			list[0] = new SphereGridTile()
			{
				QuadKey = this.QuadKey | (0x0000000000000000L >> (level * 2)),
				Level = level,
				Vertex1 = this.Vertex1,
				Vertex2 = v3,
				Vertex3 = v1
			};

			list[1] = new SphereGridTile()
			{
				QuadKey = this.QuadKey | (0x1000000000000000L >> (level * 2)),
				Level = level,
				Vertex1 = v1,
				Vertex2 = v2,
				Vertex3 = this.Vertex3
			};

			list[2] = new SphereGridTile()
			{
				QuadKey = this.QuadKey | (0x2000000000000000L >> (level * 2)),
				Level = level,
				Vertex1 = v3,
				Vertex2 = this.Vertex2,
				Vertex3 = v2
			};

			list[3] = new SphereGridTile()
			{
				QuadKey = this.QuadKey | (0x3000000000000000L >> (level * 2)),
				Level = level,
				Vertex1 = v1,
				Vertex2 = v2,
				Vertex3 = v3
			};

		}


		/// <summary>
		/// Checks if this tile is covered (fully or partially) by a circle with given center and radius represented as an angle.  
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="center">Center of the circle.</param>
		/// <param name="angle">Angle representing circle radius.</param>
		/// <returns>True if the circle covers this tile.</returns>
		public bool CoveredByCircle<T>(T center, double angle)
			where T : ICartesian
		{
			double circleCosine = Math.Cos(angle);
			double phi;

			// any vertex inside circle or circle center inside triangle

			if (SphericalExtension._Cosine(center.X, center.Y, center.Z, Vertex1.X, Vertex1.Y, Vertex1.Z) > circleCosine ||
				SphericalExtension._Cosine(center.X, center.Y, center.Z, Vertex2.X, Vertex2.Y, Vertex2.Z) > circleCosine ||
				SphericalExtension._Cosine(center.X, center.Y, center.Z, Vertex3.X, Vertex3.Y, Vertex3.Z) > circleCosine ||
				SphericalExtension._InsideTriangle(center.X, center.Y, center.Z, Vertex1, Vertex2, Vertex3)) return true;

			// check side Vertex1-Vertex2

			var (vx, vy, vz) = SphericalExtension._CrossProduct(Vertex1.X, Vertex1.Y, Vertex1.Z, Vertex2.X, Vertex2.Y, Vertex2.Z, false);
			var (cx, cy, cz) = SphericalExtension._CrossProduct(center.X, center.Y, center.Z, vx, vy, vz, false);

			if (SphericalExtension._DotProduct(cx, cy, cz, Vertex1.X, Vertex1.Y, Vertex1.Z) *
				SphericalExtension._DotProduct(cx, cy, cz, Vertex2.X, Vertex2.Y, Vertex2.Z) < 0.0)
			{
				phi = Math.Acos(SphericalExtension._Cosine(vx, vy, vz, center.X, center.Y, center.Z));
				if (phi > Math.PI / 2) phi = Math.PI - phi;
				if (angle + phi > Math.PI / 2) return true;
			}

			// check side Vertex2-Vertex3

			(vx, vy, vz) = SphericalExtension._CrossProduct(Vertex2.X, Vertex2.Y, Vertex2.Z, Vertex3.X, Vertex3.Y, Vertex3.Z, false);
			(cx, cy, cz) = SphericalExtension._CrossProduct(center.X, center.Y, center.Z, vx, vy, vz, false);

			if (SphericalExtension._DotProduct(cx, cy, cz, Vertex2.X, Vertex2.Y, Vertex2.Z) *
				SphericalExtension._DotProduct(cx, cy, cz, Vertex3.X, Vertex3.Y, Vertex3.Z) < 0.0)
			{
				phi = Math.Acos(SphericalExtension._Cosine(vx, vy, vz, center.X, center.Y, center.Z));
				if (phi > Math.PI / 2) phi = Math.PI - phi;
				if (angle + phi > Math.PI / 2) return true;
			}

			// check side Vertex3-Vertex1

			(vx, vy, vz) = SphericalExtension._CrossProduct(Vertex3.X, Vertex3.Y, Vertex3.Z, Vertex1.X, Vertex1.Y, Vertex1.Z, false);
			(cx, cy, cz) = SphericalExtension._CrossProduct(center.X, center.Y, center.Z, vx, vy, vz, false);

			if (SphericalExtension._DotProduct(cx, cy, cz, Vertex3.X, Vertex3.Y, Vertex3.Z) *
				SphericalExtension._DotProduct(cx, cy, cz, Vertex1.X, Vertex1.Y, Vertex1.Z) < 0.0)
			{
				phi = Math.Acos(SphericalExtension._Cosine(vx, vy, vz, center.X, center.Y, center.Z));
				if (phi > Math.PI / 2) phi = Math.PI - phi;
				if (angle + phi > Math.PI / 2) return true;
			}

			return false;
		}

		/// <summary>
		/// Checks if this tile is covered (fully or partially) by a polyline.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="polyline">Vertices of the polyline.</param>
		/// <param name="tolerance">Test tolerance in radians.</param>
		/// <returns>True if the polyline covers this tile.</returns>
		public bool CoveredByPolyline<T>(IEnumerable<T> polyline, double tolerance)
			where T : ICartesian
		{
			T vertex = polyline.First();

			foreach (T v in polyline.Skip(1))
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
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="polygon">Vertices of the polygon.</param>
		/// <returns>True if the polygon covers this tile.</returns>
		public bool CoveredByPolygon<T>(IEnumerable<T> polygon)
			where T : ICartesian
		{
			T vertex = polygon.Last();

			foreach (T v in polygon)
			{
				if (SphericalExtension.SectionsIntersect(Vertex1, Vertex2, vertex, v) >= 0 ||
					SphericalExtension.SectionsIntersect(Vertex2, Vertex3, vertex, v) >= 0 ||
					SphericalExtension.SectionsIntersect(Vertex3, Vertex1, vertex, v) >= 0)
					return true;

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
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="center">Center of the circle.</param>
		/// <param name="angle">Angle representing circle radius.</param>
		/// <returns>True if this tile fully encloses the circle.</returns>
		public bool EnclosesCircle<T>(T center, double angle)
			where T : ICartesian
		{
			double circleCosine = Math.Cos(angle);

			// any vertex outside circle or circle center outside triangle

			if (SphericalExtension._Cosine(center.X, center.Y, center.Z, Vertex1.X, Vertex1.Y, Vertex1.Z) > circleCosine ||
				SphericalExtension._Cosine(center.X, center.Y, center.Z, Vertex2.X, Vertex2.Y, Vertex2.Z) > circleCosine ||
				SphericalExtension._Cosine(center.X, center.Y, center.Z, Vertex3.X, Vertex3.Y, Vertex3.Z) > circleCosine ||
				!SphericalExtension._InsideTriangle(center.X, center.Y, center.Z, Vertex1, Vertex2, Vertex3))
				return false;

			// check side Vertex1-Vertex2

			var (vx, vy, vz) = SphericalExtension._CrossProduct(Vertex1.X, Vertex1.Y, Vertex1.Z, Vertex2.X, Vertex2.Y, Vertex2.Z, false);
			double phi = Math.Acos(SphericalExtension._Cosine(vx, vy, vz, center.X, center.Y, center.Z));
			if (phi > Math.PI / 2) phi = Math.PI - phi;
			if (angle + phi > Math.PI / 2) return false;

			// check side Vertex2-Vertex3
			(vx, vy, vz) = SphericalExtension._CrossProduct(Vertex2.X, Vertex2.Y, Vertex2.Z, Vertex3.X, Vertex3.Y, Vertex3.Z, false);
			phi = Math.Acos(SphericalExtension._Cosine(vx, vy, vz, center.X, center.Y, center.Z));
			if (phi > Math.PI / 2) phi = Math.PI - phi;
			if (angle + phi > Math.PI / 2) return false;

			// check side Vertex3-Vertex1
			(vx, vy, vz) = SphericalExtension._CrossProduct(Vertex3.X, Vertex3.Y, Vertex3.Z, Vertex1.X, Vertex1.Y, Vertex1.Z, false);
			phi = Math.Acos(SphericalExtension._Cosine(vx, vy, vz, center.X, center.Y, center.Z));
			if (phi > Math.PI / 2) phi = Math.PI - phi;
			return angle + phi <= Math.PI / 2; // <= because 'false' above
		}

		/// <summary>
		/// Checks if this tile fully encloses a polygon. Can be used for polylines too.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="polygon">Vertices of the polygon or polyline.</param>
		/// <returns>True if encloses.</returns>
		public bool EnclosesPolygon<T>(IEnumerable<T> polygon)
			where T : ICartesian
		{
			foreach (T vertex in polygon)
				if (!SphericalExtension._InsideTriangle(vertex.X, vertex.Y, vertex.Z, Vertex1, Vertex2, Vertex3))
					return false;

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
			if (sort)
				tiles.Sort();

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
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="center">Center of circle.</param>
		/// <param name="angle">Angle representing circle radius.</param>
		/// <param name="level">Maximum grid level of the tiles.</param>
		/// <param name="join">Join tiles if possible.</param>
		/// <returns>List of tiles covering circle.</returns>
		public static List<SphereGridTile> CoverCircleByTiles<T>(T center, double angle, int level, bool join)
			where T : ICartesian
		{
			if (level < 0 || level > SphereGridHelper.ABSOLUTE_MAX_LEVEL)
				throw new ArgumentOutOfRangeException(nameof(level));

			SphereGridTile tile = new SphereGridTile(center, 0);
			SphereGridTile t;
			List<SphereGridTile> result = new List<SphereGridTile>() { tile };
			int k = tile.LevelKey;

			if (!tile.EnclosesCircle(center, angle))
				for (int i = 0; i < SphereGridHelper.PrimaryTiles.Length; i++)
					if (i != k && (t = new SphereGridTile(SphereGridHelper.QuadKeyForPrimaryTile(i), 0)).CoveredByCircle(center, angle))
						result.Add(t);

			SphereGridTile[] tmp = new SphereGridTile[4];
			int n;

			for (int i = 0; i < level; i++)
			{
				n = result.Count;

				for (int j = 0; j < n; j++)
				{
					result[j].Split(tmp);

					if (tmp[0].CoveredByCircle(center, angle)) result.Add(tmp[0]);
					if (tmp[1].CoveredByCircle(center, angle)) result.Add(tmp[1]);
					if (tmp[2].CoveredByCircle(center, angle)) result.Add(tmp[2]);
					if (tmp[3].CoveredByCircle(center, angle)) result.Add(tmp[3]);
				}

				result.RemoveRange(0, n);
			}

			if (join && level > 0 && result.Count > 1)
				JoinTiles(result, level);

			return result;
		}

		/// <summary>
		/// Returns a list of tiles covering (fully or partially) a polyline. 
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="polyline">Vertices of the polyline.</param>
		/// <param name="level">Maximum grid level of tiles.</param>
		/// <param name="join">Join tiles if possible.</param>
		/// <returns>List of tiles covering the polyline.</returns>
		public static List<SphereGridTile> CoverPolylineByTiles<T>(IEnumerable<T> polyline, int level, double tolerance, bool join)
			where T : ICartesian
		{
			if (level < 0 || level > SphereGridHelper.ABSOLUTE_MAX_LEVEL)
				throw new ArgumentOutOfRangeException(nameof(level));

			SphereGridTile tile = new SphereGridTile(polyline.First(), 0);
			SphereGridTile t;
			List<SphereGridTile> result = new List<SphereGridTile>() { tile };
			int k = tile.LevelKey;

			for (int i = 0; i < SphereGridHelper.PrimaryTiles.Length; i++)
				if (i != k && (t = new SphereGridTile(SphereGridHelper.QuadKeyForPrimaryTile(i), 0)).CoveredByPolyline(polyline, tolerance))
					result.Add(t);

			SphereGridTile[] tmp = new SphereGridTile[4];
			int n;

			for (int i = 0; i < level; i++)
			{
				n = result.Count;

				for (int j = 0; j < n; j++)
				{
					result[j].Split(tmp);

					if (tmp[0].CoveredByPolyline(polyline, tolerance)) result.Add(tmp[0]);
					if (tmp[1].CoveredByPolyline(polyline, tolerance)) result.Add(tmp[1]);
					if (tmp[2].CoveredByPolyline(polyline, tolerance)) result.Add(tmp[2]);
					if (tmp[3].CoveredByPolyline(polyline, tolerance)) result.Add(tmp[3]);
				}

				result.RemoveRange(0, n);
			}

			if (join && level > 0 && result.Count > 1)
				JoinTiles(result, level);

			return result;
		}

		/// <summary>
		/// Returns a list of tiles covering (fully or partially) a polygon. 
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="polygon">Vertices of the polygon.</param>
		/// <param name="level">Maximum grid level of tiles.</param>
		/// <param name="join">Join tiles if possible.</param>
		/// <returns>List of tiles covering the polygon.</returns>
		public static List<SphereGridTile> CoverPolygonByTiles<T>(IEnumerable<T> polygon, int level, bool join)
			where T : ICartesian
		{
			if (level < 0 || level > SphereGridHelper.ABSOLUTE_MAX_LEVEL)
				throw new ArgumentOutOfRangeException(nameof(level));

			SphereGridTile tile = new SphereGridTile(polygon.First(), 0);
			SphereGridTile t;
			List<SphereGridTile> result = new List<SphereGridTile>() { tile };
			int k = tile.LevelKey;

			if (!tile.EnclosesPolygon(polygon))
				for (int i = 0; i < SphereGridHelper.PrimaryTiles.Length; i++)
					if (i != k && (t = new SphereGridTile(SphereGridHelper.QuadKeyForPrimaryTile(i), 0)).CoveredByPolygon(polygon))
						result.Add(t);

			SphereGridTile[] tmp = new SphereGridTile[4];
			int n;

			for (int i = 0; i < level; i++)
			{
				n = result.Count;

				for (int j = 0; j < n; j++)
				{
					result[j].Split(tmp);

					if (tmp[0].CoveredByPolygon(polygon)) result.Add(tmp[0]);
					if (tmp[1].CoveredByPolygon(polygon)) result.Add(tmp[1]);
					if (tmp[2].CoveredByPolygon(polygon)) result.Add(tmp[2]);
					if (tmp[3].CoveredByPolygon(polygon)) result.Add(tmp[3]);
				}

				result.RemoveRange(0, n);
			}

			if (join && level > 0 && result.Count > 1)
				JoinTiles(result, level);

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

		public override bool Equals(object obj)
		{
			if (obj == null || GetType() != obj.GetType())
				return false;

			return Equals((SphereGridTile)obj);
		}

		public override int GetHashCode()
		{
			return QuadKey.GetHashCode();
		}

		public override string ToString()
		{
#if NETCOREAPP2_1_OR_GREATER
			Span<int> keys = stackalloc int[SphereGridHelper.ABSOLUTE_MAX_LEVEL + 1];
#else
			int[] keys = new int[SphereGridHelper.ABSOLUTE_MAX_LEVEL + 1];
#endif

			SphereGridHelper.QuadKeyParts(QuadKey, keys);

			StringBuilder result = new StringBuilder(keys[0].ToString(), 128);

			for (int i = 1; i <= SphereGridHelper.ABSOLUTE_MAX_LEVEL; i++)
				result.Append(',').Append(keys[i]);

			return result.ToString();
		}
	}

}