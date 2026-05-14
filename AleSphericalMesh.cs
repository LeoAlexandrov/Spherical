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

namespace AleProjects.Spherical.Mesh
{
	public static class Mesh
	{
		// ABSOLUTE_MAX_LEVEL = 29 = ((64 - 1) - 5) / 2
		// "64" bits in long type,  
		// "1" is a sign bit, 
		// "5" five bits for 0..19 value indexing icosahedron face which is 1/20 of a sphere, 
		// "2" two bits to represent 0..3 value which is index of sub-tile inside parent face.

		/// <summary>
		/// Maximum subdivision level of the sphere mesh.
		/// </summary>
		public const int ABSOLUTE_MAX_LEVEL = ((64 - 1) - 5) / 2;


		#region Icosahedron-vertices

		private const double _Ring_Radius = 0.8944271909999159; // 2.0 / Math.Sqrt(5.0)

		public static readonly CartesianValue AxisPos = new CartesianValue(0.0, 0.0, 1.0);

		public static readonly CartesianValue Ring1_0 = new CartesianValue(
			_Ring_Radius * Math.Cos(0.0 * 2.0 * Math.PI / 5.0),
			_Ring_Radius * Math.Sin(0.0 * 2.0 * Math.PI / 5.0),
			1.0 / Math.Sqrt(5.0));

		public static readonly CartesianValue Ring1_1 = new CartesianValue(
			_Ring_Radius * Math.Cos(1.0 * 2.0 * Math.PI / 5.0),
			_Ring_Radius * Math.Sin(1.0 * 2.0 * Math.PI / 5.0),
			1.0 / Math.Sqrt(5.0));

		public static readonly CartesianValue Ring1_2 = new CartesianValue(
			_Ring_Radius * Math.Cos(2.0 * 2.0 * Math.PI / 5.0),
			_Ring_Radius * Math.Sin(2.0 * 2.0 * Math.PI / 5.0),
			1.0 / Math.Sqrt(5.0));

		public static readonly CartesianValue Ring1_3 = new CartesianValue(
			_Ring_Radius * Math.Cos(3.0 * 2.0 * Math.PI / 5.0),
			_Ring_Radius * Math.Sin(3.0 * 2.0 * Math.PI / 5.0),
			1.0 / Math.Sqrt(5.0));

		public static readonly CartesianValue Ring1_4 = new CartesianValue(
			_Ring_Radius * Math.Cos(4.0 * 2.0 * Math.PI / 5.0),
			_Ring_Radius * Math.Sin(4.0 * 2.0 * Math.PI / 5.0),
			1.0 / Math.Sqrt(5.0));

		public static readonly CartesianValue Ring2_0 = new CartesianValue(
			_Ring_Radius * Math.Cos(0.0 * 2.0 * Math.PI / 5.0 + Math.PI / 5.0),
			_Ring_Radius * Math.Sin(0.0 * 2.0 * Math.PI / 5.0 + Math.PI / 5.0),
			-1.0 / Math.Sqrt(5.0));

		public static readonly CartesianValue Ring2_1 = new CartesianValue(
			_Ring_Radius * Math.Cos(1.0 * 2.0 * Math.PI / 5.0 + Math.PI / 5.0),
			_Ring_Radius * Math.Sin(1.0 * 2.0 * Math.PI / 5.0 + Math.PI / 5.0),
			-1.0 / Math.Sqrt(5.0));

		public static readonly CartesianValue Ring2_2 = new CartesianValue(
			_Ring_Radius * Math.Cos(2.0 * 2.0 * Math.PI / 5.0 + Math.PI / 5.0),
			_Ring_Radius * Math.Sin(2.0 * 2.0 * Math.PI / 5.0 + Math.PI / 5.0),
			-1.0 / Math.Sqrt(5.0));

		public static readonly CartesianValue Ring2_3 = new CartesianValue(
			_Ring_Radius * Math.Cos(3.0 * 2.0 * Math.PI / 5.0 + Math.PI / 5.0),
			_Ring_Radius * Math.Sin(3.0 * 2.0 * Math.PI / 5.0 + Math.PI / 5.0),
			-1.0 / Math.Sqrt(5.0));

		public static readonly CartesianValue Ring2_4 = new CartesianValue(
			_Ring_Radius * Math.Cos(4.0 * 2.0 * Math.PI / 5.0 + Math.PI / 5.0),
			_Ring_Radius * Math.Sin(4.0 * 2.0 * Math.PI / 5.0 + Math.PI / 5.0),
			-1.0 / Math.Sqrt(5.0));

		public static readonly CartesianValue AxisNeg = new CartesianValue(0.0, 0.0, -1.0);

		#endregion

		#region Icosahedron-faces

		public static readonly CartesianValue[][] IcosahedronFaces =
		{
			new[] { AxisPos, Ring1_0, Ring1_1 },
			new[] { AxisPos, Ring1_1, Ring1_2 },
			new[] { AxisPos, Ring1_2, Ring1_3 },
			new[] { AxisPos, Ring1_3, Ring1_4 },
			new[] { AxisPos, Ring1_4, Ring1_0 },
		
			new[] { Ring1_0, Ring1_1, Ring2_0 },
			new[] { Ring1_1, Ring2_0, Ring2_1 },
			new[] { Ring1_1, Ring1_2, Ring2_1 },
			new[] { Ring1_2, Ring2_1, Ring2_2 },
			new[] { Ring1_2, Ring1_3, Ring2_2 },
			new[] { Ring1_3, Ring2_2, Ring2_3 },
			new[] { Ring1_3, Ring1_4, Ring2_3 },
			new[] { Ring1_4, Ring2_3, Ring2_4 },
			new[] { Ring1_4, Ring1_0, Ring2_4 },
			new[] { Ring1_0, Ring2_4, Ring2_0 },
		
			new[] { AxisNeg, Ring2_0, Ring2_4 },
			new[] { AxisNeg, Ring2_4, Ring2_3 },
			new[] { AxisNeg, Ring2_3, Ring2_2 },
			new[] { AxisNeg, Ring2_2, Ring2_1 },
			new[] { AxisNeg, Ring2_1, Ring2_0 }
		};

		#endregion

		/// <summary>
		/// Finds an icosahedron face containing the specified location on a sphere.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="location">Location on the sphere.</param>
		/// <returns>Index of the icosahedron face.</returns>
		public static int FindPrimaryCell<T>(T location)
			where T : ICartesian
		{
			for (int i = 0; i < IcosahedronFaces.Length; i++)
				if (SphericalExtension._InsideTriangle(location.X, location.Y, location.Z, IcosahedronFaces[i][0], IcosahedronFaces[i][1], IcosahedronFaces[i][2]))
					return i;

			throw new InvalidOperationException();
		}

		/// <summary>
		/// Finds an icosahedron face containing the specified location on a sphere.
		/// </summary>
		/// <param name="latitude">Latitude of the location in degrees.</param>
		/// <param name="longitude">Longitude of the location in degrees.</param>
		/// <returns>Index of the icosahedron face.</returns>
		public static int FindPrimaryCell(double latitude, double longitude)
		{
			var p = new CartesianValue(latitude, longitude);

			return FindPrimaryCell(p);
		}


		/// <summary>
		/// Calculates the key of the mesh cell at the specified subdivision level containing the given location.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="location">Location on the sphere.</param>
		/// <param name="level">Mesh subdivision level.</param>
		/// <returns>Key uniquely identifying mesh cell containing given location.</returns>
		public static long BuildKey<T>(T location, int level)
			where T : ICartesian
		{
			if (level < 0 || level > ABSOLUTE_MAX_LEVEL)
				throw new ArgumentOutOfRangeException(nameof(level));

			double x = location.X;
			double y = location.Y;
			double z = location.Z;

			var vertex1 = new CartesianValue();
			var vertex2 = new CartesianValue();
			var vertex3 = new CartesianValue();
			var v1 = new CartesianValue();
			var v2 = new CartesianValue();
			var v3 = new CartesianValue();

			long result = 0;

			for (int i = 0; i < IcosahedronFaces.Length; i++)
				if (SphericalExtension._InsideTriangle(x, y, z, IcosahedronFaces[i][0], IcosahedronFaces[i][1], IcosahedronFaces[i][2]))
				{
					vertex1.SetCartesian(IcosahedronFaces[i][0].X, IcosahedronFaces[i][0].Y, IcosahedronFaces[i][0].Z);
					vertex2.SetCartesian(IcosahedronFaces[i][1].X, IcosahedronFaces[i][1].Y, IcosahedronFaces[i][1].Z);
					vertex3.SetCartesian(IcosahedronFaces[i][2].X, IcosahedronFaces[i][2].Y, IcosahedronFaces[i][2].Z);
					result = (long)i << (64 - 1 - 5);
					break;
				}

			for (int i = 0; i < level; i++)
			{
				v1.SetNormalized(vertex1.X + vertex2.X, vertex1.Y + vertex2.Y, vertex1.Z + vertex2.Z);
				v2.SetNormalized(vertex2.X + vertex3.X, vertex2.Y + vertex3.Y, vertex2.Z + vertex3.Z);
				v3.SetNormalized(vertex3.X + vertex1.X, vertex3.Y + vertex1.Y, vertex3.Z + vertex1.Z);

				if (SphericalExtension._InsideTriangle(x, y, z, vertex1, v1, v3))
				{
					vertex2 = v1;
					vertex3 = v3;
					result |= (long)0x0000_0000_0000_0000 >> (i * 2);
					continue;
				}

				if (SphericalExtension._InsideTriangle(x, y, z, vertex2, v2, v1))
				{
					vertex1 = vertex2;
					vertex2 = v2;
					vertex3 = v1;
					result |= 0x0100_0000_0000_0000 >> (i * 2);
					continue;
				}

				if (SphericalExtension._InsideTriangle(x, y, z, vertex3, v3, v2))
				{
					vertex1 = vertex3;
					vertex2 = v3;
					vertex3 = v2;
					result |= 0x0200_0000_0000_0000 >> (i * 2);
					continue;
				}

				vertex1 = v1;
				vertex2 = v2;
				vertex3 = v3;
				result |= 0x0300_0000_0000_0000 >> (i * 2);
			}

			return result;
		}

		/// <summary>
		/// Calculates the key of the mesh cell at the specified subdivision level containing a location with the given latitude and longitude.
		/// </summary>
		/// <param name="latitude">Latitude of the location in degrees.</param>
		/// <param name="longitude">Longitude of the location in degrees.</param>
		/// <param name="level">Mesh subdivision level.</param>
		/// <returns>Key uniquely identifying mesh cell containing given location.</returns>
		public static long BuildKey(double latitude, double longitude, int level)
		{
			var location = new CartesianValue(latitude, longitude, true);
			return BuildKey(location, level);
		}

		/// <summary>
		/// Returns a bit mask for the key of a mesh cell at the specified subdivision level.
		/// </summary>
		/// <param name="level">Mesh subdivision level.</param>
		/// <returns>Bit mask.</returns>
		public static long KeyMask(int level)
		{
			long result = 0x7C_00_00_00_00_00_00_00;
			long mask = 0x03_00_00_00_00_00_00_00;

			while (level-- > 0)
			{
				result |= mask;
				mask >>= 2;
			}

			return result;
		}

		/// <summary>
		/// Returns a key for an icosahedron face by its index.
		/// </summary>
		/// <param name="index">Index of the icosahedron face.</param>
		/// <returns>Key for the icosahedron face.</returns>
		public static long KeyForPrimaryFace(int index)
		{
			return (long)(index & 31) << (64 - 5 - 1);
		}

		/// <summary>
		/// Decomposes a key into cell indices at each mesh subdivision level.
		/// </summary>
		/// <param name="key">Key to decompose.</param>
		/// <returns>Decomposed key.</returns>
#if NETCOREAPP2_1_OR_GREATER || NETSTANDARD2_1_OR_GREATER
		public static void KeyParts(long key, Span<int> keys)
#else
		public static void KeyParts(long key, int[] keys)
#endif
		{
			for (int i = ABSOLUTE_MAX_LEVEL; i > 0; i--)
			{
				keys[i] = (int)(key & 3);
				key >>= 2;
			}

			keys[0] = (int)(key & 31);
		}

	}



	/// <summary>
	/// Represents a triangular cell of the spherical mesh.
	/// At level 0, there are 20 cells corresponding to the icosahedron faces.
	/// Each subsequent level is obtained by splitting a cell into 4 sub-cells by connecting the midpoints of its sides.
	/// </summary>
	public struct MeshCell : IComparable<MeshCell>, IEquatable<MeshCell>
	{
		/// <summary>
		/// The unique identifier of this cell.
		/// </summary>
		public long Key { get; private set; }

		/// <summary>
		/// The mesh subdivision level for this cell.
		/// </summary>
		public int Level { get; private set; }

		/// <summary>
		/// The first vertex of the cell triangle.
		/// </summary>
		public CartesianValue Vertex1 { get; private set; }

		/// <summary>
		/// The second vertex of the cell triangle.
		/// </summary>
		public CartesianValue Vertex2 { get; private set; }

		/// <summary>
		/// The third vertex of the cell triangle.
		/// </summary>
		public CartesianValue Vertex3 { get; private set; }

		/// <summary>
		/// The index of this cell inside its parent cell at the previous level.
		/// </summary>
#if NETCOREAPP2_1_OR_GREATER || NETSTANDARD2_1_OR_GREATER
		public readonly int Index
#else
		public int Index
#endif
		{
			get
			{
				long result = Key >> ((Mesh.ABSOLUTE_MAX_LEVEL - Level) * 2);

				return Level == 0 ? (int)(result & 0x1f) : (int)(result & 0x03);
			}
		}

		/// <summary>
		///The key of the upper parent cell at the previous level.
		/// </summary>
#if NETCOREAPP2_1_OR_GREATER || NETSTANDARD2_1_OR_GREATER
		public readonly long ParentKey
#else
		public long ParentKey
#endif
		{
			get => Level > 0 ? Key & (long)(0x7f_ff_ff_ff_ff_ff_ff_ff & (0xff_ff_ff_ff_ff_ff_ff_ff << ((Mesh.ABSOLUTE_MAX_LEVEL - Level + 1) * 2))) : 0;
		}

		/// <summary>
		/// The upper value of the Key. Allows to select objects inside this cell using filter &gt;= Key and &lt;= KeyUpperValue. 
		/// </summary>
#if NETCOREAPP2_1_OR_GREATER || NETSTANDARD2_1_OR_GREATER
		public readonly long KeyUpperValue
#else
		public long KeyUpperValue
#endif
		{
			get => Key | (0x03_ff_ff_ff_ff_ff_ff_ff >> (Level * 2));
		}


		/// <summary>
		/// Initializes a new instance of the MeshCell that represents the cell containing
		/// the specified location at the given subdivision level.
		/// </summary>
		/// <remarks>The mesh cell is determined by locating the icosahedron face containing the specified
		/// location and recursively subdividing it to the requested level.</remarks>
		/// <param name="location">The location for which to find the containing mesh cell. Cannot be null.</param>
		/// <param name="level">The subdivision level. Must be between 0 and SphereMesh.ABSOLUTE_MAX_LEVEL, inclusive.</param>
		/// <exception cref="ArgumentNullException">Thrown if location is null.</exception>
		/// <exception cref="ArgumentOutOfRangeException">Thrown if the level is less than 0 or greater than SphereMesh.ABSOLUTE_MAX_LEVEL.</exception>
		public MeshCell(ICartesian location, int level)
		{
			if (location == null)
				throw new ArgumentNullException(nameof(location));

			if (level < 0 || level > Mesh.ABSOLUTE_MAX_LEVEL)
				throw new ArgumentOutOfRangeException(nameof(level));

			CartesianValue locationValue = new CartesianValue(location);
			int k = Mesh.FindPrimaryCell(locationValue);
			long key = (long)k << (64 - 1 - 5); // 1 - sign bit, 5 - bits for 0..19 values

			Level = level;

			if (level > 0)
			{
				MeshCell cell = new MeshCell()
				{
					Key = key,
					Vertex1 = new CartesianValue(Mesh.IcosahedronFaces[k][0]),
					Vertex2 = new CartesianValue(Mesh.IcosahedronFaces[k][1]),
					Vertex3 = new CartesianValue(Mesh.IcosahedronFaces[k][2]),
				};


				for (int i = 0; i < level; i++)
					cell = cell.SplitAndFind(locationValue);

				Key = cell.Key;
				Vertex1 = cell.Vertex1;
				Vertex2 = cell.Vertex2;
				Vertex3 = cell.Vertex3;
			}
			else
			{
				Key = key;
				Vertex1 = new CartesianValue(Mesh.IcosahedronFaces[k][0]);
				Vertex2 = new CartesianValue(Mesh.IcosahedronFaces[k][1]);
				Vertex3 = new CartesianValue(Mesh.IcosahedronFaces[k][2]);
			}
		}

		/// <summary>
		/// Initializes a new instance of the MeshCell class for the specified cell key and subdivision level.
		/// </summary>
		/// <remarks>This constructor determines the vertices and key for the mesh cell at the specified level. For
		/// level 0, the cell corresponds to a face of the icosahedron. For higher levels, the cell is derived by
		/// subdividing the base face according to the key and level provided.</remarks>
		/// <param name="key">The unique key identifying the mesh cell. Must be non-negative.</param>
		/// <param name="level">The subdivision level of the mesh cell. Must be between 0 and SphereMesh.ABSOLUTE_MAX_LEVEL, inclusive.</param>
		/// <exception cref="ArgumentOutOfRangeException">Thrown if key is negative, or if level is less than 0 or greater than SphereMesh.ABSOLUTE_MAX_LEVEL.</exception>
		public MeshCell(long key, int level)
		{
			if (key < 0)
				throw new ArgumentOutOfRangeException(nameof(key));

			if (level < 0 || level > Mesh.ABSOLUTE_MAX_LEVEL)
				throw new ArgumentOutOfRangeException(nameof(level));

			if (level > 0)
			{
#if NETCOREAPP2_1_OR_GREATER || NETSTANDARD2_1_OR_GREATER
				Span<int> parts = stackalloc int[Mesh.ABSOLUTE_MAX_LEVEL + 1];
#else
				int[] parts = new int[Mesh.ABSOLUTE_MAX_LEVEL + 1];
#endif
				Mesh.KeyParts(key, parts);
				int k = parts[0];

				MeshCell cell = new MeshCell()
				{
					Key = key & 0x7c00_0000_0000_0000,
					Level = 0,
					Vertex1 = new CartesianValue(Mesh.IcosahedronFaces[k][0]),
					Vertex2 = new CartesianValue(Mesh.IcosahedronFaces[k][1]),
					Vertex3 = new CartesianValue(Mesh.IcosahedronFaces[k][2]),
				};

#if NETCOREAPP2_1_OR_GREATER || NETSTANDARD2_1_OR_GREATER
				Span<MeshCell> cells = stackalloc MeshCell[4];
#else
				MeshCell[] cells = new MeshCell[4];
#endif

				for (int i = 1; i <= level; i++)
				{
					cell.Split(cells);
					cell = cells[parts[i]];
				}

				Key = cell.Key;
				Level = cell.Level;
				Vertex1 = cell.Vertex1;
				Vertex2 = cell.Vertex2;
				Vertex3 = cell.Vertex3;
			}
			else
			{
				int k = (int)(key >> (64 - 1 - 5)); // 1 - sign bit, 5 - bits for 0..19 values

				Key = key;
				Level = 0;
				Vertex1 = new CartesianValue(Mesh.IcosahedronFaces[k][0]);
				Vertex2 = new CartesianValue(Mesh.IcosahedronFaces[k][1]);
				Vertex3 = new CartesianValue(Mesh.IcosahedronFaces[k][2]);
			}
		}

		/// <summary>
		/// Splits the current mesh cell into subcells and returns the subcell that contains the specified location.
		/// </summary>
		/// <remarks>This method is typically used to refine a mesh by subdividing a cell and locating the subcell
		/// that contains a given location. The returned MeshCell will have an incremented level compared to the original
		/// cell.</remarks>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="location">The Cartesian coordinate used to determine which subcell to return.</param>
		/// <returns>A MeshCell representing the subcell that contains the specified location.</returns>
		public MeshCell SplitAndFind<T>(T location)
			where T : ICartesian
		{
			var (x, y, z) = SphericalExtension._Normalized(Vertex1.X + Vertex2.X, Vertex1.Y + Vertex2.Y, Vertex1.Z + Vertex2.Z);
			CartesianValue v1 = new CartesianValue(x, y, z);
			(x, y, z) = SphericalExtension._Normalized(Vertex2.X + Vertex3.X, Vertex2.Y + Vertex3.Y, Vertex2.Z + Vertex3.Z);
			CartesianValue v2 = new CartesianValue(x, y, z);
			(x, y, z) = SphericalExtension._Normalized(Vertex3.X + Vertex1.X, Vertex3.Y + Vertex1.Y, Vertex3.Z + Vertex1.Z);
			CartesianValue v3 = new CartesianValue(x, y, z);

			x = location.X;
			y = location.Y;
			z = location.Z;

			int shift = this.Level * 2;
			int level = this.Level + 1;

			if (SphericalExtension._InsideTriangle(x, y, z, this.Vertex1, v1, v3))
				return new MeshCell()
				{
					Key = this.Key | ((long)0x0000_0000_0000_0000 >> shift),
					Level = level,
					Vertex1 = this.Vertex1,
					Vertex2 = v1,
					Vertex3 = v3
				};

			if (SphericalExtension._InsideTriangle(x, y, z, this.Vertex2, v2, v1))
				return new MeshCell()
				{
					Key = this.Key | (0x0100_0000_0000_0000 >> shift),
					Level = level,
					Vertex1 = this.Vertex2,
					Vertex2 = v2,
					Vertex3 = v1
				};

			if (SphericalExtension._InsideTriangle(x, y, z, this.Vertex3, v3, v2))
				return new MeshCell()
				{
					Key = this.Key | (0x0200_0000_0000_0000 >> shift),
					Level = level,
					Vertex1 = this.Vertex3,
					Vertex2 = v3,
					Vertex3 = v2,
				};

			return new MeshCell()
			{
				Key = this.Key | (0x0300_0000_0000_0000 >> shift),
				Level = level,
				Vertex1 = v1,
				Vertex2 = v2,
				Vertex3 = v3,
			};
		}

		/// <summary>
		/// Splits the current mesh cell into four smaller mesh cells and writes them to the specified buffer.
		/// </summary>
		/// <remarks>This method performs a subdivision of the mesh cell, typically used for mesh refinement or
		/// tessellation. The input list is overwritten with the four resulting child cells, each representing a subdivision
		/// of the original cell. The order of the output cells is consistent and deterministic.</remarks>
		/// <param name="list">A span or array of four <see cref="MeshCell"/> instances to receive the resulting subdivided mesh cells. The list
		/// must have a length of at least four.</param>
#if NETCOREAPP2_1_OR_GREATER || NETSTANDARD2_1_OR_GREATER
		public void Split(Span<MeshCell> list)
#else
		public void Split(MeshCell[] list)
#endif
		{
			var (x, y, z) = SphericalExtension._Normalized(Vertex1.X + Vertex2.X, Vertex1.Y + Vertex2.Y, Vertex1.Z + Vertex2.Z);
			CartesianValue v1 = new CartesianValue(x, y, z);
			(x, y, z) = SphericalExtension._Normalized(Vertex2.X + Vertex3.X, Vertex2.Y + Vertex3.Y, Vertex2.Z + Vertex3.Z);
			CartesianValue v2 = new CartesianValue(x, y, z);
			(x, y, z) = SphericalExtension._Normalized(Vertex3.X + Vertex1.X, Vertex3.Y + Vertex1.Y, Vertex3.Z + Vertex1.Z);
			CartesianValue v3 = new CartesianValue(x, y, z);

			int shift = this.Level * 2;
			int level = this.Level + 1;

			list[0] = new MeshCell()
			{
				Key = this.Key | ((long)0x0000_0000_0000_0000 >> shift),
				Level = level,
				Vertex1 = this.Vertex1,
				Vertex2 = v1,
				Vertex3 = v3
			};

			list[1] = new MeshCell()
			{
				Key = this.Key | (0x0100_0000_0000_0000 >> shift),
				Level = level,
				Vertex1 = this.Vertex2,
				Vertex2 = v2,
				Vertex3 = v1
			};

			list[2] = new MeshCell()
			{
				Key = this.Key | (0x0200_0000_0000_0000 >> shift),
				Level = level,
				Vertex1 = this.Vertex3,
				Vertex2 = v3,
				Vertex3 = v2,
			};

			list[3] = new MeshCell()
			{
				Key = this.Key | (0x0300_0000_0000_0000 >> shift),
				Level = level,
				Vertex1 = v1,
				Vertex2 = v2,
				Vertex3 = v3,
			};
		}

		/// <summary>
		/// Checks if this cell is covered (fully or partially) by a circle with the given center and radius represented as an angle.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="center">Center of the circle.</param>
		/// <param name="angle">Angle representing circle radius.</param>
		/// <returns>True if the circle covers this cell.</returns>
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
		/// Checks if this cell is covered (fully or partially) by a polyline.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="polyline">Vertices of the polyline.</param>
		/// <param name="tolerance">Test tolerance in radians.</param>
		/// <returns>True if the polyline covers this cell.</returns>
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
		/// Checks if this cell is covered (fully or partially) by a polygon.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="polygon">Vertices of the polygon.</param>
		/// <returns>True if the polygon covers this cell.</returns>
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
		/// Checks if this cell fully encloses a circle with the given center and radius represented as an angle.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="center">Center of the circle.</param>
		/// <param name="angle">Angle representing circle radius.</param>
		/// <returns>True if this cell fully encloses the circle.</returns>
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
		/// Checks if this cell fully encloses a polygon. Can be used for polylines too.
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
		/// Joins cells of the same subdivision level to the upper possible level.
		/// </summary>
		/// <param name="cells">List of cells to join ordered by their keys.</param>
		/// <param name="sort">Sort if necessary.</param>
		public static void Join(List<MeshCell> cells, bool sort = true)
		{
			if (cells == null)
				return;

			if (sort)
				cells.Sort();

			int level = cells.Last().Level;

			if (level == 0 || cells.Count < 3)
				return;

			bool hasJoined = true;
			int k = 0;

			for (int i = level; i > 0 && hasJoined; i--)
			{
				hasJoined = false;

				for (int j = 0; j < cells.Count; j++)
					if (cells[j].Level == i &&
						j < cells.Count - 1 &&
						cells[j + 1].Level == i &&
						cells[j + 1].ParentKey == cells[j].ParentKey)
					{
						k++;

						if (k == 3)
						{
							j -= 2;

							cells[j] = new MeshCell()
							{
								Level = i - 1,
								Key = cells[j].ParentKey,
								Vertex1 = cells[j].Vertex1,
								Vertex2 = cells[j + 1].Vertex1,
								Vertex3 = cells[j + 2].Vertex1
							};

							cells.RemoveRange(j + 1, 3);
							k = 0;
							hasJoined = true;
						}
					}
					else k = 0;
			}
		}


		/// <summary>
		/// Fills the given list with cells at the given subdivision level covering (fully or partially) a circle with the given center and radius represented as an angle.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="center">Center of circle.</param>
		/// <param name="angle">Angle representing circle radius.</param>
		/// <param name="level">Mesh subdivision level.</param>
		/// <param name="join">Join cells if possible.</param>
		/// <param name="result">List of cells.</param>
		public static void CoverCircle<T>(T center, double angle, int level, bool join, List<MeshCell> result)
			where T : ICartesian
		{
			if (result == null)
				throw new ArgumentNullException(nameof(result));

			if (level < 0 || level > Mesh.ABSOLUTE_MAX_LEVEL)
				throw new ArgumentOutOfRangeException(nameof(level));

			if (result.Count != 0)
				result.Clear();

			MeshCell cell = new MeshCell(center, 0);
			MeshCell t;
			int k = cell.Index;

			result.Add(cell);

			if (!cell.EnclosesCircle(center, angle))
				for (int i = 0; i < Mesh.IcosahedronFaces.Length; i++)
					if (i != k && (t = new MeshCell(Mesh.KeyForPrimaryFace(i), 0)).CoveredByCircle(center, angle))
						result.Add(t);

#if NETCOREAPP2_1_OR_GREATER || NETSTANDARD2_1_OR_GREATER
			Span<MeshCell> splitTo = stackalloc MeshCell[4];
#else
			MeshCell[] splitTo = new MeshCell[4];
#endif

			int n;
			List<MeshCell> tmp = new List<MeshCell>(8);

			for (int i = 0; i < level; i++)
			{
				n = result.Count;

				for (int j = 0; j < n; j++)
				{
					result[j].Split(splitTo);

					if (splitTo[0].CoveredByCircle(center, angle)) tmp.Add(splitTo[0]);
					if (splitTo[1].CoveredByCircle(center, angle)) tmp.Add(splitTo[1]);
					if (splitTo[2].CoveredByCircle(center, angle)) tmp.Add(splitTo[2]);
					if (splitTo[3].CoveredByCircle(center, angle)) tmp.Add(splitTo[3]);
				}

				(result, tmp) = (tmp, result);
				tmp.Clear();
			}

			if (level % 2 != 0)
			{
				tmp.AddRange(result);
				result = tmp;
			}

			if (join)
				Join(result);
		}

		/// <summary>
		/// Returns a list of cells at the given subdivision level covering (fully or partially) a circle with the given center and radius represented as an angle. 
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="center">Center of circle.</param>
		/// <param name="angle">Angle representing circle radius.</param>
		/// <param name="level">Mesh subdivision level.</param>
		/// <param name="join">Join cells if possible.</param>
		/// <returns>List of cells covering circle.</returns>
		public static List<MeshCell> CoverCircle<T>(T center, double angle, int level, bool join)
			where T : ICartesian
		{
			var result = new List<MeshCell>();

			CoverCircle(center, angle, level, join, result);

			return result;
		}

		/// <summary>
		/// Fills the given list with cells at the given subdivision level covering (fully or partially) a polyline.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="polyline">Vertices of the polyline.</param>
		/// <param name="level">Mesh subdivision level.</param>
		/// <param name="join">Join cells if possible.</param>
		/// <param name="result">List of cells.</param>
		public static void CoverPolyline<T>(IEnumerable<T> polyline, int level, double tolerance, bool join, List<MeshCell> result)
			where T : ICartesian
		{
			if (result == null)
				throw new ArgumentNullException(nameof(result));

			if (level < 0 || level > Mesh.ABSOLUTE_MAX_LEVEL)
				throw new ArgumentOutOfRangeException(nameof(level));

			if (result.Count != 0)
				result.Clear();

			MeshCell cell = new MeshCell(polyline.First(), 0);
			MeshCell c;
			int k = cell.Index;

			result.Add(cell);

			for (int i = 0; i < Mesh.IcosahedronFaces.Length; i++)
				if (i != k && (c = new MeshCell(Mesh.KeyForPrimaryFace(i), 0)).CoveredByPolyline(polyline, tolerance))
					result.Add(c);

#if NETCOREAPP2_1_OR_GREATER || NETSTANDARD2_1_OR_GREATER
			Span<MeshCell> splitTo = stackalloc MeshCell[4];
#else
			MeshCell[] splitTo = new MeshCell[4];
#endif

			int n;
			List<MeshCell> tmp = new List<MeshCell>(8);

			for (int i = 0; i < level; i++)
			{
				n = result.Count;

				for (int j = 0; j < n; j++)
				{
					result[j].Split(splitTo);

					if (splitTo[0].CoveredByPolyline(polyline, tolerance)) tmp.Add(splitTo[0]);
					if (splitTo[1].CoveredByPolyline(polyline, tolerance)) tmp.Add(splitTo[1]);
					if (splitTo[2].CoveredByPolyline(polyline, tolerance)) tmp.Add(splitTo[2]);
					if (splitTo[3].CoveredByPolyline(polyline, tolerance)) tmp.Add(splitTo[3]);
				}

				(result, tmp) = (tmp, result);
				tmp.Clear();
			}

			if (level % 2 != 0)
			{
				tmp.AddRange(result);
				result = tmp;
			}

			if (join)
				Join(result);
		}

		/// <summary>
		/// Returns a list of cells at the given subdivision level covering (fully or partially) a polyline. 
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="polyline">Vertices of the polyline.</param>
		/// <param name="level">Mesh subdivision level.</param>
		/// <param name="join">Join cells if possible.</param>
		/// <returns>List of cells covering the polyline.</returns>
		public static List<MeshCell> CoverPolyline<T>(IEnumerable<T> polyline, int level, double tolerance, bool join)
			where T : ICartesian
		{
			var result = new List<MeshCell>();

			CoverPolyline(polyline, level, tolerance, join, result);

			return result;
		}

		/// <summary>
		/// Fills the given list with cells at the given subdivision level covering (fully or partially) a polygon.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="polygon">Vertices of the polygon.</param>
		/// <param name="level">Mesh subdivision level.</param>
		/// <param name="join">Join cells if possible.</param>
		/// <param name="result">List of cells.</param>
		public static void CoverPolygon<T>(IEnumerable<T> polygon, int level, bool join, List<MeshCell> result)
			where T : ICartesian
		{
			if (result == null)
				throw new ArgumentNullException(nameof(result));

			if (level < 0 || level > Mesh.ABSOLUTE_MAX_LEVEL)
				throw new ArgumentOutOfRangeException(nameof(level));

			if (result.Count != 0)
				result.Clear();

			MeshCell cell = new MeshCell(polygon.First(), 0);
			MeshCell t;
			int k = cell.Index;

			result.Add(cell);

			if (!cell.EnclosesPolygon(polygon))
				for (int i = 0; i < Mesh.IcosahedronFaces.Length; i++)
					if (i != k && (t = new MeshCell(Mesh.KeyForPrimaryFace(i), 0)).CoveredByPolygon(polygon))
						result.Add(t);

#if NETCOREAPP2_1_OR_GREATER || NETSTANDARD2_1_OR_GREATER
			Span<MeshCell> splitTo = stackalloc MeshCell[4];
#else
			MeshCell[] splitTo = new MeshCell[4];
#endif

			int n;
			List<MeshCell> tmp = new List<MeshCell>(8);

			for (int i = 0; i < level; i++)
			{
				n = result.Count;

				for (int j = 0; j < n; j++)
				{
					result[j].Split(splitTo);

					if (splitTo[0].CoveredByPolygon(polygon)) tmp.Add(splitTo[0]);
					if (splitTo[1].CoveredByPolygon(polygon)) tmp.Add(splitTo[1]);
					if (splitTo[2].CoveredByPolygon(polygon)) tmp.Add(splitTo[2]);
					if (splitTo[3].CoveredByPolygon(polygon)) tmp.Add(splitTo[3]);
				}

				(result, tmp) = (tmp, result);
				tmp.Clear();
			}

			if (level % 2 != 0)
			{
				tmp.AddRange(result);
				result = tmp;
			}

			if (join)
				Join(result);
		}

		/// <summary>
		/// Returns a list of cells at the given subdivision level covering (fully or partially) a polygon. 
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="polygon">Vertices of the polygon.</param>
		/// <param name="level">Mesh subdivision level.</param>
		/// <param name="join">Join cells if possible.</param>
		/// <returns>List of cells covering the polygon.</returns>
		public static List<MeshCell> CoverPolygon<T>(IEnumerable<T> polygon, int level, bool join)
			where T : ICartesian
		{
			var result = new List<MeshCell>();

			CoverPolygon(polygon, level, join, result);

			return result;
		}


		public int CompareTo(MeshCell other)
		{
			return Key.CompareTo(other.Key);
		}

		public bool Equals(MeshCell other)
		{
			return Key.Equals(other.Key);
		}

		public override bool Equals(object obj)
		{
			if (obj == null || GetType() != obj.GetType())
				return false;

			return Equals((MeshCell)obj);
		}

		public override int GetHashCode()
		{
			return Key.GetHashCode();
		}

		public override string ToString()
		{
#if NETCOREAPP2_1_OR_GREATER || NETSTANDARD2_1_OR_GREATER
			Span<int> keys = stackalloc int[Mesh.ABSOLUTE_MAX_LEVEL + 1];
#else
			int[] keys = new int[Mesh.ABSOLUTE_MAX_LEVEL + 1];
#endif

			Mesh.KeyParts(Key, keys);

			StringBuilder result = new StringBuilder(keys[0].ToString(), 128);

			for (int i = 1; i <= Mesh.ABSOLUTE_MAX_LEVEL; i++)
				result.Append(',').Append(keys[i]);

			return result.ToString();
		}

	}
}