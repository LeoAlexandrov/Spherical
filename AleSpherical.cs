/* (c) AleProjects.com, 2018
 * v.2.0
 * 
 * MIT License
 * 
 * Convention in variable names: 
 *      Dec, Ra, Lon, Lat (short) - radians; 
 *      Declination, RightAscension, Longitude, Latitude (full) - degrees;
 *      
 *      Headings are always in radians.
 *      
 *      
 * Install System.ValueTuple package via Nuget
 * 
 * PM> Install-Package "System.ValueTuple"
 */

using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;


namespace AleProjects.Spherical
{

	/// <summary>
	/// Represents a basic vector in the 3D Euclidean space.
	/// </summary>
	public interface ICartesian
	{
		/// <summary>
		/// X component of the vector.
		/// </summary>
		double X { get; }

		/// <summary>
		/// Y component of the vector.
		/// </summary>
		double Y { get; }

		/// <summary>
		/// Z component of the vector.
		/// </summary>
		double Z { get; }

		/// <summary>
		/// Sets x,y,z components of the vector.
		/// </summary>
		/// <param name="x">x component.</param>
		/// <param name="y">y component.</param>
		/// <param name="z">z component.</param>
		void SetCartesian(double x, double y, double z);
	}



	/// <summary>
	/// Represents a basic spherical coordinate.
	/// </summary>

	public interface IGeoCoordinate
	{
		/// <summary>
		/// Latitude in degrees.
		/// </summary>
		double Latitude { get; }

		/// <summary>
		/// Longitude in degrees.
		/// </summary>
		double Longitude { get; }
	}



	/// <summary>
	/// Basic implementation of a vector in the 3D Euclidean space.
	/// </summary>
	public class Cartesian : ICartesian
	{
		/// <summary>
		/// X component of the vector.
		/// </summary>
		public double X { get; protected set; }

		/// <summary>
		/// Y component of the vector.
		/// </summary>
		public double Y { get; protected set; }

		/// <summary>
		/// Z component of the vector.
		/// </summary>
		public double Z { get; protected set; }

		/// <summary>
		/// Vector declination in radians.
		/// </summary>
		public double Dec
		{
			get => Math.Asin(Z / Math.Sqrt(X * X + Y * Y + Z * Z));
		}

		/// <summary>
		/// Vector right ascension in radians.
		/// </summary>
		public double Ra
		{
			get => Math.Atan2(Y, X);
		}

		/// <summary>
		/// Vector declination in degrees.
		/// </summary>
		public double Declination
		{
			get => Math.Asin(Z / Math.Sqrt(X * X + Y * Y + Z * Z)) * 180.0 / Math.PI;
		}

		/// <summary>
		/// Vector right ascension in degrees.
		/// </summary>
		public double RightAscension
		{
			get => Math.Atan2(Y, X) * 180.0 / Math.PI;
		}

		/// <summary>
		/// Length of the vector. 
		/// </summary>
		public double Length
		{
			get => SphericalExtension._VectorLength(X, Y, Z);
		}


		/// <summary>
		/// Creates a vector with X=1, Y=0, Z=0 what is the same as latitude=0 and Longitude=0.
		/// </summary>
		public Cartesian()
		{
			X = 1.0;
			Y = 0.0;
			Z = 0.0;
		}

		/// <summary>
		/// Creates a vector with X,Y,Z components, optionally normalizes it.
		/// </summary>
		/// <param name="x">X component.</param>
		/// <param name="y">Y component.</param>
		/// <param name="z">Y component.</param>
		/// <param name="normalize">Normalizes the vector when true.</param>
		public Cartesian(double x, double y, double z, bool normalize)
		{
			if (normalize)
				(x, y, z) = SphericalExtension._Normalized(x, y, z);

			X = x;
			Y = y;
			Z = z;
		}

		/// <summary>
		/// Creates a normalized vector with given declination and right ascension.
		/// </summary>
		/// <param name="dec">Declination</param>
		/// <param name="ra">Right ascension</param>
		/// <param name="degrees">Considers dec and ra as given in degrees when true.</param>
		public Cartesian(double dec, double ra, bool degrees = true)
		{
			if (degrees)
			{
				dec *= Math.PI / 180.0;
				ra *= Math.PI / 180.0;
			}

			double cos_dec = Math.Cos(dec);

			X = cos_dec * Math.Cos(ra);
			Y = cos_dec * Math.Sin(ra);
			Z = Math.Sin(dec);
		}

		/// <summary>
		/// Creates a new vector from another vector.
		/// </summary>
		/// <param name="cartesian">The vector whose X,Y,Z used for the new vector.</param>
		public Cartesian(ICartesian cartesian)
		{
			X = cartesian.X;
			Y = cartesian.Y;
			Z = cartesian.Z;
		}


		/// <summary>
		/// Sets x,y,z components of this vector.
		/// </summary>
		/// <param name="x">x component.</param>
		/// <param name="y">y component.</param>
		/// <param name="z">z component.</param>
		public virtual void SetCartesian(double x, double y, double z)
		{
			X = x;
			Y = y;
			Z = z;
		}

		/// <summary>
		/// Sets declination and right ascension of this vector.
		/// </summary>
		/// <param name="dec">Declination in radians.</param>
		/// <param name="ra">Right ascension in radians.</param>
		public virtual void SetSpherical(double dec, double ra)
		{
			double cos_dec = Math.Cos(dec);

			X = cos_dec * Math.Cos(ra);
			Y = cos_dec * Math.Sin(ra);
			Z = Math.Sin(dec);
		}

		/// <summary>
		/// Returns a string representation of the vector.
		/// </summary>
		/// <returns></returns>
		public override string ToString()
		{
			return string.Format("({0}; {1}; {2})", X, Y, Z);
		}
	}



	/// <summary>
	/// Basic value implementation of a vector in the 3D Euclidean space.
	/// </summary>
	public struct CartesianValue : ICartesian
	{
		/// <summary>
		/// X component of the vector.
		/// </summary>
		public double X { get; set; }

		/// <summary>
		/// Y component of the vector.
		/// </summary>
		public double Y { get; set; }

		/// <summary>
		/// Z component of the vector.
		/// </summary>
		public double Z { get; set; }


		public CartesianValue(double x, double y, double z)
		{
			X = x;
			Y = y;
			Z = z;
		}

		public CartesianValue(ICartesian cartesian)
		{
			X = cartesian.X;
			Y = cartesian.Y;
			Z = cartesian.Z;
		}

		public CartesianValue(double dec, double ra, bool degrees = true)
		{
			if (degrees)
			{
				dec *= Math.PI / 180.0;
				ra *= Math.PI / 180.0;
			}

			double cos_dec = Math.Cos(dec);

			X = cos_dec * Math.Cos(ra);
			Y = cos_dec * Math.Sin(ra);
			Z = Math.Sin(dec);
		}

		/// <summary>
		/// Sets x,y,z components of this vector.
		/// </summary>
		/// <param name="x">x component.</param>
		/// <param name="y">y component.</param>
		/// <param name="z">z component.</param>
		public void SetCartesian(double x, double y, double z)
		{
			X = x;
			Y = y;
			Z = z;
		}

		/// <summary>
		/// Sets declination and right ascension of this vector.
		/// </summary>
		/// <param name="dec">Declination in radians.</param>
		/// <param name="ra">Right ascension in radians.</param>
		public void SetSpherical(double dec, double ra)
		{
			double cos_dec = Math.Cos(dec);

			X = cos_dec * Math.Cos(ra);
			Y = cos_dec * Math.Sin(ra);
			Z = Math.Sin(dec);
		}

		/// <summary>
		/// Sets normalized x,y,z components of this vector.
		/// </summary>
		/// <param name="x">x component.</param>
		/// <param name="y">y component.</param>
		/// <param name="z">z component.</param>
		public void SetNormalized(double x, double y, double z)
		{
			(x, y, z) = SphericalExtension._Normalized(x, y, z);

			X = x;
			Y = y;
			Z = z;
		}

		public override string ToString()
		{
			return string.Format("({0}; {1}; {2})", X, Y, Z);
		}
	}



	/// <summary>
	/// Basic value implementation of a geo point.
	/// </summary>
	public struct LatLonValue : IGeoCoordinate
	{
		/// <summary>
		/// Latitude in degrees.
		/// </summary>
		public double Latitude { get; set; }

		/// <summary>
		/// Longitude in degrees.
		/// </summary>
		public double Longitude { get; set; }

		public LatLonValue(double latitude, double longitude)
		{
			Latitude = latitude;
			Longitude = longitude;
		}

		public LatLonValue(string text)
		{
			if (string.IsNullOrEmpty(text))
			{
				Latitude = 0.0;
				Longitude = 0.0;

				return;
			}

			int commaIndex = text.IndexOf(',');

			if (commaIndex <= 0 || commaIndex == text.Length - 1)
				throw new ArgumentException("Invalid text format for LatLonValue.", nameof(text));

			var fmt = new System.Globalization.NumberFormatInfo();

			if (!double.TryParse(text.AsSpan(0, commaIndex), fmt, out double latitude) ||
				latitude > 90.0 || latitude < -90.0)
				throw new ArgumentException("Invalid latitude value in text.", nameof(text));

			if (!double.TryParse(text.AsSpan(commaIndex + 1), fmt, out double longitude) ||
				longitude > 180.0 || longitude < -180.0)
				throw new ArgumentException("Invalid longitude value in text.", nameof(text));

			Latitude = latitude;
			Longitude = longitude;
		}

		public static bool TryParse(string text, out LatLonValue val)
		{
			val = default;

			if (string.IsNullOrEmpty(text))
				return false;

			int commaIndex = text.IndexOf(',');

			if (commaIndex <= 0 || commaIndex == text.Length - 1)
				return false;

			var fmt = new System.Globalization.NumberFormatInfo();

			if (!double.TryParse(text.AsSpan(0, commaIndex), fmt, out double latitude) || latitude > 90.0 || latitude < -90.0)
				return false;

			if (!double.TryParse(text.AsSpan(commaIndex + 1), fmt, out double longitude) || longitude > 180.0 || longitude < -180.0)
				return false;

			val.Latitude = latitude;
			val.Longitude = longitude;

			return true;
		}

	}



	/// <summary>
	/// Defines delegate type for user polyline section test
	/// </summary>
	/// <param name="before"></param>
	/// <param name="start"></param>
	/// <param name="finish"></param>
	/// <param name="angle"></param>
	/// <param name="azimuth"></param>
	/// <returns>True if test succeeds.</returns>
	public delegate bool PolylineSectionTest<T>(T before, T start, T finish, double angle, double bearing)
		where T : ICartesian;



	/// <summary>
	/// Defines delegate type for user route section test
	/// </summary>
	/// <param name="before"></param>
	/// <param name="start"></param>
	/// <param name="finish"></param>
	/// <param name="angle"></param>
	/// <param name="azimuth"></param>
	/// <param name="moveAzimuth"></param>
	/// <param name="azimuthTolerance"></param>
	/// <returns>True if test succeeds.</returns>
	public delegate bool RouteSectionTest<T>(T before, T start, T finish, double angle, double bearing, double moveBearing, double bearingTolerance)
		where T : ICartesian;



	/// <summary>
	/// Represents result of a test of a point and a polyline. 
	/// </summary>
	public class PolylineTestResult
	{
		protected double _TotalDistance = -1.0;

		public int Id { get; set; }
		public int Fails { get; set; }
		public int SectionIndex { get; set; }
		public double AngleWithSectionPlane { get; set; } = -1.0;
		public double SectionDirection { get; set; } = -1.0;
		public double AngleWithNextVertex { get; set; } = -1.0;
		public double[] SectionsAngles { get; protected set; }

		public double TotalDistance(double sphereRadius)
		{
			if (_TotalDistance < 0.0 && SectionsAngles != null)
			{
				_TotalDistance = 0.0;

				for (int i = 0; i < SectionsAngles.Length; i++)
					_TotalDistance += SectionsAngles[i];
			}

			return _TotalDistance * sphereRadius;

		}

		public double DistanceToPolyline(double sphereRadius)
		{
			return AngleWithSectionPlane >= 0.0 ? AngleWithSectionPlane * sphereRadius : -1.0;
		}

		public double DistanceToEnd(double sphereRadius)
		{
			double result = AngleWithNextVertex;
			int start = SectionIndex >= 0 ? SectionIndex + 1 : -1 - SectionIndex;

			for (int i = start; i < SectionsAngles.Length; i++)
				result += SectionsAngles[i];

			return result * sphereRadius;
		}

		public static PolylineTestResult Create<T>(IEnumerable<T> polyline, bool reverse,
			int sectionIndex, double sectionDir, double angleWithSectionPlane, double angleWithNextVertex)
			where T : ICartesian

		{
			PolylineTestResult result = new PolylineTestResult()
			{
				SectionIndex = sectionIndex,
				SectionDirection = sectionDir,
				AngleWithSectionPlane = angleWithSectionPlane,
				AngleWithNextVertex = angleWithNextVertex,
				SectionsAngles = new double[polyline.Count() - 1]
			};

			IEnumerable<T> vertices = reverse ? polyline.Reverse().Skip(1) : polyline.Skip(1);
			T previous = reverse ? polyline.Last() : polyline.First();
			int i = 0;

			foreach (var vertex in vertices)
			{
				result.SectionsAngles[i++] = previous.Angle(vertex);
				previous = vertex;
			}


			return result;
		}

	}



	/// <summary>
	/// Represents result of a test of a route. 
	/// </summary>
	public class RouteTestResult
	{
		public int Id { get; set; }
		public int Fails { get; set; }
		public int SectionIndex { get; set; }
		public double AngleWithSectionPlane { get; set; } = -1.0;
		public double SectionDirection { get; set; } = -1.0;
		public double AngleWithNextVertex { get; set; } = -1.0;
	}



	/// <summary>
	/// Provides extension methods for ICartesian and IGeoCoordinate objects.
	/// </summary>
	public static class SphericalExtension
	{
		public const double EPSILON = 2.2204460492503131e-016;
		public const double EARTH_MEAN_RADIUS = 6371000;

		private const string Error_Message_Not_Polygon = "Parameter is not a polygon.";
		private const string Error_Message_Not_Convex = "Parameter is not a convex.";
		private const string Error_Message_Not_Polyline = "Parameter is not a polyline.";
		private const string Error_Message_Empty_Collection = "Collection can't be empty.";

		// properties

		/// <summary>
		/// Maximum latitude in radians a map based on Mercator projection can display.
		/// </summary>
		public static double MaxMercatorMapsLat
		{
			get => Math.Atan(Math.Sinh(Math.PI));
		}


		// Auxiliary methods

		private static T NewCartesian<T>(double x, double y, double z)
			where T : ICartesian, new()
		{
			T result = new T();
			result.SetCartesian(x, y, z);

			return result;
		}

		private static T NewCartesian<T>(double dec, double ra)
			where T : ICartesian, new()
		{
			T result = new T();

			double cos_dec = Math.Cos(dec);
			double x = cos_dec * Math.Cos(ra);
			double y = cos_dec * Math.Sin(ra);
			double z = Math.Sin(dec);

			result.SetCartesian(x, y, z);

			return result;
		}



		// Helping methods

		#region Helping-methods

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		internal static double _VectorLength(double ax, double ay, double az)
		{
			return Math.Sqrt(ax * ax + ay * ay + az * az);
		}

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		internal static (double, double, double) _Normalized(double ax, double ay, double az)
		{
			double l = Math.Sqrt(ax * ax + ay * ay + az * az);

			if (l > 0.0)
			{
				ax /= l;
				ay /= l;
				az /= l;
			}

			return (ax, ay, az);
		}

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		internal static double _DotProduct(double ax, double ay, double az, double bx, double by, double bz)
		{
			return ax * bx + ay * by + az * bz;
		}

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		internal static (double, double, double) _CrossProduct(double ax, double ay, double az, double bx, double by, double bz, bool normalize)
		{
			double x = ay * bz - az * by;
			double y = -(ax * bz - az * bx);
			double z = ax * by - ay * bx;

			if (normalize)
			{
				double l = Math.Sqrt(x * x + y * y + z * z);

				if (l > 0.0)
				{
					x /= l;
					y /= l;
					z /= l;
				}
			}

			return (x, y, z);
		}

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		internal static double _TripleProduct(double ax, double ay, double az, double bx, double by, double bz, double cx, double cy, double cz)
		{
			return ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx);
		}

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		internal static bool _InsideTriangle(double ax, double ay, double az, double bx, double by, double bz, double cx, double cy, double cz, double dx, double dy, double dz)
		{
			double x1 = ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx);
			double x2 = ax * (cy * dz - cz * dy) - ay * (cx * dz - cz * dx) + az * (cx * dy - cy * dx);

			if (x1 * x2 >= 0.0)
			{
				double x3 = ax * (dy * bz - dz * by) - ay * (dx * bz - dz * bx) + az * (dx * by - dy * bx);

				if (x2 * x3 >= 0.0 && x1 * x3 >= 0.0)
				{
					double sx = bx + cx + dx;
					double sy = by + cy + dy;
					double sz = bz + cz + dz;

					return ax * sx + ay * sy + az * sz >= 0.0;
				}
			}

			return false;
		}

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		internal static bool _InsideTriangle(double ax, double ay, double az, in CartesianValue b, in CartesianValue c, in CartesianValue d)
		{
			double x1 = ax * (b.Y * c.Z - b.Z * c.Y) - ay * (b.X * c.Z - b.Z * c.X) + az * (b.X * c.Y - b.Y * c.X);
			double x2 = ax * (c.Y * d.Z - c.Z * d.Y) - ay * (c.X * d.Z - c.Z * d.X) + az * (c.X * d.Y - c.Y * d.X);

			if (x1 * x2 >= 0.0)
			{
				double x3 = ax * (d.Y * b.Z - d.Z * b.Y) - ay * (d.X * b.Z - d.Z * b.X) + az * (d.X * b.Y - d.Y * b.X);

				if (x2 * x3 >= 0.0 && x1 * x3 >= 0.0)
				{
					double sx = b.X + c.X + d.X;
					double sy = b.Y + c.Y + d.Y;
					double sz = b.Z + c.Z + d.Z;

					return ax * sx + ay * sy + az * sz >= 0.0;
				}
			}

			return false;
		}

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		internal static double _Cosine(double ax, double ay, double az, double bx, double by, double bz)
		{
			double r = (ax * bx + ay * by + az * bz) / Math.Sqrt((ax * ax + ay * ay + az * az) * (bx * bx + by * by + bz * bz));

			if (Math.Abs(r) > 1.0)
				r = Math.Truncate(r);

			return r;
		}

		#endregion Helping-methods



		// ICartesian extension methods

		/// <summary>
		/// Extension method calculating a length of the vector.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">Vector.</param>
		/// <returns>Length of the vector.</returns>
		public static double VectorLength<T>(this T cartesian)
			where T : ICartesian
		{
			return _VectorLength(cartesian.X, cartesian.Y, cartesian.Z);
		}

		/// <summary>
		/// Extension method calculating a declination of the vector (latitude for geopoint).
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">Vector.</param>
		/// <returns>Declination in radians.</returns>
		public static double Dec<T>(this T cartesian)
			where T : ICartesian
		{
			return Math.Asin(cartesian.Z / _VectorLength(cartesian.X, cartesian.Y, cartesian.Z));
		}

		/// <summary>
		/// Extension method calculating raght ascension of the vector (longitude for geopoint).
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">Vector.</param>
		/// <returns>Raght ascension in radians.</returns>
		public static double Ra<T>(this T cartesian)
			where T : ICartesian
		{
			return Math.Atan2(cartesian.Y, cartesian.X);
		}

		/// <summary>
		/// Extension method calculating a declination of the vector (latitude for geopoint).
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">Vector.</param>
		/// <returns>Declination in degrees.</returns>
		public static double Declination<T>(this T cartesian)
			where T : ICartesian
		{
			return Math.Asin(cartesian.Z / _VectorLength(cartesian.X, cartesian.Y, cartesian.Z)) * 180.0 / Math.PI;
		}

		/// <summary>
		/// Extension method calculating raght ascension of the vector (longitude for geopoint).
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">Vector.</param>
		/// <returns>Raght ascension in degrees.</returns>
		public static double RightAscension<T>(this T cartesian)
			where T : ICartesian
		{
			return Math.Atan2(cartesian.Y, cartesian.X) * 180.0 / Math.PI;
		}

		/// <summary>
		/// Extension method calculating latitude if the vector represents geopoint. Equal to the Declination&lt;T&gt; method.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">Vector.</param>
		/// <returns>Declination in degrees.</returns>
		public static double Latitude<T>(this T cartesian)
			where T : ICartesian
		{
			return Math.Asin(cartesian.Z / _VectorLength(cartesian.X, cartesian.Y, cartesian.Z)) * 180.0 / Math.PI;
		}

		/// <summary>
		/// Extension method calculating longitude if the vector represents geopoint. Equal to the RightAscension&lt;T&gt; method.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">Vector.</param>
		/// <returns>Raght ascension in degrees.</returns>
		public static double Longitude<T>(this T cartesian)
			where T : ICartesian
		{
			return Math.Atan2(cartesian.Y, cartesian.X) * 180.0 / Math.PI;
		}

		/// <summary>
		/// Extension method calculating rounded latitude if the vector represents geopoint. Equal to the Declination&lt;T&gt; method.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">Vector.</param>
		/// <param name="digits">Number of fractional digits.</param>
		/// <returns>Declination in degrees.</returns>
		public static double Latitude<T>(this T cartesian, int digits)
			where T : ICartesian
		{
			return Math.Round(Math.Asin(cartesian.Z / _VectorLength(cartesian.X, cartesian.Y, cartesian.Z)) * 180.0 / Math.PI, digits);
		}

		/// <summary>
		/// Extension method calculating rounded longitude if the vector represents geopoint. Equal to the RightAscension&lt;T&gt; method.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">Vector.</param>
		/// <param name="digits">Number of fractional digits.</param>
		/// <returns>Raght ascension in degrees.</returns>
		public static double Longitude<T>(this T cartesian, int digits)
			where T : ICartesian
		{
			return Math.Round(Math.Atan2(cartesian.Y, cartesian.X) * 180.0 / Math.PI, digits);
		}

		/// <summary>
		/// Tests if a vector is zero-vector, i.e X,Y,Z components are equal to 0 with given precision. 
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">Vector to test.</param>
		/// <param name="precision">Precision of comparing X,Y,Z components with 0.</param>
		/// <returns></returns>
		public static bool IsZeroVector<T>(this T cartesian, double precision)
			where T : ICartesian
		{
			return Math.Abs(cartesian.X) <= precision &&
				Math.Abs(cartesian.Y) <= precision &&
				Math.Abs(cartesian.Z) <= precision;
		}

		/// <summary>
		/// Extension method normalizing the specified vector.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">Vector to normalize.</param>
		/// <returns>Previous length of the vector.</returns>
		public static void Normalize<T>(this T cartesian)
			where T : ICartesian
		{
			var (x, y, z) = _Normalized(cartesian.X, cartesian.Y, cartesian.Z);

			cartesian.SetCartesian(x, y, z);
		}


		public static T Add<T, U>(this T cartesian, U V, bool normalize)
			where T : ICartesian, new()
			where U : ICartesian
		{
			double x = cartesian.X + V.X;
			double y = cartesian.Y + V.Y;
			double z = cartesian.Z + V.Z;

			if (normalize)
				(x, y, z) = _Normalized(x, y, z);

			return NewCartesian<T>(x, y, z);
		}

		/// <summary>
		/// Extension method calculating dot product of two vectors.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface and have default constructor.</typeparam>
		/// <typeparam name="U">Type U must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">First vector.</param>
		/// <param name="V">Second vector.</param>
		/// <returns>Value of dot product.</returns>
		public static double DotProduct<T, U>(this T cartesian, U V)
			where T : ICartesian
			where U : ICartesian
		{
			return _DotProduct(cartesian.X, cartesian.Y, cartesian.Z, V.X, V.Y, V.Z);
		}

		/// <summary>
		/// Extension method calculating dot product of two vectors.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">First vector.</param>
		/// <param name="x">X component of the second vector.</param>
		/// <param name="y">Y component of the second vector.</param>
		/// <param name="z">Z component of the second vector.</param>
		/// <returns>Value of dot product.</returns>
		public static double DotProduct<T>(this T cartesian, double x, double y, double z)
			where T : ICartesian
		{
			return _DotProduct(cartesian.X, cartesian.Y, cartesian.Z, x, y, z);
		}

		/// <summary>
		/// Extension method calculating cross product of two vectors. The second vector is represented by X,Y,Z components. 
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface and have default constructor.</typeparam>
		/// <typeparam name="U">Type U must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">First vector.</param>
		/// <param name="x">X component of the second vector.</param>
		/// <param name="y">Y component of the second vector.</param>
		/// <param name="z">Z component of the second vector.</param>
		/// <param name="normalize">When true, method returns normalized vector.</param>
		/// <returns>Vector representing cross product.</returns>
		public static T CrossProduct<T, U>(this U cartesian, double x, double y, double z, bool normalize)
			where T : ICartesian, new()
			where U : ICartesian
		{
			var (x1, y1, z1) = _CrossProduct(cartesian.X, cartesian.Y, cartesian.Z, x, y, z, normalize);

			return NewCartesian<T>(x1, y1, z1);
		}

		/// <summary>
		/// Extension method calculating cross product of two vectors.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface and have default constructor.</typeparam>
		/// <typeparam name="U">Type U must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">First vector.</param>
		/// <param name="V">Second vector.</param>
		/// <param name="normalize">When true, method returns normalized vector.</param>
		/// <returns>Vector representing cross product.</returns>
		public static T CrossProduct<T, U>(this T cartesian, U V, bool normalize)
			where T : ICartesian, new()
			where U : ICartesian
		{
			var (x, y, z) = _CrossProduct(cartesian.X, cartesian.Y, cartesian.Z, V.X, V.Y, V.Z, normalize);

			return NewCartesian<T>(x, y, z);
		}

		public static T CrossProduct<T, U>(this U cartesian, in CartesianValue V, bool normalize)
			where T : ICartesian, new()
			where U : ICartesian
		{
			var (x, y, z) = _CrossProduct(cartesian.X, cartesian.Y, cartesian.Z, V.X, V.Y, V.Z, normalize);

			return NewCartesian<T>(x, y, z);
		}

		/// <summary>
		/// Extension method calculating cross product of two vectors. The second vector is represented by declination and right ascension.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface and have default constructor.</typeparam>
		/// <typeparam name="U">Type U must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian"></param>
		/// <param name="dec">Declination in radians (latitude for geopoint)</param>
		/// <param name="ra">Right ascension (longitude for geopoint)</param>
		/// <returns>Vector representing cross product.</returns>
		public static T CrossProduct<T, U>(this U cartesian, double dec, double ra)
			where T : ICartesian, new()
			where U : ICartesian
		{
			double cos_dec = Math.Cos(dec);
			var (x, y, z) = _CrossProduct(cartesian.X, cartesian.Y, cartesian.Z, cos_dec * Math.Cos(ra), cos_dec * Math.Sin(ra), Math.Sin(dec), true);

			return NewCartesian<T>(x, y, z);
		}

		/// <summary>
		/// Extension method calculating cross product of two vectors. Type T must implement ICartesian interface. 
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <typeparam name="U">Type U must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">First vector.</param>
		/// <param name="V">Second vector.</param>
		/// <returns>Cross product decomposed to x,y,z components.</returns>
		public static (double x, double y, double z) CrossProduct<T, U>(this T cartesian, U V)
			where T : ICartesian
			where U : ICartesian
		{
			return _CrossProduct(cartesian.X, cartesian.Y, cartesian.Z, V.X, V.Y, V.Z, false);
		}

		/// <summary>
		/// Extension method calculating cross product of two vectors. The second vector is represented by X,Y,Z components. Type T must implement ICartesian interface. 
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">First vector.</param>
		/// <param name="x">X component of the second vector.</param>
		/// <param name="y">Y component of the second vector.</param>
		/// <param name="z">Z component of the second vector.</param>
		/// <returns>Cross product decomposed to x,y,z components.</returns>
		public static (double x, double y, double z) CrossProduct<T>(this T cartesian, double x, double y, double z)
			where T : ICartesian
		{
			return _CrossProduct(cartesian.X, cartesian.Y, cartesian.Z, x, y, z, false);
		}

		/// <summary>
		/// Extension method calculating triple product of three vectors.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface and have default constructor.</typeparam>
		/// <typeparam name="U">Type U must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">First vector.</param>
		/// <param name="V1">Second vector.</param>
		/// <param name="V2">Third vector.</param>
		/// <returns>Value of triple product.</returns>
		public static double TripleProduct<T, U>(this T cartesian, U V1, U V2)
			where T : ICartesian
			where U : ICartesian
		{
			return _TripleProduct(cartesian.X, cartesian.Y, cartesian.Z, V1.X, V1.Y, V1.Z, V2.X, V2.Y, V2.Z);
		}

		/// <summary>
		/// Extension method rotating a vector around an axis by an angle. The result, the initial vector, and the axis are the right hand rule vectors. 
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface and have default constructor.</typeparam>
		/// <typeparam name="U">Type U must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">Vector to rotate.</param>
		/// <param name="angle">Rotation angle.</param>
		/// <param name="axisX">Rotation axis X.</param>
		/// <param name="axisY">Rotation axis Y.</param>
		/// <param name="axisZ">Rotation axis Z.</param>
		/// <returns>New vector obtained from rotation.</returns>
		public static T Rotate<T, U>(this U cartesian, double angle, double axisX, double axisY, double axisZ)
			where T : ICartesian, new()
			where U : ICartesian
		{
			var (x, y, z) = _Normalized(axisX, axisY, axisZ);

			double cosine = Math.Cos(angle);
			double sine = Math.Sin(angle);
			double one_minus_cos = 1.0 - cosine;

			// rotation matrix by lines
			// line 1
			double M11 = cosine + one_minus_cos * x * x;
			double M12 = one_minus_cos * x * y - sine * z;
			double M13 = one_minus_cos * x * z + sine * y;
			// line 2
			double M21 = one_minus_cos * y * x + sine * z;
			double M22 = cosine + one_minus_cos * y * y;
			double M23 = one_minus_cos * y * z - sine * x;
			// line 3
			double M31 = one_minus_cos * z * x - sine * y;
			double M32 = one_minus_cos * z * y + sine * x;
			double M33 = cosine + one_minus_cos * z * z;

			return NewCartesian<T>(
				_DotProduct(cartesian.X, cartesian.Y, cartesian.Z, M11, M12, M13),
				_DotProduct(cartesian.X, cartesian.Y, cartesian.Z, M21, M22, M23),
				_DotProduct(cartesian.X, cartesian.Y, cartesian.Z, M31, M32, M33));
		}

		/// <summary>
		/// This extension method rotates a sphere with given vector around some axis. 
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface and have default constructor.</typeparam>
		/// <typeparam name="U">Type U must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">Vector on the sphere.</param>
		/// <param name="angle">The angle to rotate the sphere.</param>
		/// <param name="axis">The axis to rotate the sphere around.</param>
		/// <returns></returns>
		public static T SpinAround<T, U>(this T cartesian, double angle, U axis)
			where T : ICartesian, new()
			where U : ICartesian
		{
			CartesianValue N = cartesian.CrossProduct<CartesianValue, T>(axis.X, axis.Y, axis.Z, false);
			CartesianValue V = axis.CrossProduct<CartesianValue, U>(N, false);
			CartesianValue R = V.Rotate<CartesianValue, CartesianValue>(-angle, axis.X, axis.Y, axis.Z);
			var (Nx, Ny, Nz) = R.CrossProduct(axis.X, axis.Y, axis.Z);

			return R.Rotate<T, CartesianValue>(Math.PI / 2 - cartesian.Angle(axis), Nx, Ny, Nz);
		}

		/// <summary>
		/// Extension method returning a vector representing a azimuth plane containing given location. 
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface and have default constructor.</typeparam>
		/// <typeparam name="U">Type U must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">Vector that represents location on a sphere.</param>
		/// <param name="bearing">Bearing on a sphere.</param>
		/// <returns></returns>
		public static T BearingPlane<T, U>(this U cartesian, double bearing)
			where T : ICartesian, new()
			where U : ICartesian
		{
			T result;
			double dec = cartesian.Dec();
			double ra = cartesian.Ra();

			if (Math.Abs(dec) <= EPSILON)
			{
				var (x, y, z) = SphericalToCartesian(bearing, ra - Math.PI / 2.0);
				result = NewCartesian<T>(x, y, z);
			}
			else if (Math.Abs(Math.PI / 2 - Math.Abs(dec)) <= EPSILON)
				result = default;
			else
				result = cartesian.CrossProduct<T, U>(0.0, ra - Math.Atan(Math.Sin(dec) * Math.Tan(bearing)));

			return result;
		}

		public static T Middle<T, U>(this U cartesian, U V)
			where T : ICartesian, new()
			where U : ICartesian
		{
			var (x1, y1, z1) = _Normalized(cartesian.X, cartesian.Y, cartesian.Z);
			var (x2, y2, z2) = _Normalized(V.X, V.Y, V.Z);

			return NewCartesian<T>(x1 + x2, y1 + y2, z1 + z2);
		}

		/// <summary>
		/// Extension method testing whether a vector belongs to a triangle formed by three vectors. 
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <typeparam name="U">Type U must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">Vector to test.</param>
		/// <param name="V1">First vertex of the triangle.</param>
		/// <param name="V2">Second vertex of the triangle.</param>
		/// <param name="V3">Third vertex of the triangle.</param>
		/// <returns>True if the vector belongs to the triangle.</returns>
		public static bool InsideTriangle<T, U>(this T cartesian, U V1, U V2, U V3)
			where T : ICartesian
			where U : ICartesian
		{
			return _InsideTriangle(cartesian.X, cartesian.Y, cartesian.Z, V1.X, V1.Y, V1.Z, V2.X, V2.Y, V2.Z, V3.X, V3.Y, V3.Z);
		}

		/// <summary>
		/// Extension method testing if a vector is inside a convex including its borders.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <typeparam name="U">Type U must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">Vector to test.</param>
		/// <param name="polygon">Vertices of the convex.</param>
		/// <returns>True if inside.</returns>
		public static bool InsideConvex<T, U>(this T cartesian, IEnumerable<U> polygon)
			where T : ICartesian
			where U : ICartesian
		{
			if (polygon.Count() < 3)
				throw new ArgumentException(Error_Message_Not_Polygon, nameof(polygon));

			U first = polygon.First();
			U previous = polygon.Last();
			U current = first;

			foreach (U next in polygon.Skip(1))
			{
				if (_InsideTriangle(cartesian.X, cartesian.Y, cartesian.Z,
						previous.X, previous.Y, previous.Z,
						current.X, current.Y, current.Z,
						next.X, next.Y, next.Z)) return true;

				previous = current;
				current = next;
			}

			return _InsideTriangle(cartesian.X, cartesian.Y, cartesian.Z,
				previous.X, previous.Y, previous.Z,
				current.X, current.Y, current.Z,
				first.X, first.Y, first.Z);
		}

		/// <summary>
		/// Extension method testing if a vector is inside a polygon (including borders). Uses Winding number method.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <typeparam name="U">Type U must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">Vector to test.</param>
		/// <param name="polygon">Vertices of the polygon.</param>
		/// <param name="polygonCenter">Vector inside the polygon for faster test.</param>
		/// <returns>True if inside.</returns>
		public static bool InsidePolygon<T, U>(this T cartesian, IEnumerable<U> polygon, U polygonCenter = default)
			where T : ICartesian
			where U : ICartesian
		{
			if (polygon.Count() < 3)
				throw new ArgumentException(Error_Message_Not_Polygon, nameof(polygon));

			U previous = polygon.Last();

			var (Nx, Ny, Nz) = _CrossProduct(cartesian.X, cartesian.Y, cartesian.Z, previous.X, previous.Y, previous.Z, false);

			if (Math.Abs(Nx) <= EPSILON && Math.Abs(Ny) <= EPSILON && Math.Abs(Nz) <= EPSILON)
				return true;

			double wn = 0.0;

			foreach (U current in polygon)
			{
				var (N1x, N1y, N1z) = _CrossProduct(cartesian.X, cartesian.Y, cartesian.Z, current.X, current.Y, current.Z, false);

				if (Math.Abs(N1x) <= EPSILON && Math.Abs(N1y) <= EPSILON && Math.Abs(N1z) <= EPSILON)
					return true;

				if (_TripleProduct(cartesian.X, cartesian.Y, cartesian.Z, previous.X, previous.Y, previous.Z, current.X, current.Y, current.Z) > 0.0)
					wn += Math.Acos(_Cosine(Nx, Ny, Nz, N1x, N1y, N1z));
				else
					wn -= Math.Acos(_Cosine(Nx, Ny, Nz, N1x, N1y, N1z));

				Nx = N1x;
				Ny = N1y;
				Nz = N1z;
				previous = current;
			}

			if (Math.Abs(wn) > Math.PI)
			{
				if (EqualityComparer<U>.Default.Equals(polygonCenter, default))
				{
					Nx = 0.0;
					Ny = 0.0;
					Nz = 0.0;

					foreach (U current in polygon)
					{
						Nx += current.X;
						Ny += current.Y;
						Nz += current.Z;
					}
				}
				else
				{
					Nx = polygonCenter.X;
					Ny = polygonCenter.Y;
					Nz = polygonCenter.Z;
				}

				return cartesian.X * Nx + cartesian.Y * Ny + cartesian.Z * Nz >= 0.0;
			}

			return false;
		}

		/// <summary>
		/// Inflates a convex polygon by rotating its planes by given angles.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface and have default constructor.</typeparam>
		/// <typeparam name="U">Type U must implement the ICartesian interface.</typeparam>
		/// <param name="polygon">Polygon to inflate.</param>
		/// <param name="angles">Angles of rotation in radians. Negative value deflates the polygon.</param>
		/// <param name="result">Inflated polygon.</param>
		public static void InflateConvex<T, U>(IEnumerable<U> polygon, IReadOnlyList<double> angles, List<T> result)
			where T : ICartesian, new()
			where U : ICartesian
		{
			int n = polygon.Count();

			if (n < 3)
				throw new ArgumentException(Error_Message_Not_Polygon, nameof(polygon));

			int m = angles.Count;

			if (m == 0)
				throw new ArgumentException(Error_Message_Empty_Collection, nameof(angles));

			U first = polygon.First();
			U last = polygon.Last();
			U previous = last;
			U current = first;
			bool? cw = null;
			double x;

			foreach (U next in polygon.Skip(1))
			{
				x = _TripleProduct(previous.X, previous.Y, previous.Z, current.X, current.Y, current.Z, next.X, next.Y, next.Z);

				if (Math.Abs(x) <= EPSILON)
					throw new ArgumentException(Error_Message_Not_Convex, nameof(polygon));

				if (!cw.HasValue)
					cw = Math.Sign(x) > 0;
				else if (cw.Value != (Math.Sign(x) > 0))
					throw new ArgumentException(Error_Message_Not_Convex, nameof(polygon));

				previous = current;
				current = next;
			}

			x = _TripleProduct(previous.X, previous.Y, previous.Z, current.X, current.Y, current.Z, first.X, first.Y, first.Z);

			if (Math.Abs(x) <= EPSILON)
				throw new ArgumentException(Error_Message_Not_Convex, nameof(polygon));

			if (!cw.HasValue)
				cw = Math.Sign(x) > 0.0;
			else if (cw.Value != Math.Sign(x) > 0.0)
				throw new ArgumentException(Error_Message_Not_Convex, nameof(polygon));

			bool clockwise = cw.Value;
			CartesianValue[] planes = new CartesianValue[n];
			CartesianValue N, M, axis;
			previous = last;
			int i = 0;
			int j = 0;

			foreach (U next in polygon)
			{
				N = previous.CrossProduct<CartesianValue, U>(next.X, next.Y, next.Z, false);
				M = previous.Middle<CartesianValue, U>(next);
				axis = N.CrossProduct(M, false);

				if (!clockwise)
					axis.SetCartesian(-axis.X, -axis.Y, -axis.Z);

				planes[i++] = N.Rotate<CartesianValue, CartesianValue>(angles[j == m ? m - 1 : j++], axis.X, axis.Y, axis.Z);
				previous = next;
			}

			if (result.Count != 0)
				result.Clear();

			T res;

			for (i = 0; i < n - 1; i++)
			{
				res = planes[i].CrossProduct<T, CartesianValue>(planes[i + 1], false);

				if (!clockwise)
					res.SetCartesian(-res.X, -res.Y, -res.Z);

				result.Add(res);
			}

			res = planes[n - 1].CrossProduct<T, CartesianValue>(planes[0], false);

			if (!clockwise)
				res.SetCartesian(-res.X, -res.Y, -res.Z);

			result.Add(res);
		}

		/// <summary>
		/// Inflates a convex polygon by rotating its planes by given angles.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface and have default constructor.</typeparam>
		/// <typeparam name="U">Type U must implement the ICartesian interface.</typeparam>
		/// <param name="polygon">Polygon to inflate.</param>
		/// <param name="angles">Angles of rotation in radians. Negative value deflates the polygon.</param>
		/// <returns>Inflated polygon</returns>
		public static List<T> InflateConvex<T, U>(IEnumerable<U> polygon, IReadOnlyList<double> angles)
			where T : ICartesian, new()
			where U : ICartesian
		{
			var result = new List<T>();

			InflateConvex(polygon, angles, result);

			return result;
		}

		/// <summary>
		/// Extension method calculating a cosine of an angle between two vectors.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <typeparam name="U">Type U must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">First vector</param>
		/// <param name="V">Second vector.</param>
		/// <returns>Cosine of the angle between two vectors.</returns>
		public static double Cosine<T, U>(this T cartesian, U V)
			where T : ICartesian
			where U : ICartesian
		{
			return _Cosine(cartesian.X, cartesian.Y, cartesian.Z, V.X, V.Y, V.Z);
		}

		/// <summary>
		/// Extension method calculating a cosine of an angle between two vectors, the second one is represented by X,Y,Z components.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">First vector.</param>
		/// <param name="x">X component of the second vector.</param>
		/// <param name="y">Y component of the second vector.</param>
		/// <param name="z">Z component of the second vector.</param>
		/// <returns>Cosine of the angle between two vectors.</returns>
		public static double Cosine<T>(this T cartesian, double x, double y, double z)
			where T : ICartesian
		{
			return _Cosine(cartesian.X, cartesian.Y, cartesian.Z, x, y, z);
		}

		/// <summary>
		/// Extension method calculating a cosine of anangle between two vectors, the second one is represented by declination and right ascension.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">First  vector.</param>
		/// <param name="dec">Declination of the second vector in radians (latitude for geopoint).</param>
		/// <param name="ra">Right ascension of the second vector in radians (longitude for geopoint).</param>
		/// <returns>Cosine of the angle between two vectors.</returns>
		public static double Cosine<T>(this T cartesian, double dec, double ra)
			where T : ICartesian
		{
			double cos_dec = Math.Cos(dec);

			return _Cosine(cartesian.X, cartesian.Y, cartesian.Z, cos_dec * Math.Cos(ra), cos_dec * Math.Sin(ra), Math.Sin(dec));
		}

		/// <summary>
		/// Extension method calculating an angle [0..π) between two vectors.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <typeparam name="U">Type U must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">First vector.</param>
		/// <param name="V">Second vector.</param>
		/// <returns>Angle in radians between two vectors.</returns>
		public static double Angle<T, U>(this T cartesian, U V)
			where T : ICartesian
			where U : ICartesian
		{
			return Math.Acos(_Cosine(cartesian.X, cartesian.Y, cartesian.Z, V.X, V.Y, V.Z));
		}

		/// <summary>
		/// Extension method calculating an angle [0..π) between two vectors, the second one is represented as X,Y,Z components.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">First vector.</param>
		/// <param name="x">X component of the second vector.</param>
		/// <param name="y">Y component of the second vector.</param>
		/// <param name="z">Z component of the second vector.</param>
		/// <returns>Angle in radians between two vectors.</returns>
		public static double Angle<T>(this T cartesian, double x, double y, double z)
			where T : ICartesian
		{
			return Math.Acos(_Cosine(cartesian.X, cartesian.Y, cartesian.Z, x, y, z));
		}

		/// <summary>
		/// Extension method calculating an angle [0..π) between two vectors, the second one is represented as declination and right ascension.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">First 3D vector.</param>
		/// <param name="dec">Declination of the second vector in radians (latitude for geopoint).</param>
		/// <param name="ra">Right ascension of the second vector in radians (longitude for geopoint).</param>
		/// <returns>Angle in radian between two vectorss.</returns>
		public static double Angle<T>(this T cartesian, double dec, double ra)
			where T : ICartesian
		{
			double cos_dec = Math.Cos(dec);
			double cosine = _Cosine(cartesian.X, cartesian.Y, cartesian.Z, cos_dec * Math.Cos(ra), cos_dec * Math.Sin(ra), Math.Sin(dec));

			return Math.Acos(cosine);
		}

		/// <summary>
		/// Extension method returning a bearing from one location on a sphere to another.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <typeparam name="U">Type U must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">First location.</param>
		/// <param name="V">Second location.</param>
		/// <returns>Bearing in radians.</returns>
		public static double BearingTo<T, U>(this T cartesian, U V)
			where T : ICartesian
			where U : ICartesian
		{
			// N1 - cross product of cartesian and axiz Z, N1z is always zero and not used.
			double N1x = cartesian.Y;
			double N1y = -cartesian.X;

			if (Math.Abs(N1x) <= EPSILON && Math.Abs(N1y) <= EPSILON)
				return cartesian.Z >= 0 ? Math.PI : 0.0;

			var (N2x, N2y, N2z) = _CrossProduct(cartesian.X, cartesian.Y, cartesian.Z, V.X, V.Y, V.Z, false);
			double result = Math.Acos(_Cosine(N2x, N2y, N2z, N1x, N1y, 0.0));
			double c_ra = cartesian.Ra();
			double v_ra = V.Ra();

			if (Math.Abs(v_ra - c_ra) > Math.PI)
			{
				if (v_ra > 0.0)
					result = 2.0 * Math.PI - result;
			}
			else if (v_ra < c_ra)
				result = 2.0 * Math.PI - result;

			return result;
		}

		/// <summary>
		/// Extension method returning a bearing from one location on a sphere to another.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">First location.</param>
		/// <param name="point">Second location.</param>
		/// <returns>Bearing in radians.</returns>
		public static double BearingTo<T>(this T cartesian, IGeoCoordinate point)
			where T : ICartesian
		{
			// N1 - cross product of cartesian and axiz Z, N1z is always zero and not used.
			double N1x = cartesian.Y;
			double N1y = -cartesian.X;

			if (Math.Abs(N1x) <= EPSILON && Math.Abs(N1y) <= EPSILON)
				return cartesian.Z >= 0 ? Math.PI : 0.0;

			double ra = point.Longitude * Math.PI / 180.0;
			double dec = point.Latitude * Math.PI / 180.0;
			double cos_dec = Math.Cos(point.Latitude * Math.PI / 180.0);
			var (N2x, N2y, N2z) = _CrossProduct(cartesian.X, cartesian.Y, cartesian.Z, cos_dec * Math.Cos(ra), cos_dec * Math.Sin(ra), Math.Sin(dec), false);
			double result = Math.Acos(_Cosine(N2x, N2y, N2z, N1x, N1y, 0.0));
			double c_ra = cartesian.Ra();

			if (Math.Abs(ra - c_ra) > Math.PI)
			{
				if (ra > 0.0)
					result = 2.0 * Math.PI - result;
			}
			else if (ra < c_ra)
				result = 2.0 * Math.PI - result;

			return result;
		}

		/// <summary>
		/// Extension method returning a bearing from one location on a sphere to another.
		/// </summary>
		/// <typeparam name="T">Type T must implement the IGeoCoordinate interface.</typeparam>
		/// <param name="point">First location.</param>
		/// <param name="latitude">Latitude for the second location in degrees.</param>
		/// <param name="longitude">Longitude for the second location in degrees.</param>
		/// <returns>Bearing in radians.</returns>
		public static double BearingTo<T>(this T point, double latitude, double longitude)
			where T : IGeoCoordinate
		{
			latitude *= Math.PI / 180.0;
			longitude *= Math.PI / 180.0;

			double pLat = point.Latitude * Math.PI / 180.0;
			double pLon = point.Longitude * Math.PI / 180.0;
			double cosine_c = Math.Sin(pLat);        //because angle c is pi/2-Lat
			double cosine_a = Math.Sin(latitude);

			double cosine_D = Math.Cos(pLat) * Math.Cos(latitude) * (Math.Cos(pLon - longitude) - 1.0) + Math.Cos(pLat - latitude);

			// when distance is 0 or second point is on exactly opposite side of the earth - all directions are equal and we consider heading as 0 
			// or when cosine_c = sin(Lat) == -1 it means that Latitude is -90 degrees and we are in South Pole, so every direction is to north and heading=0
			if (1 - Math.Abs(cosine_D) < EPSILON || Math.Abs(cosine_c + 1.0) < EPSILON) return 0.0;

			// when cosine_c = sin(Lat) == 1 it means that Latitude is 90 degrees and we are in North Pole, so every direction is to south and heading=180
			if (Math.Abs(cosine_c - 1.0) < EPSILON) return Math.PI;

			double cosine_alpha = (cosine_a - cosine_D * cosine_c) / (Math.Sqrt(1.0 - cosine_D * cosine_D) * Math.Sqrt(1.0 - cosine_c * cosine_c));

			// this is possible
			if (Math.Abs(cosine_alpha) > 1.0) cosine_alpha = Math.Truncate(cosine_alpha);

			double result = Math.Acos(cosine_alpha);

			if (Math.Abs(longitude - pLon) > Math.PI)
			{
				if (longitude > 0.0) result = 2.0 * Math.PI - result;
			}
			else if (longitude < pLon) result = 2.0 * Math.PI - result;

			return result;
		}

		/// <summary>
		/// Extension method calculating a distance on a sphere with the given radius between two locations.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <typeparam name="U">Type U must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">First location.</param>
		/// <param name="V">Second location.</param>
		/// <param name="sphereRadius">Radius of the sphere. Default value is the radius of the Earth.</param>
		/// <returns>Distance between two locations.</returns>
		public static double DistanceTo<T, U>(this T cartesian, U V, double sphereRadius = EARTH_MEAN_RADIUS)
			where T : ICartesian
			where U : ICartesian
		{
			return Math.Acos(_Cosine(cartesian.X, cartesian.Y, cartesian.Z, V.X, V.Y, V.Z)) * sphereRadius;
		}

		/// <summary>
		/// Extension method calculating a distance on a sphere with the given radius between two locations.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">First location.</param>
		/// <param name="V">Second location.</param>
		/// <param name="sphereRadius">Radius of the sphere. Default value is the radius of the Earth.</param>
		/// <returns>Distance between two locations.</returns>
		public static double DistanceTo<T>(this T cartesian, IGeoCoordinate V, double sphereRadius = EARTH_MEAN_RADIUS)
			where T : ICartesian
		{
			double cos_lat = Math.Cos(V.Latitude * Math.PI / 180.0);
			double x = cos_lat * Math.Cos(V.Longitude * Math.PI / 180.0);
			double y = cos_lat * Math.Sin(V.Longitude * Math.PI / 180.0);
			double z = Math.Sin(V.Latitude * Math.PI / 180.0);

			return Math.Acos(_Cosine(cartesian.X, cartesian.Y, cartesian.Z, x, y, z)) * sphereRadius;
		}

		/// <summary>
		/// Extension method calculating a distance on a sphere with the given radius between two locations.
		/// </summary>
		/// <typeparam name="T">Type T must implement the IGeoCoordinate interface.</typeparam>
		/// <param name="point">First location.</param>
		/// <param name="latitude">Latitude for the second location in degrees.</param>
		/// <param name="longitude">Longitude for the second location in degrees.</param>
		/// <param name="sphereRadius">Radius of the sphere. Default value is the radius of the Earth.</param>
		/// <returns>Distance between two locations.</returns>
		public static double DistanceTo<T>(this T point, double latitude, double longitude, double sphereRadius = EARTH_MEAN_RADIUS)
			where T : IGeoCoordinate
		{
			double pLat = point.Latitude * Math.PI / 180.0;
			double pLon = point.Longitude * Math.PI / 180.0;

			latitude *= Math.PI / 180.0;
			longitude *= Math.PI / 180.0;

			double cosine = Math.Cos(pLat) * Math.Cos(latitude) * (Math.Cos(pLon - longitude) - 1.0) + Math.Cos(pLat - latitude);
			return sphereRadius * Math.Acos(cosine);
		}

		/// <summary>
		/// Extension method returning a new location on a sphere standing at given distance with specified bearing from the given location.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface and have default constructor.</typeparam>
		/// <typeparam name="U">Type U must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">Initial location.</param>
		/// <param name="distance">Distance from the initial location.</param>
		/// <param name="bearing">Bearing from the initial location in radians.</param>
		/// <param name="sphereRadius">Radius of the sphere. Default value is the radius of the Earth.</param>
		/// <returns>Location standing at given distance and bearing.</returns>
		public static T PointTo<T, U>(this U cartesian, double distance, double bearing, double sphereRadius = EARTH_MEAN_RADIUS)
			where T : ICartesian, new()
			where U : ICartesian
		{
			double lat = cartesian.Dec();
			double lon = cartesian.Ra();

			if (distance < 0.0)
				distance = -distance;

			if (distance <= EPSILON)
				return NewCartesian<T>(cartesian.X, cartesian.Y, cartesian.Z);

			if (Math.Abs(lat - Math.PI / 2.0) <= EPSILON)
				return NewCartesian<T>(0.0, 0.0, 1.0);

			if (Math.Abs(lat + Math.PI / 2.0) <= EPSILON)
				return NewCartesian<T>(0.0, 0.0, -1.0);


			bearing = bearing.To2PIRange();

			double cosine_c = Math.Sin(lat);  //because angle c is pi/2-lat
			double cosine_D = Math.Cos(distance / sphereRadius);

			double cosine_a = cosine_D * cosine_c + Math.Sqrt((1.0 - cosine_D * cosine_D) * (1.0 - cosine_c * cosine_c)) * Math.Cos(bearing);
			lat = Math.PI / 2.0 - Math.Acos(cosine_a);

			double cosine_beta = (cosine_D - cosine_c * cosine_a) / Math.Sqrt((1.0 - cosine_c * cosine_c) * (1.0 - cosine_a * cosine_a));

			if (Math.Abs(cosine_beta) > 1.0)
				cosine_beta = Math.Truncate(cosine_beta);

			double beta = Math.Acos(cosine_beta);

			if (bearing < Math.PI)
				lon += beta;
			else
				lon -= beta;

			if (Math.Abs(lon) > Math.PI)
				if (lon > 0.0)
					lon -= 2.0 * Math.PI;
				else
					lon += 2.0 * Math.PI;

			return NewCartesian<T>(lat, lon);
		}

		/// <summary>
		/// Returns a new location on a sphere standing at given distance with specified bearing from the given location.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface and have default constructor.</typeparam>
		/// <typeparam name="U">Type U must implement the IGeoCoordinate interface.</typeparam>
		/// <param name="location">Initial location.</param>
		/// <param name="distance">Distance from the initial location.</param>
		/// <param name="bearing">Bearing from the initial location in radians.</param>
		/// <param name="sphereRadius">Radius of the sphere. Default value is the radius of the Earth.</param>
		/// <returns>Location standing at given distance and bearing.</returns>
		public static T PointFrom<T, U>(this U location, double distance, double bearing, double sphereRadius = EARTH_MEAN_RADIUS)
			where T : ICartesian, new()
			where U : IGeoCoordinate
		{
			CartesianValue loc = new CartesianValue(location.Latitude, location.Longitude);

			return loc.PointTo<T, CartesianValue>(distance, bearing, sphereRadius);
		}

		/// <summary>
		/// Tests two sections for intersection.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <typeparam name="U">Type U must implement the ICartesian interface.</typeparam>
		/// <param name="V1S1">First section start.</param>
		/// <param name="V2S1">First section end.</param>
		/// <param name="V1S2">Second section start.</param>
		/// <param name="V2S2">Second section end.</param>
		/// <returns>1 if intersect, -1 if not, 0 if one section contain start or end of another section.</returns>
		public static int SectionsIntersect<T, U>(T V1S1, T V2S1, U V1S2, U V2S2)
			where T : ICartesian
			where U : ICartesian
		{
			var (Nx, Ny, Nz) = _CrossProduct(V1S1.X, V1S1.Y, V1S1.Z, V2S1.X, V2S1.Y, V2S1.Z, false);
			double tp1 = _DotProduct(Nx, Ny, Nz, V1S2.X, V1S2.Y, V1S2.Z);
			double tp2 = _DotProduct(Nx, Ny, Nz, V2S2.X, V2S2.Y, V2S2.Z);

			if (tp1 * tp2 <= 0.0)
			{
				(Nx, Ny, Nz) = _CrossProduct(V1S2.X, V1S2.Y, V1S2.Z, V2S2.X, V2S2.Y, V2S2.Z, false);

				tp1 = _DotProduct(Nx, Ny, Nz, V1S1.X, V1S1.Y, V1S1.Z);
				tp2 = _DotProduct(Nx, Ny, Nz, V2S1.X, V2S1.Y, V2S1.Z);

				if (tp1 * tp2 <= 0.0)
				{
					double x1 = V1S1.X + V2S1.X;
					double y1 = V1S1.Y + V2S1.Y;
					double z1 = V1S1.Z + V2S1.Z;
					double x2 = V1S2.X + V2S2.X;
					double y2 = V1S2.Y + V2S2.Y;
					double z2 = V1S2.Z + V2S2.Z;

					if (_Cosine(x1, y1, z1, x2, y2, z2) >= 0.0)
						return 1;
				}
			}

			return -1;
		}

		/// <summary>
		/// Tests two sections for intersection with given tolerance. Sections may not intersect, but to be close to each other with some tolerance.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <typeparam name="U">Type U must implement the ICartesian interface.</typeparam>
		/// <param name="V1S1">First section start.</param>
		/// <param name="V2S1">First section end.</param>
		/// <param name="V1S2">Second section start.</param>
		/// <param name="V2S2">Second section end.</param>
		/// <param name="tolerance">Tolerance given as an angle at the origin.</param>
		/// <returns></returns>
		public static int SectionsIntersect<T, U>(T V1S1, T V2S1, U V1S2, U V2S2, double tolerance)
			where T : ICartesian
			where U : ICartesian
		{
			int sign = SectionsIntersect(V1S1, V2S1, V1S2, V2S2);

			if (sign >= 0 || tolerance <= 0.0)
				return sign;

			double cos_t = Math.Cos(tolerance);
			double res;

			if (((res = V1S1.TestSection(V1S2, V2S2)) >= 0.0 && res < tolerance) ||
				((res = V2S1.TestSection(V1S2, V2S2)) >= 0.0 && res < tolerance) ||
				((res = V1S2.TestSection(V1S1, V2S1)) >= 0.0 && res < tolerance) ||
				((res = V2S2.TestSection(V1S1, V2S1)) >= 0.0 && res < tolerance) ||
				_Cosine(V1S1.X, V1S1.Y, V1S1.Z, V1S2.X, V1S2.Y, V1S2.Z) >= cos_t ||
				_Cosine(V1S1.X, V1S1.Y, V1S1.Z, V2S2.X, V2S2.Y, V2S2.Z) >= cos_t ||
				_Cosine(V2S1.X, V2S1.Y, V2S1.Z, V1S2.X, V1S2.Y, V1S2.Z) >= cos_t ||
				_Cosine(V2S1.X, V2S1.Y, V2S1.Z, V2S2.X, V2S2.Y, V2S2.Z) >= cos_t) sign = 1;

			return sign;
		}

		/// <summary>
		/// Extension method testing if a vector is between two another vectors. 
		/// Returned positive value indicates success and represents a minimum angle from the vector to a geodesical arc on a sphere connecting two another vectors. 
		/// Negative value indicates test failure and its absolute value represents a minimum angle to the geodesical arc.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <typeparam name="U">Type U must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">Vector to test.</param>
		/// <param name="V1">First vector.</param>
		/// <param name="V2">Second vector.</param>
		/// <returns>Minimum angle from the vector to the geodesical arc. Positive when succeeds.</returns>
		public static double TestSection<T, U>(this T cartesian, U V1, U V2)
			where T : ICartesian
			where U : ICartesian
		{
			var (Nx, Ny, Nz) = _CrossProduct(V1.X, V1.Y, V1.Z, V2.X, V2.Y, V2.Z, false);
			double alpha = Math.Acos(_Cosine(cartesian.X, cartesian.Y, cartesian.Z, Nx, Ny, Nz));

			if (alpha > Math.PI / 2.0)
				alpha -= Math.PI / 2.0;
			else
				alpha = Math.PI / 2.0 - alpha;

			var (N1x, N1y, N1z) = _CrossProduct(cartesian.X, cartesian.Y, cartesian.Z, Nx, Ny, Nz, false);
			double l1 = _DotProduct(N1x, N1y, N1z, V1.X, V1.Y, V1.Z);
			double l2 = _DotProduct(N1x, N1y, N1z, V2.X, V2.Y, V2.Z);

			if (l1 * l2 > 0.0 ||
				_DotProduct(cartesian.X, cartesian.Y, cartesian.Z, V1.X, V1.Y, V1.Z) < 0.0 ||
				_DotProduct(cartesian.X, cartesian.Y, cartesian.Z, V2.X, V2.Y, V2.Z) < 0.0) alpha = -alpha - EPSILON;

			return alpha;
		}

		/// <summary>
		/// Returns nearest point on a section.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface and have default constructor.</typeparam>
		/// <typeparam name="U">Type U must implement the IGeoCoordinate interface.</typeparam>
		/// <param name="cartesian">Location.</param>
		/// <param name="V1">Section start.</param>
		/// <param name="V2">Second end.</param>
		/// <returns>Azimuth.</returns>
		public static T SnapTo<T, U>(this T cartesian, U V1, U V2)
			where T : ICartesian, new()
			where U : ICartesian
		{
			var (N1x, N1y, N1z) = _CrossProduct(V1.X, V1.Y, V1.Z, V2.X, V2.Y, V2.Z, false);
			var (N12, N2y, N2z) = _CrossProduct(cartesian.X, cartesian.Y, cartesian.Z, N1x, N1y, N1z, false);
			var (resX, resY, resZ) = _CrossProduct(N1x, N1y, N1z, N12, N2y, N2z, true);

			T result = new T();
			result.SetCartesian(resX, resY, resZ);

			return result;
		}

		/// <summary>
		/// Extension method testing if a point represented by a vector belongs to a polyline. 
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <typeparam name="U">Type U must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">Point to test.</param>
		/// <param name="polyline">Vertices of the polyline.</param>
		/// <param name="tolerance">Tolerance given as an angle at the origin.</param>
		/// <param name="reverse">Iterate vertices in reverse order if true.</param>
		/// <param name="userTest">Additional test delegate. If the point moves and has direction, this is the place to test it.</param>
		/// <param name="angleFromPolyline">Angle representing minimum distance to polyline if test succeeds</param>
		/// <returns>Index of a vertex which is a start of closest polyline segment. -1 if test fails.</returns>
		public static int TestPolyline<T, U>(this T cartesian, IEnumerable<U> polyline,
			double tolerance, bool reverse,
			PolylineSectionTest<U> userTest, out double angleFromPolyline)
			where T : ICartesian
			where U : ICartesian
		{
			if (polyline.Count() < 2)
				throw new ArgumentException(Error_Message_Not_Polyline, nameof(polyline));

			double res, bearing, angle;

			IEnumerable<U> vertices = reverse ? polyline.Reverse().Skip(1) : polyline.Skip(1);
			U previous = reverse ? polyline.Last() : polyline.First();
			U before = default;
			int i = 0;

			int closestIndex = -1;
			double closestAngle = double.MaxValue;

			foreach (var vertex in vertices)
			{
				res = cartesian.TestSection(previous, vertex);

				if (res >= 0.0)
				{
					if (res <= tolerance && res < closestAngle)
					{
						bearing = 0.0;

						if (userTest == null || userTest(default, previous, vertex, res, bearing = previous.BearingTo(vertex)))
						{
							closestIndex = i;
							closestAngle = res;
						}
					}

					before = default;
				}
				else
				{
					if (!EqualityComparer<U>.Default.Equals(before, default) &&
						(angle = cartesian.Angle(previous)) <= tolerance &&
						angle < cartesian.Angle(before) &&
						angle < cartesian.Angle(vertex))
					{
						bearing = 0.0;

						if (userTest == null || userTest(before, previous, vertex, res, bearing = previous.BearingTo(vertex)))
						{
							closestIndex = i;
							closestAngle = angle;
						}
					}

					before = previous;
				}

				previous = vertex;
				i++;
			}

			angleFromPolyline = closestIndex >= 0 ? Math.Abs(closestAngle) : -1.0;

			return closestIndex;
		}

		/// <summary>
		/// Extension method testing if a point represented by a vector belongs to a polyline using results of previous test for better perfomance. 
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <typeparam name="U">Type U must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">Point to test.</param>
		/// <param name="polyline">Vertices of the polyline.</param>
		/// <param name="tolerance">Tolerance given as an angle at the origin.</param>
		/// <param name="reverse">Iterate vertices in reverse order if true.</param>
		/// <param name="userTest">Additional test delegate. If the point moves and has direction, this is the place to test it.</param>
		/// <param name="result">Result of the previous test.</param>
		/// <returns>Returns new or modified previous PolylineTestResult object. Null if first test fails or increased property Fails of returned PolylineTestResult object</returns>
		public static PolylineTestResult TestPolyline<T, U>(this T cartesian, IEnumerable<U> polyline,
			double tolerance, bool reverse,
			PolylineSectionTest<U> userTest,
			PolylineTestResult result)
			where T : ICartesian
			where U : ICartesian
		{
			if (polyline.Count() < 2)
				throw new ArgumentException(Error_Message_Not_Polyline, nameof(polyline));

			double res, bearing, angle;

			IEnumerable<U> vertices = reverse ? polyline.Reverse().Skip(1) : polyline.Skip(1);
			U previous = reverse ? polyline.Last() : polyline.First();
			U before = default;
			int i = 0;

			if (result != null)
			{
				int sectionIndex = result.SectionIndex;
				if (sectionIndex < 0) sectionIndex = 1 - sectionIndex;

				foreach (var vertex in vertices)
				{
					if (i >= sectionIndex)
					{
						res = cartesian.TestSection(previous, vertex);

						if (res >= 0.0)
						{
							if (res <= tolerance)
							{
								bearing = i != sectionIndex ? previous.BearingTo(vertex) : result.SectionDirection;

								if (userTest == null || userTest(default, previous, vertex, res, bearing))
								{
									if (sectionIndex != i)
									{
										result.SectionIndex = reverse ? -1 - i : i;
										result.SectionDirection = bearing;
									}

									result.AngleWithSectionPlane = res;
									result.AngleWithNextVertex = cartesian.Angle(vertex);
									result.Fails = 0;
									return result;
								}
							}

							before = default;
						}
						else
						{
							if (!EqualityComparer<U>.Default.Equals(before, default) &&
								(angle = cartesian.Angle(previous)) <= tolerance &&
								angle < cartesian.Angle(before) &&
								angle < cartesian.Angle(vertex))
							{
								bearing = previous.BearingTo(vertex);

								if (userTest == null || userTest(before, previous, vertex, res, bearing))
								{
									result.AngleWithSectionPlane = angle;
									result.AngleWithNextVertex = cartesian.Angle(vertex); ;
									result.Fails = 0;
									return result;
								}
							}

							before = previous;
						}
					}

					previous = vertex;
					i++;
				}

				result.Fails++;
			}
			else
			{
				U closestVertexAhead = default;
				int closestIndex = -1;
				double closestAngle = double.MaxValue;
				double closestBearing = double.MinValue;

				foreach (var vertex in vertices)
				{
					res = cartesian.TestSection(previous, vertex);

					if (res >= 0.0)
					{
						if (res <= tolerance && res < closestAngle)
						{
							bearing = previous.BearingTo(vertex);

							if (userTest == null || userTest(default, previous, vertex, res, bearing))
							{
								closestIndex = i;
								closestAngle = res;
								closestVertexAhead = vertex;
								closestBearing = bearing;
							}
						}

						before = default;
					}
					else
					{
						if (!EqualityComparer<U>.Default.Equals(before, default) &&
							(angle = cartesian.Angle(previous)) <= tolerance &&
							angle < cartesian.Angle(before) &&
							angle < cartesian.Angle(vertex))
						{
							bearing = previous.BearingTo(vertex);

							if (userTest == null || userTest(before, previous, vertex, res, bearing))
							{
								closestIndex = i;
								closestAngle = angle;
								closestVertexAhead = vertex;
								closestBearing = bearing;
							}
						}

						before = previous;
					}

					previous = vertex;
					i++;
				}


				if (closestIndex >= 0)
				{
					result = PolylineTestResult.Create(polyline,
						reverse,
						reverse ? -1 - closestIndex : closestIndex,
						closestBearing,
						Math.Abs(closestAngle),
						cartesian.Angle(closestVertexAhead)
					);
				}
			}

			return result;
		}

		private static bool RouteSectionTest<T>(T before, T start, T finish, double angle, double bearing, double moveBearing, double bearingTolerance)
			where T : ICartesian
		{
			double bearing2;

			return Math.Abs(moveBearing.BearingsDiff(bearing)) < bearingTolerance ||
				(!EqualityComparer<T>.Default.Equals(before, default) &&
				 (Math.Abs(moveBearing.BearingsDiff(bearing2 = before.BearingTo(start))) < bearingTolerance ||
				 moveBearing.BearingInRange(bearing, bearing2)));
		}

		/// <summary>
		/// 
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <typeparam name="U">Type U must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian"></param>
		/// <param name="route"></param>
		/// <param name="moveABearing"></param>
		/// <param name="tolerance"></param>
		/// <param name="bearingTolerance"></param>
		/// <param name="userTest"></param>
		/// <param name="angleFromRoute"></param>
		/// <returns></returns>
		public static int TestRoute<T, U>(this T cartesian, IEnumerable<U> route,
			double moveABearing, double tolerance, double bearingTolerance,
			RouteSectionTest<U> userTest,
			out double angleFromRoute)
			where T : ICartesian
			where U : ICartesian
		{
			if (route.Count() < 2)
				throw new ArgumentException(Error_Message_Not_Polyline, nameof(route));

			double res, bearing, angle;

			IEnumerable<U> vertices = route.Skip(1);
			U previous = route.First();
			U before = default;
			int i = 0;

			int closestIndex = -1;
			double closestAngle = double.MaxValue;

			if (userTest == null)
				userTest = RouteSectionTest;

			foreach (var vertex in vertices)
			{
				res = cartesian.TestSection(previous, vertex);

				if (res >= 0.0)
				{
					if (res <= tolerance && res < closestAngle)
					{
						bearing = 0.0;

						if (moveABearing < 0.0 ||
							userTest(default, previous, vertex, res, bearing = previous.BearingTo(vertex), moveABearing, bearingTolerance))
						{
							closestIndex = i;
							closestAngle = res;
						}
					}

					before = default;
				}
				else
				{
					if (!EqualityComparer<U>.Default.Equals(before, default) &&
						(angle = cartesian.Angle(previous)) <= tolerance &&
						angle < cartesian.Angle(before) &&
						angle < cartesian.Angle(vertex))
					{
						bearing = 0.0;

						if (moveABearing < 0.0 ||
							userTest(before, previous, vertex, res, bearing = previous.BearingTo(vertex), moveABearing, bearingTolerance))
						{
							closestIndex = i;
							closestAngle = angle;
						}
					}

					before = previous;
				}

				previous = vertex;
				i++;
			}

			angleFromRoute = closestIndex >= 0 ? Math.Abs(closestAngle) : -1.0;

			return closestIndex;
		}

		/// <summary>
		/// 
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <typeparam name="U">Type U must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian"></param>
		/// <param name="route"></param>
		/// <param name="moveBearing"></param>
		/// <param name="tolerance"></param>
		/// <param name="bearingTolerance"></param>
		/// <param name="userTest"></param>
		/// <param name="result"></param>
		/// <returns></returns>
		public static RouteTestResult TestRoute<T, U>(this T cartesian, IEnumerable<U> route,
			double moveBearing, double tolerance, double bearingTolerance,
			RouteSectionTest<U> userTest,
			RouteTestResult result)
			where T : ICartesian
			where U : ICartesian
		{
			if (route.Count() < 2)
				throw new ArgumentException(Error_Message_Not_Polyline, nameof(route));

			double res, bearing, angle;

			IEnumerable<U> vertices = route.Skip(1);
			U previous = route.First();
			U before = default;
			int i = 0;

			if (userTest == null)
				userTest = RouteSectionTest;

			if (result != null)
			{
				int sectionIndex = result.SectionIndex;

				if (sectionIndex < 0)
					sectionIndex = 1 - sectionIndex;

				foreach (var vertex in vertices)
				{
					if (i >= sectionIndex)
					{
						res = cartesian.TestSection(previous, vertex);

						if (res >= 0.0)
						{
							if (res <= tolerance)
							{
								bearing = i != sectionIndex ? previous.BearingTo(vertex) : result.SectionDirection;

								if (moveBearing < 0.0 || userTest(default, previous, vertex, res, bearing, moveBearing, bearingTolerance))
								{
									if (sectionIndex != i)
									{
										result.SectionIndex = i;
										result.SectionDirection = bearing;
									}

									result.AngleWithSectionPlane = res;
									result.AngleWithNextVertex = cartesian.Angle(vertex);
									result.Fails = 0;
									return result;
								}
							}

							before = default;
						}
						else
						{
							if (!EqualityComparer<U>.Default.Equals(before, default) &&
								(angle = cartesian.Angle(previous)) <= tolerance &&
								angle < cartesian.Angle(before) &&
								angle < cartesian.Angle(vertex))
							{
								bearing = previous.BearingTo(vertex);

								if (moveBearing < 0.0 || userTest(before, previous, vertex, res, bearing, moveBearing, bearingTolerance))
								{
									result.AngleWithSectionPlane = angle;
									result.AngleWithNextVertex = cartesian.Angle(vertex); ;
									result.Fails = 0;
									return result;
								}
							}

							before = previous;
						}
					}

					previous = vertex;
					i++;
				}

				result.Fails++;
			}
			else
			{
				U closestVertexAhead = default;
				int closestIndex = -1;
				double closestAngle = double.MaxValue;
				double closestBearing = double.MinValue;

				foreach (var vertex in vertices)
				{
					res = cartesian.TestSection(previous, vertex);

					if (res >= 0.0)
					{
						if (res <= tolerance && res < closestAngle)
						{
							bearing = previous.BearingTo(vertex);

							if (moveBearing < 0.0 || userTest(default, previous, vertex, res, bearing, moveBearing, bearingTolerance))
							{
								closestIndex = i;
								closestAngle = res;
								closestVertexAhead = vertex;
								closestBearing = bearing;
							}
						}

						before = default;
					}
					else
					{
						if (!EqualityComparer<U>.Default.Equals(before, default) &&
							(angle = cartesian.Angle(previous)) <= tolerance &&
							angle < cartesian.Angle(before) &&
							angle < cartesian.Angle(vertex))
						{
							bearing = previous.BearingTo(vertex);

							if (moveBearing < 0.0 || userTest(before, previous, vertex, res, bearing, moveBearing, bearingTolerance))
							{
								closestIndex = i;
								closestAngle = angle;
								closestVertexAhead = vertex;
								closestBearing = bearing;
							}
						}

						before = previous;
					}

					previous = vertex;
					i++;
				}


				if (closestIndex >= 0)
				{
					result = new RouteTestResult()
					{
						SectionIndex = closestIndex,
						SectionDirection = closestBearing,
						AngleWithSectionPlane = Math.Abs(closestAngle),
						AngleWithNextVertex = cartesian.Angle(closestVertexAhead)
					};
				}
			}

			return result;
		}


		// IGeoCoordinate extension methods


		/// <summary>
		/// Extension method checking if a location with given latitude and longitude is inside an area with given South-West and North-East corners.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="lat">Latitude of the location to test in radians.</param>
		/// <param name="lon">Longitude of the location to test in radians.</param>
		/// <param name="SW">Area South-West corner.</param>
		/// <param name="NE">Area North-East corner.</param>
		/// <param name="tolerance">Tolerance in radians representing distance on the sphere.</param>
		/// <returns>True when location is inside.</returns>
		public static bool WithinBounds<T>(double lat, double lon, T SW, T NE, double tolerance)
			where T : ICartesian
		{
			double swLat = SW.Dec();
			double swLon = SW.Ra();
			double neLat = NE.Dec();
			double neLon = NE.Ra();

			double toleranceRad = Math.Abs(tolerance) / EARTH_MEAN_RADIUS;
			double nLat = neLat + toleranceRad;

			if (nLat > Math.PI / 2.0)
				nLat = Math.PI / 2.0;

			if (lat > nLat)
				return false;

			double sLat = swLat - toleranceRad;

			if (sLat < -Math.PI / 2.0)
				sLat = -Math.PI / 2.0;

			if (lat < sLat)
				return false;

			if (neLon - swLon >= 2.0 * Math.PI - 2.0 * toleranceRad)
				return true;

			double eLon = (neLon + toleranceRad).ToLongitudeRange();
			double wLon = (swLon - toleranceRad).ToLongitudeRange();

			return eLon < wLon ?
				(lon >= wLon && lon <= Math.PI) || (lon >= -Math.PI && lon <= eLon) :
				lon >= wLon && lon <= eLon;
		}

		/// <summary>
		/// Extension method checking if a location is inside an area with given South-West and North-East corners.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <typeparam name="U">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="location">Location to test.</param>
		/// <param name="SW">Area South-West corner.</param>
		/// <param name="NE">Area North-East corner.</param>
		/// <param name="tolerance">Tolerance in radians representing distance on the sphere.</param>
		/// <returns>True when location is inside.</returns>
		public static bool WithinBounds<T, U>(this T location, U SW, U NE, double tolerance)
			where T : ICartesian
			where U : ICartesian
		{
			return WithinBounds(location.Dec(), location.Ra(), SW, NE, tolerance);
		}

		/// <summary>
		/// Extension method returning a key of a mercator map tile at given zoom contaning given location.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="coordinate">Location.</param>
		/// <param name="zoom">Map zoom level.</param>
		/// <returns>Tile key.</returns>
		public static long MercatorMapTileKey<T>(this T coordinate, int zoom)
			where T : ICartesian
		{
			long res;
			double cLat = coordinate.Dec();

			if (Math.Abs(cLat) < MaxMercatorMapsLat)
			{
				int N = 1 << zoom;
				double sinY = Math.Sin(cLat);

				int tileX = (int)Math.Floor(N * (0.5 + coordinate.Ra().ToLongitudeRange() / 2.0 / Math.PI));
				int tileY = (int)Math.Floor(N * (0.5 - Math.Log((1.0 + sinY) / (1.0 - sinY)) / 4.0 / Math.PI));

				res = ((long)tileX << 32) + tileY;

				if (res < 0)
					throw new InvalidOperationException();
			}
			else res = -1;

			return res;
		}


		// General purpose methods

		public static (double x, double y, double z) SphericalToCartesian(double dec, double ra)
		{
			double cos_dec = Math.Cos(dec);
			return (cos_dec * Math.Cos(ra), cos_dec * Math.Sin(ra), Math.Sin(dec));
		}

		public static (double dec, double ra) CartesianToSpherical(double x, double y, double z)
		{
			return (Math.Asin(z / Math.Sqrt(x * x + y * y + z * z)), Math.Atan2(y, x));
		}

		public static CartesianValue AsCartesian<T>(T location)
			where T : IGeoCoordinate
		{
			return new CartesianValue(location.Latitude, location.Longitude);
		}

		/// <summary>
		/// Returns random unit-vector.
		/// </summary>
		/// <typeparam name="T">Must implement ICartesian interface.</typeparam>
		/// <returns>Random unit-vector.</returns>
		public static T RandomVector<T>()
			where T : ICartesian, new()
		{
			System.Random rnd = new System.Random();

			double x = 1.0 - rnd.NextDouble();
			double y = 1.0 - rnd.NextDouble();
			double z = 1.0 - rnd.NextDouble();

			(x, y, z) = _Normalized(x, y, z);

			T result = new T();
			result.SetCartesian(x, y, z);

			return result;
		}

		/// <summary>
		/// Extension method "normalizing" a bearing to the [0..360) range.
		/// </summary>
		/// <param name="bearing">Bearing in degrees.</param>
		/// <returns>Bearing within [0..360) range.</returns>
		public static double To360Range(this double bearing)
		{
			if (bearing >= 360.0 || bearing < 0.0)
				bearing -= Math.Floor(bearing / 360.0) * 360.0;

			return bearing;
		}

		/// <summary>
		/// Extension method "normalizing" a bearing to the [0..2π) range.
		/// </summary>
		/// <param name="bearing">Bearing in radians.</param>
		/// <returns>Bearing in the [0..2π) range.</returns>
		public static double To2PIRange(this double bearing)
		{
			if (bearing >= 2.0 * Math.PI || bearing < 0.0)
				bearing -= Math.Floor(bearing / (2.0 * Math.PI)) * 2.0 * Math.PI;

			return bearing;
		}

		/// <summary>
		/// Extension method "normalizing" a longitude in radians to the [-π..π] range. 
		/// For example, longitude 182 degrees is -178 degrees.
		/// </summary>
		/// <param name="lon">Longitude in radians.</param>
		/// <returns>Longitude in the [-π..π] range.</returns>
		public static double ToLongitudeRange(this double lon)
		{
			if (lon > Math.PI || lon < -Math.PI)
				lon -= Math.Floor((lon + Math.PI) / (2.0 * Math.PI)) * 2.0 * Math.PI;

			return lon;
		}

		/// <summary>
		/// Extension method "normalizing" a longitude in degrees to the [-180..180] range. 
		/// For example, longitude 182 degrees is -178 degrees.
		/// </summary>
		/// <param name="longitude">Longitude in degrees.</param>
		/// <returns>Longitude in the [-180..180] range.</returns>
		public static double ToLongitudeRangeDegrees(this double longitude)
		{
			if (longitude > 180.0 || longitude < -180.0)
				longitude -= Math.Floor((longitude + 180.0) / (2.0 * 180.0)) * 2.0 * 180.0;

			return longitude;
		}

		/// <summary>
		/// Extension method "normalizing" a latitude in radians to the [-π/2..π/2] range. 
		/// For example, latitude 92 degrees is 88 degrees.
		/// </summary>
		/// <param name="lat">Latitude in radians.</param>
		/// <returns>Latitude in the [-π/2..π/2] range.</returns>
		public static double ToLatitudeRange(this double lat)
		{
			double res = lat;

			if (res < -Math.PI / 2)
				do { res = -Math.PI - lat; } while (res < -Math.PI / 2);
			else if (res > Math.PI / 2)
				do { res = Math.PI - lat; } while (res >= Math.PI / 2);

			return res;
		}

		/// <summary>
		/// Extension method "normalizing" a latitude in degrees to the [-90..90] range. 
		/// For example, latitude 92 degrees is 88 degrees.
		/// </summary>
		/// <param name="latitude">Latitude in degrees.</param>
		/// <returns>Latitude in the [-90..90] range.</returns>
		public static double ToLatitudeRangeDegrees(this double latitude)
		{
			double res = latitude;

			if (res < -90.0)
				do { res = -180.0 - latitude; } while (res < -90.0);
			else if (res > 90.0)
				do { res = 180.0 - latitude; } while (res >= 90.0);

			return res;
		}

		/// <summary>
		/// Extension method returning an opposite bearing for a given bearing.
		/// </summary>
		/// <param name="bearing">Bearing in radians or degrees depending on the degrees argument.</param>
		/// <param name="degrees">Bearing in degrees when true.</param>
		/// <returns>Opposite bearing.</returns>
		public static double OppositeBearing(this double bearing, bool degrees = false)
		{
			if (degrees)
				if ((bearing = bearing.To360Range()) < 180.0)
					bearing += 180.0;
				else
					bearing -= 180.0;
			else if ((bearing = bearing.To2PIRange()) < Math.PI)
				bearing += Math.PI;
			else
				bearing -= Math.PI;

			return bearing;
		}

		/// <summary>
		/// Extension method returning a difference of two bearings.
		/// </summary>
		/// <param name="bearing">First bearing in radians or degrees depending on the degrees argument.</param>
		/// <param name="otherBearing">Second bearing in radians or degrees depending on the degrees argument.</param>
		/// <param name="degrees">Bearings in degrees when true.</param>
		/// <returns>Difference of two bearings.</returns>
		public static double BearingsDiff(this double bearing, double otherBearing, bool degrees = false)
		{
			otherBearing -= bearing;

			if (degrees)
			{
				if (otherBearing > 180.0)
					otherBearing -= 360.0;
				else if (otherBearing < -180.0)
					otherBearing += 360.0;
			}
			else
			{
				if (otherBearing > Math.PI)
					otherBearing -= 2 * Math.PI;
				else if (otherBearing < -Math.PI)
					otherBearing += 2 * Math.PI;
			}

			return otherBearing;
		}

		/// <summary>
		/// Extension method returning a difference of two bearings.
		/// </summary>
		/// <param name="bearing">First bearing in radians or degrees depending on the degrees argument.</param>
		/// <param name="otherBearing">Second bearing in radians or degrees depending on the degrees argument.</param>
		/// <param name="degrees">Bearings in degrees when true.</param>
		/// <returns>Difference of two bearings.</returns>
		public static float BearingsDiff(this float bearing, float otherBearing, bool degrees = false)
		{
			otherBearing -= bearing;

			if (degrees)
			{
				if (otherBearing > 180.0f)
					otherBearing -= 360.0f;
				else if (otherBearing < -180.0f)
					otherBearing += 360.0f;
			}
			else
			{
				if (otherBearing > Math.PI)
					otherBearing -= (float)(2.0 * Math.PI);
				else if (otherBearing < -Math.PI)
					otherBearing += (float)(2.0 * Math.PI);
			}

			return otherBearing;
		}

		/// <summary>
		/// Extension method testing if a bearing is within given range. The range can't be more than π. The end of the range can be larger than the start. 
		/// </summary>
		/// <param name="bearing">Bearing to test.</param>
		/// <param name="h1">Start value of the range.</param>
		/// <param name="h2">End value of the range.</param>
		/// <returns>True if within the range.</returns>
		public static bool BearingInRange(this double bearing, double h1, double h2)
		{
			if (Math.Abs(h1 - h2) <= Math.PI)
				return h2 > h1 ? bearing >= h1 && bearing <= h2 : bearing >= h2 && bearing <= h1;

			if (h2 > h1)
				return (bearing >= h2 && bearing < 2.0 * Math.PI) || (bearing >= 0.0 && bearing <= h1);

			return (bearing >= h1 && bearing < 2.0 * Math.PI) || (bearing >= 0.0 && bearing <= h2);
		}

		/// <summary>
		/// Returns South-West and North-East corners of an area that encloses given locations.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface and have default constructor.</typeparam>
		/// <typeparam name="U">Type U must implement the ICartesian interface.</typeparam>
		/// <param name="locations">Locations.</param>
		/// <param name="SW">Out parameter with the area South-West corner.</param>
		/// <param name="NE">Out parameter with the area North-East corner.</param>
		public static void GetBounds<T, U>(IEnumerable<U> locations, out T SW, out T NE)
			where T : ICartesian, new()
			where U : ICartesian
		{
			if (!locations.Any())
				throw new ArgumentException(Error_Message_Empty_Collection, nameof(locations));

			double minLat = double.MaxValue;
			double minLon = double.MaxValue;
			double maxLat = double.MinValue;
			double maxLon = double.MinValue;
			double lat;
			double lon;

			foreach (U location in locations)
			{
				lat = location.Dec();
				lon = location.Ra();

				if (lat > maxLat) maxLat = lat;
				if (lat < minLat) minLat = lat;
				if (lon > maxLon) maxLon = lon;
				if (lon < minLon) minLon = lon;
			}

			if (maxLon - minLon > Math.PI * 2)
			{
				SW = NewCartesian<T>(minLat, maxLon);
				NE = NewCartesian<T>(maxLat, minLon);
			}
			else
			{
				SW = NewCartesian<T>(minLat, minLon);
				NE = NewCartesian<T>(maxLat, maxLon);
			}
		}

		/// <summary>
		/// Returns length of a polyline on a sphere with given radius.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="locations">Polyline points.</param>
		/// <param name="sphereRadius">Radius of the sphere in meters. Default value is the radius of the Earth.</param>
		/// <returns>Length of the polyline in meters.</returns>
		public static double GetLength<T>(IEnumerable<T> locations, double sphereRadius = EARTH_MEAN_RADIUS)
			where T : ICartesian
		{
			if (locations.Count() < 2)
				throw new ArgumentException(Error_Message_Not_Polyline, nameof(locations));

			double result = 0.0;
			T previous = locations.First();

			foreach (T loc in locations.Skip(1))
			{
				result += previous.DistanceTo(loc, sphereRadius);
				previous = loc;
			}

			return result;
		}

		/// <summary>
		/// Returns a longitude of a mercator map tile by its horizontal position.
		/// </summary>
		/// <param name="tileX">Horizontal position of the tile.</param>
		/// <param name="zoom">Zoom level.</param>
		/// <returns>Latitude of tile in radians.</returns>
		public static double MercatorTileLon(int tileX, int zoom)
		{
			return 2.0 * Math.PI * tileX / Math.Pow(2, zoom) - Math.PI;
		}

		/// <summary>
		/// Returns a latitude of a mercator map tile by its vertical position.
		/// </summary>
		/// <param name="tileY">Vertical position of the tile.</param>
		/// <param name="zoom">Zoom level.</param>
		/// <returns>Longitude of the tile in radians.</returns>
		public static double MercatorTileLat(int tileY, int zoom)
		{
			return Math.Atan(Math.Sinh(Math.PI - 2.0 * Math.PI * tileY / Math.Pow(2, zoom)));
		}

		/// <summary>
		/// Returns bounds of a mercator map tile by its horizontal and vertical positions.
		/// </summary>
		/// <param name="tileX">Horizontal position of the tile.</param>
		/// <param name="tileY">Vertical position of the tile.</param>
		/// <param name="zoom">Zoom level.</param>
		/// <returns>Bounds of the tile as an array {lon1, lat1, lon2, lat2}. All values in radians.</returns>
		public static double[] MercatorTileBounds(int tileX, int tileY, int zoom)
		{
			double N = Math.Pow(2, zoom);
			double lon1 = 2.0 * Math.PI * tileX / N - Math.PI;
			double lat1 = Math.Atan(Math.Sinh(Math.PI - 2.0 * Math.PI * tileY / N));

			int maxTileXY = (1 << zoom) - 1;

			double lon2 = tileX < maxTileXY ?
				2.0 * Math.PI * (tileX + 1) / N - Math.PI :
				Math.PI;

			double lat2 = tileY < maxTileXY ?
				Math.Atan(Math.Sinh(Math.PI - 2.0 * Math.PI * (tileY + 1) / N)) :
				-MaxMercatorMapsLat;

			return new double[] { lon1, lat1, lon2, lat2 };
		}


		// Google polyline

		/// <summary>
		/// Auxiliary method encoding the next point in the polyline. Thank you, stackoverflow.com
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="location">Current point to encode.</param>
		/// <param name="previous">Previous encoded point.</param>
		/// <param name="escape">Add second back-slash when true.</param>
		/// <param name="result">Where result is accumulated.</param>
		/// <returns>Encoded point in the polyline.</returns>
		private static void EncodeToGooglePolyline<T>(this T location, T previous, bool escape, StringBuilder result)
			where T : ICartesian
		{
			char c;

			bool prevIsDefault = EqualityComparer<T>.Default.Equals(previous, default);

			int lastLat = !prevIsDefault ? (int)(previous.Latitude() * 100000.0) : 0;
			int lastLon = !prevIsDefault ? (int)(previous.Longitude() * 100000.0) : 0;

			int diff = (int)(location.Latitude() * 100000.0) - lastLat;
			int shifted = diff << 1;

			if (diff < 0)
				shifted = ~shifted;

			int rem = shifted;

			while (rem >= 0x20)
			{
				c = (char)((0x20 | (rem & 0x1f)) + 63);
				result.Append(c);

				if (c == '\\' && escape)
					result.Append('\\');

				rem >>= 5;
			}

			c = (char)(rem + 63);
			result.Append(c);

			if (c == '\\' && escape)
				result.Append('\\');

			diff = (int)(location.Longitude() * 100000.0) - lastLon;
			shifted = diff << 1;

			if (diff < 0) shifted = ~shifted;

			rem = shifted;

			while (rem >= 0x20)
			{
				c = (char)((0x20 | (rem & 0x1f)) + 63);
				result.Append(c);

				if (c == '\\' && escape)
					result.Append('\\');

				rem >>= 5;
			}

			c = (char)(rem + 63);
			result.Append(c);

			if (c == '\\' && escape)
				result.Append('\\');
		}

		/// <summary>
		/// Extension method encoding a polyline into the string.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="locations">Polyline points.</param>
		/// <param name="escape">Add second back-slash when true.</param>
		/// <returns>Encoded polyline.</returns>
		public static string EncodeToGooglePolyline<T>(IEnumerable<T> locations, bool escape)
			where T : ICartesian
		{
			if (locations == null || !locations.Any())
				return null;

			StringBuilder result = new StringBuilder();
			T previous = default;

			foreach (T location in locations)
			{
				location.EncodeToGooglePolyline(previous, escape, result);
				previous = location;
			}

			return result.ToString();
		}

		/// <summary>
		/// Decodes a polyline from the string.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface and have default constructor.</typeparam>
		/// <param name="encodedPoints">Encoded polyline.</param>
		/// <returns>Points of the polyline.</returns>
		public static List<T> DecodeGooglePolyline<T>(string encodedPoints)
			where T : ICartesian, new()
		{
			if (string.IsNullOrEmpty(encodedPoints))
				return null;

			List<T> result = new List<T>();

			return DecodeGooglePolyline<T>(encodedPoints, result);
		}

		/// <summary>
		/// Decodes a polyline from the string.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface and have default constructor.</typeparam>
		/// <param name="encodedPoints">Encoded polyline.</param>
		/// <param name="result">Contains points of the polyline.</param>
		/// <returns>Points of the polyline.</returns>
		public static List<T> DecodeGooglePolyline<T>(string encodedPoints, List<T> result)
			where T : ICartesian, new()
		{
			if (string.IsNullOrEmpty(encodedPoints))
				return null;

			if (result == null)
				result = new List<T>();
			else
				result.Clear();

			T vertex;
			char[] polylineChars = encodedPoints.ToCharArray();

			int index = 0;
			int currentLat = 0;
			int currentLng = 0;
			int next5bits;
			int sum;
			int shifter;

			while (index < polylineChars.Length)
			{
				// calculate next latitude
				sum = 0;
				shifter = 0;

				do
				{
					next5bits = (int)polylineChars[index++] - 63;
					sum |= (next5bits & 31) << shifter;
					shifter += 5;
				} while (next5bits >= 32 && index < polylineChars.Length);

				if (index >= polylineChars.Length)
					break;

				currentLat += (sum & 1) == 1 ? ~(sum >> 1) : (sum >> 1);

				//calculate next longitude
				sum = 0;
				shifter = 0;

				do
				{
					next5bits = (int)polylineChars[index++] - 63;
					sum |= (next5bits & 31) << shifter;
					shifter += 5;
				} while (next5bits >= 32 && index < polylineChars.Length);

				if (index >= polylineChars.Length && next5bits >= 32)
					break;

				currentLng += (sum & 1) == 1 ? ~(sum >> 1) : (sum >> 1);

				var (x, y, z) = SphericalToCartesian((Convert.ToDouble(currentLat) / 100000.0) * Math.PI / 180.0, (Convert.ToDouble(currentLng) / 100000.0) * Math.PI / 180.0);
				vertex = new T();
				vertex.SetCartesian(x, y, z);
				result.Add(vertex);
			}

			return result;
		}


		// Sunrise & sunset

		/// <summary>
		/// Extension method returning a Julian day since 2000 for a given date. 
		/// </summary>
		/// <param name="date">Date.</param>
		/// <returns>Julian day since 2000.</returns>
		public static double DaysSinceJ2000(this DateTime date)
		{
			return 367.0 * date.Year -
				Math.Truncate(7.0 * (date.Year + Math.Truncate((date.Month + 9.0) / 12.0)) / 4.0) +
				Math.Truncate(275.0 * date.Month / 9.0) +
				date.Day - 730531.5;
		}

		/// <summary>
		/// Extension method returning a Julian day since 2000 for a given date. 
		/// </summary>
		/// <param name="date">Date.</param>
		/// <returns>Julian day since 2000.</returns>
		public static double DaysSinceJ2000(this DateTimeOffset date)
		{
			return 367.0 * date.Year -
				Math.Truncate(7.0 * (date.Year + Math.Truncate((date.Month + 9.0) / 12.0)) / 4.0) +
				Math.Truncate(275.0 * date.Month / 9.0) +
				date.Day - 730531.5;
		}

		/// <summary>
		/// Extension method returning sunrise time for a given location and date with offset. 
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="location">Location to calculate sunrise time for.</param>
		/// <param name="date">Date and offset.</param>
		/// <returns>Sunrise time.</returns>
		public static DateTime Sunrise<T>(this T location, DateTimeOffset date)
			where T : ICartesian
		{
			double lat = location.Dec();
			double JDays = date.DaysSinceJ2000();
			double JCenturies = JDays / 36525.0;
			double elevation = -5.0 / 6.0; // - 6 degrees for civil twilight, 12 - nautical, 18 - astronomical

			double L = (4.8949504201433 + 628.331969753199 * JCenturies) % (2.0 * Math.PI);
			double G = (6.2400408 + 628.3019501 * JCenturies) % (2.0 * Math.PI);
			double ec = 0.033423 * Math.Sin(G) + 0.00034907 * Math.Sin(2.0 * G);
			double lambda = L + ec;
			double GHA = 0.0430398 * Math.Sin(2 * lambda) - 0.00092502 * Math.Sin(4 * lambda) - ec;
			double obl = 0.409093 - 0.0002269 * JCenturies;
			double delta = Math.Asin(Math.Sin(obl) * Math.Sin(lambda));
			double cosc = (Math.Sin(elevation * Math.PI / 180.0) - Math.Sin(lat) * Math.Sin(delta)) / (Math.Cos(lat) * Math.Cos(delta));

			if (cosc > 1.0)
				return DateTime.MaxValue; // never rises 

			if (cosc < -1.0)
				return DateTime.MinValue; // never sets

			double correction = Math.Acos(cosc);
			double utnew = Math.PI - (GHA + location.Ra() + correction);
			double eventUT = utnew * 57.29577951 / 15.0;

			if (eventUT >= 24.0)
				eventUT %= 24.0;

			DateTime result = new DateTime(date.Year, date.Month, date.Day);
			return result.Add(date.Offset).AddHours(eventUT);
		}

		/// <summary>
		/// Extension method returning sunset time for a given location and date with offset.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="location">Location to calculate sunset time for.</param>
		/// <param name="date">Date and offset.</param>
		/// <returns>Sunset time.</returns>
		public static DateTime Sunset<T>(this T location, DateTimeOffset date)
			where T : ICartesian
		{
			double lat = location.Dec();
			double JDays = date.DaysSinceJ2000();
			double JCenturies = JDays / 36525.0;
			double elevation = -5.0 / 6.0; // - 6 degrees for civil twilight, 12 - nautical, 18 - astronomical

			double L = (4.8949504201433 + 628.331969753199 * JCenturies) % (2.0 * Math.PI);
			double G = (6.2400408 + 628.3019501 * JCenturies) % (2.0 * Math.PI);
			double ec = 0.033423 * Math.Sin(G) + 0.00034907 * Math.Sin(2.0 * G);
			double lambda = L + ec;
			double GHA = 0.0430398 * Math.Sin(2 * lambda) - 0.00092502 * Math.Sin(4 * lambda) - ec;
			double obl = 0.409093 - 0.0002269 * JCenturies;
			double delta = Math.Asin(Math.Sin(obl) * Math.Sin(lambda));
			double cosc = (Math.Sin(elevation * Math.PI / 180.0) - Math.Sin(lat) * Math.Sin(delta)) / (Math.Cos(lat) * Math.Cos(delta));  // B25

			if (cosc > 1.0)
				return DateTime.MaxValue; // never rises 

			if (cosc < -1.0)
				return DateTime.MinValue; // never sets

			double correction = Math.Acos(cosc);
			double utnew = Math.PI - (GHA + location.Ra() - correction);
			double eventUT = utnew * 57.29577951 / 15.0;

			if (eventUT >= 24.0)
				eventUT %= 24.0;

			DateTime result = new DateTime(date.Year, date.Month, date.Day);
			return result.Add(date.Offset).AddHours(eventUT);
		}

		/// <summary>
		/// Returns times of nearest sunrise and sunset for a given location and date/time with the offset. 
		/// Nearest sunrise can be for the next day, nearest sunset can be for the previous day.
		/// Polar night at given location at given date, if method returns (DateTime.MaxValue, DateTime.MinValue).
		/// Polar day, if method returns (DateTime.MinValue, DateTime.MaxValue).
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="location">Location.</param>
		/// <param name="date">Date, time and offset.</param>
		/// <returns>Times of sunrise and sunset.</returns>
		public static (DateTime sunrise, DateTime sunset) NearestSunriseSunset<T>(this T location, DateTimeOffset date)
			where T : ICartesian
		{
			// Install System.ValueTuple package via Nuget if this fails to compile.
			// PM > Install - Package "System.ValueTuple"
			DateTime sunrise, sunset;

			if (location != null)
			{
				sunrise = location.Sunrise(date);

				if (sunrise == DateTime.MaxValue)
					return (DateTime.MaxValue, DateTime.MinValue); // never rises

				if (sunrise == DateTime.MinValue)
					return (DateTime.MinValue, DateTime.MaxValue); // never sets

				sunset = location.Sunset(date);

				if (date < sunrise)
				{
					DateTimeOffset dayBefore = date.AddDays(-1);
					sunset = location.Sunset(dayBefore);
				}
				else if (date > sunset)
				{
					DateTimeOffset dayAfter = date.AddDays(1);
					sunrise = location.Sunrise(dayAfter);
				}
			}
			else
			{
				sunrise = new DateTime(date.Year, date.Month, date.Day, 7, 0, 0);
				sunset = new DateTime(date.Year, date.Month, date.Day, 19, 0, 0);

				if (date < sunrise)
				{
					DateTimeOffset dayBefore = date.AddDays(-1);
					sunset = new DateTime(dayBefore.Year, dayBefore.Month, dayBefore.Day, 19, 0, 0);
				}
				else if (date > sunset)
				{
					DateTimeOffset dayAfter = date.AddDays(1);
					sunrise = new DateTime(dayAfter.Year, dayAfter.Month, dayAfter.Day, 7, 0, 0);
				}
			}

			return (sunrise, sunset);
		}

		/// <summary>
		/// Extension method testing if the Sun is in the sky for a given location and date/time.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="location">Location to test.</param>
		/// <param name="date">Date, time, and offset for the test.</param>
		/// <returns>True if the Sun is in the sky.</returns>
		public static bool DayTime<T>(this T location, DateTimeOffset date)
			where T : ICartesian
		{
			var (sunrise, sunset) = location.NearestSunriseSunset(date);

			return date.DateTime > sunrise && date.DateTime < sunset;
		}
	}

}