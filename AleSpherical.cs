/* (c) AleProjects.com, 2018
 * v.1.0
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
using System.Globalization;
using System.Linq;
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
		/// Latitude in radians.
		/// </summary>
		double Lat { get; }

		/// <summary>
		/// Longitude in radians.
		/// </summary>
		double Lon { get; }

		/// <summary>
		/// Sets latitude and longitude in radians for the coordinate.
		/// </summary>
		/// <param name="lat">Latitude in radians.</param>
		/// <param name="lon">Longitude in radians.</param>
		void SetSpherical(double lat, double lon);
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
			get => Math.Asin(Z / Length);
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
			get => Math.Asin(Z / this.VectorLength()) * 180.0 / Math.PI; 
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
			get => Math.Sqrt(X * X + Y * Y + Z * Z);
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
			double len;

			if (normalize && (len = Math.Sqrt(x * x + y * y + z * z)) > 0.0)
			{
				X = x / len;
				Y = y / len;
				Z = z / len;
			}
			else
			{
				X = x;
				Y = y;
				Z = z;
			}
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
		public Cartesian(ICartesian cartesian) : this(cartesian.X, cartesian.Y, cartesian.Z, false) { }

		/// <summary>
		/// Creates a new vector from coordinates on a sphere.
		/// </summary>
		/// <param name="geoCoordinate">Coordinates on a sphere</param>
		public Cartesian(IGeoCoordinate geoCoordinate) : this(geoCoordinate.Lat, geoCoordinate.Lon, false) { }


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
		/// Sets declination and right ascension for this coordinate.
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
	/// Implementation of a basic location on a sphere with Latitude and Longitude.
	/// </summary>
	public class GeoCoordinate : Cartesian, IGeoCoordinate
	{
		private double _Lat, _Lon;


		/// <summary>
		/// Latitude in radians.
		/// </summary>
		public double Lat
		{
			get => _Lat;
			protected set => SetSpherical(value, _Lon);
		}

		/// <summary>
		/// Longitude in radians.
		/// </summary>
		public double Lon
		{
			get => _Lon;
			protected set => SetSpherical(_Lat, value);
		}

		/// <summary>
		/// Latitude in degrees.
		/// </summary>
		public double Latitude
		{
			get => _Lat * 180.0 / Math.PI;
		}

		/// <summary>
		/// Longitude in degrees.
		/// </summary>
		public double Longitude
		{
			get => _Lon * 180.0 / Math.PI;
		}

		/// <summary>
		/// Altitude in meters.
		/// </summary>
		public double Altitude { get; set; }


		/// <summary>
		/// Creates a new geocoordinate with Latitude=0 Longitude=0.
		/// </summary>
		public GeoCoordinate() : base()
		{
			_Lat = 0.0;
			_Lon = 0.0;
		}

		/// <summary>
		/// Creates a new geocoordinate with the specified latitude and longitude.
		/// </summary>
		/// <param name="lat">Latitude.</param>
		/// <param name="lon">Longitude.</param>
		/// <param name="degrees">When true, indicates that the lon and lat parameters are in degrees.</param>
		public GeoCoordinate(double lat, double lon, bool degrees = true) : base(lat, lon, degrees)
		{
			if (degrees)
			{
				lat *= Math.PI / 180.0;
				lon *= Math.PI / 180.0;
			}

			_Lat = lat;
			_Lon = lon;
		}

		/// <summary>
		/// Creates a new geocoordinate as a copy of the given location.
		/// </summary>
		/// <param name="location">Geocoordinate whose Latitude and Longitude used for a new geocoordinate.</param>
		public GeoCoordinate(IGeoCoordinate location) : base(location.Lat, location.Lon, false)
		{
			_Lat = location.Lat;
			_Lon = location.Lon;
		}

		/// <summary>
		/// Creates a new geocoordinate from the vector.
		/// </summary>
		/// <param name="cartesian">The vector whose Dec and Ra used as Latitude and Longitude for a new geocoordinate.</param>
		public GeoCoordinate(ICartesian cartesian) : base(cartesian)
		{
			_Lat = cartesian.Dec();
			_Lon = cartesian.Ra();
		}


		/// <summary>
		/// Sets x,y,z components for the vector.
		/// </summary>
		/// <param name="x">x component.</param>
		/// <param name="y">y component.</param>
		/// <param name="z">z component.</param>
		public override void SetCartesian(double x, double y, double z)
		{
			double L = Math.Sqrt(x * x + y * y + z * z);
			base.SetCartesian(x / L, y / L, z / L);
			_Lat = Math.Asin(Z);
			_Lon = Math.Atan2(Y, X);
		}

		/// <summary>
		/// Sets latitude and longitude given in radians for the coordinate.
		/// </summary>
		/// <param name="dec">Latitude in radians.</param>
		/// <param name="ra">Longitude in radians.</param>
		public override void SetSpherical(double lat, double lon)
		{
			base.SetSpherical(lat, lon);
			_Lat = lat;
			_Lon = lon;
		}


		/// <summary>
		/// Returns a string representation of the geocoordinate.
		/// </summary>
		/// <returns></returns>
		public override string ToString()
		{
			IFormatProvider fmt = new NumberFormatInfo();
			return Latitude.ToString(fmt) + "," + Longitude.ToString(fmt);
		}
	}



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

		public PolylineTestResult(IEnumerable<ICartesian> polyline, bool reverse)
		{
			if (polyline == null) throw new ArgumentNullException();

			SectionsAngles = new double[polyline.Count() - 1];
			IEnumerable<ICartesian> vertices = reverse ? polyline.Reverse().Skip(1) : polyline.Skip(1);
			ICartesian previous = reverse ? polyline.Last() : polyline.First();
			int i = 0;

			foreach (var vertex in vertices)
			{
				SectionsAngles[i++] = previous.Angle(vertex);
				previous = vertex;
			}
		}
	}



	/// <summary>
	/// Provides extension methods for ICartesian and IGeoCoordinate objects.
	/// </summary>
	public static class SphericalExtension
	{
		public const double EARTH_MEAN_RADIUS = 6371000;
		public const double EPSILON = 1.0e-15;
		public const double EPSILON_GEO = 1.0e-8;


		private struct CartesianStruct : ICartesian
		{
			public double X { get; set; }
			public double Y { get; set; }
			public double Z { get; set; }

			public void SetCartesian(double x, double y, double z)
			{
				X = x;
				Y = y;
				Z = z;
			}
		}

		// properties

		//public static double Precision { get; set; } = EPSILON;

		/// <summary>
		/// Maximum latitude in radians a map based on Mercator projection can display.
		/// </summary>
		public static double MaxMercatorMapsLat { get => Math.Atan(Math.Sinh(Math.PI)); }


		// Auxiliary methods

		private static T NewCartesianDerivative<T>(double x, double y, double z) where T : ICartesian, new()
		{
			T result = new T();
			result.SetCartesian(x, y, z);

			return result;
		}

		private static T NewGeoCoordinateDerivative<T>(double dec, double ra) where T : IGeoCoordinate, new()
		{
			T result = new T();
			result.SetSpherical(dec, ra);

			return result;
		}


		// ICartesian extension methods

		/// <summary>
		/// Extension method calculating a length of the vector.
		/// </summary>
		/// <param name="cartesian">Vector.</param>
		/// <returns>Length of the vector.</returns>
		public static double VectorLength(this ICartesian cartesian)
		{
			return Math.Sqrt(cartesian.X * cartesian.X + cartesian.Y * cartesian.Y + cartesian.Z * cartesian.Z);
		}

		/// <summary>
		/// Extension method calculating a declination of the vector.
		/// </summary>
		/// <param name="cartesian">Vector.</param>
		/// <returns>Declination in radians.</returns>
		public static double Dec(this ICartesian cartesian)
		{
			return Math.Asin(cartesian.Z / cartesian.VectorLength());
		}

		/// <summary>
		/// Extension method calculating raght ascension of the vector.
		/// </summary>
		/// <param name="cartesian">Vector.</param>
		/// <returns>Raght ascension in radians.</returns>
		public static double Ra(this ICartesian cartesian)
		{
			return Math.Atan2(cartesian.Y, cartesian.X);
		}

		/// <summary>
		/// Extension method calculating a declination of the vector.
		/// </summary>
		/// <param name="cartesian">Vector.</param>
		/// <returns>Declination in degrees.</returns>
		public static double Declination(this ICartesian cartesian)
		{
			return Math.Asin(cartesian.Z / cartesian.VectorLength()) * 180.0 / Math.PI;
		}

		/// <summary>
		/// Extension method calculating raght ascension of the vector.
		/// </summary>
		/// <param name="cartesian">Vector.</param>
		/// <returns>Raght ascension in degrees.</returns>
		public static double RightAscension(this ICartesian cartesian)
		{
			return Math.Atan2(cartesian.Y, cartesian.X) * 180.0 / Math.PI;
		}


		/// <summary>
		/// Extension method normalizing the specified vector.
		/// </summary>
		/// <param name="cartesian">Vector to normalize.</param>
		/// <returns>Previous length of the vector.</returns>
		public static double Normalize(this ICartesian cartesian)
		{
			double len = cartesian.VectorLength();

			if (len > 0) cartesian.SetCartesian(cartesian.X / len, cartesian.Y / len, cartesian.Z / len);

			return len;
		}

		public static T Add<T>(this ICartesian cartesian, ICartesian V, bool normalize) where T : ICartesian, new()
		{
			double x = cartesian.X + V.X;
			double y = cartesian.Y + V.Y;
			double z = cartesian.Z + V.Z;

			if (normalize)
			{
				double l = Math.Sqrt(x * x + y * y + z * z);
				x /= l;
				y /= l;
				z /= l;
			}

			return NewCartesianDerivative<T>(x, y, z);
		}

		/// <summary>
		/// Extension method calculating cross product of two vectors. The second vector is represented by X,Y,Z components. 
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">First vector.</param>
		/// <param name="x">X component of the second vector.</param>
		/// <param name="y">Y component of the second vector.</param>
		/// <param name="z">Z component of the second vector.</param>
		/// <param name="normalize">When true, method returns normalized vector.</param>
		/// <returns>Vector representing cross product.</returns>
		public static T CrossProduct<T>(this ICartesian cartesian, double x, double y, double z, bool normalize) where T : ICartesian, new()
		{
			double x1 = cartesian.Y * z - cartesian.Z * y;
			double y1 = -(cartesian.X * z - cartesian.Z * x);
			double z1 = cartesian.X * y - cartesian.Y * x;

			if (normalize)
			{
				double l = Math.Sqrt(x1 * x1 + y1 * y1 + z1 * z1);
				x1 /= l;
				y1 /= l;
				z1 /= l;
			}

			return NewCartesianDerivative<T>(x1, y1, z1);
		}

		/// <summary>
		/// Extension method calculating cross product of two vectors.
		/// </summary>
		/// <typeparam name="T">Type T must implement ICartesian interface.</typeparam>
		/// <param name="cartesian">First vector.</param>
		/// <param name="V">Second vector.</param>
		/// <param name="normalize">When true, method returns normalized vector.</param>
		/// <returns>Vector representing cross product.</returns>
		public static T CrossProduct<T>(this ICartesian cartesian, ICartesian V, bool normalize) where T : ICartesian, new()
		{
			return cartesian.CrossProduct<T>(V.X, V.Y, V.Z, normalize);
		}

		/// <summary>
		/// Extension method calculating cross product of two vectors. The second vector is represented by declination and right ascension.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian"></param>
		/// <param name="dec">Declination in radians</param>
		/// <param name="ra">Right ascension</param>
		/// <returns>Vector representing cross product.</returns>
		public static T CrossProduct<T>(this ICartesian cartesian, double dec, double ra) where T : ICartesian, new()
		{
			double cos_dec = Math.Cos(dec);

			return cartesian.CrossProduct<T>(cos_dec * Math.Cos(ra), cos_dec * Math.Sin(ra), Math.Sin(dec), false);
		}

		/// <summary>
		/// Extension method calculating cross product of two vectors. Type T must implement ICartesian interface. 
		/// </summary>
		/// <param name="cartesian">First vector.</param>
		/// <param name="V">Second vector.</param>
		/// <returns>Cross product decomposed to x,y,z components.</returns>
		public static (double x, double y, double z) CrossProduct(this ICartesian cartesian, ICartesian V)
		{
			return (cartesian.Y * V.Z - cartesian.Z * V.Y, -(cartesian.X * V.Z - cartesian.Z * V.X), cartesian.X * V.Y - cartesian.Y * V.X);
		}

		/// <summary>
		/// Extension method calculating cross product of two vectors. The second vector is represented by X,Y,Z components. Type T must implement ICartesian interface. 
		/// </summary>
		/// <param name="cartesian">First vector.</param>
		/// <param name="x">X component of the second vector.</param>
		/// <param name="y">Y component of the second vector.</param>
		/// <param name="z">Z component of the second vector.</param>
		/// <returns>Cross product decomposed to x,y,z components.</returns>
		public static (double x, double y, double z) CrossProduct(this ICartesian cartesian, double x, double y, double z)
		{
			return (cartesian.Y * z - cartesian.Z * y, -(cartesian.X * z - cartesian.Z * x), cartesian.X * y - cartesian.Y * x);
		}

		/// <summary>
		/// Extension method calculating dot product of two vectors.
		/// </summary>
		/// <param name="cartesian">First vector.</param>
		/// <param name="V">Second vector.</param>
		/// <returns>Value of dot product.</returns>
		public static double DotProduct(this ICartesian cartesian, ICartesian V)
		{
			return cartesian.X * V.X + cartesian.Y * V.Y + cartesian.Z * V.Z;
		}

		/// <summary>
		/// Extension method calculating triple product of three vectors.
		/// </summary>
		/// <param name="cartesian">First vector.</param>
		/// <param name="V1">Second vector.</param>
		/// <param name="V2">Third vector.</param>
		/// <returns>Value of triple product.</returns>
		public static double TripleProduct(this ICartesian cartesian, ICartesian V1, ICartesian V2)
		{
			return cartesian.X * (V1.Y * V2.Z - V1.Z * V2.Y) -
				cartesian.Y * (V1.X * V2.Z - V1.Z * V2.X) +
				cartesian.Z * (V1.X * V2.Y - V1.Y * V2.X);
		}

		/// <summary>
		/// Extension method calculating triple product of three vectors. Two another vectors are represented as X,Y,Z components.
		/// </summary>
		/// <param name="cartesian">First vector.</param>
		/// <param name="x1">X component of the second vector.</param>
		/// <param name="y1">Y component of the second vector.</param>
		/// <param name="z1">Z component of the second vector.</param>
		/// <param name="x2">X component of the third vector.</param>
		/// <param name="y2">Y component of the third vector.</param>
		/// <param name="z2">Z component of the third vector.</param>
		/// <returns>Value of triple product.</returns>
		public static double TripleProduct(this ICartesian cartesian, double x1, double y1, double z1, double x2, double y2, double z2)
		{
			return cartesian.X * (y1 * z2 - z1 * y2) -
				cartesian.Y * (x1 * z2 - z1 * x2) +
				cartesian.Z * (x1 * y2 - y1 * x2);
		}

		/// <summary>
		/// Extension method rotating a vector around an axis by an angle. The result, the initial vector, and the axis are the right hand rule vectors. 
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">Vector to rotate.</param>
		/// <param name="angle">Rotation angle.</param>
		/// <param name="axis">Rotation axis.</param>
		/// <returns>New vector obtained from rotation.</returns>
		public static T Rotate<T>(this ICartesian cartesian, double angle, ICartesian axis) where T : ICartesian, new()
		{
			double L = axis.VectorLength();

			if (Math.Abs(L) < EPSILON) throw new ArgumentException();

			double cosine = Math.Cos(angle);
			double sine = Math.Sin(angle);

			double x = axis.X / L;
			double y = axis.Y / L;
			double z = axis.Z / L;

			// rotation matrix by lines
			// line 1
			double M11 = cosine + (1.0 - cosine) * x * x;
			double M12 = (1.0 - cosine) * x * y - sine * z;
			double M13 = (1.0 - cosine) * x * z + sine * y;
			// line 2
			double M21 = (1.0 - cosine) * y * x + sine * z;
			double M22 = cosine + (1.0 - cosine) * y * y;
			double M23 = (1.0 - cosine) * y * z - sine * x;
			// line 3
			double M31 = (1.0 - cosine) * z * x - sine * y;
			double M32 = (1.0 - cosine) * z * y + sine * x;
			double M33 = cosine + (1.0 - cosine) * z * z;

			return NewCartesianDerivative<T>(
				M11 * cartesian.X + M12 * cartesian.Y + M13 * cartesian.Z,
				M21 * cartesian.X + M22 * cartesian.Y + M23 * cartesian.Z,
				M31 * cartesian.X + M32 * cartesian.Y + M33 * cartesian.Z);
		}

		/// <summary>
		/// Extension method returning a vector representing a azimuth plane containing given location. 
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">Vector that represents location on a sphere.</param>
		/// <param name="azimuthOnSphere">Azimuth on a sphere.</param>
		/// <returns></returns>
		public static T AzimuthPlane<T>(this ICartesian cartesian, double azimuthOnSphere) where T : ICartesian, new()
		{
			T result;
			double dec = cartesian.Dec();
			double ra = cartesian.Ra();

			if (Math.Abs(dec) < EPSILON)
			{
				var (x, y, z) = SphericalToCartesian(azimuthOnSphere, ra - Math.PI / 2.0);
				result = NewCartesianDerivative<T>(x, y, z);
			}
			else if (Math.Abs(Math.PI / 2 - Math.Abs(dec)) < EPSILON)
				result = new T();
			else
				result = cartesian.CrossProduct<T>(0.0, ra - Math.Atan(Math.Sin(dec) * Math.Tan(azimuthOnSphere)));

			return result;
		}

		/// <summary>
		/// Extension method dividing a line (not an arc on a sphere!) connecting two vectors in the specified ratio.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="cartesian">First vector.</param>
		/// <param name="lineEnd">Second vector.</param>
		/// <param name="ratio">Ratio in the interval 0..1</param>
		/// <param name="normalize">When true, method returns a normalized vector.</param>
		/// <returns>Vector dividing a line in the specified ratio.</returns>
		public static T DivideLineInRatio<T>(this ICartesian cartesian, ICartesian lineEnd, double ratio, bool normalize) where T : ICartesian, new()
		{
			if (ratio < 0.0 || ratio > 1.0) ratio = 0.5;
			double r = 1.0 - ratio;

			double x = r * cartesian.X + ratio * lineEnd.X;
			double y = r * cartesian.Y + ratio * lineEnd.Y;
			double z = r * cartesian.Z + ratio * lineEnd.Z;

			if (normalize)
			{
				double l = Math.Sqrt(x * x + y * y + z * z);
				x /= l;
				y /= l;
				z /= l;
			}

			return NewCartesianDerivative<T>(x, y, z);
		}

		/// <summary>
		/// Extension method testing whether a vector belongs to a triangle formed by three vectors. 
		/// </summary>
		/// <param name="cartesian">Vector to test.</param>
		/// <param name="vertex1">First vertex of the triangle.</param>
		/// <param name="vertex2">Second vertex of the triangle.</param>
		/// <param name="vertex3">Third vertex of the triangle.</param>
		/// <returns>True if the vector belongs to the triangle.</returns>
		public static bool InsideTriangle(this ICartesian cartesian, ICartesian vertex1, ICartesian vertex2, ICartesian vertex3)
		{
			double x1 = cartesian.TripleProduct(vertex1, vertex2);
			double x2 = cartesian.TripleProduct(vertex2, vertex3);
			double x3;

			if (x1 * x2 >= 0.0 &&
				x2 * (x3 = cartesian.TripleProduct(vertex3, vertex1)) >= 0.0 &&
				x1 * x3 >= 0.0)
			{
				x1 = vertex2.X - vertex1.X;
				x2 = vertex3.X - vertex1.X;
				double y1 = vertex2.Y - vertex1.Y;
				double y2 = vertex3.Y - vertex1.Y;
				double z1 = vertex2.Z - vertex1.Z;
				double z2 = vertex3.Z - vertex1.Z;
				double l1 = Math.Sqrt(x1 * x1 + y1 * y1 + z1 * z1);
				double l2 = Math.Sqrt(x2 * x2 + y2 * y2 + z2 * z2);
				x1 /= l1;
				y1 /= l1;
				z1 /= l1;
				x2 /= l2;
				y2 /= l2;
				z2 /= l2;

				double A = y1 * z2 - z1 * y2;
				double B = -(x1 * z2 - x2 * z1);
				double C = x1 * y2 - x2 * y1;

				l1 = cartesian.X * A + cartesian.Y * B + cartesian.Z * C;
				l2 = vertex1.X * A + vertex1.Y * B + vertex1.Z * C;

				return l1 * l2 >= 0.0;
			}

			return false;
		}

		/// <summary>
		/// Extension method testing whether a vector belongs to a triangle formed by three another vectors represented by X,Y,Z components. 
		/// </summary>
		/// <param name="cartesian">Vector to test.</param>
		/// <param name="x1">X component of the first vertex of the triangle</param>
		/// <param name="y1">Y component of the first vertex of the triangle</param>
		/// <param name="z1">Z component of the first vertex of the triangle</param>
		/// <param name="x2">X component of the second vertex of the triangle</param>
		/// <param name="y2">Y component of the second vertex of the triangle</param>
		/// <param name="z2">Z component of the second vertex of the triangle</param>
		/// <param name="x3">X component of the third vertex of the triangle</param>
		/// <param name="y3">Y component of the third vertex of the triangle</param>
		/// <param name="z3">Z component of the third vertex of the triangle</param>
		/// <returns>True if the  vector belongs to the triangle.</returns>
		public static bool InsideTriangle(this ICartesian cartesian, double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3)
		{
			double p1 = cartesian.TripleProduct(x1, y1, z1, x2, y2, z2);
			double p2 = cartesian.TripleProduct(x2, y2, z2, x3, y3, z3);
			double p3;

			if (p1 * p2 >= 0.0 &&
				p2 * (p3 = cartesian.TripleProduct(x3, y3, z3, x1, y1, z1)) >= 0.0 &&
				p1 * p3 >= 0.0)
			{
				x2 -= x1;
				x3 -= x1;
				y2 -= y1;
				y3 -= y1;
				z2 -= z1;
				z3 -= z1;
				double l1 = Math.Sqrt(x2 * x2 + y2 * y2 + z2 * z2);
				double l2 = Math.Sqrt(x3 * x3 + y3 * y3 + z3 * z3);
				x2 /= l1;
				y2 /= l1;
				z2 /= l1;
				x3 /= l2;
				y3 /= l2;
				z3 /= l2;

				double A = y2 * z3 - z2 * y3;
				double B = -(x2 * z3 - x3 * z2);
				double C = x2 * y3 - x3 * y2;

				l1 = cartesian.X * A + cartesian.Y * B + cartesian.Z * C;
				l2 = x1 * A + y1 * B + z1 * C;

				return l1 * l2 >= 0.0;
			}

			return false;
		}

		/// <summary>
		/// Extension method testing whether a polygon is a convex.
		/// </summary>
		/// <param name="polygon">Polygon to test.</param>
		/// <returns>True if a convex.</returns>
		public static bool IsConvex(IList<ICartesian> polygon)
		{
			int n = polygon.Count - 2;
			int sign;
			int k = 0;

			do
			{
				sign = Math.Sign(polygon[k].TripleProduct(polygon[k + 1], polygon[k + 2]));
				k++;
			}
			while (k < n && sign == 0);

			int tsign;
			bool result = true;

			for (int i = k; i < n; i++)
				if ((tsign = Math.Sign(polygon[i].TripleProduct(polygon[i + 1], polygon[i + 2]))) != 0 && tsign * sign < 0)
				{
					result = false;
					break;
				}

			if (result &&
				(tsign = Math.Sign(polygon[n].TripleProduct(polygon[n + 1], polygon[0]))) != 0) result = tsign * sign > 0;

			if (result &&
				(tsign = Math.Sign(polygon[n + 1].TripleProduct(polygon[0], polygon[1]))) != 0) result = tsign * sign > 0;

			return result;
		}

		/// <summary>
		/// Returns an array of angles between planes formed by vertices of a polygon. Any negative angle in the returned array indicates that the polygon is not a convex.
		/// </summary>
		/// <param name="polygon">Polygon to examine.</param>
		/// <returns>Array of double with angles between planes.</returns>
		public static double[] PolygonInfo(IEnumerable<ICartesian> polygon)
		{
			double[] result = new double[polygon.Count()];

			ICartesian first = polygon.First();
			ICartesian previous = polygon.Last();
			ICartesian current = first;
			CartesianStruct N1, N2;
			int i = 0;

			foreach (ICartesian next in polygon.Skip(1))
			{
				N1 = previous.CrossProduct<CartesianStruct>(current, false);
				N2 = current.CrossProduct<CartesianStruct>(next, false);
				result[i++] = previous.TripleProduct(current, next) > 0.0 ? N1.Angle(N2) : -N1.Angle(N2);

				previous = current;
				current = next;
			}

			N1 = previous.CrossProduct<CartesianStruct>(current, false);
			N2 = current.CrossProduct<CartesianStruct>(first, false);
			result[i] = previous.TripleProduct(current, first) > 0.0 ? N1.Angle(N2) : -N1.Angle(N2);

			i = result.Count(a => a < 0.0);

			if (i > result.Length - i)
				for (i = 0; i < result.Length; i++) result[i] = -result[i];

			return result;
		}

		/// <summary>
		/// Extension method testing if a vector is inside a convex.
		/// </summary>
		/// <param name="cartesian">Vector to test.</param>
		/// <param name="vertices">Vertices of the convex.</param>
		/// <returns>True if inside.</returns>
		public static bool InsideConvex(this ICartesian cartesian, IList<ICartesian> polygon)
		{
			int n = polygon.Count - 1;
			bool result = false;

			for (int i = 1; i < n & !result; i++)
				result |= cartesian.InsideTriangle(polygon[0], polygon[i], polygon[i + 1]);

			return result;
		}

		/// <summary>
		/// Extension method testing if a vector is inside a polygon. Uses Winding number method.
		/// </summary>
		/// <param name="cartesian">Vector to test.</param>
		/// <param name="polygon">Vertices of the polygon.</param>
		/// <param name="includeBorders">Include borders of the polygon when true.</param>
		/// <returns>True if inside.</returns>
		public static bool InsidePolygon(this ICartesian cartesian, IEnumerable<ICartesian> polygon, bool includeBorders)
		{
			CartesianStruct N1;
			double wn = 0.0;
			ICartesian previous = polygon.Last();

			CartesianStruct N = cartesian.CrossProduct<CartesianStruct>(previous, false);

			foreach (ICartesian current in polygon)
			{
				N1 = cartesian.CrossProduct<CartesianStruct>(current, false);

				if (cartesian.TripleProduct(previous, current) > 0.0)
					wn += N.Angle(N1);
				else
					wn -= N.Angle(N1);

				N = N1;
				previous = current;
			}

			bool result = Math.Abs(wn) >= 2.0 * Math.PI;

			if (!result && includeBorders)
			{
				previous = polygon.Last();

				foreach (ICartesian current in polygon)
				{
					if (Math.Abs(cartesian.TripleProduct(previous, current)) < EPSILON)
					{
						result = true;
						break;
					}

					previous = current;
				}
			}

			return result;
		}

		/// <summary>
		/// Inflates a convex polygon by rotating its planes by given angles.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="polygon">Polygon to inflate.</param>
		/// <param name="angles">Angles of rotation in radians. Negative value deflates the polygon.</param>
		/// <returns>Inflated polygon</returns>
		public static List<T> InflateConvex<T>(IEnumerable<ICartesian> polygon, IList<double> angles) where T : ICartesian, new()
		{
			int n = polygon.Count();
			int m = angles.Count;

			if (n < 3 || m == 0) throw new ArgumentException();

			bool clockwise = false;
			int sign;
			Cartesian center = null;
			ICartesian first = polygon.First();
			ICartesian last = polygon.Last();
			ICartesian previous = last;
			ICartesian current = first;

			foreach (ICartesian next in polygon.Skip(1))
			{
				if ((sign = Math.Sign(previous.TripleProduct(current, next))) != 0)
				{
					clockwise = sign > 0;
					center = previous.DivideLineInRatio<Cartesian>(next, 0.5, false);
					break;
				}

				previous = current;
				current = next;
			}

			if (center == null)
			{
				ICartesian next = first;

				if ((sign = Math.Sign(previous.TripleProduct(current, next))) != 0)
				{
					clockwise = sign > 0;
					center = previous.DivideLineInRatio<Cartesian>(next, 0.5, false);
				}
				else throw new ArgumentException();
			}

			CartesianStruct[] planes = new CartesianStruct[n];
			CartesianStruct N, axis;
			previous = last;
			int i = 0;
			int j = 0;

			foreach (ICartesian next in polygon)
			{
				N = previous.CrossProduct<CartesianStruct>(next, false);
				axis = N.CrossProduct<CartesianStruct>(previous.DivideLineInRatio<CartesianStruct>(next, 0.5, false), false);

				if (!clockwise) axis.SetCartesian(-axis.X, -axis.Y, -axis.Z);

				planes[i++] = N.Rotate<CartesianStruct>(angles[j == m ? m - 1 : j++], axis);
				previous = next;
			}

			List<T> result = new List<T>(n);
			T res;

			for (i = 0; i < n - 1; i++)
			{
				res = planes[i].CrossProduct<T>(planes[i + 1], true);
				if (!clockwise) res.SetCartesian(-res.X, -res.Y, -res.Z);
				result.Add(res);
			}

			res = planes[n - 1].CrossProduct<T>(planes[0], true);
			if (!clockwise) res.SetCartesian(-res.X, -res.Y, -res.Z);
			result.Add(res);

			if (result.Any(p => double.IsNaN(p.X) || double.IsNaN(p.Y) || double.IsNaN(p.Z) || p.VectorLength() < EPSILON)) result = null;

			return result;
		}

		/// <summary>
		/// Extension method calculating a cosine of an angle between two vectors.
		/// </summary>
		/// <param name="cartesian">First vector</param>
		/// <param name="V">Second vector.</param>
		/// <returns>Cosine of the angle between two vectors.</returns>
		public static double Cosine(this ICartesian cartesian, ICartesian V)
		{
			double result = (cartesian.X * V.X + cartesian.Y * V.Y + cartesian.Z * V.Z) / (cartesian.VectorLength() * V.VectorLength());
			if (Math.Abs(result) > 1.0) result = Math.Truncate(result);

			return result;
		}

		/// <summary>
		/// Extension method calculating a cosine of an angle between two vectors, the second one is represented by X,Y,Z components.
		/// </summary>
		/// <param name="cartesian">First vector.</param>
		/// <param name="x">X component of the second vector.</param>
		/// <param name="y">Y component of the second vector.</param>
		/// <param name="z">Z component of the second vector.</param>
		/// <returns>Cosine of the angle between two vectors.</returns>
		public static double Cosine(this ICartesian cartesian, double x, double y, double z)
		{
			double result = (cartesian.X * x + cartesian.Y * y + cartesian.Z * z) / (cartesian.VectorLength() * Math.Sqrt(x * x + y * y + z * z));
			if (Math.Abs(result) > 1.0) result = Math.Truncate(result);

			return result;
		}

		/// <summary>
		/// Extension method calculating a cosine of anangle between two vectors, the second one is represented by declination and right ascension.
		/// </summary>
		/// <param name="cartesian">First  vector.</param>
		/// <param name="dec">Declination of the second vector in radians.</param>
		/// <param name="ra">Right ascension of the second vector in radians.</param>
		/// <returns>Cosine of the angle between two vectors.</returns>
		public static double Cosine(this ICartesian cartesian, double dec, double ra)
		{
			double c_dec = cartesian.Dec();
			double result = Math.Cos(c_dec) * Math.Cos(dec) * (Math.Cos(cartesian.Ra() - ra) - 1.0) + Math.Cos(c_dec - dec);
			if (Math.Abs(result) > 1.0) result = Math.Truncate(result);

			return result;
		}

		/// <summary>
		/// Extension method calculating an angle [0..π) between two vectors.
		/// </summary>
		/// <param name="cartesian">First vector.</param>
		/// <param name="V">Second vector.</param>
		/// <returns>Angle in radians between two vectors.</returns>
		public static double Angle(this ICartesian cartesian, ICartesian V)
		{
			return Math.Acos(cartesian.Cosine(V));
		}

		/// <summary>
		/// Extension method calculating an angle [0..π) between two vectors, the second one is represented as X,Y,Z components.
		/// </summary>
		/// <param name="cartesian">First vector.</param>
		/// <param name="x">X component of the second vector.</param>
		/// <param name="y">Y component of the second vector.</param>
		/// <param name="z">Z component of the second vector.</param>
		/// <returns>Angle in radians between two vectors.</returns>
		public static double Angle(this ICartesian cartesian, double x, double y, double z)
		{
			return Math.Acos(cartesian.Cosine(x, y, z));
		}

		/// <summary>
		/// Extension method calculating an angle [0..π) between two vectors, the second one is represented as declination and right ascension.
		/// </summary>
		/// <param name="cartesian">First 3D vector.</param>
		/// <param name="dec">Declination ofthe second vector in radians.</param>
		/// <param name="ra">Right ascension of the second vector in radians.</param>
		/// <returns>Angle in radian between two vectorss.</returns>
		public static double Angle(this ICartesian cartesian, double dec, double ra)
		{
			double c_dec = cartesian.Dec();
			double result = Math.Cos(c_dec) * Math.Cos(dec) * (Math.Cos(cartesian.Ra() - ra) - 1.0) + Math.Cos(c_dec - dec);
			if (Math.Abs(result) > 1.0) result = Math.Truncate(result);

			return Math.Acos(result);
		}

		/// <summary>
		/// Returns initial azimuth from one location on a sphere to another.
		/// </summary>
		/// <param name="cartesian">First location.</param>
		/// <param name="V">Second location.</param>
		/// <returns>Azimuth.</returns>
		public static double AzimuthOnSphere(this ICartesian cartesian, ICartesian V)
		{
			// N1 - cross product of cartesian and axiz Z, N1z is always zero and not used.
			double N1x = cartesian.Y;
			double N1y = -cartesian.X;

			if (Math.Abs(N1x) < EPSILON && Math.Abs(N1y) < EPSILON)
				return cartesian.Z >= 0 ? Math.PI : 0.0;

			// N2 - cross product of cartesian and V
			double N2x = cartesian.Y * V.Z - cartesian.Z * V.Y;
			double N2y = -(cartesian.X * V.Z - cartesian.Z * V.X);
			double N2z = cartesian.X * V.Y - cartesian.Y * V.X;

			double cosine = (N1x * N2x + N1y * N2y) / (Math.Sqrt(N1x * N1x + N1y * N1y) * Math.Sqrt(N2x * N2x + N2y * N2y + N2z * N2z));
			if (Math.Abs(cosine) > 1.0) cosine = Math.Truncate(cosine);
			double result = Math.Acos(cosine);

			double ra = cartesian.Ra();
			double v_ra = V.Ra();

			if (Math.Abs(v_ra - ra) > Math.PI)
			{
				if (v_ra > 0.0) result = 2.0 * Math.PI - result;
			}
			else if (v_ra < ra) result = 2.0 * Math.PI - result;

			return result;
		}

		/// <summary>
		/// Tests two sections for intersection.
		/// </summary>
		/// <param name="vertex1S1">First section start.</param>
		/// <param name="vertex2S1">First section end.</param>
		/// <param name="vertex1S2">Second section start.</param>
		/// <param name="vertex2S2">Second section end.</param>
		/// <returns>True if intersect.</returns>
		public static int SectionsIntersect(ICartesian vertex1S1, ICartesian vertex2S1, ICartesian vertex1S2, ICartesian vertex2S2)
		{
			var (x1, y1, z1) = vertex1S1.CrossProduct(vertex2S1);
			double length = Math.Sqrt(x1 * x1 + y1 * y1 + z1 * z1);
			x1 /= length;
			y1 /= length;
			z1 /= length;

			double l1 = x1 * vertex1S2.X + y1 * vertex1S2.Y + z1 * vertex1S2.Z;
			double l2 = x1 * vertex2S2.X + y1 * vertex2S2.Y + z1 * vertex2S2.Z;

			if (l1 * l2 < 0)
			{
				var (x2, y2, z2) = vertex1S2.CrossProduct(vertex2S2);
				length = Math.Sqrt(x2 * x2 + y2 * y2 + z2 * z2);
				x2 /= length;
				y2 /= length;
				z2 /= length;

				l1 = x2 * vertex1S1.X + y2 * vertex1S1.Y + z2 * vertex1S1.Z;
				l2 = x2 * vertex2S1.X + y2 * vertex2S1.Y + z2 * vertex2S1.Z;
			}

			return -Math.Sign(l1 * l2);
		}

		/// <summary>
		/// Tests two sections for intersection with given tolerance. Sections may not intersect, but to be close to each other with some tolerance.
		/// </summary>
		/// <param name="vertex1S1">First section start.</param>
		/// <param name="vertex2S1">First section end.</param>
		/// <param name="vertex1S2">Second section start.</param>
		/// <param name="vertex2S2">Second section end.</param>
		/// <param name="tolerance">Tolerance given as an angle at the origin.</param>
		/// <returns></returns>
		public static int SectionsIntersect(ICartesian vertex1S1, ICartesian vertex2S1, ICartesian vertex1S2, ICartesian vertex2S2, double tolerance)
		{
			var (x1, y1, z1) = vertex1S1.CrossProduct(vertex2S1);
			double length = Math.Sqrt(x1 * x1 + y1 * y1 + z1 * z1);
			x1 /= length;
			y1 /= length;
			z1 /= length;

			double l1 = x1 * vertex1S2.X + y1 * vertex1S2.Y + z1 * vertex1S2.Z;
			double l2 = x1 * vertex2S2.X + y1 * vertex2S2.Y + z1 * vertex2S2.Z;

			if (l1 * l2 < 0.0)
			{
				var (x2, y2, z2) = vertex1S2.CrossProduct(vertex2S2);
				length = Math.Sqrt(x2 * x2 + y2 * y2 + z2 * z2);
				x2 /= length;
				y2 /= length;
				z2 /= length;

				l1 = x2 * vertex1S1.X + y2 * vertex1S1.Y + z2 * vertex1S1.Z;
				l2 = x2 * vertex2S1.X + y2 * vertex2S1.Y + z2 * vertex2S1.Z;
			}

			int sign = -Math.Sign(l1 * l2);

			if (sign >= 0 || tolerance <= 0.0) return sign;

			double cos_t = Math.Cos(tolerance);

			if (((l1 = vertex1S1.TestSection(vertex1S2, vertex2S2)) >= 0.0 && l1 < tolerance) ||
				((l1 = vertex2S1.TestSection(vertex1S2, vertex2S2)) >= 0.0 && l1 < tolerance) ||
				((l1 = vertex1S2.TestSection(vertex1S1, vertex2S1)) >= 0.0 && l1 < tolerance) ||
				((l1 = vertex2S2.TestSection(vertex1S1, vertex2S1)) >= 0.0 && l1 < tolerance) ||
				vertex1S1.Cosine(vertex1S2) >= cos_t ||
				vertex1S1.Cosine(vertex2S2) >= cos_t ||
				vertex2S1.Cosine(vertex1S2) >= cos_t ||
				vertex2S1.Cosine(vertex2S2) >= cos_t) sign = 1;

			return sign;
		}

		/// <summary>
		/// Extension method testing if a vector is between two another vectors. 
		/// Returned positive value indicates success and represents a minimum angle from the vector to a geodesical arc on a sphere connecting two another vectors. 
		/// Negative value indicates test failure and its absolute value represents a minimum angle to the geodesical arc.
		/// </summary>
		/// <param name="cartesian">Vector to test.</param>
		/// <param name="vertex1">First vector.</param>
		/// <param name="vertex2">Second vector.</param>
		/// <returns>Minimum angle from the vector to the geodesical arc. Positive when succeeds.</returns>
		public static double TestSection(this ICartesian cartesian, ICartesian vertex1, ICartesian vertex2)
		{
			var N = vertex1.CrossProduct(vertex2);
			double alpha = cartesian.Angle(N.x, N.y, N.z);

			if (alpha > Math.PI / 2.0) alpha -= Math.PI / 2.0;
			else alpha = Math.PI / 2.0 - alpha;

			var N1 = cartesian.CrossProduct(N.x, N.y, N.z);
			double l1 = N1.x * vertex1.X + N1.y * vertex1.Y + N1.z * vertex1.Z;
			double l2 = N1.x * vertex2.X + N1.y * vertex2.Y + N1.z * vertex2.Z;

			if (l1 * l2 > 0.0 || cartesian.DotProduct(vertex1) < 0.0 || cartesian.DotProduct(vertex2) < 0.0) alpha = -alpha;

			return alpha;
		}

		/// <summary>
		/// Extension method testing if a point represented by a vector belongs to a polyline. 
		/// </summary>
		/// <param name="cartesian">Point to test.</param>
		/// <param name="polyline">Vertices of the polyline.</param>
		/// <param name="tolerance">Tolerance given as an angle at the origin.</param>
		/// <param name="reverse">Iterate vertices in reverse order if true.</param>
		/// <param name="userTest">Additional test delegate. If the point moves and has direction, this is the place to test it.</param>
		/// <param name="angleFromPolyline">Angle representing minimum distance to polyline if test succeeds</param>
		/// <returns>Index of a vertex which is a start of closest polyline segment. -1 if test fails.</returns>
		public static int TestPolyline(this ICartesian cartesian, IEnumerable<ICartesian> polyline, double tolerance, bool reverse, Func<ICartesian, ICartesian, ICartesian, double, double, bool> userTest, out double angleFromPolyline)
		{
			double res, azimuth, angle;

			IEnumerable<ICartesian> vertices = reverse ? polyline.Reverse().Skip(1) : polyline.Skip(1);
			ICartesian previous = reverse ? polyline.Last() : polyline.First();
			ICartesian before = null;
			int i = 0;

			ICartesian closestVertexAhead = null;
			int closestIndex = -1;
			double closestAngle = double.MaxValue;
			double closestAzimuth = double.MinValue;

			foreach (var vertex in vertices)
			{
				res = cartesian.TestSection(previous, vertex);

				if (res >= 0.0)
				{
					if (res <= tolerance && res < closestAngle)
					{
						azimuth = 0;

						if (userTest == null || userTest(null, previous, vertex, res, azimuth = previous.AzimuthOnSphere(vertex)))
						{
							closestIndex = i;
							closestAngle = res;
							closestVertexAhead = vertex;
							closestAzimuth = azimuth;
						}
					}

					before = null;
				}
				else
				{
					if (before != null &&
						(angle = cartesian.Angle(previous)) <= tolerance &&
						angle < cartesian.Angle(before) &&
						angle < cartesian.Angle(vertex))
					{
						azimuth = 0;

						if (userTest == null || userTest(before, previous, vertex, res, azimuth = previous.AzimuthOnSphere(vertex)))
						{
							closestIndex = i;
							closestAngle = angle;
							closestVertexAhead = vertex;
							closestAzimuth = azimuth;
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
		/// <param name="cartesian">Point to test.</param>
		/// <param name="polyline">Vertices of the polyline.</param>
		/// <param name="tolerance">Tolerance given as an angle at the origin.</param>
		/// <param name="reverse">Iterate vertices in reverse order if true.</param>
		/// <param name="userTest">Additional test delegate. If the point moves and has direction, this is the place to test it.</param>
		/// <param name="result">Result of the previous test.</param>
		/// <returns>Returns new or modified previous PolylineTestResult object. Null if first test fails or increased property Fails of returned PolylineTestResult object</returns>
		public static PolylineTestResult TestPolyline(this ICartesian cartesian, IEnumerable<ICartesian> polyline, double tolerance, bool reverse, Func<ICartesian, ICartesian, ICartesian, double, double, bool> userTest, PolylineTestResult result)
		{
			double res, azimuth, angle;

			IEnumerable<ICartesian> vertices = reverse ? polyline.Reverse().Skip(1) : polyline.Skip(1);
			ICartesian previous = reverse ? polyline.Last() : polyline.First();
			ICartesian before = null;
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
								azimuth = i != sectionIndex ? previous.AzimuthOnSphere(vertex) : result.SectionDirection;

								if (userTest == null || userTest(null, previous, vertex, res, azimuth))
								{
									if (sectionIndex != i)
									{
										result.SectionIndex = reverse ? -1 - i : i;
										result.SectionDirection = azimuth;
									}

									result.AngleWithSectionPlane = res;
									result.AngleWithNextVertex = cartesian.Angle(vertex);
									result.Fails = 0;
									return result;
								}
							}

							before = null;
						}
						else
						{
							if (before != null &&
								(angle = cartesian.Angle(previous)) <= tolerance &&
								angle < cartesian.Angle(before) &&
								angle < cartesian.Angle(vertex))
							{
								azimuth = previous.AzimuthOnSphere(vertex);

								if (userTest == null || userTest(before, previous, vertex, res, azimuth))
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
				ICartesian closestVertexAhead = null;
				int closestIndex = -1;
				double closestAngle = double.MaxValue;
				double closestAzimuth = double.MinValue;

				foreach (var vertex in vertices)
				{
					res = cartesian.TestSection(previous, vertex);

					if (res >= 0.0)
					{
						if (res <= tolerance && res < closestAngle)
						{
							azimuth = previous.AzimuthOnSphere(vertex);

							if (userTest == null || userTest(null, previous, vertex, res, azimuth))
							{
								closestIndex = i;
								closestAngle = res;
								closestVertexAhead = vertex;
								closestAzimuth = azimuth;
							}
						}

						before = null;
					}
					else
					{
						if (before != null &&
							(angle = cartesian.Angle(previous)) <= tolerance &&
							angle < cartesian.Angle(before) &&
							angle < cartesian.Angle(vertex))
						{
							azimuth = previous.AzimuthOnSphere(vertex);

							if (userTest == null || userTest(before, previous, vertex, res, azimuth))
							{
								closestIndex = i;
								closestAngle = angle;
								closestVertexAhead = vertex;
								closestAzimuth = azimuth;
							}
						}

						before = previous;
					}

					previous = vertex;
					i++;
				}


				if (closestIndex >= 0)
				{
					result = new PolylineTestResult(polyline, reverse)
					{
						SectionIndex = reverse ? -1 - closestIndex : closestIndex,
						SectionDirection = closestAzimuth,
						AngleWithSectionPlane = Math.Abs(closestAngle),
						AngleWithNextVertex = cartesian.Angle(closestVertexAhead)
					};
				}
			}

			return result;
		}


		// IGeoCoordinate extension methods

		/// <summary>
		/// Extension method calculating a distance on a sphere with the given radius between two locations with given geocoordinates.
		/// </summary>
		/// <param name="from">First location.</param>
		/// <param name="point">Second location.</param>
		/// <param name="sphereRadius">Radius of the sphere. Default value is the radius of the Earth.</param>
		/// <returns>Distance between two locations.</returns>
		public static double DistanceTo(this IGeoCoordinate from, IGeoCoordinate point, double sphereRadius = EARTH_MEAN_RADIUS)
		{
			double cosine = Math.Cos(from.Lat) * Math.Cos(point.Lat) * (Math.Cos(from.Lon - point.Lon) - 1.0) + Math.Cos(from.Lat - point.Lat);
			return sphereRadius * Math.Acos(cosine);
		}

		/// <summary>
		/// Extension method calculating a distance on a sphere with the given radius between two locations with given geocoordinates. The second location is represented by Latitude and Logitude in radians. 
		/// </summary>
		/// <param name="from">First location.</param>
		/// <param name="lat">Latitude of the second location in radians.</param>
		/// <param name="lon">Longitude of the second location in radians.</param>
		/// <param name="sphereRadius">Radius of the sphere. Default value is the radius of the Earth.</param>
		/// <returns>Distance between two locations.</returns>
		public static double DistanceTo(this IGeoCoordinate from, double lat, double lon, double sphereRadius = EARTH_MEAN_RADIUS)
		{
			double cosine = Math.Cos(from.Lat) * Math.Cos(lat) * (Math.Cos(from.Lon - lon) - 1.0) + Math.Cos(from.Lat - lat);
			return sphereRadius * Math.Acos(cosine);
		}

		/// <summary>
		/// Extension method calculating an initial azimuth from one location on a sphere to another represented by Latitude and Longitude in radians. Azimuth to the North is 0, to the East is π/2, etc.
		/// </summary>
		/// <param name="from">First location.</param>
		/// <param name="lat">Latitude of the second location in radians.</param>
		/// <param name="lon">Longitude of the second location in radians.</param>
		/// <returns>Azimuth to another location in radians.</returns>
		public static double AzimuthTo(this IGeoCoordinate from, double lat, double lon)
		{
			double cosine_c = Math.Sin(from.Lat);        //because angle c is pi/2-Lat
			double cosine_a = Math.Sin(lat);

			double cosine_D = Math.Cos(from.Lat) * Math.Cos(lat) * (Math.Cos(from.Lon - lon) - 1.0) + Math.Cos(from.Lat - lat);

			// when distance is 0 or second point is on exactly opposite side of the earth - all directions are equal and we consider heading as 0 
			// or when cosine_c = sin(Lat) == -1 it means that Latitude is -90 degrees and we are in South Pole, so every direction is to north and heading=0
			if (1 - Math.Abs(cosine_D) < EPSILON || Math.Abs(cosine_c + 1.0) < EPSILON) return 0.0;

			// when cosine_c = sin(Lat) == 1 it means that Latitude is 90 degrees and we are in North Pole, so every direction is to south and heading=180
			if (Math.Abs(cosine_c - 1.0) < EPSILON) return Math.PI;

			double cosine_alpha = (cosine_a - cosine_D * cosine_c) / (Math.Sqrt(1.0 - cosine_D * cosine_D) * Math.Sqrt(1.0 - cosine_c * cosine_c));

			// this is possible
			if (Math.Abs(cosine_alpha) > 1.0) cosine_alpha = Math.Truncate(cosine_alpha);

			double heading = Math.Acos(cosine_alpha);

			if (Math.Abs(lon - from.Lon) > Math.PI)
			{
				if (lon > 0.0) heading = 2.0 * Math.PI - heading;
			}
			else if (lon < from.Lon) heading = 2.0 * Math.PI - heading;

			return heading;
		}

		/// <summary>
		/// Extension method calculating an initial azimuth from one location on a sphere to another. Azimuth to the North is 0, to the East is π/2, etc.
		/// </summary>
		/// <param name="from">First location.</param>
		/// <param name="point">Second location.</param>
		/// <returns>Azimuth to another location in radians.</returns>
		public static double AzimuthTo(this IGeoCoordinate from, IGeoCoordinate point)
		{
			return from.AzimuthTo(point.Lat, point.Lon);
		}

		/// <summary>
		/// Extension method calculating a distance and an initial heading from one location on a sphere with given radius to another represented by Latitude and Longitude in radians. Heading to the North is 0, to the East is π/2, etc.
		/// This combined method is little faster than DistanceTo() together with HeadingTo().
		/// </summary>
		/// <param name="from">First location.</param>
		/// <param name="lat">Latitude of the second location in radians.</param>
		/// <param name="lon">Longitude of the second location in radians.</param>
		/// <param name="heading">Out parameter with initial heading in radians.</param>
		/// <param name="sphereRadius">Radius of sphere. Default value is radius of the Earth.</param>
		/// <returns>Distance between two locations.</returns>
		public static double DistanceAndAzimuthTo(this IGeoCoordinate from, double lat, double lon, out double azimuth, double sphereRadius = EARTH_MEAN_RADIUS)
		{
			double cosine_c = Math.Sin(from.Lat);        //because angle c is pi/2-Lat
			double cosine_a = Math.Sin(lat);

			double cosine_D = Math.Cos(from.Lat) * Math.Cos(lat) * (Math.Cos(from.Lon - lon) - 1.0) + Math.Cos(from.Lat - lat);
			double distance = Math.Acos(cosine_D) * sphereRadius;

			if (1 - Math.Abs(cosine_D) < EPSILON || Math.Abs(cosine_c + 1.0) < EPSILON)
			{
				// when distance is 0 or second point is on exactly opposite side of the earth - all directions are equal and we consider heading as 0 
				// or when cosine_c = sin(Lat) == -1 it means that Latitude is -90 degrees and we are at South Pole, so every direction is to north and heading=0
				azimuth = 0;
				return distance;
			}

			if (Math.Abs(cosine_c - 1.0) < EPSILON)
			{
				// when cosine_c = sin(Lat) == 1 it means that Latitude is 90 degrees and we are at North Pole, so every direction is to south and heading=180
				azimuth = Math.PI;
				return distance;
			}

			double cosine_alpha = (cosine_a - cosine_D * cosine_c) / (Math.Sqrt(1 - cosine_D * cosine_D) * Math.Sqrt(1.0 - cosine_c * cosine_c));

			// this is possible
			if (Math.Abs(cosine_alpha) > 1.0) cosine_alpha = Math.Truncate(cosine_alpha);

			azimuth = Math.Acos(cosine_alpha);

			if (Math.Abs(lon - from.Lon) > Math.PI)
			{
				if (lon > 0) azimuth = 2.0 * Math.PI - azimuth;
			}
			else if (lon < from.Lon) azimuth = 2.0 * Math.PI - azimuth;

			return distance;
		}

		/// <summary>
		/// Extension method calculating distance and initial heading from one location on a sphere with given radius to another. Heading to the North is 0, to the East is π/2, etc.
		/// This combined method is little faster than the DistanceTo() together with the HeadingTo().
		/// </summary>
		/// <param name="from">First location.</param>
		/// <param name="point">Second location.</param>
		/// <param name="heading">Out parameter with initial heading in radians.</param>
		/// <param name="sphereRadius">Radius of the sphere. Default value is the radius of the Earth.</param>
		/// <returns>Distance between two locations.</returns>
		public static double DistanceAndAzimuthTo(this IGeoCoordinate from, IGeoCoordinate point, out double azimuth, double sphereRadius = EARTH_MEAN_RADIUS)
		{
			return from.DistanceAndAzimuthTo(point.Lat, point.Lon, out azimuth, sphereRadius);
		}

		/// <summary>
		/// Extension method returning a new location on a sphere standing at given distance and initial azimuth from the given location.
		/// </summary>
		/// <typeparam name="T">Type T must implement the IGeoCoordinate interface.</typeparam>
		/// <param name="from">Initial location.</param>
		/// <param name="distance">Distance from the initial location.</param>
		/// <param name="azimuth">Azimuth in radians.</param>
		/// <param name="sphereRadius">Radius of the sphere. Default value is the radius of the Earth.</param>
		/// <returns>Location standing at given distance and initial azimuth.</returns>
		public static T PointTo<T>(this IGeoCoordinate from, double distance, double azimuth, double sphereRadius = EARTH_MEAN_RADIUS) where T : IGeoCoordinate, new()
		{
			double lat = from.Lat;
			double lon = from.Lon;

			if (distance < 0.0) distance = -distance;

			if (distance < EPSILON)
				return NewGeoCoordinateDerivative<T>(lat, lon);

			if (Math.Abs(lat - Math.PI / 2.0) < EPSILON)
				return NewGeoCoordinateDerivative<T>(Math.PI / 2.0 - distance / sphereRadius, 0.0);

			if (Math.Abs(lat + Math.PI / 2.0) < EPSILON)
				return NewGeoCoordinateDerivative<T>(distance / sphereRadius - Math.PI / 2.0, 0.0);


			azimuth = azimuth.To2PIRange();

			double cosine_c = Math.Sin(lat);  //because angle c is pi/2-lat
			double cosine_D = Math.Cos(distance / sphereRadius);

			double cosine_a = cosine_D * cosine_c + Math.Sqrt(1.0 - cosine_D * cosine_D) * Math.Sqrt(1.0 - cosine_c * cosine_c) * Math.Cos(azimuth);
			lat = Math.PI / 2.0 - Math.Acos(cosine_a);

			double cosine_beta = (cosine_D - cosine_c * cosine_a) / (Math.Sqrt(1.0 - cosine_c * cosine_c) * Math.Sqrt(1.0 - cosine_a * cosine_a));

			if (Math.Abs(cosine_beta) > 1.0) cosine_beta = Math.Truncate(cosine_beta);

			double beta = Math.Acos(cosine_beta);

			if (azimuth < Math.PI) lon += beta; else lon -= beta;

			if (Math.Abs(lon) > Math.PI)
				if (lon > 0.0) lon -= 2.0 * Math.PI; else lon += 2.0 * Math.PI;

			return NewGeoCoordinateDerivative<T>(lat, lon);
		}

		/// <summary>
		/// Extension method returning the location in the middle between two locations. This method uses DivideLineInRatio() with ratio 0.5.
		/// </summary>
		/// <typeparam name="T">Type T must implement the IGeoCoordinate interface.</typeparam>
		/// <param name="from">First location.</param>
		/// <param name="point">Second location.</param>
		/// <returns>Location in the middle between two locations.</returns>
		public static T PointToMiddle<T>(this IGeoCoordinate from, IGeoCoordinate point) where T : IGeoCoordinate, new()
		{
			Cartesian v1 = new Cartesian(from.Lat, from.Lon, false);
			Cartesian result = v1.DivideLineInRatio<Cartesian>(new Cartesian(point.Lat, point.Lon, false), 0.5, true);

			return NewGeoCoordinateDerivative<T>(result.Dec, result.Ra);
		}

		/// <summary>
		/// Extension method checking if a location with given latitude and longitude is inside an area with given South-West and North-East corners.
		/// </summary>
		/// <param name="location">Location to test.</param>
		/// <param name="SW">Area South-West corner.</param>
		/// <param name="NE">Area North-East corner.</param>
		/// <param name="tolerance">Tolerance in radians representing distance on the sphere.</param>
		/// <returns>True when location is inside.</returns>
		public static bool WithinBounds(double lat, double lon, IGeoCoordinate SW, IGeoCoordinate NE, double tolerance)
		{
			double toleranceRad = Math.Abs(tolerance) / EARTH_MEAN_RADIUS;

			double nLat = NE.Lat + toleranceRad;
			if (nLat > Math.PI / 2.0) nLat = Math.PI / 2.0;

			if (lat > nLat) return false;

			double sLat = SW.Lat - toleranceRad;
			if (sLat < -Math.PI / 2.0) sLat = -Math.PI / 2.0;

			if (lat < sLat) return false;

			if (NE.Lon - SW.Lon >= 2.0 * Math.PI - 2.0 * toleranceRad) return true;

			double eLon = (NE.Lon + toleranceRad).ToLongitudeRange();
			double wLon = (SW.Lon - toleranceRad).ToLongitudeRange();

			return eLon < wLon ?
				(lon >= wLon && lon <= Math.PI) || (lon >= -Math.PI && lon <= eLon) :
				lon >= wLon && lon <= eLon;
		}

		/// <summary>
		/// Extension method checking if a location is inside an area with given South-West and North-East corners.
		/// </summary>
		/// <param name="location">Location to test.</param>
		/// <param name="SW">Area South-West corner.</param>
		/// <param name="NE">Area North-East corner.</param>
		/// <param name="tolerance">Tolerance in radians representing distance on the sphere.</param>
		/// <returns>True when location is inside.</returns>
		public static bool WithinBounds(this IGeoCoordinate location, IGeoCoordinate SW, IGeoCoordinate NE, double tolerance)
		{
			return WithinBounds(location.Lat, location.Lon, SW, NE, tolerance);
		}

		/// <summary>
		/// Extension method returning a key of a mercator map tile at given zoom contaning given location.
		/// </summary>
		/// <param name="coordinate">Location.</param>
		/// <param name="zoom">Map zoom level.</param>
		/// <returns>Tile key.</returns>
		public static long MercatorMapTileKey(this IGeoCoordinate coordinate, int zoom)
		{
			long res;

			if (Math.Abs(coordinate.Lat) < MaxMercatorMapsLat)
			{
				int N = 1 << zoom;
				double sinY = Math.Sin(coordinate.Lat);

				int tileX = (int)Math.Floor(N * (0.5 + coordinate.Lon.ToLongitudeRange() / 2.0 / Math.PI));
				int tileY = (int)Math.Floor(N * (0.5 - Math.Log((1.0 + sinY) / (1.0 - sinY)) / 4.0 / Math.PI));

				res = ((long)tileX << 32) + tileY;

				if (res < 0) throw new InvalidOperationException();
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

		/// <summary>
		/// Extension method "normalizing" the heading to the [0..360) range.
		/// </summary>
		/// <param name="heading">Heading in degrees.</param>
		/// <returns>Heading within [0..360) range.</returns>
		public static double To360Range(this double heading)
		{
			if (heading >= 360.0 ||
				heading < 0) heading -= Math.Floor(heading / 360.0) * 360.0;

			return heading;
		}

		/// <summary>
		/// Extension method "normalizing" a heading to the [0..2π) range.
		/// </summary>
		/// <param name="heading">Heading in radians.</param>
		/// <returns>Heading in the [0..2π) range.</returns>
		public static double To2PIRange(this double heading)
		{
			if (heading >= 2.0 * Math.PI ||
				heading < 0) heading -= Math.Floor(heading / (2.0 * Math.PI)) * 2.0 * Math.PI;

			return heading;
		}

		/// <summary>
		/// Extension method "normalizing" a longitude in radians to the [-π..π] range. 
		/// For example, longitude 182 degrees is -178 degrees.
		/// </summary>
		/// <param name="lon">Longitude in radians.</param>
		/// <returns>Longitude in the [-π..π] range.</returns>
		public static double ToLongitudeRange(this double lon)
		{
			if (lon > Math.PI ||
				lon < -Math.PI) lon -= Math.Floor((lon + Math.PI) / (2.0 * Math.PI)) * 2.0 * Math.PI;

			return lon;
		}

		/// <summary>
		/// Extension method "normalizing" a longitude in degrees to the [-180..180] range. 
		/// For example, longitude 182 degrees is -178 degrees.
		/// </summary>
		/// <param name="lon">Longitude in degrees.</param>
		/// <returns>Longitude in the [-180..180] range.</returns>
		public static double ToLongitudeRangeDegrees(this double longitude)
		{
			if (longitude > 180.0 ||
				longitude < -180.0) longitude -= Math.Floor((longitude + 180.0) / (2.0 * 180.0)) * 2.0 * 180.0;

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
		/// Extension method returning an opposite heading for a given heading.
		/// </summary>
		/// <param name="heading">Heading in radians or degrees depending on the degrees argument.</param>
		/// <param name="degrees">Heading in degrees when true.</param>
		/// <returns>Opposite heading.</returns>
		public static double OppositeHeading(this double heading, bool degrees = false)
		{
			if (degrees)
				if ((heading = heading.To360Range()) < 180.0) heading += 180.0; else heading -= 180.0;
			else
				if ((heading = heading.To2PIRange()) < Math.PI) heading += Math.PI; else heading -= Math.PI;

			return heading;
		}

		/// <summary>
		/// Extension method returning a difference of two headings.
		/// </summary>
		/// <param name="heading">First heading in radians or degrees depending on the degrees argument.</param>
		/// <param name="otherHeading">Second heading in radians or degrees depending on the degrees argument.</param>
		/// <param name="degrees">Headings in degrees when true.</param>
		/// <returns>Difference of two headings.</returns>
		public static double HeadingsDiff(this double heading, double otherHeading, bool degrees = false)
		{
			otherHeading -= heading;

			if (degrees)
			{
				if (otherHeading > 180.0) otherHeading -= 360.0;
				else if (otherHeading < -180.0) otherHeading += 360.0;
			}
			else
			{
				if (otherHeading > Math.PI) otherHeading -= 2 * Math.PI;
				else if (otherHeading < -Math.PI) otherHeading += 2 * Math.PI;
			}

			return otherHeading;
		}

		/// <summary>
		/// Extension method returning a difference of two headings.
		/// </summary>
		/// <param name="heading">First heading in radians or degrees depending on the degrees argument.</param>
		/// <param name="otherHeading">Second heading in radians or degrees depending on the degrees argument.</param>
		/// <param name="degrees">Headings in degrees when true.</param>
		/// <returns>Difference of two headings.</returns>
		public static float HeadingsDiff(this float heading, float otherHeading, bool degrees = false)
		{
			otherHeading -= heading;

			if (degrees)
			{
				if (otherHeading > 180.0f) otherHeading -= 360.0f;
				else if (otherHeading < -180.0f) otherHeading += 360.0f;
			}
			else
			{
				if (otherHeading > Math.PI) otherHeading -= (float)(2.0 * Math.PI);
				else if (otherHeading < -Math.PI) otherHeading += (float)(2.0 * Math.PI);
			}

			return otherHeading;
		}

		/// <summary>
		/// Extension method testing if a heading is within given range. The range can't be more than π. The end of the range can be larger than the start. 
		/// </summary>
		/// <param name="heading">Heading to test.</param>
		/// <param name="h1">Start value of the range.</param>
		/// <param name="h2">End value of the range.</param>
		/// <returns>True if within the range.</returns>
		public static bool HeadingInRange(this double heading, double h1, double h2)
		{
			if (Math.Abs(h1 - h2) <= Math.PI)
				return h2 > h1 ? heading >= h1 && heading <= h2 : heading >= h2 && heading <= h1;

			if (h2 > h1)
				return (heading >= h2 && heading < 2.0 * Math.PI) || (heading >= 0.0 && heading <= h1);

			return (heading >= h1 && heading < 2.0 * Math.PI) || (heading >= 0.0 && heading <= h2);
		}

		/// <summary>
		/// Returns South-West and North-East corners of an area that encloses given locations.
		/// </summary>
		/// <param name="locations">Locations.</param>
		/// <param name="SW">Out parameter with the area South-West corner.</param>
		/// <param name="NE">Out parameter with the area North-East corner.</param>
		public static void GetBounds<T>(IEnumerable<IGeoCoordinate> locations, out T SW, out T NE) where T : IGeoCoordinate, new()
		{
			if (locations == null || 
				!locations.Any()) throw new ArgumentException();

			double minLat = double.MaxValue;
			double minLon = double.MaxValue;
			double maxLat = double.MinValue;
			double maxLon = double.MinValue;

			foreach (IGeoCoordinate location in locations)
			{
				if (location.Lat > maxLat) maxLat = location.Lat;
				if (location.Lat < minLat) minLat = location.Lat;
				if (location.Lon > maxLon) maxLon = location.Lon;
				if (location.Lon < minLon) minLon = location.Lon;
			}

			if (maxLon - minLon > Math.PI * 2)
			{
				SW = NewGeoCoordinateDerivative<T>(minLat, maxLon);
				NE = NewGeoCoordinateDerivative<T>(maxLat, minLon); 
			}
			else
			{
				SW = NewGeoCoordinateDerivative<T>(minLat, minLon);
				NE = NewGeoCoordinateDerivative<T>(maxLat, maxLon);
			}
		}

		/// <summary>
		/// Returns length of a polyline on a sphere with given radius.
		/// </summary>
		/// <param name="locations">Polyline points.</param>
		/// <param name="sphereRadius">Radius of the sphere in meters. Default value is the radius of the Earth.</param>
		/// <returns>Length of the polyline in meters.</returns>
		public static double GetLength(IEnumerable<IGeoCoordinate> locations, double sphereRadius = EARTH_MEAN_RADIUS)
		{
			if (locations == null || 
				!locations.Any()) throw new ArgumentException();

			double result = 0.0;
			IGeoCoordinate previous = locations.First();

			foreach (IGeoCoordinate loc in locations.Skip(1))
			{
				result += previous.DistanceTo(loc);
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
		/// Auxiliary extension method encoding the next point in the polyline. Thank you, stackoverflow.com
		/// </summary>
		/// <param name="location">Current point to encode.</param>
		/// <param name="previous">Previous encoded point.</param>
		/// <param name="escape">Add second back-slash when true.</param>
		/// <returns>Encoded point in the polyline.</returns>
		public static string EncodeToGooglePolyline(this IGeoCoordinate location, IGeoCoordinate previous, bool escape)
		{
			StringBuilder result = new StringBuilder();
			char c;

			int lastLat = previous != null ? (int)(previous.Lat * 100000.0 * 180.0 / Math.PI) : 0;
			int lastLon = previous != null ? (int)(previous.Lon * 100000.0 * 180.0 / Math.PI) : 0;

			int diff = (int)(location.Lat * 100000.0 * 180.0 / Math.PI) - lastLat;
			int shifted = diff << 1;

			if (diff < 0) shifted = ~shifted;

			int rem = shifted;

			while (rem >= 0x20)
			{
				c = (char)((0x20 | (rem & 0x1f)) + 63);
				result.Append(c);
				if (c == '\\' && escape) result.Append('\\');

				rem >>= 5;
			}

			c = (char)(rem + 63);
			result.Append(c);
			if (c == '\\' && escape) result.Append('\\');

			diff = (int)(location.Lon * 100000.0 * 180.0 / Math.PI) - lastLon;
			shifted = diff << 1;

			if (diff < 0) shifted = ~shifted;

			rem = shifted;

			while (rem >= 0x20)
			{
				c = (char)((0x20 | (rem & 0x1f)) + 63);
				result.Append(c);
				if (c == '\\' && escape) result.Append('\\');

				rem >>= 5;
			}

			c = (char)(rem + 63);
			result.Append(c);
			if (c == '\\' && escape) result.Append('\\');

			return result.ToString();
		}

		/// <summary>
		/// Extension method encoding a polyline into the string.
		/// </summary>
		/// <param name="locations">Polyline points.</param>
		/// <param name="escape">Add second back-slash when true.</param>
		/// <returns>Encoded polyline.</returns>
		public static string EncodeToGooglePolyline(IEnumerable<IGeoCoordinate> locations, bool escape)
		{
			if (locations == null) return null;

			StringBuilder result = new StringBuilder();
			IGeoCoordinate previous = null;

			foreach (IGeoCoordinate location in locations)
			{
				result.Append(location.EncodeToGooglePolyline(previous, escape));
				previous = location;
			}

			return result.ToString();
		}

		/// <summary>
		/// Decodes a polyline from the string.
		/// </summary>
		/// <typeparam name="T">Type T must implement the ICartesian interface.</typeparam>
		/// <param name="encodedPoints">Encoded polyline.</param>
		/// <returns>Points of the polyline.</returns>
		public static List<T> DecodeGooglePolyline<T>(string encodedPoints) where T : ICartesian, new()
		{
			if (string.IsNullOrEmpty(encodedPoints)) return null;

			List<T> result = new List<T>();
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

				if (index >= polylineChars.Length) break;

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

				if (index >= polylineChars.Length && next5bits >= 32) break;

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
		/// <param name="location">Location to calculate sunrise time for.</param>
		/// <param name="date">Date and offset.</param>
		/// <returns>Sunrise time.</returns>
		public static DateTime Sunrise(this IGeoCoordinate location, DateTimeOffset date)
		{


			double JDays = date.DaysSinceJ2000();
			double JCenturies = JDays / 36525.0;
			double elevation = -5.0 / 6.0; // - 6 degrees for civil twilight

			double L = (4.8949504201433 + 628.331969753199 * JCenturies) % (2.0 * Math.PI);
			double G = (6.2400408 + 628.3019501 * JCenturies) % (2.0 * Math.PI);
			double ec = 0.033423 * Math.Sin(G) + 0.00034907 * Math.Sin(2.0 * G);
			double lambda = L + ec;
			double GHA = 0.0430398 * Math.Sin(2 * lambda) - 0.00092502 * Math.Sin(4 * lambda) - ec;
			double obl = 0.409093 - 0.0002269 * JCenturies;
			double delta = Math.Asin(Math.Sin(obl) * Math.Sin(lambda));
			double cosc = (Math.Sin(elevation * Math.PI / 180.0) - Math.Sin(location.Lat) * Math.Sin(delta)) / (Math.Cos(location.Lat) * Math.Cos(delta));

			if (cosc > 1.0) return DateTime.MaxValue; // never rises 
			if (cosc < -1.0) return DateTime.MinValue; // never sets

			double correction = Math.Acos(cosc);
			double utnew = Math.PI - (GHA + location.Lon + correction);
			double eventUT = utnew * 57.29577951 / 15.0;

			if (eventUT >= 24.0) eventUT %= 24.0;

			DateTime result = new DateTime(date.Year, date.Month, date.Day);
			return result.Add(date.Offset).AddHours(eventUT);
		}

		/// <summary>
		/// Extension method returning sunset time for a given location and date with offset.
		/// </summary>
		/// <param name="location">Location to calculate sunset time for.</param>
		/// <param name="date">Date and offset.</param>
		/// <returns>Sunset time.</returns>
		public static DateTime Sunset(this IGeoCoordinate location, DateTimeOffset date)
		{
			double JDays = date.DaysSinceJ2000();
			double JCenturies = JDays / 36525.0;
			double elevation = -5.0 / 6.0; // - 6 degrees for civil twilight

			double L = (4.8949504201433 + 628.331969753199 * JCenturies) % (2.0 * Math.PI);
			double G = (6.2400408 + 628.3019501 * JCenturies) % (2.0 * Math.PI);
			double ec = 0.033423 * Math.Sin(G) + 0.00034907 * Math.Sin(2.0 * G);
			double lambda = L + ec;
			double GHA = 0.0430398 * Math.Sin(2 * lambda) - 0.00092502 * Math.Sin(4 * lambda) - ec;
			double obl = 0.409093 - 0.0002269 * JCenturies;
			double delta = Math.Asin(Math.Sin(obl) * Math.Sin(lambda));
			double cosc = (Math.Sin(elevation * Math.PI / 180.0) - Math.Sin(location.Lat) * Math.Sin(delta)) / (Math.Cos(location.Lat) * Math.Cos(delta));  // B25

			if (cosc > 1.0) return DateTime.MaxValue; // never rises 
			if (cosc < -1.0) return DateTime.MinValue; // never sets

			double correction = Math.Acos(cosc);
			double utnew = Math.PI - (GHA + location.Lon - correction);
			double eventUT = utnew * 57.29577951 / 15.0;

			if (eventUT >= 24.0) eventUT %= 24.0;

			DateTime result = new DateTime(date.Year, date.Month, date.Day);
			return result.Add(date.Offset).AddHours(eventUT);
		}

		/// <summary>
		/// Returns times of nearest sunrise and sunset for a given location and date/time with the offset. 
		/// Nearest sunrise can be for the next day, nearest sunset can be for the previous day.
		/// Polar night at given location at given date, if method returns (DateTime.MaxValue, DateTime.MinValue).
		/// Polar day, if method returns (DateTime.MinValue, DateTime.MaxValue).
		/// </summary>
		/// <param name="location">Location.</param>
		/// <param name="date">Date, time and offset.</param>
		/// <returns>Times of sunrise and sunset.</returns>
		public static (DateTime sunrise, DateTime sunset) NearestSunriseSunset(IGeoCoordinate location, DateTimeOffset date)
		{
			// Install System.ValueTuple package via Nuget if this fails to compile.
			// PM > Install - Package "System.ValueTuple"
			DateTime sunrise, sunset;

			if (location != null)
			{
				sunrise = location.Sunrise(date);
				if (sunrise == DateTime.MaxValue) return (DateTime.MaxValue, DateTime.MinValue); // never rises
				if (sunrise == DateTime.MinValue) return (DateTime.MinValue, DateTime.MaxValue); // never sets

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
		/// <param name="location">Location to test.</param>
		/// <param name="date">Date, time, and offset for the test.</param>
		/// <returns>True if the Sun is in the sky.</returns>
		public static bool DayTime(IGeoCoordinate location, DateTimeOffset date)
		{
			var (sunrise, sunset) = NearestSunriseSunset(location, date);
			return date.DateTime > sunrise && date.DateTime < sunset;
		}
	}

}