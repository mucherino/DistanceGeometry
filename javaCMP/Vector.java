
/* javaCMP project
 *
 * Vector class
 *
 * This class collects some static methods for manipulating vectors (arrays of double).
 *
 * last update: April 23rd, 2023
 *
 * AM
*/

import java.util.Random;

public class Vector
{
   // fixed internal tolerance error
   final static double eps = 1.e-6;

   // generating a random vector of a given length, with all elements in the range [0,1].
   public static double[] random(int n,Random R)
   {
      try
      {
         if (n <= 0) throw new IllegalArgumentException("The given length for the vector is nonpositive");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }
      if (R == null)  R = new Random();

      double[] v = new double[n];
      for (int i = 0; i < n; i++)  v[i] = R.nextDouble();
      return v;
   }

   // computing the linear combination of two vectors with two coefficients alpha and beta
   public static double[] linear(double alpha,double[] v,double beta,double[] w)
   {
      try
      {
         if (v == null) throw new IllegalArgumentException("The first vector is null");
         if (w == null) throw new IllegalArgumentException("The second vector is null");
         if (v.length != w.length) throw new IllegalArgumentException("The input vectors do not have the same length");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      double[] x = new double[v.length];
      for (int i = 0; i < v.length; i++)  x[i] = alpha*v[i] + beta*w[i];
      return x;
   }

   // performing the sum of two vectors with the same length
   public static double[] sum(double[] v,double[] w)
   {
      return linear(1.0,v,1.0,w);
   }

   // performing the vectorial difference of two vectors with the same length
   public static double[] diff(double[] v,double[] w)
   {
      return linear(1.0,v,-1.0,w);
   }

   // computing the cross product of two three-dimensional vectors
   public static double[] crossProduct(double[] v,double[] w)
   {
      try
      {
         if (v == null) throw new IllegalArgumentException("The array of double representing the first vector is null");
         if (w == null) throw new IllegalArgumentException("The array of double representing the second vector is null");
         if (v.length != 3) throw new IllegalArgumentException("The length of the first vector is not 3");
         if (w.length != 3) throw new IllegalArgumentException("The length of the second vector is not 3");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      double[] result = new double[3];
      result[0] = v[1] * w[2] - v[2] * w[1];
      result[1] = v[2] * w[0] - v[0] * w[2];
      result[2] = v[0] * w[1] - v[1] * w[0];
      return result;
   }

   // computing the dot product of two vectors with the same length
   public static double dotProduct(double[] v,double[] w)
   {
      try
      {
         if (v == null) throw new IllegalArgumentException("The array of double representing the first vector is null");
         if (w == null) throw new IllegalArgumentException("The array of double representing the second vector is null");
         if (v.length != w.length) throw new IllegalArgumentException("The two vectors do not have the same dimensions");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      double sum = 0;
      for (int i = 0; i < v.length; i++)  sum = sum + v[i]*w[i];
      return sum;
   }

   // computing the norm of a vector
   public static double norm(double[] v)
   {
      try
      {
         if (v == null) throw new IllegalArgumentException("The array of double representing the vector is null");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      return Math.sqrt(dotProduct(v,v));
   }

   // computing the unit vector related to the input vector
   public static double[] normalize(double[] v)
   {
      double l = norm(v);
      for(int i = 0; i < v.length; i++)  v[i] = v[i]/l;
      return v;
   }

   // computing the unit vector normal to the plane defined by three vectors
   public static double[] normal2plane(double[] v,double[] w,double[] u)
   {
      double[] cross = null;
      try
      {
         if (v == null) throw new IllegalArgumentException("The array of double representing the vector v is null");
         if (w == null) throw new IllegalArgumentException("The array of double representing the vector w is null");
         if (u == null) throw new IllegalArgumentException("The array of double representing the vector u is null");
         if (v.length != 3) throw new IllegalArgumentException("The length of vector v is not 3");
         if (w.length != 3) throw new IllegalArgumentException("The length of vector w is not 3");
         if (u.length != 3) throw new IllegalArgumentException("The length of vector u is not 3");
         double[] vw = diff(w,v);
         double[] vu = diff(u,v);
         cross = crossProduct(vw,vu);
         if (cross[0] == 0.0 && cross[1] == 0.0 && cross[2] == 0.0) throw new IllegalArgumentException("The three given vectors are aligned");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      return normalize(cross);
   }

   // computing the Euclidean distance between two vectors
   public static double EuclideanDistance(double[] v,double[] w)
   {
      try
      {
         if (v == null) throw new IllegalArgumentException("The array of double representing the vector is null");
         if (w == null) throw new IllegalArgumentException("The array of double representing the vector is null");
         if (v.length != w.length) throw new IllegalArgumentException("The two vectors differ in size");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      double sqd = 0.0;
      for (int i = 0; i < v.length; i++)
      {
         double term = w[i] - v[i];
         sqd = sqd + term*term;
      }
      return Math.sqrt(sqd);
   }

   // verifying whether the triangular inequality is satisfied for the triangle (A,B,C)
   // the triangle is representated through the three vector lengths AB, BC and AC
   public static boolean triangularInequality(double dAB,double dBC,double dAC)
   {
      try
      {
         if (dAB <= 0.0) throw new IllegalArgumentException("The distance between A and B cannot be nonpositive");
         if (dBC <= 0.0) throw new IllegalArgumentException("The distance between B and C cannot be nonpositive");
         if (dAC <= 0.0) throw new IllegalArgumentException("The distance between A and C cannot be nonpositive");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      if (dAB > dBC && dAB > dAC)
         return dAB <= dBC + dAC + eps;
      else if (dBC > dAB && dBC > dAC)
         return dBC <= dAB + dAC + eps;
      return dAC <= dAB + dBC + eps;
   }

   // verifying whether the triangular inequality is satisfied for the triangle (A,B,C)
   // the triangle is represented through the 3D coordinates of the points A, B and C
   public static boolean triangularInequality(double[] A,double[] B,double[] C)
   {
      try
      {
         if (A == null) throw new IllegalArgumentException("The array of double for the coordinates of A is null");
         if (B == null) throw new IllegalArgumentException("The array of double for the coordinates of B is null");
         if (C == null) throw new IllegalArgumentException("The array of double for the coordinates of C is null");
         if (A.length != B.length || A.length != C.length) throw new IllegalArgumentException("The three arrays of double have different lengths");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      return triangularInequality(EuclideanDistance(A,B),EuclideanDistance(B,C),EuclideanDistance(A,C));
   }

   // computing the cosine of the vector angle in B in the triangle (A,B,C)
   // the triangle is representated through the three vector lengths AB, BC and AC
   public static double cos(double dAB,double dBC,double dAC)
   {
      try
      {
         if (!triangularInequality(dAB,dBC,dAC))
            throw new IllegalArgumentException("The given distances do not form a triangle (triangular inequality not satisfied)");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      double val = dAB*dAB + dBC*dBC - dAC*dAC;
      val = val / (2.0*dAB*dBC);
      if (val < -1.0)  val = -1.0;
      if (val >  1.0)  val =  1.0;
      return val;
   }

   // computing the cosine of the vector angle in B in the triangle (A,B,C)
   // the triangle is represented through the 3D coordinates of the points A, B and C
   public static double cos(double[] A,double[] B,double[] C)
   {
      try
      {
         if (A == null) throw new IllegalArgumentException("The array of double for the coordinates of A is null");
         if (B == null) throw new IllegalArgumentException("The array of double for the coordinates of B is null");
         if (C == null) throw new IllegalArgumentException("The array of double for the coordinates of C is null");
         if (A.length != B.length || A.length != C.length) throw new IllegalArgumentException("The three arrays of double have different lengths");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      return cos(EuclideanDistance(A,B),EuclideanDistance(B,C),EuclideanDistance(A,C));
   }

   // embedding a triangle (A,B,C) in the three-dimensional Euclidean space
   public static double[][] embedTriangle(double dAB,double dBC,double dAC)
   {
      double cosABC = cos(dAB,dBC,dAC);
      double sinABC = Math.sqrt(1.0 - cosABC*cosABC);
      double[][] X = new double[3][];
      X[0] = new double[] {0.0,0.0,0.0};
      X[1] = new double[] {-dAB,0.0,0.0};
      X[2] = new double[] {-dAB + dBC*cosABC,dBC*sinABC,0.0};
      return X;
   }

   // generating a String representation for a vector
   public static String toString(double[] v)
   {
      try
      {
         if (v == null) throw new IllegalArgumentException("The array of double representing the vector is null");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      String print = "(" + v[0];
      for (int i = 1; i < v.length; i++)  print = print + "," + v[i];
      return print + ")";
   }

   // main method -- performing some basic tests
   public static void main(String[] args) throws Exception
   {
      System.out.print("Vector class ");

      // seed for random numbers
      long seed;
      Random R = new Random();
      if (args != null && args.length > 0)
         seed = Long.parseLong(args[0]);
      else
         seed = Math.abs(R.nextLong());
      System.out.print("(seed is " + seed + ")");
      R = new Random(seed);

      System.out.print(" ... ");
      int NTESTS = 100000;

      for (int itests = 0; itests < NTESTS; itests++)
      {
         int n = 10 + R.nextInt(90);
         int k = R.nextInt(n);

         // random
         Exception E = new Exception("public static double[] random(int,Random)");
         double[] a = random(n,R);
         if (a == null) throw E;
         if (a.length != n) throw E;
         if (a[k] < 0.0 || a[k] > 1.0) throw E;

         // sum
         E = new Exception("public static double[] sum(double[],double[]");
         double[] b = random(n,R);
         double[] c = sum(a,b);
         if (c == null) throw E;
         k = R.nextInt(n);
         if (c[k] != a[k] + b[k]) throw E;

         // diff
         E = new Exception("public static double[] diff(double[],double[])");
         double[] A = diff(c,b);
         if (A == null) throw E;
         k = R.nextInt(n);
         if (Math.abs(A[k] - a[k]) > eps) throw E;
         double[] B = diff(c,a);
         if (B == null) throw E;
         k = R.nextInt(n);
         if (Math.abs(B[k] - b[k]) > eps) throw E;

         // linear
         E = new Exception("public static double[] linear(double,double[],double,double[])");
         c = linear(2.0,a,3.0,b);
         k = R.nextInt(n);
         if (Math.abs(c[k] - 2.0*a[k] - 3.0*b[k]) > eps) throw E;

         // crossProduct
         E = new Exception("public static double[] crossProduct(double[],double[])");
         a = random(3,R);
         b = random(3,R);
         k = R.nextInt(3);
         a[k] = 0.0;
         b[k] = 0.0;
         c = crossProduct(a,b);
         k++;
         if (k == 3)  k = R.nextInt(2);
         if (c[k] != 0.0) throw E;

         // dotProduct
         E = new Exception("public static double dotProduct(double[])");
         a = random(3,R);
         b = random(3,R);
         c = crossProduct(a,b);
         if (Math.abs(dotProduct(a,c)) > eps) throw E;
         if (Math.abs(dotProduct(b,c)) > eps) throw E;

         // norm and normalize
         E = new Exception("public static double norm(double[])");
         a = new double[n];
         for (int i = 0; i < n; i++)  a[i] = 0.0;
         if (Math.abs(norm(a)) > eps) throw E;
         for (int i = 0; i < n; i++)  a[i] = 1.0/Math.sqrt(n);
         if (Math.abs(norm(a) - 1.0) > eps) throw E;
         E = new Exception("public static double normalize(double[])");
         a = random(n,R);
         b = normalize(a);
         if (Math.abs(norm(a) - 1.0) > eps) throw E;
         if (Math.abs(norm(b) - 1.0) > eps) throw E;
         if (a != b) throw E;

         // normal2plane
         E = new Exception("public static double[] normal2plane(double[],double[],double[])");
         a = random(3,R);
         b = random(3,R);
         c = random(3,R);
         double[] normal = normal2plane(a,b,c);
         double[] orth = null;
         int s = R.nextInt(3);
         if (s == 0)
            orth = diff(a,b);
         else if (s == 1)
            orth = diff(a,c);
         else
            orth = diff(b,c);
         if (Math.abs(dotProduct(normal,orth)) > eps) throw E;

         // Euclidean distance
         E = new Exception("public static double EuclideanDistance(double[],double[])");
         a = random(n,R);
         double dist = EuclideanDistance(a,a);
         if (Math.abs(dist) > eps) throw E;
         b = sum(a,a);
         dist = EuclideanDistance(a,b);
         if (Math.abs(dist - norm(a)) > eps) throw E;

         // triangularInequality
         b = random(n,R);
         c = random(n,R);
         if (!triangularInequality(a,b,c)) throw new Exception("public static boolean triangularInequality(double[],double[],double[])");
         double dAB = EuclideanDistance(a,b);
         double dBC = EuclideanDistance(b,c);
         double dAC = dAB + 2.1*dBC;
         if (triangularInequality(dAB,dBC,dAC)) throw new Exception("public static boolean triangularInequality(double,double,double)");

         // cos
         dAC = EuclideanDistance(a,c);
         double cABC = cos(dAB,dBC,dAC);
         double cBCA = cos(dBC,dAC,dAB);
         double cCAB = cos(dAC,dAB,dBC);
         if (Math.abs(Math.acos(cABC) + Math.acos(cBCA) + Math.acos(cCAB) - Math.PI) > eps)
            throw new Exception("public static double cos_vector_angle(double,double,double)");
         cABC = cos(a,b,c);
         if (Math.abs(Math.acos(cABC) + Math.acos(cBCA) + Math.acos(cCAB) - Math.PI) > eps)
            throw new Exception("public static double cos_vector_angle(double[],double[],double[])");

         // embedTriangle
         E = new Exception("public static double[][] embedTriangle(double,double,double)");
         double[][] X = embedTriangle(dAB,dBC,dAC);
         if (X == null) throw E;
         a = X[0];
         b = X[1];
         c = X[2];
         if (Math.abs(dAB - EuclideanDistance(a,b)) > eps) throw E;
         if (Math.abs(dBC - EuclideanDistance(b,c)) > eps) throw E;
         if (Math.abs(dAC - EuclideanDistance(a,c)) > eps) throw E;
      }
      System.out.println("OK");
   }
}

