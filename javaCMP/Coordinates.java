
/* javaCMP project
 *
 * Coordinates class
 *
 * An instance of this class is able to simultaneously hold various types of coordinates.
 * Notice that every Coordinates instance is only aware of its own references;
 * there is for example no mechanism to avoid the generation of cycles in the implied graph structure.
 * This initial implementation is for dimension 3 only.
 *
 * last update: April 23rd, 2023
 *
 * AM
*/

import java.util.Set;
import java.util.HashSet;
import java.util.List;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Map;
import java.util.HashMap;
import java.util.Random;

public class Coordinates
{
   /* Data structures for the different coordinate types */

   // Cartesian type
   private class Cartesians
   {
      // Cartesians hold the reference of their Coordinates class
      Coordinates mainInstance;

      // internal coordinates
      BoundedDouble[] coords;

      // Cartesian constructor
      public Cartesians(Coordinates X)
      {
         this.mainInstance = X;
         this.coords = null;
      }

      // verifying whether the coordinates are defined
      public boolean isDefined()
      {
         if (this.coords == null)  return false;
         for (int i = 0; i < 3; i++)  if (this.coords[i] == null)  return false;
         return true;
      }

      // aligning to other Cartesians types (even child types)
      public void alignTo(Cartesians C)
      {
         try
         {
            if (C == null) throw new IllegalArgumentException("Cartesians type is null");
            if (!C.isDefined()) throw new IllegalArgumentException("Cartesians type coordinates not defined");
         }
         catch (Exception e)
         {
            e.printStackTrace();
            System.exit(1);
         }

         // this is not a simple copy, bounds may be enforced here
         if (this.coords == null)  this.coords = BoundedDouble.allocateArray(3);
         this.coords[0].setValue(C.coords[0].getValue());
         this.coords[1].setValue(C.coords[1].getValue());
         this.coords[2].setValue(C.coords[2].getValue());
      }

      // destroying this instance
      public void destroy()
      {
         this.coords = null;
      }
   }

   // Spherical type
   private class Spherical extends Cartesians
   {
      // additional internal attributes
      Coordinates R;  // reference for spherical coordinates
      BoundedDouble distance;  // distance to reference
      BoundedDouble theta;  // theta angle
      BoundedDouble phi;  // phi angle

      // Spherical constructor
      public Spherical(Coordinates X,Coordinates R)
      {
         super(X);
         try
         {
            if (R == null) throw new IllegalArgumentException("The Coordinates reference for spherical type is null");
            if (!R.areCartesiansDefined()) throw new IllegalArgumentException("The Coordinates reference has no Cartesian type");
            this.R = R;
            this.distance = null;
            this.theta = null;
            this.phi = null;

            // verifying whether there exists a path from R already leading to this
            List<Coordinates> successors = new LinkedList<Coordinates> ();
            successors.add(R);
            while (!successors.isEmpty())
            {
               Coordinates current = successors.get(0);
               List<Coordinates> tmp = current.listOfCoordinatesDependingOnThis();
               if (tmp.contains(this))  break;
               successors.remove(current);
               for (Coordinates S : tmp)  if (!successors.contains(S))  successors.add(S);
            }
            if (successors.isEmpty())  R.referenceFor(this);
         }
         catch (Exception e)
         {
            e.printStackTrace();
            System.exit(1);
         }
      }

      // verifying whether the coordinates are defined
      @Override
      public boolean isDefined()
      {
         if (this.distance == null)  return false;
         if (this.theta == null)  return false;
         if (this.phi == null)  return false;
         return true;
      }

      // computing the internal Cartesian coordinates
      public void compute(boolean fromSpherical2Cartesians)
      {
         if (fromSpherical2Cartesians)
         {
            try
            {
               if (!this.isDefined()) throw new IllegalStateException("This Cartesians type has no defined coordinates");
            }
            catch (Exception e)
            {
               e.printStackTrace();
               System.exit(1);
            }

            // computing the Cartesian coordinates
            if (this.coords == null)  this.coords = BoundedDouble.allocateArray(3);
            double[] ref = this.R.getCartesians();
            double distance = this.distance.getValue();
            double theta = this.theta.getValue();
            double phi = this.phi.getValue();
            this.coords[0].setValue(ref[0] + distance*Math.sin(theta)*Math.cos(phi));
            this.coords[1].setValue(ref[1] + distance*Math.sin(theta)*Math.sin(phi));
            this.coords[2].setValue(ref[2] + distance*Math.cos(theta));
         }
         else
         {
            double x = 0.0;
            double y = 0.0;
            double z = 0.0;
            double distance = 0.0;
            double[] ref = this.R.getCartesians();
            try
            {
               if (this.coords == null) throw new IllegalArgumentException("Internal Cartesians coordinates not defined");
               x = this.coords[0].getValue() - ref[0];
               y = this.coords[1].getValue() - ref[1];
               z = this.coords[2].getValue() - ref[2];
               distance = x*x + y*y + z*z;
               if (distance == 0.0) throw new IllegalStateException("Error in Spherical type: distance to the reference is 0");
               distance = Math.sqrt(distance);
            }
            catch (Exception e)
            {
               e.printStackTrace();
               System.exit(1);
            }

            // computing the spherical coordinates
            if (this.distance == null)  this.distance = new BoundedDouble();
            this.distance.setValue(distance);
            if (this.theta == null)  this.theta = new BoundedDouble();
            this.theta.setValue(Math.acos(z/distance));
            double absphi = 0.0;
            if (this.phi == null)  this.phi = new BoundedDouble();
            double zplane = x*x + y*y;
            if (zplane != 0.0)
            {
               absphi = Math.acos(x/Math.sqrt(zplane));
               if (y < 0)  absphi = -absphi;
            }
            this.phi.setValue(absphi);
         }
      }

      // destroying this instance
      @Override
      public void destroy()
      {
         this.coords = null;
         this.distance = null;
         this.theta = null;
         this.phi = null;
         this.R.forget(this);
         this.R = null;
      }
   }

   // Torsion type
   private class Torsion extends Cartesians
   {
      // additional internal attributes
      Coordinates A;  // first reference for quadruplet
      Coordinates B;  // second reference for quadruplet
      Coordinates C;  // third reference for quadruplet
      double cosBCD;  // cosine of vector angle BCD
      double sinBCD;  // sine of vector angle BCD
      double dCD;     // distance between the third reference and the variable coordinates (D)
      BoundedDouble torsion;  // the torsion angle value (possibly with bounds)

      // Torsion constructor
      public Torsion(Coordinates X,Coordinates A,Coordinates B,Coordinates C,double cosBCD,double dCD)
      {
         super(X);
         try
         {
            if (A == null) throw new IllegalArgumentException("The Coordinates reference A is null");
            this.A = A;
            if (B == null) throw new IllegalArgumentException("The Coordinates reference B is null");
            this.B = B;
            if (C == null) throw new IllegalArgumentException("The Coordinates reference C is null");
            this.C = C;
            if (cosBCD < -1.0 || cosBCD > 1.0) throw new IllegalArgumentException("cosBCD is not a cosine");
            this.cosBCD = cosBCD;
            this.sinBCD = Math.sqrt(1.0 - cosBCD*cosBCD);
            if (dCD <= 0.0) throw new IllegalArgumentException("dCD is not a distance value (is either 0 or negative)");
            this.dCD = dCD;
            this.torsion = null;

            // verifying whether there exists a path from A already leading to B, C or this
            List<Coordinates> successors = new LinkedList<Coordinates> ();
            successors.add(A);
            while (!successors.isEmpty())
            {
               Coordinates current = successors.get(0);
               List<Coordinates> tmp = current.listOfCoordinatesDependingOnThis();
               if (tmp.contains(B) || tmp.contains(C) || tmp.contains(this))  break;
               successors.remove(current);
               for (Coordinates S : tmp)  if (!successors.contains(S))  successors.add(S);
            }
            if (successors.isEmpty())  A.referenceFor(this);

            // verifying whether there exists a path from B already leading to A, C or this
            successors = new LinkedList<Coordinates> ();
            successors.add(B);
            while (!successors.isEmpty())
            {  
               Coordinates current = successors.get(0);
               List<Coordinates> tmp = current.listOfCoordinatesDependingOnThis();
               if (tmp.contains(A) || tmp.contains(C) || tmp.contains(this))  break;
               successors.remove(current);
               for (Coordinates S : tmp)  if (!successors.contains(S))  successors.add(S);
            }
            if (successors.isEmpty())  B.referenceFor(this);

            // verifying whether there exists a path from C already leading to A, B or this
            successors = new LinkedList<Coordinates> ();
            successors.add(C);
            while (!successors.isEmpty())
            {  
               Coordinates current = successors.get(0);
               List<Coordinates> tmp = current.listOfCoordinatesDependingOnThis();
               if (tmp.contains(A) || tmp.contains(B) || tmp.contains(this))  break;
               successors.remove(current);
               for (Coordinates S : tmp)  if (!successors.contains(S))  successors.add(S);
            }
            if (successors.isEmpty())  C.referenceFor(this);
         }
         catch (Exception e)
         {
            e.printStackTrace();
            System.exit(1);
         }
      }

      // verifying whether the coordinates are defined
      @Override
      public boolean isDefined()
      {
         return this.torsion != null;
      }

      // computing the internal Cartesian coordinates
      public void compute(boolean fromTorsion2Cartesians)
      {
         // extracting the Cartesian coordinates from the reference Coordinates instances
         double[] a = this.A.getCartesians();
         double[] b = this.B.getCartesians();
         double[] c = this.C.getCartesians();

         if (fromTorsion2Cartesians)
         {
            try
            {
               if (!this.isDefined()) throw new IllegalArgumentException("Torsion type coordinates not defined");
            }
            catch (Exception e)
            {
               e.printStackTrace();
               System.exit(1);
            }

            // x axis (first column of U matrix)
            double[] v1 = Vector.diff(c,b);

            // z axis (third column)
            double[] v2 = Vector.diff(a,b);
            double[] v3 = Vector.crossProduct(v1,v2);
            Vector.normalize(v1);
            Vector.normalize(v3);

            // y axis (second column)
            v2 = Vector.crossProduct(v3,v1);
            Vector.normalize(v2);

            // computing vector v (depends on the vector and torsion angles)
            double[] v = new double[3]; 
            v[0] = -this.dCD*this.cosBCD;
            v[1] =  this.dCD*this.sinBCD*Math.cos(this.torsion.getValue());
            v[2] =  this.dCD*this.sinBCD*Math.sin(this.torsion.getValue());

            // generation of the Cartesian coordinates
            if (this.coords == null)  this.coords = BoundedDouble.allocateArray(3);
            this.coords[0].setValue(c[0] + v[0]*v1[0] + v[1]*v2[0] + v[2]*v3[0]);
            this.coords[1].setValue(c[1] + v[0]*v1[1] + v[1]*v2[1] + v[2]*v3[1]);
            this.coords[2].setValue(c[2] + v[0]*v1[2] + v[1]*v2[2] + v[2]*v3[2]);
         }
         else
         {
            try
            {  
               if (this.coords == null) throw new IllegalArgumentException("Internal Cartesians coordinates not defined");
            }  
            catch (Exception e)
            {
               e.printStackTrace();
               System.exit(1);
            }

            // extracting the Cartesian coordinates of D
            double[] d = new double[3];
            d[0] = this.coords[0].getValue();
            d[1] = this.coords[1].getValue();
            d[2] = this.coords[2].getValue();

            try
            {
               if (Math.abs(Vector.EuclideanDistance(c,d) - this.dCD) > epsilon)
                  throw new IllegalStateException("The original distance dCD for the torsion angle is not satisfied");
               if (Math.abs(Vector.cos(b,c,d) - cosBCD) > epsilon)
                  throw new IllegalStateException("The original cosine BCD for the torsion angle is not satisfited");
            }
            catch (Exception e)
            {
               e.printStackTrace();
               System.exit(1);
            }

            // computing the torsion angle by simply applying the definition
            double[] n1 = Vector.normal2plane(a,b,c);
            double[] n2 = Vector.normal2plane(b,c,d);
            double[] bc = Vector.diff(c,b);
            Vector.normalize(bc);
            double[] n1bc = Vector.crossProduct(n1,bc);
            double x = Vector.dotProduct(n1,n2);
            double y = Vector.dotProduct(n1bc,n2);
            if (this.torsion == null)  this.torsion = new BoundedDouble();
            this.torsion.setValue(-Math.atan2(y,x));
         }
      }

      // destroying this instance
      @Override
      public void destroy()
      {
         this.coords = null;
         this.torsion = null;
         this.A.forget(this);
         this.B.forget(this);
         this.C.forget(this);
         this.A = null;
         this.B = null;
         this.C = null;
      }
   }

   /* Coordinates class */

   // static attribute
   static double epsilon = 1.e-6;

   // Coordinates attributes
   private int counter;
   private HashMap<String,Cartesians> representations;
   private HashSet<Cartesians> successors;

   // Coordinates constructor
   public Coordinates()
   {
      this.counter = 1;
      this.representations = new HashMap<String,Cartesians> ();
      this.successors = new HashSet<Cartesians> ();
   }

   // verifying whether the 'pure' Cartesians type is in use
   public boolean containsCartesians()
   {
      return this.representations.containsKey("Cartesians0");
   }

   // verifying whether the Spherical type is in use
   public boolean containsSpherical()
   {
      for (String name : this.representations.keySet())
      {
         if (name.contains("Spherical"))  return true;
      }
      return false;
   }

   // verifying whether the Torsion type is in use
   public boolean containsTorsion()
   {
      for (String name : this.representations.keySet())
      {
         if (name.contains("Torsion"))  return true;
      }
      return false;
   }

   // counting the total number of types that are in use
   public int numberOfTypes()
   {
      return this.representations.size();
   }

   // counting the number of Spherical types that are in use
   public int numberOfSpherical()
   {
      int n = 0;
      for (String name : this.representations.keySet())
      {
         if (name.contains("Spherical"))  n++;
      }
      return n;
   }

   // counting the number of Torsion types that are in use
   public int numberOfTorsion()
   {
      int n = 0;
      for (String name : this.representations.keySet())
      {
         if (name.contains("Torsion"))  n++;
      }
      return n;
   }

   // selecting a random in-use type
   private Cartesians getAnyType()
   {
      try
      {
         if (this.numberOfTypes() == 0) throw new IllegalStateException("No types are currently in use");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      Map.Entry<String,Cartesians> anyEntry = this.representations.entrySet().iterator().next();
      String key = anyEntry.getKey();
      Cartesians value = anyEntry.getValue();
      return value;
   }

   // verifying whether the type whose name is given by the String is defined
   public boolean isDefined(String typeName)
   {
      this.checkTypeName(typeName);
      Cartesians C = this.representations.get(typeName);
      return C.isDefined();
   }

   // verifying whether the 'pure' Cartesians type is defined
   public boolean areCartesiansDefined()
   {
      this.checkTypeName("Cartesians0");
      Cartesians C = this.representations.get("Cartesians0");
      return C.coords != null;
   }

   // verifying whether all currently used types are defined
   public boolean allTypesDefined()
   {
      for (String name : this.representations.keySet())
      {
         Cartesians C = this.representations.get(name);
         if (!C.isDefined())  return false;
      }
      return true;
   }

   // adding the Cartesian representation
   public String addCartesians()
   {
      try
      {
         if (this.containsCartesians()) throw new IllegalStateException("This Coordinates instance can only contain one Cartesians type");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      String typeName = "Cartesians0";
      Cartesians C = new Cartesians(this);
      if (this.numberOfTypes() > 0 && this.allTypesDefined() && this.isCoherent())
      {
         Cartesians refC = this.getAnyType();
         Coordinates.update(refC,C);
      }
      this.representations.put(typeName,C);
      return typeName;
   }

   // adding a Spherical representation
   public String addSpherical(Coordinates R)
   {
      String typeName = "Spherical" + this.counter;
      Spherical S = new Spherical(this,R);
      if (this.numberOfTypes() > 0 && this.allTypesDefined() && this.isCoherent())
      {
         Cartesians refC = this.getAnyType();
         Coordinates.update(refC,S);
      }
      this.representations.put(typeName,S);
      this.counter++;
      return typeName;
   }

   // adding a Torsion representation
   public String addTorsion(Coordinates A,Coordinates B,Coordinates C,double cosBCD,double dCD)
   {
      String typeName = "Torsion" + this.counter;
      Torsion T = new Torsion(this,A,B,C,cosBCD,dCD);
      if (this.numberOfTypes() > 0 && this.allTypesDefined() && this.isCoherent())
      {
         Cartesians refC = this.getAnyType();
         Coordinates.update(refC,T);
      }
      this.representations.put(typeName,T);
      this.counter++;
      return typeName;
   }

   // adding the Cartesian representation with bounds
   public String addCartesians(double[] lb,double[] ub)
   {
      try
      {
         if (lb == null) throw new IllegalArgumentException("The array of double containing the lower bounds is null");
         if (ub == null) throw new IllegalArgumentException("The array of double containing the upper bounds is null");
         if (lb.length != 3) throw new IllegalArgumentException("The array of double containing the lower bounds has size different from 3");
         if (ub.length != 3) throw new IllegalArgumentException("The array of double containing the upper bounds has size different from 3");
         for (int i = 0; i < 3; i++)  if (lb[i] > ub[i])
             throw new IllegalArgumentException("Some given lower bounds are strictly greater than the corresponding upper bounds");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      String typeName = "Cartesians0";
      Cartesians C = new Cartesians(this);
      C.coords = new BoundedDouble[3];
      C.coords[0] = new BoundedDouble(lb[0],ub[0]);
      C.coords[1] = new BoundedDouble(lb[1],ub[1]);
      C.coords[2] = new BoundedDouble(lb[2],ub[2]);
      if (this.numberOfTypes() > 0 && this.allTypesDefined() && this.isCoherent())
      {
         Cartesians refC = this.getAnyType();
         Coordinates.update(refC,C);
      }
      else
      {
         C.coords[0].recenter();
         C.coords[1].recenter();
         C.coords[2].recenter();
      }
      this.representations.put(typeName,C);
      return typeName;
   }

   // adding a Spherical representation with bounds on the distance
   public String addSpherical(Coordinates R,double lb,double ub)
   {
      try
      {
         if (lb > ub) throw new IllegalArgumentException("The lower bound is strictly greater than the upper bound");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      String typeName = "Spherical" + this.counter;
      Spherical S = new Spherical(this,R);
      S.distance = new BoundedDouble(lb,ub);
      if (this.numberOfTypes() > 0 && this.allTypesDefined() && this.isCoherent())
      {
         Cartesians refC = this.getAnyType();
         Coordinates.update(refC,S);
      }
      else
      {
         S.distance.recenter();
      }
      this.representations.put(typeName,S);
      this.counter++;
      return typeName;
   }

   // adding a Torsion representation with bounds on the torsion angle
   public String addTorsion(Coordinates A,Coordinates B,Coordinates C,double cosBCD,double dCD,double lb,double ub)
   {

      try
      {
         if (lb > ub) throw new IllegalArgumentException("The lower bound is strictly greater than the upper bound");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      String typeName = "Torsion" + this.counter;
      Torsion T = new Torsion(this,A,B,C,cosBCD,dCD);
      T.torsion = new BoundedDouble(lb,ub);
      if (this.numberOfTypes() > 0 && this.allTypesDefined() && this.isCoherent())
      {
         Cartesians refC = this.getAnyType();
         Coordinates.update(refC,T);
      }
      else
      {
         T.torsion.recenter();
      }
      this.representations.put(typeName,T);
      this.counter++;
      return typeName;
   }

   // setting the Cartesian representation
   public void setCartesians(double[] coords)
   {
      try
      {
         if (!this.containsCartesians()) throw new IllegalStateException("This Coordinates instance does not contain any Cartesians type");
         if (coords == null) throw new IllegalArgumentException("The array of double containing the cartesian coordinates is null");
         if (coords.length != 3) throw new IllegalArgumentException("The array of double has size different from 3");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      Cartesians C = this.representations.get("Cartesians0");
      if (C.coords == null)  C.coords = BoundedDouble.allocateArray(3);
      for (int i = 0; i < 3; i++)  C.coords[i].setValue(coords[i]);
      this.updateFrom(C);
      this.notification();
   }

   // setting the distance in the Spherical representation
   public void setSphericalDistance(String typeName,double distance)
   {
      try
      {
         this.checkTypeName(typeName);
         if (distance <= 0) throw new IllegalArgumentException("Input distance value is nonpositive");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      Spherical S = (Spherical) this.representations.get(typeName);
      if (S.distance == null)  S.distance = new BoundedDouble();
      S.distance.setValue(distance);
      if (S.isDefined())
      {
         S.compute(true);
         this.updateFrom(S);
         this.notification();
      }
   }

   // setting the angles in the Spherical representation
   public void setSphericalAngles(String typeName,double theta,double phi)
   {
      try
      {
         this.checkTypeName(typeName);
         if (!isAngle(theta)) throw new IllegalArgumentException("Input theta value is not an angle");
         if (!isAngle(phi)) throw new IllegalArgumentException("Input phi value is not an angle");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      Spherical S = (Spherical) this.representations.get(typeName);
      if (S.theta == null)  S.theta = new BoundedDouble();
      if (S.phi == null)  S.phi = new BoundedDouble();
      S.theta.setValue(theta);
      S.phi.setValue(phi);
      if (S.isDefined())
      {
         S.compute(true);
         this.updateFrom(S);
         this.notification();
      }
   }

   // setting the torsion angle
   public void setTorsion(String typeName,double torsion)
   {
      try
      {
         this.checkTypeName(typeName);
         if (!isAngle(torsion)) throw new IllegalArgumentException("Input torsion value is not an angle");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      Torsion T = (Torsion) this.representations.get(typeName);
      if (T.torsion == null)  T.torsion = new BoundedDouble();
      T.torsion.setValue(torsion);
      T.compute(true);
      this.updateFrom(T);
      this.notification();
   }

   // getting the Cartesian coordinates
   public double[] getCartesians()
   {
      try
      {
         if (!this.containsCartesians()) throw new IllegalStateException("This Coordinates instance does not contain the Cartesians type");
         if (!this.areCartesiansDefined()) throw new IllegalStateException("The Cartesian coordinates in this instance are not defined");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      double[] coords = new double[3];
      Cartesians C = this.representations.get("Cartesians0");
      coords[0] = C.coords[0].getValue();
      coords[1] = C.coords[1].getValue();
      coords[2] = C.coords[2].getValue();
      return coords;
   }

   // getting the distance value in the Spherical type
   public double getSphericalDistance(String typeName)
   {
      Spherical S = this.isSphericalType(typeName);
      try
      {
         if (S.distance == null) throw new IllegalStateException("The distance in the Spherical type is not defined");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      return S.distance.getValue();
   }

   // getting the two angles in the Spherical type
   public double[] getSphericalAngles(String typeName)
   {
      Spherical S = this.isSphericalType(typeName);
      try
      {
         if (S.theta == null || S.phi == null)
            throw new IllegalArgumentException("At least one of the two angles in the Spherical type is not defined");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      double[] angles = new double[2];
      angles[0] = S.theta.getValue();
      angles[1] = S.phi.getValue();
      return angles;
   }

   // getting the torsion angle in the Torsion type
   public double getTorsion(String typeName)
   {
      Torsion T = this.isTorsionType(typeName);
      try
      {
         if (T.torsion == null) throw new IllegalArgumentException("The torsion angle in the Torsion type is not defined");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      return T.torsion.getValue();
   }

   // removing the Cartesian representation
   public void removeCartesians()
   {
      try
      {
         if (!this.containsCartesians()) throw new IllegalStateException("This Coordinates instance does not contain any Cartesians type");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      Cartesians C = this.representations.get("Cartesians0");
      C.destroy();
      this.representations.remove("Cartesians0");
   }

   // removing a Spherical representation
   public void removeSpherical(String typeName)
   {
      Spherical S = this.isSphericalType(typeName);
      S.destroy();
      this.representations.remove(typeName);
   }

   // removing a Torsion representation
   public void removeTorsion(String typeName)
   {
      Torsion T = this.isTorsionType(typeName);
      T.destroy();
      this.representations.remove(typeName);
   }

   // getting the "expected" Cartesian coordinates
   public double[] getExpectedCartesians()
   {
      try
      {
         if (!this.allTypesDefined()) throw new IllegalStateException("Not all types are defined");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      // initialization
      double[] coords = new double[3];
      for (int i = 0; i < 3; i++)  coords[i] = 0.0;

      // the computation depends on the coherence state
      if (this.isCoherent())
      {
         Cartesians C = this.getAnyType();
         for (int i = 0; i < 3; i++)  coords[i] = C.coords[i].getValue();
      }
      else
      {
         int r = 0;
         for (String name : this.representations.keySet())
         {
            Cartesians C = this.representations.get(name);
            for (int i = 0; i < 3; i++)  coords[i] = coords[i] + C.coords[i].getValue();
            r++;
         }
         for (int i = 0; i < 3; i++)   coords[i] = coords[i]/r;
      }
      return coords;
   }

   // verifying whether all in-use and defined coordinate types are coherent
   public boolean isCoherent(double eps)
   {
      try
      {
         if (eps < 0.0) throw new IllegalArgumentException("The tolerance epsilon is negative");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      // checking if all types are defined
      Set<String> typeNames = this.representations.keySet();
      if (!this.allTypesDefined()) throw new IllegalStateException("Not all coordinates types are defined");

      // each type against all others
      for (String type1 : typeNames)
      {
         Cartesians C1 = this.representations.get(type1);
         for (String type2 : typeNames)
         {
            Cartesians C2 = this.representations.get(type2);
            if (C1 != C2)
            {
               if (Math.abs(C2.coords[0].getValue() - C1.coords[0].getValue()) > eps)  return false;
               if (Math.abs(C2.coords[1].getValue() - C1.coords[1].getValue()) > eps)  return false;
               if (Math.abs(C2.coords[2].getValue() - C1.coords[2].getValue()) > eps)  return false;
            }
         }
      }

      return true;
   }

   // verifying whether the various coordinates types are coherent (with fixed precision)
   public boolean isCoherent()
   {
      return this.isCoherent(epsilon);
   }

   // associating this Coordinates instance to an internal Cartesians type, which is supposed to use this Coordinates instance as a reference
   private void referenceFor(Cartesians C)
   {
      try
      {
         if (C == null) throw new IllegalArgumentException("The input Cartesians type is null");
         if (this.representations.containsValue(C)) throw new IllegalArgumentException("The input Cartesians type belongs to this Coordinates instance");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      this.successors.add(C);
   }

   // verifying whether this Coordinates instance has an associated Cartesian type from the input instance
   private boolean isReferenceFor(Coordinates other)
   {
      for (Cartesians C : this.successors)  if (C.mainInstance == other)  return true;
      return false;
   }

   // forgetting the input Coordinates instance successor
   private void forget(Cartesians C)
   {
      this.successors.remove(C);
   }

   // constructing a list of Coordinates instances this instance is reference for
   public List<Coordinates> listOfCoordinatesDependingOnThis()
   {
      List<Coordinates> list = new ArrayList<Coordinates> ();
      for (Cartesians C : this.successors)  list.add(C.mainInstance);
      return list;
   }

   // updating a given representation from another one
   private static void update(Cartesians reference,Cartesians C)
   {
      try
      {
         if (reference == null) throw new IllegalArgumentException("The reference Cartesians instance is null");
         if (C == null) throw new IllegalArgumentException("The Cartesians instance to update is null");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      // copying the Cartesian coordinates from the reference
      C.alignTo(reference);

      // updating by using the appropriate conversion method
      if (C instanceof Spherical)
      {
         Spherical S = (Spherical) C;
         S.compute(false);
         S.compute(true);
      }
      else if (C instanceof Torsion)
      {
         Torsion T = (Torsion) C;
         T.compute(false);
         T.compute(true);
      }
      // pure Cartesian types do not need to invoke 'compute'
   }

   // updating the various representations from a given one
   private void updateFrom(String typeName)
   {
      // checking typeName
      this.checkTypeName(typeName);

      // updating
      Cartesians refC = this.representations.get(typeName);
      this.updateFrom(refC);
   }

   // updating the various representations from a given one
   private void updateFrom(Cartesians refC)
   {
      try
      {
         if (refC == null) throw new IllegalArgumentException("The reference Cartesians instance is null");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      // updating
      for (String name : this.representations.keySet())
      {
         Cartesians C = this.representations.get(name);
         if (refC != C)
         {
            Coordinates.update(refC,C);
         }
      }
   }

   // notifying the successors (they should recompute their inner variables)
   private void notification()
   {
      for (Cartesians C : this.successors)
      {
         if (C != null)
         {
            Coordinates Main = C.mainInstance;
            if (Main != null)
            {
               if (C instanceof Spherical)
               {
                  Spherical S = (Spherical) C;
                  S.compute(true);
               }
               else if (C instanceof Torsion)
               {
                  Torsion T = (Torsion) C;
                  T.compute(true);
               }
               Main.updateFrom(C);
               Main.notification();
            }
         }
      }
   }

   // checking the input names of Coordinates types
   private void checkTypeName(String typeName)
   {
      try
      {
         if (typeName == null) throw new IllegalArgumentException("The String supposed to contain the type name is null");
         if (typeName.isEmpty()) throw new IllegalArgumentException("The String supposed to contain the type name is empty");
         if (!typeName.contains("Cartesians") && !typeName.contains("Spherical") && !typeName.contains("Torsion"))
            throw new IllegalArgumentException("The given type name '" + typeName + "' is not valid");
         if (!this.representations.containsKey(typeName))
            throw new IllegalArgumentException("The given type name '" + typeName + "' is unknown");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }
   }

   // checking whether a given type belongs to the Spherical group
   private Spherical isSphericalType(String typeName)
   {
      Spherical S = null;
      try
      {
         this.checkTypeName(typeName);
         Cartesians C = this.representations.get(typeName);
         boolean isSpherical = (C instanceof Spherical);
         if (!isSpherical) throw new IllegalArgumentException("The type name '" + typeName + "' does not correspond to Spherical");
         S = (Spherical) C;
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      return S;
   }

   // checking whether a given type belongs to the Torsion group
   private Torsion isTorsionType(String typeName)
   {
      Torsion T = null;
      try
      {
         this.checkTypeName(typeName);
         Cartesians C = this.representations.get(typeName);
         boolean isTorsion = (C instanceof Torsion);
         if (!isTorsion) throw new IllegalArgumentException("The type name '" + typeName + "' does not correspond to Torsion");
         T = (Torsion) C;
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      return T;
   }

   // toString
   public String toString()
   {
      String info = "Cartesians ";
      if (this.containsCartesians())
         info = info + "1; ";
      else
         info = info + "0; ";
      info = info + "Spherical " + this.numberOfSpherical() + "; ";
      info = info + "Torsion " + this.numberOfTorsion();
      try
      {
         if (this.isCoherent())  info = info + "; coherent";
      }
      catch (IllegalStateException E)
      {
         info = info + "; with undefined types";
      }
      return "Coordinates [" + info + "]";
   }

   /* static methods */

   // verifying whether a given real value lies in the interval [-2pi,2pi]
   public static boolean isAngle(double value)
   {
      if (!Double.isFinite(value))  return false;
      double twopi = 2.0*Math.PI;
      return value >= -twopi && value <= twopi;
   }

   // projecting a given double value on the internal [-2pi,pi]
   public static double closestAngle(double value)
   {
      double twopi = 2.0*Math.PI;
      if (value < -twopi)  return -twopi;
      if (value >  twopi)  return  twopi;
      return value;
   }

   // main -- performing some basic tests
   public static void main(String[] args) throws Exception
   {
      System.out.print("Coordinates class ");

      // seed for random numbers
      long seed;
      Random R = new Random();
      if (args != null && args.length > 0)
         seed = Long.parseLong(args[0]);
      else
         seed = Math.abs(R.nextLong());
      System.out.print("(seed is " + seed + ")");
      R = new Random(seed);

      // isAngle (static method)
      Exception E = new Exception("public static boolean isAngle(double)");
      if (!isAngle(0.0)) throw E;
      if (!isAngle(3.0)) throw E;
      if (!isAngle(-3.0 + R.nextDouble())) throw E;
      if (!isAngle(R.nextDouble())) throw E;
      if (isAngle(100.0 - R.nextDouble())) throw E;
      if (isAngle(Double.MAX_VALUE)) throw E;
      if (isAngle(Double.POSITIVE_INFINITY)) throw E;

      // closestAngle (static method)
      E = new Exception("public static double closestAngle(double)");
      if (!isAngle(closestAngle(-100.0))) throw E;
      if (!isAngle(closestAngle(200.0))) throw E;
      if (!isAngle(closestAngle(-Double.POSITIVE_INFINITY))) throw E;

      // getting ready to test the other methods
      System.out.print(" ... ");
      int NTESTS = 100000;
      for (int itests = 0; itests < NTESTS; itests++)
      {
         // Coordinates B
         // Basic example, only with Cartesian type, no references
         Coordinates B = new Coordinates();
         if (B == null) throw new Exception("public Coordinates()");
         if (B.counter != 1) throw new Exception("public Coordinates()");
         if (B.representations == null) throw new Exception("public Coordinates()");
         if (B.successors == null) throw new Exception("public Coordinates()");
         if (B.numberOfTypes() != 0) throw new Exception("public int numberOfTypes()");
         if (B.containsCartesians()) throw new Exception("public boolean containsCartesians()");
         String CartesianName = B.addCartesians();
         if (CartesianName != "Cartesians0") throw new Exception("public String addCartesians()");
         if (!B.containsCartesians()) throw new Exception("public boolean containsCartesians()");
         if (B.numberOfTypes() != 1) throw new Exception("public int numberOfTypes()");
         if (B.areCartesiansDefined()) throw new Exception("public boolean areCartesiansDefined()");
         if (B.isDefined(CartesianName)) throw new Exception("public boolean isDefined(String)");
         double[] Bcoords = Vector.random(3,R);
         B.setCartesians(Bcoords);
         if (!B.isDefined(CartesianName)) throw new Exception("public boolean isDefined(String)");
         if (!B.areCartesiansDefined()) throw new Exception("public boolean areCartesiansDefined()");

         // Coordinates C
         // Cartesian and Spherical types, where B is the reference for Spherical type
         Coordinates C = new Coordinates();
         if (C == null) throw new Exception("public Coordinates()");
         CartesianName = C.addCartesians();
         if (!C.containsCartesians()) throw new Exception("public boolean containsCartesians()");
         double distance = 0.1 + 0.8*R.nextDouble();
         double[] Ccoords = Vector.random(3,R);
         Vector.normalize(Ccoords);
         C.setCartesians(Vector.linear(1.0,Bcoords,distance,Ccoords));
         if (!C.isDefined(CartesianName)) throw new Exception("public boolean isDefined(String)");
         if (C.containsSpherical()) throw new Exception("public boolean containsSpherical()");
         if (C.numberOfSpherical() > 0) throw new Exception("public int numberOfSpherical()");
         String SphericalName = C.addSpherical(B);
         if (SphericalName == null) throw new Exception("public String addSpherical(Coordinates)");
         if (!SphericalName.contains("Spherical")) throw new Exception("public String addSpherical(Coordinates)");
         if (C.numberOfTypes() != 2) throw new Exception("public int numberOfTypes()");
         if (!C.containsSpherical()) throw new Exception("public boolean containsSpherical()");
         if (C.numberOfSpherical() != 1) throw new Exception("public int numberOfSpherical()");
         if (!C.isDefined(SphericalName)) throw new Exception("public boolean isDefined(Coordinates)");
         if (Math.abs(distance - C.getSphericalDistance(SphericalName)) > epsilon) throw new Exception("public double getSphericalDistance(String)");
         C.setSphericalDistance(SphericalName,1.0);
         if (!C.isDefined(SphericalName)) throw new Exception("public boolean isDefined(Coordinates)");
         double[] obtained = C.getCartesians();
         double[] expected = Vector.linear(1.0,Bcoords,1.0,Ccoords);
         if (Vector.EuclideanDistance(obtained,expected) > epsilon) throw new Exception("Internal updating system");
         Bcoords = Vector.sum(Bcoords,new double[]{0.1*R.nextDouble(),0.1*R.nextDouble(),0.1*R.nextDouble()});
         B.setCartesians(Bcoords);
         if (!B.areCartesiansDefined()) throw new Exception("public boolean areCartesiansDefined()");
         if (!C.areCartesiansDefined()) throw new Exception("public boolean areCartesiansDefined()");
         if (!C.allTypesDefined()) throw new Exception("Internal notification system");
         if (!C.isCoherent()) throw new Exception("Internal notification system");
         distance = 1.0;
         obtained = C.getCartesians();
         expected = Vector.linear(1.0,Bcoords,1.0,Ccoords);
         if (Vector.EuclideanDistance(obtained,expected) > epsilon)  throw new Exception("Internal notification system");
         String OtherSpherical = C.addSpherical(B,distance - 0.1*R.nextDouble(),distance + 0.1*R.nextDouble());
         if (C.numberOfSpherical() < 2) throw new Exception("public String addSpherical(Coordinates,double,double)");
         if (!C.isDefined(OtherSpherical)) throw new Exception("public String addSpherical(Coordinates,double,double)");
         obtained = C.getCartesians();
         if (!C.isCoherent()) throw new Exception("public String addSpherical(Coordinates,double,double)");
         if (Vector.EuclideanDistance(obtained,expected) > epsilon) throw new Exception("public String addSpherical(Coordinates,double,double)");
         C.removeSpherical(OtherSpherical);
         if (C.numberOfSpherical() >= 2) throw new Exception("public void removeSpherical(String)");
         OtherSpherical = C.addSpherical(B,distance - 0.20,distance - 0.10 - 0.09*R.nextDouble());
         if (!C.isDefined(OtherSpherical)) throw new Exception("public String addSpherical(Coordinates,double,double)");
         if (C.isCoherent()) throw new Exception("public String addSpherical(Coordinates,double,double)");
         C.removeSpherical(OtherSpherical);
         Ccoords = expected;

         // Coordinates A
         // Cartesians to be used as a reference for B, which will have an impact on C as well
         Coordinates A = new Coordinates();
         if (A == null) throw new Exception("public Coordinates()");
         CartesianName = A.addCartesians();
         if (A.areCartesiansDefined()) throw new Exception("public boolean areCartesiansDefined()");
         double[] Acoords = new double[]{0.0,0.0,0.0};
         A.setCartesians(Acoords);
         if (!A.areCartesiansDefined()) throw new Exception("public boolean areCartesiansDefined()");
         if (!A.allTypesDefined()) throw new Exception("public boolean allTypesDefined()");
         if (!A.isCoherent()) throw new Exception("public boolean isCoherent()");
         SphericalName = B.addSpherical(A);
         if (B.numberOfTypes() != 2) throw new Exception("public int numberOfTypes()");
         if (!B.allTypesDefined()) throw new Exception("public boolean allTypesDefined()");
         if (!B.isCoherent()) throw new Exception("public boolean isCoherent()");
         double[] translation = new double[]{R.nextDouble(),R.nextDouble(),R.nextDouble()};
         Acoords = Vector.sum(Acoords,translation);
         Bcoords = Vector.sum(Bcoords,translation);
         Ccoords = Vector.sum(Ccoords,translation);
         A.setCartesians(Acoords);
         if (!A.isCoherent()) throw new Exception("public boolean isCoherent()");
         if (!B.isCoherent()) throw new Exception("public boolean isCoherent()");
         if (Vector.EuclideanDistance(Bcoords,B.getCartesians()) > epsilon) throw new Exception("Internal notification system");
         if (!C.isCoherent()) throw new Exception("public boolean isCoherent()");
         if (Vector.EuclideanDistance(Ccoords,C.getCartesians()) > epsilon) throw new Exception("Internal notification system");
         if (Vector.EuclideanDistance(C.getCartesians(),C.getExpectedCartesians()) > 0.0)  throw new Exception("public double[] getExpectedCartesians()");
         double[] angles = B.getSphericalAngles(SphericalName);
         int k = R.nextInt(2);
         double rotation = 0.1 - 0.2*R.nextDouble();
         double newangle = closestAngle(angles[k] + rotation);
         rotation = newangle - angles[k];
         angles[k] = newangle;
         B.setSphericalAngles(SphericalName,angles[0],angles[1]);
         if (!B.allTypesDefined()) throw new Exception("public boolean allTypesDefined()");
         if (!B.isCoherent()) throw new Exception("public boolean isCoherent()");
         double[] newBcoords = B.getCartesians();
         if (Math.abs(Math.cos(rotation) - Vector.cos(Bcoords,Acoords,newBcoords)) > 0.1) throw new Exception("Internal notification system");
         Bcoords = newBcoords;
         Ccoords = C.getCartesians();

         // Coordinates D
         // Cartesian and Torsion angle types
         Coordinates D = new Coordinates();
         if (D.containsTorsion()) throw new Exception("public boolean containsTorsion()");
         if (D.numberOfTorsion() > 0) throw new Exception("public int numberOfTorsion()");
         double cosBCD = R.nextDouble();
         double dCD = 0.5 + R.nextDouble();
         String Torsion1 = D.addTorsion(A,B,C,cosBCD,dCD);
         if (!D.containsTorsion()) throw new Exception("public boolean containsTorsion()");
         if (D.numberOfTorsion() != 1) throw new Exception("public int numberOfTorsion()");
         if (D.isDefined(Torsion1)) throw new Exception("public boolean isDefined()");
         if (D.containsCartesians()) throw new Exception("public boolean containsCartesians()");
         CartesianName = D.addCartesians();
         if (!D.containsCartesians()) throw new Exception("public boolean containsCartesians()");
         if (D.areCartesiansDefined()) throw new Exception("public boolean areCartesiansDefined()");
         double torAngle = 3.0*R.nextDouble();
         D.setTorsion(Torsion1,torAngle);
         if (!D.isDefined(Torsion1)) throw new Exception("public boolean isDefined()");
         if (!D.isDefined(CartesianName)) throw new Exception("public boolean isDefined()");
         if (!D.areCartesiansDefined()) throw new Exception("public boolean areCartesiansDefined()");
         if (!D.allTypesDefined()) throw new Exception("public boolean allTypesDefined()");
         if (!D.isCoherent()) throw new Exception("public boolean isCoherent()");
         double[] Dcoords = D.getCartesians();
         String Torsion2 = D.addTorsion(A,B,C,cosBCD,dCD);
         if (D.numberOfTorsion() != 2) throw new Exception("public int numberOfTorsion()");
         if (D.numberOfTypes() != 3) throw new Exception("public int numberOfTypes()");
         if (!D.isDefined(Torsion2)) throw new Exception("public boolean isDefined()");
         if (!D.allTypesDefined()) throw new Exception("public boolean allTypesDefined()");
         if (!D.isCoherent()) throw new Exception("public boolean isCoherent()");
         if (Math.abs(D.getTorsion(Torsion2) - torAngle) > epsilon) throw new Exception("Internal updating system");
         D.removeTorsion(Torsion2);
         if (D.numberOfTorsion() > 1) throw new Exception("public int numberOfTorsion()");
         if (!D.isCoherent()) throw new Exception("public boolean isCoherent()");
         if (torAngle > 0.8 && torAngle < 1.6)
         {
            String OtherTorsion = D.addTorsion(A,B,C,cosBCD,dCD,torAngle - 0.7,torAngle - 0.5);
            if (D.numberOfTorsion() < 2) throw new Exception("public String addTorsion(Coordinates,Coordinates,Coordinates,double,double,double,double)");
            boolean coherence1 = D.isCoherent();
            D.removeTorsion(OtherTorsion);
            if (D.numberOfTorsion() > 2) throw new Exception("public void removeTorsion(String)");
            OtherTorsion = D.addTorsion(A,B,C,cosBCD,dCD,torAngle + 0.5,torAngle + 0.7);
            boolean coherence2 = D.isCoherent();
            if (coherence1 || coherence2) throw new Exception("public String addTorsion(Coordinates,Coordinates,Coordinates,double,double,double,double)");
            if (Vector.EuclideanDistance(D.getCartesians(),D.getExpectedCartesians()) < epsilon)
               throw new Exception("public String addTorsion(Coordinates,Coordinates,Coordinates,double,double,double,double)");
            D.removeTorsion(OtherTorsion);
            if (!D.isCoherent()) throw new Exception("public void removeTorsion(String)");
         }

         // repositioning A in {0,0,0}
         A.setCartesians(new double[]{0.0,0.0,0.0});
         Bcoords = Vector.diff(Bcoords,Acoords);
         if (Vector.EuclideanDistance(B.getCartesians(),Bcoords) > epsilon) throw new Exception("Internal notification system");
         Ccoords = Vector.diff(Ccoords,Acoords);
         if (Vector.EuclideanDistance(C.getCartesians(),Ccoords) > epsilon) throw new Exception("Internal notification system");
         Dcoords = Vector.diff(Dcoords,Acoords);
         if (Vector.EuclideanDistance(D.getCartesians(),Dcoords) > epsilon) throw new Exception("Internal notification system");

         // massively destroying the Coordinates types
         Dcoords = D.getCartesians();
         D.removeTorsion(Torsion1);
         if (D.containsTorsion()) throw new Exception("public boolean containsTorsion()");
         if (!D.isCoherent()) throw new Exception("public boolean isCoherent()");
         C.removeSpherical(SphericalName);
         if (C.containsSpherical()) throw new Exception("public boolean containsSpherical()");
         if (!C.allTypesDefined()) throw new Exception("public boolean allTypesDefined()");
         B.removeSpherical(SphericalName);
         if (B.containsSpherical()) throw new Exception("public boolean containsSpherical()");
         if (B.numberOfSpherical() != 0) throw new Exception("public int numberOfSpherical()");
         if (!B.allTypesDefined()) throw new Exception("public boolean allTypesDefined()");
         if (!B.isCoherent()) throw new Exception("public boolean isCoherent()");
         A.setCartesians(new double[]{R.nextDouble(),0.0,0.0});
         if (Vector.EuclideanDistance(D.getCartesians(),Dcoords) > 0.0) throw new Exception("Internal notification system");
      }
      System.out.println("OK");
   }
}

