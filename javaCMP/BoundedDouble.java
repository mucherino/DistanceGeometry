
/* javaCMP project
 *
 * BoundedDouble class
 *
 * last update: April 23rd, 2023
 *
 * AM
*/

import java.util.Random;

public class BoundedDouble
{
   private double lb;  // the lower bound
   private double ub;  // the upper bound
   private double value;  // the actual value

   // BoundedDouble constructor from the lower and upper bounds
   public BoundedDouble(double lb,double ub)
   {
      try
      {
         if (lb > ub) throw new IllegalArgumentException("Lower bound cannot be strictly greater than upper bound");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      this.lb = lb;
      this.ub = ub;
      this.recenter();
   }

   // BoundedDouble constructor from one bound (the other is set to infinity, choice of infinity bound via the boolean variable)
   public BoundedDouble(double bound,boolean infinityAtUpperBound)
   {
      if (infinityAtUpperBound)
      {
         this.lb = bound;
         this.ub = Double.POSITIVE_INFINITY;
      }
      else
      {
         this.lb = -Double.POSITIVE_INFINITY;
         this.ub = bound;
      }
      this.value = bound;
   }

   // BoundedDouble constructor without bounds (both are set to infinity)
   public BoundedDouble()
   {
      this.lb = -Double.POSITIVE_INFINITY;
      this.ub = Double.POSITIVE_INFINITY;
      this.value = 0.0;
   }

   // BoundedDouble constructor cloning from an existing BoundedDouble instance
   public BoundedDouble(BoundedDouble BD)
   {
      try
      {
         if (BD == null) throw new IllegalArgumentException("The input BoundedDouble instance is null");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      this.lb = BD.lb;
      this.ub = BD.ub;
      this.value = BD.value;
   }

   // setter for the actual double value
   // when it is outside the given bounds (notice use of epsilon tolerance), the given value is projected on the bounds;
   // the stored value (projected or not) is the returning value
   public double setValue(double value,double epsilon)
   {
      try
      {  
         if (epsilon < 0.0) throw new IllegalArgumentException("The tolerance epsilon cannot be negative");
      }
      catch (Exception e)
      {  
         e.printStackTrace();
         System.exit(1);
      }

      if (value < this.lb - epsilon)
      {
         this.value = this.lb;
      }
      else if (value > this.ub + epsilon)
      {
         this.value = this.ub;
      }
      else
      {
         this.value = value;
      }
      
      return this.value;
   }

   // setter for the actual double value
   // when it is outside the given bounds (strict verification), the given value is projected on the bounds;
   // the stored value (projected or not) is the returning value
   public double setValue(double value)
   {
      return this.setValue(value,0.0);
   }

   // getter for lower bound
   public double getLowerBound()
   {
      return this.lb;
   }

   // getter for upper bound
   public double getUpperBound()
   {
      return this.ub;
   }

   // getter for actual value
   public double getValue()
   {
      return this.value;
   }

   // verifying whether this instance is degenerate
   public boolean isDegenerate()
   {
      return this.lb == this.ub;
   }

   // verifying whether this instance admits a given value (not to be confused with this.value)
   // notice the use of the tolerance epsilon
   public boolean admits(double value,double epsilon)
   {
      try
      {
         if (epsilon < 0.0) throw new IllegalArgumentException("The tolerance epsilon cannot be negative");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      if (this.isDegenerate())  return Math.abs(value - this.lb) <= epsilon;
      return this.lb - epsilon < value && value < this.ub + epsilon;
   }

   // verifying whether this instance admits a given value (strict verification)
   public boolean admits(double value)
   {
      return this.admits(value,0.0);
   }

   // positioning the value at the center of the bounds (same rules as in constructors apply)
   public void recenter()
   {
      if (this.isDegenerate())
      {
         this.value = this.lb;
      }
      else if (Double.isFinite(this.lb) && Double.isFinite(this.ub))
      {
         this.value = 0.5*(this.lb + this.ub);
      }
      else if (Double.isFinite(this.lb))
      {
         this.value = this.lb;
      }
      else if (Double.isFinite(this.ub))
      {
         this.value = this.ub;
      }
      else
      {
         this.value = 0.0;
      }
   }

   // computing the distance on the real line from a given real number to this BoundedDouble instance
   public double distanceFrom(double value)
   {
      if (this.admits(value))  return 0.0;
      double dl = Double.MAX_VALUE;
      double du = Double.MAX_VALUE;
      if (value < this.lb)  dl = this.lb - value;
      if (this.ub < value)  du = value - this.ub;
      return Math.min(dl,du);
   }

   // equals method 
   public boolean equals(Object o)
   {
      boolean isBoundedDouble = (o instanceof BoundedDouble);
      if (!isBoundedDouble)  return false;
      BoundedDouble BD = (BoundedDouble) o;
      if (this.lb != BD.lb)  return false;
      if (this.ub != BD.ub)  return false;
      return true;
   }

   // hashCode method
   public int hashCode()
   {
      int hash = 0;
      if (Double.isFinite(this.lb))  hash = hash + Double.valueOf(this.lb).intValue();
      if (Double.isFinite(this.ub))  hash = hash + Double.valueOf(this.ub).intValue();
      return hash;
   }

   // toString method
   public String toString()
   {
      return "[" + this.lb + "(" + this.value + ")" + this.ub + "]";
   }

   // allocating memory for an array of BoundedDouble types (all with infinity bounds)
   public static BoundedDouble[] allocateArray(int n)
   {
      try
      {
         if (n <= 0) throw new IllegalArgumentException("Specified size for array is nonpositive");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      BoundedDouble[] BD = new BoundedDouble[n];
      for (int i = 0; i < n; i++)  BD[i] = new BoundedDouble();
      return BD;
   }

   // printing an array of BoundedDouble
   public static String printArray(BoundedDouble[] BD)
   {
      try
      {
         if (BD == null) throw new IllegalArgumentException("The input array of BoundedDouble is null");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      String print = "BoundedDouble\n";
      for (int i = 0; i < BD.length; i++)
      {
         if (BD[i] == null)
            print = print + "{NULL}\n";
         else
            print = print + BD[i].toString() + "\n";
      }
      return print;
   }

   // main method -- performing some basic tests
   public static void main(String[] args) throws Exception
   {
      System.out.print("BoundedDouble class ");

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
      double eps = 1.e-6;
      boolean print = false;

      for (int itests = 0; itests < NTESTS; itests++)
      {
         int r = R.nextInt(6);
         double lb = 2.0*R.nextDouble() - 1.0;
         double ub = lb;
         while (Math.abs(lb - ub) < eps)  ub = lb + R.nextDouble();
         Exception E = new Exception("public BoundedDouble(double,double)");
         BoundedDouble[] BD = new BoundedDouble[6];
         BD[0] = new BoundedDouble(lb,ub);
         if (BD[0] == null) throw E;
         if (BD[0].lb != lb || BD[0].ub != ub) throw E;
         if (Math.abs(BD[0].value - 0.5*(lb+ub)) > eps) throw E;
         BD[1] = new BoundedDouble(lb,lb);
         if (BD[1] == null) throw E;
         if (BD[1].lb != lb || BD[1].ub != lb) throw E;
         if (BD[1].value != lb) throw E;
         E = new Exception("public BoundedDouble(double,boolean)");
         BD[2] = new BoundedDouble(lb,true);
         if (BD[2] == null) throw E;
         if (BD[2].lb != lb) throw E;
         if (BD[2].ub < ub) throw E;
         if (BD[2].value != lb) throw E;
         BoundedDouble tmp = new BoundedDouble(BD[2].lb,BD[2].ub);
         if (BD[2].value != tmp.value) throw E;
         BD[3] = new BoundedDouble(ub,false);
         if (BD[3] == null) throw E;
         if (BD[3].lb > lb) throw E;
         if (BD[3].ub != ub) throw E;
         if (BD[3].value != ub) throw E;
         tmp = new BoundedDouble(BD[3].lb,BD[3].ub);
         if (BD[3].value != tmp.value) throw E;
         E = new Exception("public BoundedDouble()");
         BD[4] = new BoundedDouble();
         if (BD[4] == null) throw E;
         if (BD[4].lb > lb || BD[4].ub < ub) throw E;
         if (BD[4].value != 0.0) throw E;
         tmp = new BoundedDouble(BD[4].lb,BD[4].ub);
         if (BD[4].value != tmp.value) throw E;
         E = new Exception("public BoundedDouble(BoundedDouble)");
         BD[5] = new BoundedDouble(BD[0]);
         if (BD[5] == null) throw E;
         if (BD[5].lb != BD[0].lb) throw E;
         if (BD[5].ub != BD[0].ub) throw E;
         if (BD[5].value != BD[0].value) throw E;
         if (print && itests == 0)  for (int i = 0; i < 6; i++)  System.out.println(BD[i]);
         double rnd = 0.0;
         while (rnd < eps)  rnd = R.nextDouble();
         E = new Exception("public void setValue(double,double)");
         BD[0].setValue(lb - rnd,eps);
         if (Math.abs(BD[0].value - lb) > eps) throw E;
         BD[1].setValue(ub,eps);
         if (Math.abs(BD[1].value - lb) > eps) throw E;
         BD[2].setValue(ub,eps);
         if (Math.abs(BD[2].value - ub) > eps) throw E;
         BD[3].setValue(lb,eps);
         if (Math.abs(BD[3].value - lb) > eps) throw E;
         BD[4].setValue(rnd,eps);
         if (Math.abs(BD[4].value - rnd) > eps) throw E;
         BD[5].setValue(ub + rnd,eps);
         if (Math.abs(BD[5].value - ub) > eps) throw E;
         E = new Exception("public double getLowerBound()");
         if (BD[0].getLowerBound() != BD[0].lb) throw E;
         if (BD[2].getLowerBound() != BD[2].lb) throw E;
         if (BD[4].getLowerBound() != BD[4].lb) throw E;
         E = new Exception("public double getUpperBound()");
         if (BD[0].getUpperBound() != BD[0].ub) throw E;
         if (BD[3].getUpperBound() != BD[3].ub) throw E;
         if (BD[4].getUpperBound() != BD[4].ub) throw E;
         E = new Exception("public double getValue()");
         if (BD[r].getValue() != BD[r].value) throw E;
         E = new Exception("public boolean isDegenerate()");
         if (BD[0].isDegenerate()) throw E;
         if (!BD[1].isDegenerate()) throw E;
         if (r > 1)  if (BD[r].isDegenerate()) throw E;
         E = new Exception("public boolean admits(double,double)");
         if (BD[0].admits(lb - rnd,eps)) throw E;
         if (lb + eps < ub)  if (!BD[0].admits(lb + eps,eps)) throw E;
         if (BD[0].admits(ub + rnd,eps)) throw E;
         if (!BD[1].admits(lb,eps)) throw E;
         if (BD[1].admits(ub,eps)) throw E;
         if (lb < 0.0)  if (!BD[2].admits(0.0)) throw E;
         if (ub > 0.0)  if (!BD[3].admits(0.0)) throw E;
         if (!BD[4].admits(0.0)) throw E;
         if (lb < ub - eps)  if (!BD[5].admits(ub - eps,eps)) throw E;
         E = new Exception("public void recenter()");
         BD[0].setValue(R.nextDouble());
         BD[5].setValue(R.nextDouble());
         BD[0].recenter();
         BD[5].recenter();
         if (BD[0].getValue() != BD[5].getValue()) throw E;
         E = new Exception("public double distanceFrom(double)");
         if (Double.isFinite(BD[r].lb))  if (BD[r].distanceFrom(BD[r].lb - rnd) > rnd + eps) throw E;
         if (Double.isFinite(BD[r].ub))  if (BD[r].distanceFrom(BD[r].ub + rnd) > rnd + eps) throw E; 
         E = new Exception("public boolean equals(Object)");
         if (BD[r].equals(E)) throw E;
         if (!BD[0].equals(BD[5])) throw E;
         BD[0].setValue(BD[0].ub);
         if (!BD[0].equals(BD[5])) throw E;
         BD[5].setValue(BD[5].lb);
         if (!BD[0].equals(BD[5])) throw E;
         if (r > 1)  if (BD[1].equals(BD[r])) throw E;
      }
      System.out.println("OK");
   }
}

