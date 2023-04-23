
/* Java project
 *
 * Experiments class
 *
 * last update: April 23rd, 2023
 *
 * AM
*/

public class Experiments
{
   // This is the experiment reported in the LNBI conference paper (IWBBIO23, Gran Canaria)
   public static void main(String[] args)
   {
      // creating an alpha-helix with given torsion angles
      System.out.println("Experiment reported on LNBI conference paper (IWBBIO23, Gran Canaria)");
      int nAtoms = 80;
      double bondLength = 1.45;
      double bondAngleDistance = 2.46;
      double cosVectorAngle = Vector.cos(bondLength,bondLength,bondAngleDistance);
      double[] torsionAngles = new double[] {-1.0472,-0.8727,3.14159};
      double[][] clique = Vector.embedTriangle(bondLength,bondLength,bondAngleDistance);
      Coordinates[] C = new Coordinates[nAtoms];
      C[0] = new Coordinates();
      C[0].addCartesians();
      C[0].setCartesians(clique[0]);
      C[1] = new Coordinates();
      C[1].addCartesians();
      C[1].setCartesians(clique[1]);
      C[2] = new Coordinates();
      C[2].addCartesians();
      C[2].setCartesians(clique[2]);
      int next = 3;
      String cut = null;
      while (next < nAtoms)
      {
         C[next] = new Coordinates();
         String typeName = C[next].addTorsion(C[next-3],C[next-2],C[next-1],cosVectorAngle,bondLength);
         if (next == nAtoms/2)  cut = typeName;
         C[next].setTorsion(typeName,torsionAngles[next%3]);
         C[next].addCartesians();
         next++;
      }

      // printing the original and modified helix
      System.out.println("This is the original helix");
      System.out.println(printPDB(C));
      System.out.println("This is the modified helix (change in only one torsion angle)");
      C[60].setTorsion("Torsion1",Math.PI);
      System.out.println(printPDB(C));
   }

   // printing in PDB format (for the visualization of the molecules)
   public static String printPDB(Coordinates[] C)
   {
      String PDB = "HEADER    JAVA MULTIREPRESENTATION PROJECT COMMIT ONE\n";
      PDB = PDB + "REMARK 1   CODED BY AM\n";
      try
      {
         if (C == null) throw new IllegalArgumentException("The array of Coordinates instances is null");
         for (int i = 0; i < C.length; i++)
         {
            if (C[i] == null) throw new IllegalArgumentException("The Coordinates instance with index " + i + " is null");
            if (!C[i].allTypesDefined()) throw new IllegalArgumentException("The Coordinates instance with index " + i + " has no all types defined");
            if (!C[i].isCoherent()) throw new IllegalArgumentException("The Coordinates instance with index " + i + " is not coherent");
            if (!C[i].areCartesiansDefined())  C[i].addCartesians();
            double[] coords = C[i].getCartesians();
            PDB = PDB + String.format("ATOM%7d  C   XXX A   1    %8.3f%8.3f%8.3f\n",(i+1),coords[0],coords[1],coords[2]);
         }
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }
      return PDB;
   }
}

