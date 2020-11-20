import java.util.Scanner;

public class Section1Problem3 {

    // Performs Gaussian Elimination using matrix A and solution matrix B
    public static void pivotingGE(double[][] A, double[] B) {

        int N = B.length;

        // For every row
        for (int k = 0; k < N; k++) {

            // Find pivot row and swap
            int max = k;
            for (int i = k + 1; i < N; i++)
                if (Math.abs(A[i][k]) > Math.abs(A[max][k]))
                    max = i;
            double[] temp = A[k];
            A[k] = A[max];
            A[max] = temp;
            double t = B[k];
            B[k] = B[max];
            B[max] = t;

            // Identify row below
            for (int i = k + 1; i < N; i++) {

                // Use the multiplier to row reduce
                double multiplier = A[i][k] / A[k][k];
                B[i] -= multiplier * B[k];
                for (int j = k; j < N; j++)
                    A[i][j] -= multiplier * A[k][j];
            }
        }

        // Solve using back substitution
        double[] solution = new double[N];
        for (int i = N - 1; i >= 0; i--) {
            double sum = 0.0;
            for (int j = i + 1; j < N; j++)
                sum += A[i][j] * solution[j];
            solution[i] = (B[i] - sum) / A[i][i];
        }

        // Output the solution
        printSolution(solution);
    }

    // Prints the solution
    public static void printSolution(double[] sol) {
        int N = sol.length;
        System.out.println("\nSolutions: ");
        for (int i = 0; i < N; i++)
            System.out.printf("A%d = %.3f\n", i+1, sol[i]);
        System.out.println();
    }

    // Main function
    public static void main(String[] args) {

        Scanner scan = new Scanner(System.in);

        System.out.println("\nEnter number of equations:");
        int N = scan.nextInt();

        double[] B = new double[N];
        double[][] A = new double[N][N];

        System.out.println("\nEnter "+ N +" equations coefficients:");
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                A[i][j] = scan.nextDouble();

        System.out.println("\nEnter "+ N +" solutions:");
        for (int i = 0; i < N; i++)
            B[i] = scan.nextDouble();

        /**For Section1, enter input as follow (You can copy-paste into as input):

         Enter number of equations:
         6

         Enter 6  equations coefficients:
         1 1 0 1 1 0
         -8 -7 1 -6 -9 1
         22 16 -3 12 29 -9
         -26 -16 3 -12 -39 29
         21 15 -3 11 18 -39
         -18 -9 2 -6 0 18

         Enter 6 solutions:
         0
         0
         0
         1
         1
         1
        */
        pivotingGE(A,B);
    }
}
