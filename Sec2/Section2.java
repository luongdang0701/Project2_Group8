public class Section2 {

    // Performs Gaussian Elimination with no pivoting
    public static void noPivotingGE(double[][] A, double[] B) {

        int N = B.length;

        // For every row
        for (int k = 0; k < N; k++) {

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

        //Find Error with x_i=[1]
        double [] e= new double [N];
        for(int i=0; i<N; i++){
            e[i]=Math.abs(solution[i]-1);
        }

        // Output the solution
        printSolution(solution,e,N);
    }

    // Performs Gaussian Elimination using partial pivoting
    public static void pivotingGE(double[][] A, double[] B) {
        int N = B.length;
        for (int k = 0; k < N; k++) { // For every row
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

        //Find Error with x_i=[1]
        double [] e= new double [N];
        for(int i=0; i<N; i++){
            e[i]=Math.abs(solution[i]-1);
        }
        printSolution(solution,e,N); // Output the solution
    }

    // Performs Gaussian Elimination using scaled partial pivoting
    // Algorithm from: https://www.youtube.com/watch?v=4YzIfcSFVCU&ab_channel=ThomasBingham
    public static void scaledPivotingGE(double[][] A, double[] B) {

        int N = B.length;

        double[] S = new double[N];
        for (int i = 0; i < N; i++) {
            S[i] = arrayMax(A[i],false);
        }

        // For every row
        for (int k = 0; k < N; k++) {

            // Scale the rows using highest magnitude elements and find the max row
            int max = k;
            double[] RV = initRV(N);
            for (int i = k; i < N; i++) {
                RV[i] = Math.abs(A[i][k])/S[i];
            }
            max = (int)arrayMax(RV, true);

            // Pivot the rows
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

        //Find Error with x_i=[1]
        double [] e= new double [N];
        for(int i=0; i<N; i++){
            e[i]=Math.abs(solution[i]-1);
        }

        // Output the solution
        printSolution(solution,e,N);

    }

    // Returns the max from an array
    public static double arrayMax(double[] A, boolean index) {
        int maxInd = 0;
        for (int i = 0; i < A.length; i++) {
            if (Math.abs(A[i]) > Math.abs(A[maxInd])) {
                maxInd = i;
            }
        }

        return ((index) ? maxInd : Math.abs(A[maxInd]));
    }

    // Initialize the ratio vector for scaled partial pivoting
    public static double [] initRV (int N) {
        double[] RV = new double [N];
        for (int i = 0; i < N; i++)
            RV[i] = 0;
        return RV;
    }

    // Initialize matrix A
    public static double [][] initA (int N) {
        double[][] A = new double [N][N];
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                A[i][j] = 1/(i+j+1.0);
            }
        }
        return A;
    }

    // Initialize matrix B
    public static double [] initB (double[][] A) {
        int N = A.length;
        double[] B = new double [N];
        for (int i = 0; i < N; i++) {
            B[i] = 0;
            for (int j = 0; j < N; j++) {
                B[i] += A[i][j];
            }
        }
        return B;
    }

    // Prints the solution with error
    public static void printSolution(double[] sol, double []e,int n ){
        for (int i = 0; i < n; i++) {
            if (i < 9) {
                System.out.printf("%s%d = %-20.6f" + "%s%d = %.6f\n", "X", i + 1, sol[i], "E", i + 1, e[i]);
            } else {
                System.out.printf("%s%d = %-19.6f" + "%s%d = %.6f\n", "X", i + 1, sol[i], "E", i + 1, e[i]);
            }
        }
        System.out.println("||x-x~|| = " + arrayMax(e, false));
        System.out.println();
    }

    // Prints a vector
    public static void printVector(double[] sol, String letter) {
        int N = sol.length;
        for (int i = 0; i < N; i++)
            System.out.printf("%s%d = %.6f\n", letter, i+1, sol[i]);
        System.out.println();
    }

    // Print matrix
    public static void printMatrix(double[][] A) {
        int N = A.length;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                System.out.printf("%.3f ", A[i][j]);
            }
            System.out.println();
        }
        System.out.println();
    }



    //#####################################PART 2 ###########################


    public static void LUDecomposition(double[][] A) {
        int n = A.length;
        double[][] upper= new double[n][n];
        double[][] lower= new double[n][n];

        //LU Decomposition --------------------------------------------------------------------
        double[][] temp = A;
        int N = A.length;
        // For every row
        for (int k = 0; k < N; k++) {
            // Identify row below
            for (int i = k + 1; i < N; i++) {
                // Use the multiplier to row reduce
                double multiplier = temp[i][k] / temp[k][k];
                lower[i][k] = multiplier; // <= Lower --------------------
                for (int j = k; j < N; j++) {
                    temp[i][j] -= multiplier * temp[k][j];
                }
            }
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i > j) {
                } else if (i == j) {
                    lower[i][j] = 1;
                } else {
                    lower[i][j] = 0;
                }
            }
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i <= j) {
                    upper[i][j] = temp[i][j];
                } else {
                    upper[i][j] = 0.0;
                }
            }
        }

        // -------------------------------------------------------------------

        System.out.println("-----------------LU Factorial --------------");
        System.out.println("Upper - U");
        printMatrix(upper);
        double[][] inverseU = findInverseUpper(upper);

        System.out.println("Lower - L");
        printMatrix(lower);
        double[][] inverseL = findInverseLower(lower);

        System.out.println("---------Inverse of U and L----------------");
        System.out.println("InverseUpper - U^-1");
        printMatrix(inverseU);
        System.out.println("InverseLower - L^-1");
        printMatrix(inverseL);

        System.out.println("-----------------Inverse of A---------------------");
        System.out.println("InverseA - A^-1");
        double[][] inverseA = multiplyMatrices(inverseU,inverseL);
        printMatrix(inverseA);


        System.out.println("------- Estimated X ( x~ = A^-1 * b) -----------");
        double[] B = initB(A);
        double[] xsenor = multiplyMatrixWithVector(inverseA,B);
        printVector(xsenor,"x~ ");
    }

    static double[] multiplyMatrixWithVector(double A[][], double B[]) {
        int n = A.length;
        double C[] = new double[n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                C[i] += A[i][j] * B[j];
            }
        }
        return C;
    }

    static double[][] multiplyMatrices(double A[][], double B[][]) {
        int n = A.length;
        double C[][] = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                for (int k = 0; k < n; k++) {
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        return C;

    }

    public static double[][] findInverseLower(double[][] lower){
        int n = lower.length;
        double[][] inverseL = new double[n][n];
        for(int u = 0; u < n; u ++) {
            double[] B = new double[n];

            for (int i = 0; i < n; i++) {
                if(i==u){
                    B[i]=1;
                }else{
                    B[i] = 0;}
            }

            double[] solution = new double[n];
            for (int i = 0; i <n; i++) {
                double sum = 0.0;
                for (int j = 0; j < n; j++)
                    sum += lower[i][j] * solution[j];
                solution[i] = (B[i] - sum) / lower[i][i];
            }
            for(int a = 0; a < n; a++){
                inverseL[a][u] = solution[a];
            }
        }

        return inverseL;
    }




    public static double[][] findInverseUpper(double[][] upper){
        int n = upper.length;
        double[][] inverseU = new double[n][n];
        for(int u = 0; u < n; u ++) {
            double[] B = new double[n];

            for (int i = 0; i < n; i++) {
                if(i==u){
                    B[i]=1;
                }else{
                    B[i] = 0;}
            }

            double[] solution = new double[n];
            for (int i = n - 1; i >= 0; i--) {
                double sum = 0.0;
                for (int j = i + 1; j < n; j++)
                    sum += upper[i][j] * solution[j];
                solution[i] = (B[i] - sum) / upper[i][i];
            }
            for(int a = 0; a < n; a++){
                inverseU[a][u] = solution[a];
            }
        }
        return inverseU;
    }

    // Main function
    public static void main(String[] args) {

        int N = 12;
        double[][] A;
        double[] B;

        A = initA(N);
        B = initB(A);

        System.out.println("Matrix A: ");
        printMatrix(A);

        System.out.println("Solution Vector B: ");
        printVector(B,"B");

        A = initA(N);
        B = initB(A);

        System.out.println("Solved using no pivoting:");
        noPivotingGE(A,B);

        A = initA(N);
        B = initB(A);

        System.out.println("Solved using partial pivoting:");
        pivotingGE(A,B);

        A = initA(N);
        B = initB(A);

        System.out.println("Solved using scaled partial pivoting:");
        scaledPivotingGE(A,B);

        // ##########################  Part2 #########################

        A = initA(N);
        B = initB(A);
        LUDecomposition(A);
    }
}
