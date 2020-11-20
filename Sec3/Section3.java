public class Section3 {
    public static void main(String[] args) throws Exception {

        // First, define the tolerance
        double epsilon = 5E-6;

        // Define maximum number of iterations
        int maxIters = 500;

        // Define the initial values/seed (x0^0)
        double[] x0 = {0,0,0,0,0,0};

        // Define the matrix of coefficients
        double[][] M = {
                {4, -1, 0, -2, 0, 0}, // coefficients for first equation
                {-1, 4, -1, 0, -2, 0}, // second equation
                {0, -1, 4, 0, 0, -2}, // third
                {-1, 0, 0, 4, -1, 0}, // fourth
                {0, -1, 0, -1, 4, -1}, // ...
                {0, 0, -1, 0, -1, 4} // last equation
        };

        // Define vector b (equalities)
        double[] b = {-1, 0, 1, -2, 1, 2};

        // Call method
        Jacobi(M, b, x0, epsilon, maxIters);
        System.out.println("");
        GaussSeidel(M, b, x0, epsilon, maxIters);
    }

    public static void Jacobi(double[][] M, double[] b, double[] x0, double epsilon, int maxIters)
    {
        System.out.println(" *** JACOBI METHOD *** ");
        System.out.println("Solving the following system:");
        System.out.println("Number of variables: " + b.length);
        System.out.println("Toleance: " + epsilon);
        System.out.println("Max. number of iterations: " + maxIters);
        System.out.println("The system is:\n");

        // Define initial error
        double err = 1E+10; // a rally big number

        // Define number of equations/variables
        int N = b.length;

        // Print equations so user can see them
        for(int i = 0; i < N; i++)
        {
            for(int j = 0; j < N; j++)
            {
                System.out.print(M[i][j]+"*x" + (j+1) + " ");
            }
            System.out.print("= " + b[i] + "\n");
        }
        System.out.println("");


        System.out.println("Starting iterations... ");
        // Start iterations
        double[] xold = x0; // varray to store old values (from previous iteration)
        double[] xnew = xold.clone(); // array to store the result of the new iteration
        int iterCounter = 0; // variable to count the iterations

        double xi, erri; // helper variables to be used in the method

        while(err > epsilon)
        {
            iterCounter++;

            // Jacobi Method
            for(int i = 0; i < N ;i++)
            {
                xi = b[i];
                for(int j = 0; j < N; j++)
                {
                    if(i != j)
                        xi -= M[i][j] * xold[j];
                }
                xi = xi/M[i][i];
                xnew[i] = xi;
            }

            // Calculate errors and pick maximum

            /*
                For this method the error is calculated as the difference between the old solution and
                the new one.
            */
            double maxErr = -1;
            for(int i = 0; i < N; i++)
            {
                erri = Math.abs((xnew[i]-xold[i]));
                if(erri > maxErr)
                    maxErr = erri;
            }
            xold = xnew.clone(); // For the next iteration, the old_values are the new_values in this one
            err = maxErr;
            System.out.println("Iteration " + iterCounter + ", Err: " + err);
            if(iterCounter == maxIters) // reached maximum number of iterations
            {
                System.out.println("Maximum number of iterations reached. Last error: " + err);
                break;
            }
        }

        // Display results only if the system converged
        if(iterCounter < maxIters && err <= epsilon)
        {
            System.out.println("");
            System.out.println("Convergence achieved in " + iterCounter + " iterations. Min error was: " + err);
            System.out.println("The solution is:");
            for(int i = 0; i < N; i++)
            {
                System.out.println("x[" + (i+1) + "] = " + xnew[i]);
            }
        }
    }

    public static void GaussSeidel(double[][] M, double[] b, double[] x0, double epsilon, int maxIters)
    {
        System.out.println(" *** GAUSS-SEIDEL METHOD *** ");
        System.out.println("Solving the following system:");
        System.out.println("Number of variables: " + b.length);
        System.out.println("Toleance: " + epsilon);
        System.out.println("Max. number of iterations: " + maxIters);
        System.out.println("The system is:\n");

        // Define initial error
        double err = 1E+10; // a rally big number

        // Define number of equations/variables
        int N = b.length;

        // Display equations
        for(int i = 0; i < N; i++)
        {
            for(int j = 0; j < N; j++)
            {
                System.out.print(M[i][j]+"*x" + (j+1) + " ");
            }
            System.out.print("= " + b[i] + "\n");
        }
        System.out.println("");
        System.out.println("Starting iterations... ");

        // Start iterations
        double[] x = x0; // vector of solutions is setted to the initial values
        int iterCounter = 0; // variable to count interations

        double xi, erri; // helper variables
        double sigma; // helper variable for gauss-Seidel method

        // Start iterations
        while(err > epsilon)
        {
            iterCounter++;

            // Gauss-Seidel method
            for(int i = 0; i < N ;i++)
            {
                sigma = 0;
                for(int j = 0; j < N; j++)
                {
                    if(i != j)
                        sigma += M[i][j]*x[j];
                }
                xi = (b[i]-sigma)/M[i][i];
                x[i] = xi;
            }

            // Calculate errors and pick maximum
            /*
                For this method the error is calculated as the value of the functions for the current solution.
                The method converges when the value of the function f(xi)-bi = 0 is less than epsilon
            */
            double maxErr = -1;
            for(int i = 0; i < N; i++)
            {
                erri = 0;
                for(int j = 0; j < N; j++)
                {
                    erri += M[i][j]*x[j];
                }
                erri -= b[i];
                erri = Math.abs(erri);
                if(erri > maxErr)
                    maxErr = erri;
            }
            //xold = xnew.clone(); // the results of this iterations are the old for the next iteration
            err = maxErr;
            System.out.println("Iteration " + iterCounter + ", Err: " + err);
            if(iterCounter == maxIters) // reached max number of iters
            {
                System.out.println("Maximum number of iterations reached. Last error: " + err);
                break;
            }
        }

        if(iterCounter < maxIters && err <= epsilon) // display results only if system converged
        {
            System.out.println("");
            System.out.println("Convergence achieved in " + iterCounter + " iterations. Min error was: " + err);
            System.out.println("The solution is:");
            for(int i = 0; i < N; i++)
            {
                System.out.println("x[" + (i+1) + "] = " + x[i]);
            }
        }
    }
}
