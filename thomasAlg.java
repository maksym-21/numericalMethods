public class thomasAlg {

    public static double[] function(double[][] mainMatrix, double[] bFreeValues){
        //A * x = B

        double[] a = new double[mainMatrix.length]; //our alpha coefficients
        double[] b = new double[mainMatrix.length]; //our beta coefficients
        double[] results = new double[mainMatrix.length];

        int N1 = mainMatrix.length - 1;

        //forward substitution
        double y = mainMatrix[0][0]; //here we are getting our coefficients
        a[0] = (-1) * mainMatrix[0][1] / y;
        b[0] = bFreeValues[0] / y;

        for (int i = 1; i < N1 ; i++) {
            y = mainMatrix[i][i] + mainMatrix[i][i-1] * a[i-1];
            a[i] = - mainMatrix[i][i+1] / y;
            b[i] = (bFreeValues[i] - mainMatrix[i][i-1] * b[i-1] ) / y;
        }

        y = mainMatrix[N1][N1] + mainMatrix[N1][N1 - 1] * a[N1 - 1];
        b[N1] = ( bFreeValues[N1] - mainMatrix[N1][N1 - 1] * b[N1 - 1] ) / y;
        //end of forward substitution


        //start of backward substitution
        results[N1] = b[N1];

        for (int i = N1 - 1; i >= 0 ; i--) {
            results[i] = a[i] * results[i+1] + b[i];
        }
        //end of backward substitution

        return results;
    }


    public static void main(String[] args) {
        double[][] aMatrix = new double[][] {
                {4,1,0,0,0,0,0},
                {1,4,1,0,0,0,0},
                {0,1,4,1,0,0,0},
                {0,0,1,4,1,0,0},
                {0,0,0,1,4,1,0},
                {0,0,0,0,1,4,1},
                {0,0,0,0,0,1,4}
        };

        double[] bValues = new double[]{ 1, 2, 3, 4, 5, 6, 7};

        long ms1 = System.nanoTime();
        double[] myArray = function(aMatrix,bValues);
        long ms2 = System.nanoTime();

        for (int i = 0; i < myArray.length; i++) {
            System.out.println("x" + i + " -> " + myArray[i]);
        }

        System.out.println("\n" + (ms2-ms1) + " nanoseconds");
    }
}
