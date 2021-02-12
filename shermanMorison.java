public class shermanMorison {

    public static double[] thomasAlg(double[][] mainMatrix, double[] bFreeValues){
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

    public static double roundTo(double input){
        double scale = Math.pow(10,14);

        return  Math.ceil(input * scale) / scale;
    }

    public static double scalar(double[] array,double[] array2){
        double result = 0d;
        int i = 0;

        for (double n : array) {
            result += n*array2[i];
            i++;
        }

        return roundTo(result);
    }

    public static double[] function(double[][] mainMatrix,double[] bFreeValues){
        //firstly we need implement our vectors u[] and v[]
        double[] u = new double[mainMatrix.length];
        double[] v = new double[mainMatrix.length];

        int N = mainMatrix.length - 1;

        double gammaU = 1; // our coefficient can be set freely
        u[0] = gammaU; // first element in u[] is gamma coefficient
        u[u.length - 1] = mainMatrix[N][0]; // last element lowest-left element of A matrix

        v[0] = 1;
        v[v.length - 1] = mainMatrix[0][N] / gammaU; // last element in v[] is upper-right element of A

        mainMatrix[0][0] -= gammaU;
        mainMatrix[N][N] -= mainMatrix[N][0]*mainMatrix[0][N] / gammaU;

        mainMatrix[0][N] = 0;
        mainMatrix[N][0] = 0;


        //after that we have modified 3diagonal matrix

        double[] y = thomasAlg(mainMatrix,bFreeValues);
        double[] z = thomasAlg(mainMatrix,u);

        for (int i = 0; i < y.length; i++) {
            y[i] = roundTo(y[i]);
            z[i] = roundTo(z[i]);
        }

        double[] result = new double[mainMatrix.length]; //


        for (int i = 0; i < result.length; i++) {
            result[i] =  y[i] - ( scalar(v,y) ) / ( (1 + scalar(v,z)) )   * z[i] ;
        }
        
        return result;
    }

    public static void main(String[] args) {

        double[][] A = new double[][]{
                {4 , 1 , 0 , 0 , 0 , 0 , 1},
                {1 , 4 , 1 , 0 , 0 , 0 , 0},
                {0 , 1 , 4 , 1 , 0 , 0 , 0},
                {0 , 0 , 1 , 4 , 1 , 0 , 0},
                {0 , 0 , 0 , 1 , 4 , 1 , 0},
                {0 , 0 , 0 , 0 , 1 , 4 , 1},
                {1 , 0 , 0 , 0 , 0 , 1 , 4}
        };

        double[] bValues = new double[]{1 , 2 , 3 , 4 , 5 , 6 , 7};


        long nanos1 = System.nanoTime();
        double[] x = function(A,bValues);
        long nanos2 = System.nanoTime();

        for (int i = 0; i < x.length; i++) {
            System.out.println("x"+ (i+1) + " -> " + x[i]);
        }

        System.out.println();
        System.out.println("time -> " + (nanos2 - nanos1) + " nanoseconds");

    }
}
