public class gaussElimination {

    public static double[] gaussElimanation(double[][] matrix,double[] bFreeValues){
        double[][] matrixx =new double[matrix.length][matrix.length];

        for (int i = 0; i < matrixx.length; i++) {
            System.arraycopy(matrix[i], 0, matrixx[i], 0, matrixx.length);
        }

        double[] b = new double[bFreeValues.length];

        System.arraycopy(bFreeValues, 0, b, 0, bFreeValues.length);

        double[] x = new double[matrix.length];

        double var,var2;
        int LENGTH = matrix.length - 1;

        for (int k = 0; k < matrixx.length; k++) { //forward substitution
            for (int j = k + 1; j < matrixx.length; j++) {
                var = matrixx[j][k] / matrixx[k][k];
                for (int i = k; i < matrixx.length; i++) {
                    matrixx[j][i] = matrixx[j][i] - var * matrixx[k][i];
                }
                b[j] = b[j] - var * b[k];
            }
        }


        x[LENGTH] = b[LENGTH] / matrixx[LENGTH][LENGTH]; // x(n) value

        for (int j = matrixx.length - 2; j >= 0 ; j--) {
            var = 0;
            for (int k = j + 1; k < matrixx.length; k++) {
                var2 = matrixx[j][k] * x [k];
                var = var + var2;
            }
            x[j] = (b[j] - var) / matrixx[j][j];
        }

        return x;
    }

    public static double[] subtraction(double[] a,double[] b){
        double[] c = new double[a.length];

        for (int i = 0; i < a.length; i++) {
            c[i] = a[i] - b[i];
        }

        return c;
    }

    public static double getNorma(double[] vector){
        double var = 0;
        for (double a : vector) {
            var += Math.pow(a,2);
        }
        return Math.sqrt(var);
    }

    public static void main(String[] args) {
        final double[][] A = new double[][]{
            { -116.66654, 583.33346 , -333.33308 , 100.00012 , 100.00012},
            { 583.33346 , -116.66654 , -333.33308, 100.00012 , 100.00012},
            { -333.33308, -333.33308 , 133.33383 , 200.00025 , 200.00025},
            { 100.00012 , 100.00012  , 200.00025 , 50.000125 , -649.99988},
            { 100.00012 , 100.00012 , 200.00025 , -649.99988 , 50.000125}};

        double[] b1 = new double[]{-0.33388066 , 1.08033290 ,-0.98559856,1.31947922,-0.09473435};
        double[] b2 = new double[]{-0.33388066,1.08033290,-0.98559855,1.32655028,-0.10180541};
        double[] b3 = new double[]{0.72677951,0.72677951,-0.27849178,0.96592583,0.96592583};
        double[] b4 = new double[]{0.73031505,0.73031505,-0.27142071,0.96946136,0.96946136};


        double[] z1 = gaussElimanation(A,b1);
        double[] z2 = gaussElimanation(A,b2);
        double[] z3 = gaussElimanation(A,b3);
        double[] z4 = gaussElimanation(A,b4);

        long ns1 = System.nanoTime();
        System.out.println(getNorma(subtraction(b1,b2))); // ||b1 - b2||
        long ns2 = System.nanoTime();
        System.out.println(getNorma(subtraction(b3,b4))); // ||b3 - b4||
        long ns3 = System.nanoTime();
        // ||z1 - z2|| / ||b1 - b2||
        System.out.println(getNorma(subtraction(z1,z2)) / getNorma(subtraction(b1,b2)));
        long ns4 = System.nanoTime();
        // ||z3 - z4|| / ||b3 - b4||
        System.out.println(getNorma(subtraction(z3,z4)) / getNorma(subtraction(b3,b4)));
        long ns5 = System.nanoTime();

        System.out.println();
        System.out.println((ns2-ns1) + "\n" + (ns3-ns2) + "\n" + (ns4-ns3) + "\n" + (ns5-ns4) + "\n");
    }
}
