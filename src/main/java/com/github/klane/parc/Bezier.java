package com.github.klane.parc;

import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;

public class Bezier extends BasisCurve {

    public Bezier(final DoubleMatrix P) {
        super(P);
    }

    public static int binomial(final int n, final int k) {
        double b = 1;

        for (double i=1; i<=k; i++) {
            b *= (n-i+1) / i;
        }

        return (int) b;
    }

    @Override
    public DoubleMatrix matrix(final DoubleMatrix t) {
        DoubleMatrix N = DoubleMatrix.zeros(t.length, super.n+1);

        for (int i=0; i<=super.n; i++) {
            N.putColumn(i, MatrixFunctions.pow(t, i)
                    .mul(MatrixFunctions.pow(t.neg().add(1), super.n-i))
                    .mul(Bezier.binomial(super.n, i)));
        }

        return N;
    }
}
