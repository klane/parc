package com.github.klane.parc;

import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;

public class Bezier extends BasisCurve {

    public Bezier(final DoubleMatrix P) {
        super(P);
    }

    public static long binomial(final int n, final int k) {
        if (n < 0) {
            throw new IllegalArgumentException("n < 0");
        }

        if (k < 0) {
            throw new IllegalArgumentException("k < 0");
        }

        if (n < k) {
            throw new IllegalArgumentException("n < k");
        }

        long b = 1;

        for (int i=1, m=n; i<=Math.min(k, n-k); i++, m--) {
            b = b*m/i;
        }

        return b;
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
