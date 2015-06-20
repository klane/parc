package com.github.klane.parc;

import org.jblas.DoubleMatrix;

public class BSpline extends BasisCurve {

    private final int order;
    private final DoubleMatrix knots;

    public BSpline(final DoubleMatrix P, final int order, final DoubleMatrix knots) {
        super(P);
        this.order = order;
        this.knots = knots.dup();
    }

    @Override
    public DoubleMatrix matrix(final DoubleMatrix t) {
        DoubleMatrix N = DoubleMatrix.zeros(t.length, super.n+1);
        DoubleMatrix matrix;
        double denominator;

        for (int i=0; i<=super.n; i++) {
            if (this.knots.get(i+1) < this.knots.get(this.knots.length-1)) {
                N.put(t.ge(this.knots.get(i)).and(t.lt(this.knots.get(i+1))), i, 1);
            } else {
                N.put(t.ge(this.knots.get(i)), i, 1);
            }
        }

        for (int k=2; k<=this.order; k++) {
            for (int i=0; i<=super.n; i++) {
                matrix = DoubleMatrix.zeros(t.length, 1);
                denominator = this.knots.get(i+k-1) - this.knots.get(i);

                if (denominator != 0) {
                    matrix.addi(t.sub(this.knots.get(i)).mulColumnVector(N.getColumn(i)).div(denominator));
                }

                denominator = this.knots.get(i+k) - this.knots.get(i+1);

                if (denominator != 0) {
                    matrix.addi(t.neg().add(this.knots.get(i+k)).mulColumnVector(N.getColumn(i+1)).div(denominator));
                }

                N.putColumn(i, matrix);
            }
        }

        return N;
    }
}
