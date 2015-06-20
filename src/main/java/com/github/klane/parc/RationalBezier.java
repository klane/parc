package com.github.klane.parc;

import org.jblas.DoubleMatrix;

public final class RationalBezier extends Bezier {

    private final DoubleMatrix weights;

    public RationalBezier(final DoubleMatrix P) {
        super(P);
        this.weights = DoubleMatrix.ones(P.rows);
    }

    public RationalBezier(final DoubleMatrix P, final DoubleMatrix weights) {
        super(P);
        this.weights = weights.dup();
    }

    @Override
    public DoubleMatrix matrix(final DoubleMatrix t) {
        DoubleMatrix N = super.matrix(t);
        return N.mulRowVector(this.weights).divColumnVector(N.mmul(this.weights));
    }
}
