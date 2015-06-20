package com.github.klane.parc;

import org.jblas.DoubleMatrix;

public final class NURBS extends BSpline {

    private final DoubleMatrix weights;

    public NURBS(final DoubleMatrix P, final int order, final DoubleMatrix knots) {
        super(P, order, knots);
        this.weights = DoubleMatrix.ones(P.rows);
    }

    public NURBS(final DoubleMatrix P, final int order, final DoubleMatrix knots, final DoubleMatrix weights) {
        super(P, order, knots);
        this.weights = weights.dup();
    }

    @Override
    public DoubleMatrix matrix(final DoubleMatrix t) {
        DoubleMatrix N = super.matrix(t);
        return N.mulRowVector(this.weights).divColumnVector(N.mmul(this.weights));
    }
}
