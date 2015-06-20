package com.github.klane.parc;

import org.jblas.DoubleMatrix;

public abstract class BasisCurve implements ParametricCurve {

    protected final int n;
    private DoubleMatrix P;

    public BasisCurve(final int n) {
        this.n = n;
    }

    public BasisCurve(final DoubleMatrix P) {
        this(P.rows-1);
        this.P = P.dup();
    }

    @Override
    public DoubleMatrix eval(final DoubleMatrix t) {
        return this.matrix(t).mmul(this.P);
    }

    public abstract DoubleMatrix matrix(final DoubleMatrix t);
}
