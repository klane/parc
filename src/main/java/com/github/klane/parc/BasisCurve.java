package com.github.klane.parc;

import org.jblas.DoubleMatrix;

public abstract class BasisCurve implements ParametricCurve {

    protected final int n;
    protected final int d;
    protected DoubleMatrix P;

    public BasisCurve(final int n, final int d) {
        this.n = n;
        this.d = d;
    }

    public BasisCurve(final DoubleMatrix P) {
        this(P.rows-1, P.columns);
        this.P = P.dup();
    }

    @Override
    public DoubleMatrix eval(final DoubleMatrix t) {
        return this.matrix(t).mmul(this.P);
    }

    public abstract DoubleMatrix matrix(final DoubleMatrix t);
}
