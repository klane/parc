package com.github.klane.parc;

import org.jblas.DoubleMatrix;
import org.jblas.ranges.RangeUtils;

public final class HomogeneousBezier extends Bezier {

    public HomogeneousBezier(final DoubleMatrix P) {
        super(DoubleMatrix.concatHorizontally(P, DoubleMatrix.ones(P.rows)));
    }

    public HomogeneousBezier(final DoubleMatrix P, final DoubleMatrix weights) {
        super(DoubleMatrix.concatHorizontally(P.mulColumnVector(weights), weights));
    }

    @Override
    public DoubleMatrix eval(final DoubleMatrix t) {
        DoubleMatrix x = super.eval(t);
        return x.getColumns(RangeUtils.interval(0, x.columns-1)).divColumnVector(x.getColumn(x.columns-1));
    }
}
