package com.github.klane.parc;

import org.jblas.DoubleMatrix;
import org.jblas.ranges.RangeUtils;

public final class HomogeneousBSpline extends BSpline {

    public HomogeneousBSpline(final DoubleMatrix P, final int order, final DoubleMatrix knots) {
        super(DoubleMatrix.concatHorizontally(P, DoubleMatrix.ones(P.rows)), order, knots);
    }

    public HomogeneousBSpline(final DoubleMatrix P, final int order, final DoubleMatrix knots, final DoubleMatrix weights) {
        super(DoubleMatrix.concatHorizontally(P.mulColumnVector(weights), weights), order, knots);
    }

    @Override
    public DoubleMatrix eval(final DoubleMatrix t) {
        DoubleMatrix x = super.eval(t);
        return x.getColumns(RangeUtils.interval(0, x.columns-1)).divColumnVector(x.getColumn(x.columns-1));
    }
}
