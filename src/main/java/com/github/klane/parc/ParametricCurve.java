package com.github.klane.parc;

import org.jblas.DoubleMatrix;

public interface ParametricCurve {

    DoubleMatrix eval(final DoubleMatrix t);
}
