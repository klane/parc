package com.github.klane.parc;

import static org.junit.Assert.assertEquals;

import org.jblas.DoubleMatrix;
import org.junit.Before;
import org.junit.Test;

public final class BezierTest {

    private static final DoubleMatrix T = DoubleMatrix.linspace(0, 1, 100);
    private static final double TOLERANCE = 1E-10;
    private DoubleMatrix points;
    private DoubleMatrix weights;

    @Before
    public void setUp() {
        points = new DoubleMatrix(2, 9, 1, 0, 1, 1, 0, 1, -1, 1, -1, 0, -1, -1, 0, -1, 1, -1, 1, 0).transpose();
        weights = new DoubleMatrix(9, 1, 1, Math.sqrt(2)/2, 1, Math.sqrt(2)/2, 1, Math.sqrt(2)/2, 1, Math.sqrt(2)/2, 1);
    }

    @Test
    public void binomial() {
        assertEquals(Bezier.binomial(8, 3), 56);
        assertEquals(Bezier.binomial(50, 5), 2118760);
    }

    @Test
    public void homogeneousEqualsRational() {
        DoubleMatrix r = new RationalBezier(points, weights).eval(T);
        DoubleMatrix h = new HomogeneousBezier(points, weights).eval(T);

        for (int i=0; i<r.rows; i++) {
            assertEquals(r.get(i,0), h.get(i,0), TOLERANCE);
            assertEquals(r.get(i,1), h.get(i,1), TOLERANCE);
        }
    }

    @Test
    public void unitWeightedEqualsBezier() {
        DoubleMatrix b = new Bezier(points).eval(T);
        DoubleMatrix r = new RationalBezier(points).eval(T);
        DoubleMatrix h = new HomogeneousBezier(points).eval(T);

        for (int i=0; i<b.rows; i++) {
            assertEquals(b.get(i,0), r.get(i,0), TOLERANCE);
            assertEquals(b.get(i,1), r.get(i,1), TOLERANCE);
            assertEquals(b.get(i,0), h.get(i,0), TOLERANCE);
            assertEquals(b.get(i,1), h.get(i,1), TOLERANCE);
        }
    }
}
