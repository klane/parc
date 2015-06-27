package com.github.klane.parc;

import static org.junit.Assert.assertEquals;

import org.jblas.DoubleMatrix;
import org.jblas.util.Logger;
import org.junit.Before;
import org.junit.Test;

public final class BSplineTest {

    private static final DoubleMatrix T = DoubleMatrix.linspace(0, 1, 100);
    private static final double TOLERANCE = 1E-10;
    private static final double RADIUS = 1.0;
    private static final int ORDER = 3;
    private DoubleMatrix points;
    private DoubleMatrix weights;
    private DoubleMatrix knots;

    @Before
    public void setUp() {
        Logger.getLogger().setLevel(Logger.WARNING);
        points = new DoubleMatrix(2, 9, 1, 0, 1, 1, 0, 1, -1, 1, -1, 0, -1, -1, 0, -1, 1, -1, 1, 0).transpose();
        weights = new DoubleMatrix(9, 1, 1, Math.sqrt(2)/2, 1, Math.sqrt(2)/2, 1, Math.sqrt(2)/2, 1, Math.sqrt(2)/2, 1);
        knots = new DoubleMatrix(12, 1, 0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1);
    }

    @Test
    public void circle() {
        DoubleMatrix c = new NURBS(points, ORDER, knots, weights).eval(T);

        for (int i=0; i<c.rows; i++) {
            assertEquals(Math.sqrt(Math.pow(c.get(i,0), 2) + Math.pow(c.get(i,1), 2)), RADIUS, TOLERANCE);
        }
    }

    @Test
    public void homogeneousCircle() {
        DoubleMatrix c = new HomogeneousBSpline(points, ORDER, knots, weights).eval(T);

        for (int i=0; i<c.rows; i++) {
            assertEquals(Math.sqrt(Math.pow(c.get(i,0), 2) + Math.pow(c.get(i,1), 2)), RADIUS, TOLERANCE);
        }
    }

    @Test
    public void homogeneousEqualsNURBS() {
        DoubleMatrix n = new NURBS(points, ORDER, knots, weights).eval(T);
        DoubleMatrix h = new HomogeneousBSpline(points, ORDER, knots, weights).eval(T);

        for (int i=0; i<n.rows; i++) {
            assertEquals(n.get(i,0), h.get(i,0), TOLERANCE);
            assertEquals(n.get(i,1), h.get(i,1), TOLERANCE);
        }
    }

    @Test
    public void unitWeightedEqualsBSpline() {
        DoubleMatrix b = new BSpline(points, ORDER, knots).eval(T);
        DoubleMatrix n = new NURBS(points, ORDER, knots).eval(T);
        DoubleMatrix h = new HomogeneousBSpline(points, ORDER, knots).eval(T);

        for (int i=0; i<b.rows; i++) {
            assertEquals(b.get(i,0), n.get(i,0), TOLERANCE);
            assertEquals(b.get(i,1), n.get(i,1), TOLERANCE);
            assertEquals(b.get(i,0), h.get(i,0), TOLERANCE);
            assertEquals(b.get(i,1), h.get(i,1), TOLERANCE);
        }
    }
}
