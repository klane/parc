package com.github.klane.parc;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import org.jblas.DoubleMatrix;
import org.jblas.util.Logger;
import org.junit.Before;
import org.junit.Test;

public final class BezierTest {

    private static final DoubleMatrix T = DoubleMatrix.linspace(0, 1, 100);
    private static final double TOLERANCE = 1E-10;
    private DoubleMatrix points;
    private DoubleMatrix weights;

    @Before
    public void setUp() {
        Logger.getLogger().setLevel(Logger.WARNING);
        points = new DoubleMatrix(2, 9, 1, 0, 1, 1, 0, 1, -1, 1, -1, 0, -1, -1, 0, -1, 1, -1, 1, 0).transpose();
        weights = new DoubleMatrix(9, 1, 1, Math.sqrt(2)/2, 1, Math.sqrt(2)/2, 1, Math.sqrt(2)/2, 1, Math.sqrt(2)/2, 1);
    }

    @Test
    public void binomial() {
        assertEquals(Bezier.binomial(0, 0), 1);
        assertEquals(Bezier.binomial(1, 0), 1);
        assertEquals(Bezier.binomial(1, 1), 1);
        assertEquals(Bezier.binomial(2, 2), 1);
        assertEquals(Bezier.binomial(8, 3), 56);
        assertEquals(Bezier.binomial(8, 5), 56);
        assertEquals(Bezier.binomial(11, 5), 462);
        assertEquals(Bezier.binomial(11, 6), 462);
        assertEquals(Bezier.binomial(50, 5), 2118760);
        assertEquals(Bezier.binomial(50, 45), 2118760);

        try {
            Bezier.binomial(-1, 0);
            fail();
        } catch (IllegalArgumentException e) {}

        try {
            Bezier.binomial(0, -1);
            fail();
        } catch (IllegalArgumentException e) {}

        try {
            Bezier.binomial(0, 1);
            fail();
        } catch (IllegalArgumentException e) {}
    }

    @Test
    public void elevate() {
        Bezier b = new Bezier(points);
        Bezier e = b.elevate(0);
        DoubleMatrix x = b.eval(T);
        DoubleMatrix xi;

        assertEquals(b.n, e.n);

        for (int i=0; i<b.n; i++) {
            assertEquals(b.P.get(i,0), e.P.get(i,0), 0);
            assertEquals(b.P.get(i,1), e.P.get(i,1), 0);
        }

        for (int i=1; i<=20; i++) {
            e = b.elevate(i);
            xi = e.eval(T);
            assertEquals(b.n+i, e.n);

            for (int j=0; j<x.rows; j++) {
                assertEquals(x.get(j,0), xi.get(j,0), TOLERANCE);
                assertEquals(x.get(j,1), xi.get(j,1), TOLERANCE);
            }
        }
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
