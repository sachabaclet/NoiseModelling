/**
 * NoiseModelling is a library capable of producing noise maps. It can be freely used either for research and education, as well as by experts in a professional use.
 * <p>
 * NoiseModelling is distributed under GPL 3 license. You can read a copy of this License in the file LICENCE provided with this software.
 * <p>
 * Official webpage : http://noise-planet.org/noisemodelling.html
 * Contact: contact@noise-planet.org
 */

package org.noise_planet.noisemodelling.pathfinder.path;

import org.junit.jupiter.api.Test;
import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineSegment;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.prep.PreparedGeometry;
import org.locationtech.jts.geom.prep.PreparedGeometryFactory;

import java.util.Random;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * The shortcuts of {@link MirrorReceiversCompute} decide like the JTS calls they replace, on random inputs at
 * projected coordinates and on degenerate ones (collinear points, empty segments, distances at the limit)
 */
public class TestMirrorReceiversComputeShortcuts {

    private static Coordinate[] randomPoints(Random random, int count) {
        Coordinate[] points = new Coordinate[count];
        double x = 670000 + random.nextDouble() * 1000;
        double y = 6580000 + random.nextDouble() * 1000;
        double scale = random.nextBoolean() ? 1 : 100;
        for (int i = 0; i < count; i++) {
            points[i] = new Coordinate(x + random.nextGaussian() * scale, y + random.nextGaussian() * scale);
            if (random.nextInt(8) == 0) {
                // on a grid, to get exactly collinear points
                points[i].x = Math.rint(points[i].x);
                points[i].y = Math.rint(points[i].y);
            }
        }
        if (random.nextInt(50) == 0) {
            points[1] = points[0].copy();
        }
        return points;
    }

    @Test
    public void testDistanceShortcuts() {
        Random random = new Random(2);
        int decided = 0;
        for (int i = 0; i < 200000; i++) {
            Coordinate[] p = randomPoints(random, 4);
            LineSegment wall = new LineSegment(p[0], p[1]);
            double distance = wall.distance(p[2]);
            // a limit far from the distance, at the distance, or a few micrometers from it
            double limit = random.nextInt(3) == 0 ? random.nextDouble() * 200 :
                    random.nextBoolean() ? distance : distance + (random.nextDouble() - 0.5) * 4e-6;
            double squaredDistance = MirrorReceiversCompute.squaredDistance(wall, p[2]);
            if (MirrorReceiversCompute.isClearlyBelow(squaredDistance, limit)) {
                assertTrue(distance < limit);
                decided++;
            }
            if (MirrorReceiversCompute.isClearlyAbove(squaredDistance, limit)) {
                assertTrue(distance > limit);
                decided++;
            }
            LineSegment segment = new LineSegment(p[2], p[3]);
            double segmentDistance = wall.distance(segment);
            double segmentLimit = random.nextBoolean() ? random.nextDouble() * 200 :
                    segmentDistance + (random.nextDouble() - 0.5) * 4e-6;
            if (MirrorReceiversCompute.isClearlyCloser(wall, segment, segmentLimit)) {
                assertTrue(segmentDistance < segmentLimit);
                decided++;
            }
        }
        assertTrue(decided > 100000);
    }

    @Test
    public void testFirstPointInCone() {
        Random random = new Random(3);
        GeometryFactory factory = new GeometryFactory();
        int inCone = 0;
        for (int i = 0; i < 20000; i++) {
            Coordinate[] p = randomPoints(random, 3);
            Polygon cone = MirrorReceiversCompute.createWallReflectionVisibilityCone(p[2], new LineSegment(p[0], p[1]),
                    random.nextDouble() * 1000, 100);
            if (cone.isEmpty()) {
                continue;
            }
            PreparedGeometry prepared = PreparedGeometryFactory.prepare(cone);
            Coordinate[] ring = cone.getExteriorRing().getCoordinates();
            Envelope envelope = cone.getEnvelopeInternal();
            for (int j = 0; j < 10; j++) {
                Coordinate[] wall = randomPoints(random, 2);
                if (j % 3 == 0) {
                    // on a vertex or on an edge of the cone
                    int edge = random.nextInt(ring.length - 1);
                    double f = random.nextBoolean() ? 0 : random.nextDouble();
                    wall[0] = new Coordinate(ring[edge].x + f * (ring[edge + 1].x - ring[edge].x),
                            ring[edge].y + f * (ring[edge + 1].y - ring[edge].y));
                } else if (j % 3 == 1) {
                    wall[0] = new Coordinate(envelope.getMinX() + random.nextDouble() * envelope.getWidth(),
                            envelope.getMinY() + random.nextDouble() * envelope.getHeight());
                }
                LineString line = factory.createLineString(wall);
                if (MirrorReceiversCompute.isFirstPointInCone(cone, line)) {
                    assertTrue(prepared.intersects(line));
                    inCone++;
                }
            }
        }
        assertTrue(inCone > 10000);
    }

    @Test
    public void testWallPointTest() {
        Random random = new Random(1);
        for (int i = 0; i < 200000; i++) {
            Coordinate[] p = randomPoints(random, 3);
            if (random.nextInt(4) == 0) {
                // on the line of the wall, or next to it
                double f = random.nextDouble() * 3 - 1;
                p[2] = new Coordinate(p[0].x + f * (p[1].x - p[0].x), p[0].y + f * (p[1].y - p[0].y));
                if (random.nextBoolean()) {
                    p[2].x = Math.nextUp(p[2].x);
                }
            }
            if (random.nextInt(20) == 0) {
                p[1].y = p[0].y;
            }
            assertEquals(Orientation.isCCW(new Coordinate[]{p[0], p[1], p[2], p[0]}),
                    MirrorReceiversCompute.wallPointTest(new LineSegment(p[0], p[1]), p[2]));
        }
    }
}
