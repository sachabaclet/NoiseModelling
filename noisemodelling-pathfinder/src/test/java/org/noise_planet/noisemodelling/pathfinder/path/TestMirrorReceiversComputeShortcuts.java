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
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.LineSegment;

import java.util.Random;

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
}
