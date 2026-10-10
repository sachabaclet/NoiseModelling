/**
 * NoiseModelling is a library capable of producing noise maps. It can be freely used either for research and education, as well as by experts in a professional use.
 * <p>
 * NoiseModelling is distributed under GPL 3 license. You can read a copy of this License in the file LICENCE provided with this software.
 * <p>
 * Official webpage : http://noise-planet.org/noisemodelling.html
 * Contact: contact@noise-planet.org
 */

package org.noise_planet.noisemodelling.pathfinder.profilebuilder;

import org.junit.jupiter.api.Test;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.LineSegment;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Compare the single pass upper hull of {@link CutProfile#getConvexHullIndices(List, boolean)} with the JTS path
 */
public class CutProfileConvexHullTest {

    /**
     * @param types S source, R receiver, T topography, W top of a wall, w bottom of a wall, G ground effect
     * @param xy distance from the source and height of each point
     */
    private static CutProfile profile(String types, double... xy) {
        CutProfile profile = new CutProfile();
        for (int i = 0; i < types.length(); i++) {
            Coordinate c = new Coordinate(xy[2 * i], 0, xy[2 * i + 1]);
            switch (types.charAt(i)) {
                case 'S': profile.cutPoints.add(new CutPointSource(c)); break;
                case 'R': profile.cutPoints.add(new CutPointReceiver(c)); break;
                case 'T': profile.cutPoints.add(new CutPointTopography(c)); break;
                case 'G': profile.cutPoints.add(new CutPointGroundEffect(0, c, 0.5)); break;
                default:
                    CutPointWall wall = new CutPointWall(0, c, new LineSegment(new Coordinate(0, 0, c.z),
                            new Coordinate(1, 0, c.z)), Collections.emptyList());
                    wall.setZGround(types.charAt(i) == 'W' ? c.z - 10 : c.z);
                    profile.cutPoints.add(wall);
            }
        }
        return profile;
    }

    private static List<Coordinate> pts2D(CutProfile profile) {
        List<Coordinate> pts = new ArrayList<>();
        for (CutPoint cutPoint : profile.cutPoints) {
            pts.add(new Coordinate(cutPoint.getCoordinate().x, cutPoint.getCoordinate().z));
        }
        return pts;
    }

    private static void assertUpperHull(List<Integer> expected, CutProfile profile, boolean ignoreWall) {
        List<Coordinate> pts = pts2D(profile);
        assertEquals(expected, profile.getUpperHullIndices(pts, ignoreWall));
        assertEquals(expected, profile.getConvexHullIndicesJts(pts, ignoreWall));
    }

    /** The upper hull falls back to JTS, so the result (or the exception) is the one of JTS */
    private static void assertFallback(CutProfile profile) {
        List<Coordinate> pts = pts2D(profile);
        assertNull(profile.getUpperHullIndices(pts, false));
        List<Integer> jts;
        try {
            jts = profile.getConvexHullIndicesJts(pts, false);
        } catch (IllegalArgumentException ex) {
            assertThrows(IllegalArgumentException.class, () -> profile.getConvexHullIndices(pts, false));
            return;
        }
        assertEquals(jts, profile.getConvexHullIndices(pts, false));
    }

    @Test
    public void testSameHullAsJts() {
        // drone 120 m above the ground, receiver 4 m above the ground
        assertUpperHull(List.of(0, 5), profile("STGTTR",
                0, 130, 0, 10, 40, 10, 80, 11, 150, 12, 200, 14), false);
        // building with a 50 m roof between 150 m and 170 m: diffraction on its last wall
        CutProfile building = profile("STTwWWwTR",
                0, 130, 0, 10, 60, 11, 150, 12, 150, 50, 170, 50, 170, 12, 190, 12, 200, 14);
        assertUpperHull(List.of(0, 5, 8), building, false);
        assertUpperHull(List.of(0, 8), building, true);
        // hill and thin wall
        assertUpperHull(List.of(0, 3, 10), profile("STTTGTwWwTR",
                0, 20, 0, 10, 30, 18, 60, 25, 70, 22, 90, 15, 120, 14, 120, 19, 120, 14, 150, 13, 160, 15), false);
    }

    @Test
    public void testFallbackToJts() {
        // single point
        assertFallback(profile("S", 0, 12));
        // not sorted by distance
        assertFallback(profile("STR", 2, 6, 0, 1, 3, 11));
        // ordinate too large
        assertFallback(profile("STR", 0, 11, 2, 3, 4, Double.MAX_VALUE));
        // negative zero, JTS keeps the second point
        assertFallback(profile("STTR", 0, -0.0, -0.0, -0.0, 2, 2, 3, 7));
        // a point above the source at the same distance
        assertFallback(profile("SWWR", 0, 6, 0, 14, 0, 12, 5, 25));
        // receiver between two points at the same distance, JTS removes it
        assertFallback(profile("STTR", 0, 11, 1, 0, 1, 6, 1, 3));
        // all points aligned, JTS gives a line
        assertFallback(profile("STR", 0, 8, 2, 6.5, 4, 5));
        // more than 50 points with a flat inner octagon in JTS: duplicated points down to the receiver
        StringBuilder types = new StringBuilder("S");
        double[] xy = new double[124];
        xy[1] = 200;
        for (int i = 1; i <= 30; i++) {
            types.append("TT");
            xy[4 * i - 2] = xy[4 * i] = i;
            xy[4 * i - 1] = xy[4 * i + 1] = 200 - 2 * i + i % 3;
        }
        types.append("R");
        xy[122] = 31;
        xy[123] = 100;
        assertFallback(profile(types.toString(), xy));
    }
}
