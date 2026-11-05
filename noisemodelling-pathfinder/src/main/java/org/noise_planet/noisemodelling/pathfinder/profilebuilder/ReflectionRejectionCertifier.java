/**
 * NoiseModelling is a library capable of producing noise maps of cities and of evaluating noise levels
 * at specific points. It is distributed under the GPL 3 license.
 * <p>
 * Official webpage : http://noise-planet.org/noisemodelling.html
 * Contact: contact@noise-planet.org
 */

package org.noise_planet.noisemodelling.pathfinder.profilebuilder;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.LineSegment;
import org.locationtech.jts.triangulate.quadedge.Vertex;
import org.noise_planet.noisemodelling.pathfinder.utils.geometry.JTSUtility;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Decide, without computing the topography, that the homogeneous path of a reflection profile is rejected.
 * <p>
 * CnossosPathBuilder.computeCnossosPathFromCutProfile moves each reflection point up to the convex hull of the profile
 * (when the hull has diffraction points) and rejects the path if a moved point is above its wall. This is the case of
 * most reflection paths of road traffic maps, and their topography profile is the main part of their cost.
 * <p>
 * The same test is done on the profile without topography, with a margin of 1e-6 m. If the path is rejected, it is
 * also rejected with the topography when no topography point can change the hull: every point where the propagation
 * line crosses a DEM triangle side must be below the hull by 1e-4 m. This is proven with the bounds of
 * {@link TerrainMaxGrid} along each leg, and near the ends of a leg (where the hull is low, at the source for
 * example) with the exact topography points of {@link ProfileBuilder#fetchTopographicProfile}. When a test is not
 * conclusive the path is not rejected here, it is computed as before. So the results are the same.
 * <p>
 * The tests on the hull depend on the code of CnossosPathBuilder (lift of the reflection points, hull inputs):
 * both must be kept the same.
 */
final class ReflectionRejectionCertifier {
    /** Margin of the reflection height test */
    private static final double HEIGHT_MARGIN = 1e-6;
    /** Same as CnossosPathBuilder.EPSILON */
    private static final double CNOSSOS_EPSILON = 1e-7;
    /** Margin between the topography points and the hull */
    private static final double HULL_MARGIN = 1e-4;
    /** Length of the leg ends where the topography points are computed exactly */
    private static final double EXACT_LEG_END_LENGTH = 20;
    /** Triangle of a leg end not looked up yet */
    private static final int UNKNOWN_TRIANGLE = -2;

    private final ProfileBuilder profileBuilder;
    private final TerrainMaxGrid terrainMaxGrid;

    ReflectionRejectionCertifier(ProfileBuilder profileBuilder, TerrainMaxGrid terrainMaxGrid) {
        this.profileBuilder = profileBuilder;
        this.terrainMaxGrid = terrainMaxGrid;
    }

    /**
     * The reflection profile without topography, assembled as PathFinder.computeReflexion assembles the profile
     */
    private static final class ProfileWithoutTopography {
        final List<Coordinate> legPoints;
        final List<CutPoint> points = new ArrayList<>();
        /** Index in points of the first point of each leg, and of the receiver */
        final int[] legStart;
        /** Triangle of each leg end, as found by fetchTopographicProfile, looked up when needed */
        final int[] triangles;
        List<Coordinate> pts2D;
        List<Integer> hull;
        double[] hullX;
        double[] hullY;

        ProfileWithoutTopography(List<Coordinate> legPoints) {
            this.legPoints = legPoints;
            legStart = new int[legPoints.size()];
            triangles = new int[legPoints.size()];
            Arrays.fill(triangles, UNKNOWN_TRIANGLE);
        }

        /** Profile distance of a leg end */
        double getLegStartDistance(int leg) {
            return pts2D.get(legStart[leg]).x;
        }

        /** Hull height at a profile distance (the hull vertices are in increasing distance) */
        double getHullHeight(double x) {
            return interpolateHull(hullX, hullY, x);
        }
    }

    /**
     * @param legPoints Source, reflection points, receiver
     * @param obstructionCuts Building, wall and ground cut points of each leg
     * @param reflectionWalls Wall of each reflection point
     * @return True if the homogeneous path of the reflection profile is rejected by CnossosPathBuilder, false if it
     * is not or if it is unknown
     */
    boolean isCertainlyRejected(List<Coordinate> legPoints, ProfileBuilder.ObstructionCut[] obstructionCuts,
                                List<LineSegment> reflectionWalls) {
        ProfileWithoutTopography profile = assembleProfile(legPoints, obstructionCuts, reflectionWalls);
        if (profile == null || !isRejectedByReflectionHeight(profile)) {
            return false;
        }
        for (int leg = 0; leg < obstructionCuts.length; leg++) {
            if (!isLegTopographyBelowHull(profile, leg)) {
                return false;
            }
        }
        return true;
    }

    /**
     * @return The profile with the obstruction cut points, its 2D coordinates and its upper hull, null if the hull is
     * not computed by the upper hull method (degenerate cases) or a triangle is missing
     */
    private ProfileWithoutTopography assembleProfile(List<Coordinate> legPoints,
                                                     ProfileBuilder.ObstructionCut[] obstructionCuts,
                                                     List<LineSegment> reflectionWalls) {
        ProfileWithoutTopography profile = new ProfileWithoutTopography(legPoints);
        int legCount = obstructionCuts.length;
        CutPointSource source = null;
        CutPointReceiver receiver = null;
        for (int leg = 0; leg < legCount; leg++) {
            CutProfile legProfile = new CutProfile(new CutPointSource(legPoints.get(leg)),
                    new CutPointReceiver(legPoints.get(leg + 1)));
            legProfile.insertCutPoint(true, obstructionCuts[leg].cutPoints);
            List<CutPoint> legCutPoints = legProfile.cutPoints;
            profile.legStart[leg] = profile.points.size();
            if (leg == 0) {
                source = (CutPointSource) legCutPoints.get(0);
                profile.points.add(source);
            } else {
                CutPointReflection reflection = new CutPointReflection(legCutPoints.get(0),
                        reflectionWalls.get(leg - 1), Collections.emptyList());
                if (!setReflectionZGround(profile, leg, reflection)) {
                    return null;
                }
                profile.points.add(reflection);
            }
            profile.points.addAll(legCutPoints.subList(1, legCutPoints.size() - 1));
        }
        receiver = new CutPointReceiver(legPoints.get(legCount));
        profile.legStart[legCount] = profile.points.size();
        profile.points.add(receiver);
        CutProfile cutProfile = new CutProfile(source, receiver);
        cutProfile.insertCutPoint(false,
                profile.points.subList(1, profile.points.size() - 1).toArray(CutPoint[]::new));
        profile.pts2D = cutProfile.computePts2D();
        profile.hull = cutProfile.getUpperHullIndices(profile.pts2D, false);
        if (profile.hull == null) {
            return null;
        }
        profile.hullX = new double[profile.hull.size()];
        profile.hullY = new double[profile.hull.size()];
        for (int i = 0; i < profile.hull.size(); i++) {
            profile.hullX[i] = profile.pts2D.get(profile.hull.get(i)).x;
            profile.hullY[i] = profile.pts2D.get(profile.hull.get(i)).y;
        }
        return profile;
    }

    /**
     * The height test only checks that the z of a reflection point is not its z ground: when the point is above the
     * grid bound the z ground stays NaN (the test gives the same answer), else it is computed as
     * fetchTopographicProfile computes it.
     * @return False if the triangle of the point is not found
     */
    private boolean setReflectionZGround(ProfileWithoutTopography profile, int leg, CutPointReflection reflection) {
        Coordinate position = profile.legPoints.get(leg);
        if (position.z > terrainMaxGrid.getMaxAt(position.x, position.y)) {
            return true;
        }
        if (!lookupTriangle(profile, leg)) {
            return false;
        }
        reflection.zGround = getZGround(profile, leg);
        return true;
    }

    /**
     * Height test of CnossosPathBuilder.computeCnossosPathFromCutProfile: with diffraction points, each reflection
     * point between two hull vertices is moved to the hull segment, and the path is rejected if it is then above its
     * wall.
     * @return True if a moved reflection point is above its wall by more than the margin
     */
    private static boolean isRejectedByReflectionHeight(ProfileWithoutTopography profile) {
        List<Integer> hull = profile.hull;
        if (hull.size() <= 2) {
            return false;
        }
        for (int i = 1; i < hull.size(); i++) {
            int i0 = hull.get(i - 1);
            int i1 = hull.get(i);
            LineSegment hullSegment = new LineSegment(profile.pts2D.get(i0), profile.pts2D.get(i1));
            for (int pointIndex = i0 + 1; pointIndex < i1; pointIndex++) {
                CutPoint point = profile.points.get(pointIndex);
                if (point instanceof CutPointReflection && Double.compare(point.getCoordinate().z, point.getzGround()) != 0) {
                    LineSegment wall = ((CutPointReflection) point).wall;
                    double movedHeight = hullSegment.closestPoint(profile.pts2D.get(pointIndex)).y;
                    double wallHeight = Vertex.interpolateZ(point.coordinate, wall.p0, wall.p1);
                    if (movedHeight > wallHeight + CNOSSOS_EPSILON + HEIGHT_MARGIN) {
                        return true;
                    }
                }
            }
        }
        return false;
    }

    /**
     * Test the topography of a leg against the hull, with the grid bounds. When it fails, the start of the leg, then
     * its end, are computed exactly and the test is done again. The start of the first leg is the source, near the
     * ground, so it is always computed exactly.
     * @return True if no topography point of the leg can change the hull
     */
    private boolean isLegTopographyBelowHull(ProfileWithoutTopography profile, int leg) {
        boolean exactStart = leg == 0;
        boolean exactEnd = false;
        while (!isLegTopographyBelowHull(profile, leg, exactStart, exactEnd)) {
            if (!exactStart) {
                exactStart = true;
            } else if (!exactEnd) {
                exactEnd = true;
            } else {
                return false;
            }
        }
        return true;
    }

    /**
     * The topography points of a leg are the points where it crosses a triangle side. They must be below the hull by
     * the margin. The z ground of the wall points is interpolated from them and from the z ground of the leg ends: the
     * wall points that are hull inputs because their z is not their z ground must stay above it.
     * @param exactStart Compute the topography points of the start of the leg exactly
     * @param exactEnd Compute the topography points of the end of the leg exactly
     */
    private boolean isLegTopographyBelowHull(ProfileWithoutTopography profile, int leg, boolean exactStart,
                                             boolean exactEnd) {
        Coordinate start = profile.legPoints.get(leg);
        Coordinate end = profile.legPoints.get(leg + 1);
        double length = JTSUtility.dist2D(start, end);
        double startDistance = profile.getLegStartDistance(leg);
        if ((exactStart && !lookupTriangle(profile, leg)) || (exactEnd && !lookupTriangle(profile, leg + 1))) {
            return false;
        }
        double maxZGround = Math.max(getZGroundBound(profile, leg), getZGroundBound(profile, leg + 1));
        // part of the leg left to the grid bounds
        double gridStart = 0;
        double gridEnd = length;
        if (exactStart) {
            double[] exact = testExactLegEnd(profile, start, end, length, startDistance, 1);
            if (exact == null) {
                return false;
            }
            gridStart = exact[0];
            maxZGround = Math.max(maxZGround, exact[1]);
        }
        if (exactEnd) {
            double[] exact = testExactLegEnd(profile, end, start, length, startDistance + length, -1);
            if (exact == null) {
                return false;
            }
            gridEnd = length - exact[0];
            maxZGround = Math.max(maxZGround, exact[1]);
        }
        double gridMax = testGridCells(profile, start, end, length, startDistance, gridStart, gridEnd);
        if (Double.isNaN(gridMax)) {
            return false;
        }
        maxZGround = Math.max(maxZGround, gridMax);
        return areWallsAboveGround(profile, leg, maxZGround);
    }

    /**
     * Compute the topography points of the first {@link #EXACT_LEG_END_LENGTH} meters of the leg from p0 to p1: they
     * are the crossings of the triangle sides found by fetchTopographicProfile on the whole leg (in the reverse
     * direction for the end of a leg). Test them against the hull.
     * @param x0 Profile distance of p0
     * @param direction 1 if the profile distance increases from p0 to p1, -1 otherwise
     * @return The distance from p0 after which the grid bounds are used (+infinity if the whole leg is computed), and
     * the highest topography point, null if a point is not below the hull by the margin
     */
    private double[] testExactLegEnd(ProfileWithoutTopography profile, Coordinate p0, Coordinate p1, double length,
                                     double x0, int direction) {
        double exactLength = Math.min(length, EXACT_LEG_END_LENGTH);
        Coordinate exactEnd = length > 0 ? new LineSegment(p0, p1).pointAlong(exactLength / length) : new Coordinate(p1);
        exactEnd.z = p1.z;
        List<Coordinate> topography = new ArrayList<>();
        profileBuilder.fetchTopographicProfile(topography, p0, exactEnd, false);
        double maxZ = Double.NEGATIVE_INFINITY;
        // the first and last points are the ends, not crossings of triangle sides
        for (int i = 1; i < topography.size() - 1; i++) {
            Coordinate point = topography.get(i);
            double x = x0 + direction * JTSUtility.dist2D(p0, point);
            if (!(point.z + HULL_MARGIN < profile.getHullHeight(x))) {
                return null;
            }
            maxZ = Math.max(maxZ, point.z);
        }
        // the crossings of the whole leg within one centimeter of the end of this part are left to the grid bounds
        double gridFrom = exactLength >= length ? Double.POSITIVE_INFINITY : exactLength - 0.01;
        return new double[]{gridFrom, maxZ};
    }

    /**
     * Walk the grid cells crossed by the part [gridStart, gridEnd] of a leg (distances from its start). The bound of a
     * cell must be below the hull by the margin over the part of the leg in the cell: as the hull is concave, it is
     * above its lowest end.
     * @return The highest bound of the walked cells, NaN if a cell is not below the hull
     */
    private double testGridCells(ProfileWithoutTopography profile, Coordinate start, Coordinate end, double length,
                                 double startDistance, double gridStart, double gridEnd) {
        gridStart = Math.max(0, gridStart);
        gridEnd = Math.min(length, gridEnd);
        if (!(gridStart <= gridEnd && length > 0)) {
            return Double.NEGATIVE_INFINITY;
        }
        double ux = (end.x - start.x) / length;
        double uy = (end.y - start.y) / length;
        double startX = start.x + ux * gridStart;
        double startY = start.y + uy * gridStart;
        if (!terrainMaxGrid.isInside(startX, startY)
                || !terrainMaxGrid.isInside(start.x + ux * gridEnd, start.y + uy * gridEnd)) {
            return Double.NaN;
        }
        final double cellSize = TerrainMaxGrid.CELL_SIZE;
        int column = terrainMaxGrid.column(startX);
        int row = terrainMaxGrid.row(startY);
        // leg distance of the next vertical and horizontal grid lines, and between two of them
        double nextColumnDistance = Double.POSITIVE_INFINITY;
        double nextRowDistance = Double.POSITIVE_INFINITY;
        double columnStep = Double.POSITIVE_INFINITY;
        double rowStep = Double.POSITIVE_INFINITY;
        if (ux != 0) {
            double lineX = terrainMaxGrid.originX + (ux > 0 ? column + 1 : column) * cellSize;
            nextColumnDistance = gridStart + (lineX - startX) / ux;
            columnStep = cellSize / Math.abs(ux);
        }
        if (uy != 0) {
            double lineY = terrainMaxGrid.originY + (uy > 0 ? row + 1 : row) * cellSize;
            nextRowDistance = gridStart + (lineY - startY) / uy;
            rowStep = cellSize / Math.abs(uy);
        }
        HullCursor hull = new HullCursor(profile.hullX, profile.hullY);
        double hullAtCellStart = hull.getHeight(startDistance + gridStart);
        double maxBound = Double.NEGATIVE_INFINITY;
        while (true) {
            double cellEnd = Math.min(Math.min(nextColumnDistance, nextRowDistance), gridEnd);
            double hullAtCellEnd = hull.getHeight(startDistance + cellEnd);
            double bound = terrainMaxGrid.getCellMax(column, row);
            if (!(bound + HULL_MARGIN < Math.min(hullAtCellStart, hullAtCellEnd))) {
                return Double.NaN;
            }
            maxBound = Math.max(maxBound, bound);
            if (cellEnd >= gridEnd) {
                return maxBound;
            }
            if (nextColumnDistance < nextRowDistance) {
                column += ux > 0 ? 1 : -1;
                nextColumnDistance += columnStep;
            } else {
                row += uy > 0 ? 1 : -1;
                nextRowDistance += rowStep;
            }
            hullAtCellStart = hullAtCellEnd;
        }
    }

    /**
     * A wall point is a hull input when its z is not its z ground. A wall point below the hull by the margin is not a
     * hull vertex in both cases: only the other wall points must be above the highest possible z ground.
     */
    private static boolean areWallsAboveGround(ProfileWithoutTopography profile, int leg, double maxZGround) {
        for (int pointIndex = profile.legStart[leg] + 1; pointIndex < profile.legStart[leg + 1]; pointIndex++) {
            CutPoint point = profile.points.get(pointIndex);
            if (point instanceof CutPointWall && !(point.getCoordinate().z > maxZGround + HULL_MARGIN)) {
                Coordinate point2D = profile.pts2D.get(pointIndex);
                if (!(point2D.y + HULL_MARGIN < profile.getHullHeight(point2D.x))) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * @return False if the leg end is outside the DEM triangles
     */
    private boolean lookupTriangle(ProfileWithoutTopography profile, int index) {
        if (profile.triangles[index] == UNKNOWN_TRIANGLE) {
            profile.triangles[index] = profileBuilder.getTriangleIdByCoordinate(profile.legPoints.get(index));
        }
        return profile.triangles[index] >= 0;
    }

    /** @return The z ground of a leg end, computed as fetchTopographicProfile computes it */
    private double getZGround(ProfileWithoutTopography profile, int index) {
        Coordinate[] vertices = profileBuilder.getTriangleVertices(profile.triangles[index]);
        return Vertex.interpolateZ(profile.legPoints.get(index), vertices[0], vertices[1], vertices[2]);
    }

    /** @return The z ground of a leg end if its triangle is known, else its grid bound */
    private double getZGroundBound(ProfileWithoutTopography profile, int index) {
        if (profile.triangles[index] >= 0) {
            return getZGround(profile, index);
        }
        Coordinate position = profile.legPoints.get(index);
        return terrainMaxGrid.getMaxAt(position.x, position.y);
    }

    /** Hull height at a distance, the hull vertices being in increasing distance */
    private static double interpolateHull(double[] hullX, double[] hullY, double x) {
        int last = hullX.length - 1;
        if (x <= hullX[0]) {
            return hullY[0];
        }
        if (x >= hullX[last]) {
            return hullY[last];
        }
        int i = Arrays.binarySearch(hullX, x);
        if (i >= 0) {
            return hullY[i];
        }
        i = -i - 1;
        return hullY[i - 1] + (x - hullX[i - 1]) / (hullX[i] - hullX[i - 1]) * (hullY[i] - hullY[i - 1]);
    }

    /** Hull height at increasing distances, without search */
    private static final class HullCursor {
        private final double[] hullX;
        private final double[] hullY;
        private int segment = 1;

        HullCursor(double[] hullX, double[] hullY) {
            this.hullX = hullX;
            this.hullY = hullY;
        }

        double getHeight(double x) {
            int last = hullX.length - 1;
            if (x <= hullX[0]) {
                return hullY[0];
            }
            if (x >= hullX[last]) {
                return hullY[last];
            }
            while (hullX[segment] < x) {
                segment++;
            }
            return hullY[segment - 1] + (x - hullX[segment - 1]) / (hullX[segment] - hullX[segment - 1])
                    * (hullY[segment] - hullY[segment - 1]);
        }
    }
}
