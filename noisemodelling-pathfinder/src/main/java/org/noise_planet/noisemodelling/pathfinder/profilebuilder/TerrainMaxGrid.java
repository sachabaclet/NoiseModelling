/**
 * NoiseModelling is a library capable of producing noise maps of cities and of evaluating noise levels
 * at specific points. It is distributed under the GPL 3 license.
 * <p>
 * Official webpage : http://noise-planet.org/noisemodelling.html
 * Contact: contact@noise-planet.org
 */

package org.noise_planet.noisemodelling.pathfinder.profilebuilder;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.noise_planet.noisemodelling.pathfinder.delaunay.Triangle;

import java.util.Arrays;
import java.util.List;

/**
 * Upper bound of the height of the DEM over each cell of a regular grid.
 * <p>
 * The bound of a cell is the maximum, over the triangles whose envelope (plus one millimeter) intersects the cell, of
 * the lower of two bounds of the triangle: its highest vertex, and the maximum of its plane over the cell (plus one
 * centimeter), which is at a corner of the cell. So any point of the DEM surface closer than one millimeter to a cell,
 * in particular any point where a propagation line crosses a triangle side, is below the bound of the cell.
 */
final class TerrainMaxGrid {
    static final double CELL_SIZE = 10;
    private final double[] maxZ;
    final double originX;
    final double originY;
    final int columns;
    final int rows;

    private TerrainMaxGrid(double originX, double originY, int columns, int rows) {
        this.originX = originX;
        this.originY = originY;
        this.columns = columns;
        this.rows = rows;
        maxZ = new double[columns * rows];
        Arrays.fill(maxZ, Double.NEGATIVE_INFINITY);
    }

    /**
     * @return The grid of the triangles, null if there is no triangle
     */
    static TerrainMaxGrid build(List<Triangle> triangles, List<Coordinate> vertices) {
        Envelope bounds = new Envelope();
        for (Triangle triangle : triangles) {
            bounds.expandToInclude(vertices.get(triangle.getA()));
            bounds.expandToInclude(vertices.get(triangle.getB()));
            bounds.expandToInclude(vertices.get(triangle.getC()));
        }
        if (bounds.isNull()) {
            return null;
        }
        bounds.expandBy(1);
        TerrainMaxGrid grid = new TerrainMaxGrid(bounds.getMinX(), bounds.getMinY(),
                (int) (bounds.getWidth() / CELL_SIZE) + 1, (int) (bounds.getHeight() / CELL_SIZE) + 1);
        for (Triangle triangle : triangles) {
            grid.addTriangle(vertices.get(triangle.getA()), vertices.get(triangle.getB()), vertices.get(triangle.getC()));
        }
        return grid;
    }

    private void addTriangle(Coordinate a, Coordinate b, Coordinate c) {
        double vertexMax = Math.max(a.z, Math.max(b.z, c.z));
        if (Double.isNaN(a.z) || Double.isNaN(b.z) || Double.isNaN(c.z)) {
            vertexMax = Double.POSITIVE_INFINITY;
        }
        // plane z = a.z + slopeX * (x - a.x) + slopeY * (y - a.y), from the normal of the triangle
        double ux = b.x - a.x, uy = b.y - a.y, uz = b.z - a.z;
        double vx = c.x - a.x, vy = c.y - a.y, vz = c.z - a.z;
        double normalZ = ux * vy - uy * vx;
        boolean hasPlane = Math.abs(normalZ) > 1e-9;
        double slopeX = hasPlane ? -(uy * vz - uz * vy) / normalZ : 0;
        double slopeY = hasPlane ? -(uz * vx - ux * vz) / normalZ : 0;
        Envelope envelope = new Envelope(a, b);
        envelope.expandToInclude(c);
        envelope.expandBy(ProfileBuilder.MILLIMETER);
        for (int j = row(envelope.getMinY()); j <= row(envelope.getMaxY()); j++) {
            for (int i = column(envelope.getMinX()); i <= column(envelope.getMaxX()); i++) {
                double bound = vertexMax;
                if (hasPlane) {
                    // corners of the cell plus one centimeter, relative to a
                    double x0 = originX + i * CELL_SIZE - 0.01 - a.x;
                    double x1 = x0 + CELL_SIZE + 0.02;
                    double y0 = originY + j * CELL_SIZE - 0.01 - a.y;
                    double y1 = y0 + CELL_SIZE + 0.02;
                    double planeMax = a.z + Math.max(slopeX * x0, slopeX * x1) + Math.max(slopeY * y0, slopeY * y1);
                    // margin for the rounding of the plane computation
                    planeMax += 1e-6 * (1 + Math.abs(planeMax));
                    bound = Math.min(bound, planeMax);
                }
                int cell = j * columns + i;
                maxZ[cell] = Math.max(maxZ[cell], bound);
            }
        }
    }

    int column(double x) {
        return Math.max(0, Math.min(columns - 1, (int) Math.floor((x - originX) / CELL_SIZE)));
    }

    int row(double y) {
        return Math.max(0, Math.min(rows - 1, (int) Math.floor((y - originY) / CELL_SIZE)));
    }

    /**
     * @return The bound of the cell, +infinity outside the grid
     */
    double getCellMax(int column, int row) {
        if (column < 0 || row < 0 || column >= columns || row >= rows) {
            return Double.POSITIVE_INFINITY;
        }
        return maxZ[row * columns + column];
    }

    /**
     * @return An upper bound of the DEM height at this position, +infinity outside the grid
     */
    double getMaxAt(double x, double y) {
        if (!isInside(x, y)) {
            return Double.POSITIVE_INFINITY;
        }
        return maxZ[row(y) * columns + column(x)];
    }

    /**
     * @return True if the position is inside the grid, farther than one millimeter from its border
     */
    boolean isInside(double x, double y) {
        double gx = x - originX;
        double gy = y - originY;
        return gx > ProfileBuilder.MILLIMETER && gy > ProfileBuilder.MILLIMETER
                && gx < columns * CELL_SIZE - ProfileBuilder.MILLIMETER && gy < rows * CELL_SIZE - ProfileBuilder.MILLIMETER;
    }
}
