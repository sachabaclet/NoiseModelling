/**
 * NoiseModelling is a library capable of producing noise maps. It can be freely used either for research and education, as well as by experts in a professional use.
 * <p>
 * NoiseModelling is distributed under GPL 3 license. You can read a copy of this License in the file LICENCE provided with this software.
 * <p>
 * Official webpage : http://noise-planet.org/noisemodelling.html
 * Contact: contact@noise-planet.org
 */


package org.noise_planet.noisemodelling.pathfinder.delaunay;

import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.geom.*;
import org.locationtech.jts.geom.prep.PreparedPolygon;
import org.locationtech.jts.index.quadtree.Quadtree;
import org.locationtech.jts.index.strtree.STRtree;
import org.locationtech.jts.io.ParseException;
import org.locationtech.jts.io.WKTReader;
import org.locationtech.jts.io.WKTWriter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.tinfour.common.IConstraint;
import org.tinfour.common.LinearConstraint;
import org.tinfour.common.PolygonConstraint;
import org.tinfour.common.SimpleTriangle;
import org.tinfour.common.Vertex;
import org.tinfour.standard.IncrementalTin;

import java.io.*;
import java.util.*;
import java.util.stream.IntStream;
import java.util.stream.StreamSupport;

public class LayerTinfour implements LayerDelaunay {
    private double epsilon = 0.001; // merge of Vertex instances below this distance
    private static final Logger LOGGER = LoggerFactory.getLogger(LayerTinfour.class);
    public String dumpFolder = "";
    List<IConstraint> constraints = new ArrayList<>();
    Quadtree ptsIndex = new Quadtree();
    /** Cell size of {@link #vertexGrid} (m) */
    private static final double VERTEX_GRID_CELL_SIZE = 1.0;
    /** The vertices of {@link #ptsIndex} by cell of a regular grid: a cell holds a Vertex or a list of vertices */
    private final Map<Long, Object> vertexGrid = new HashMap<>();
    private boolean computeNeighbors = false;
    private double maxArea = 0;
    private boolean verbose = true;
    // Output data
    private List<Coordinate> vertices = new ArrayList<>();
    private final List<Triangle> triangles = new ArrayList<Triangle>();
    private final List<Triangle> neighbors = new ArrayList<Triangle>(); // The first neighbor triangle is opposite the first corner of triangle  i
    /**
     * RTree for polygon spatial index
     */
    private final STRtree polygonRtree = new STRtree();
    /**
     * Use PreparedPolygon for faster point-in-polygon tests
     */
    private final Map<Integer, PreparedPolygon> polygonMap = new HashMap<>();

    public boolean isVerbose() {
        return verbose;
    }

    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }

    /**
     * Add a vertex, or reuse the vertex closer than epsilon to the coordinate
     * @param coordinate Vertex position
     * @param index Attribute of a new vertex
     * @return The vertex at this position
     */
    private Vertex addCoordinate(Coordinate coordinate, int index) {
        Vertex found = mayHaveCloseVertex(coordinate) ? findCloseVertex(coordinate) : null;
        if(found == null) {
            found = new Vertex(coordinate.x, coordinate.y, Double.isNaN(coordinate.z) ? 0 : coordinate.z, index);
            ptsIndex.insert(new Envelope(coordinate),  found);
            addToVertexGrid(found);
        }
        return found;
    }

    /**
     * @return The first vertex of the quadtree query closer than epsilon to the coordinate, or null
     */
    private Vertex findCloseVertex(Coordinate coordinate) {
        final Envelope env = new Envelope(coordinate);
        env.expandBy(epsilon);
        List result = ptsIndex.query(env);
        for(Object vertex : result) {
            if(vertex instanceof Vertex) {
                if(((Vertex) vertex).getDistance(coordinate.x, coordinate.y) < epsilon) {
                    return (Vertex) vertex;
                }
            }
        }
        return null;
    }

    /**
     * Cheap test before {@link #findCloseVertex(Coordinate)}, which queries the quadtree: most coordinates have no
     * close vertex, and the vertex grid shows it reading a few cells. When there is a close vertex, the quadtree query
     * still chooses the vertex, so the vertices are the same.
     * @return False if no vertex is closer than epsilon to the coordinate (same distance test as findCloseVertex)
     */
    @SuppressWarnings("unchecked")
    private boolean mayHaveCloseVertex(Coordinate coordinate) {
        if (epsilon > VERTEX_GRID_CELL_SIZE) {
            return true;
        }
        // cells of the envelope widened by 1.5 epsilon, so that the rounding cannot miss a cell
        double reach = 1.5 * epsilon;
        long column0 = getVertexGridIndex(coordinate.x - reach);
        long column1 = getVertexGridIndex(coordinate.x + reach);
        long row0 = getVertexGridIndex(coordinate.y - reach);
        long row1 = getVertexGridIndex(coordinate.y + reach);
        for (long column = column0; column <= column1; column++) {
            for (long row = row0; row <= row1; row++) {
                Object cell = vertexGrid.get(getVertexGridKey(column, row));
                if (cell instanceof Vertex) {
                    if (((Vertex) cell).getDistance(coordinate.x, coordinate.y) < epsilon) {
                        return true;
                    }
                } else if (cell != null) {
                    for (Vertex vertex : (List<Vertex>) cell) {
                        if (vertex.getDistance(coordinate.x, coordinate.y) < epsilon) {
                            return true;
                        }
                    }
                }
            }
        }
        return false;
    }

    @SuppressWarnings("unchecked")
    private void addToVertexGrid(Vertex vertex) {
        long key = getVertexGridKey(getVertexGridIndex(vertex.getX()), getVertexGridIndex(vertex.getY()));
        Object cell = vertexGrid.get(key);
        if (cell == null) {
            vertexGrid.put(key, vertex);
        } else if (cell instanceof Vertex) {
            List<Vertex> cellVertices = new ArrayList<>(2);
            cellVertices.add((Vertex) cell);
            cellVertices.add(vertex);
            vertexGrid.put(key, cellVertices);
        } else {
            ((List<Vertex>) cell).add(vertex);
        }
    }

    private static long getVertexGridIndex(double ordinate) {
        return (long) Math.floor(ordinate / VERTEX_GRID_CELL_SIZE);
    }

    private static long getVertexGridKey(long column, long row) {
        return (column << 32) ^ (row & 0xFFFFFFFFL);
    }


    /**
     * @return When an exception occur, this folder with receiver the input data
     */
    public String getDumpFolder() {
        return dumpFolder;
    }

    /**
     * @param dumpFolder When an exception occur, this folder with receiver the input data
     */
    public void setDumpFolder(String dumpFolder) {
        this.dumpFolder = dumpFolder;
    }

    /**
     * @return Merge vertices closer than specified epsilon
     */
    public double getEpsilon() {
        return epsilon;
    }

    /**
     * @param epsilon Merge vertices closer than specified epsilon
     */
    public void setEpsilon(double epsilon) {
        this.epsilon = epsilon;
    }


    /**
     *
     * @param triangle
     * @return
     */
    private static Coordinate getCentroid(SimpleTriangle triangle) {
        Vertex va = triangle.getVertexA();
        Vertex vb = triangle.getVertexB();
        Vertex vc = triangle.getVertexC();
        double cx = ( va.getX() + vb.getX() + vc.getX() ) / 3d;
        double cy = ( va.getY() + vb.getY() + vc.getY() ) / 3d;
        double cz = ( va.getZ() + vb.getZ() + vc.getZ() ) / 3d;
        return new Coordinate( cx, cy, cz);
    }

    public void dumpData(File destinationCsv) {
        GeometryFactory factory = new GeometryFactory();
        WKTWriter wktWriter = new WKTWriter(3);
        try {
            try (BufferedWriter writer = new BufferedWriter(new FileWriter(destinationCsv))) {
                for(Object vObj : ptsIndex.queryAll()) {
                    if(vObj instanceof Vertex) {
                        final Vertex v = (Vertex)vObj;
                        Point p = factory.createPoint(toCoordinate(v));
                        writer.write(wktWriter.write(p));
                        writer.write("\n");
                    }
                }
                for (IConstraint constraint : constraints) {
                    if (constraint instanceof LinearConstraint) {
                        List<Vertex> vertices = constraint.getVertices();
                        Coordinate[] coordinates = new Coordinate[vertices.size()];
                        for (int i = 0; i < vertices.size(); i++) {
                            final Vertex v = vertices.get(i);
                            coordinates[i] = new Coordinate(v.getX(), v.getY(), v.getZ());
                        }
                        LineString l = factory.createLineString(coordinates);
                        writer.write(wktWriter.write(l));
                        writer.write("\n");
                    } else if (constraint instanceof PolygonConstraint) {
                        List<Vertex> vertices = constraint.getVertices();
                        if(vertices != null && vertices.size() >= 3) {
                            Coordinate[] coordinates = new Coordinate[vertices.size() + 1];
                            for (int i = 0; i < vertices.size() ; i++) {
                                final Vertex v = vertices.get(i);
                                coordinates[i] = new Coordinate(v.getX(), v.getY(), v.getZ());
                            }
                            coordinates[coordinates.length - 1] = coordinates[0];
                            Polygon l = factory.createPolygon(coordinates);
                            writer.write(wktWriter.write(l));
                            writer.write("\n");
                        } else {
                            LOGGER.info("Weird null polygon " + constraint);
                        }
                    }
                }
            }
        }  catch (IOException ioEx) {
            // ignore
        }
    }

    /**
     * Find polygon index by point
     * @param point Point to search
     * @return Polygon index or -1 if not found
     */
    public int findPolygonIndexByPoint(Point point) {
        for(Integer polygonIndex : (List<Integer>)polygonRtree.query(point.getEnvelopeInternal())) {
            PreparedPolygon polygon = polygonMap.get(polygonIndex);
            if(polygon.contains(point)) {
                return polygonIndex;
            }
        }
        return -1;
    }

    /**
     * Read dump file
     * @param dumpPath Path to dump file
     * @throws IOException Read error
     * @throws LayerDelaunayError Error in delaunay processing
     * @throws ParseException Error in WKT parsing
     */
    public void readDump(File dumpPath) throws IOException, LayerDelaunayError, ParseException {
        WKTReader wktReader = new WKTReader();
        try(BufferedReader reader = new BufferedReader(new FileReader(dumpPath))) {
            String line;
            int index = 0;
            while ((line = reader.readLine()) != null) {
                Geometry obj = wktReader.read(line);
                if(obj instanceof Point) {
                    addVertex(obj.getCoordinate());
                } else if(obj instanceof Polygon) {
                    addPolygon((Polygon)obj, index++);
                } else if (obj instanceof LineString) {
                    addLineString((LineString)obj, index++);
                }
            }
        }
    }

    /**
     * Launch delaunay process
     */
    @Override
    public void processDelaunay() throws LayerDelaunayError {
        polygonRtree.build(); // build the rtree as it must not be built when using a multi-thread query
        triangles.clear();
        vertices.clear();

        IncrementalTin tin;
        boolean refine;
        // Triangulate
        tin = new IncrementalTin(epsilon);
        // Add points
        tin.add(ptsIndex.queryAll(), null);
        // Add constraints
        try {
            tin.addConstraints(constraints, false);
        }catch (IllegalStateException ex) {
            // Got error
            // Dump input data
            if(!dumpFolder.isEmpty()) {
                dumpData(new File(dumpFolder, "tinfour_dump.csv"));
            }
            throw new LayerDelaunayError(ex);
        }

        do {
            refine = false;
            // Will triangulate multiple time if refinement is necessary
            if (maxArea > 0) {
                ArrayList<Vertex> newSteinerPoints = new ArrayList<>();
                GeometryFactory gf = new GeometryFactory();
                for (SimpleTriangle triangle : tin.triangles()) {
                    if (triangle.getArea() > maxArea) {
                        Vertex centroid = triangle.getCentroid();
                        // Check if the refinement point does not lie into a polygon
                        if (findPolygonIndexByPoint(gf.createPoint(toCoordinate(centroid))) != -1) {
                            continue;
                        }
                        newSteinerPoints.add(centroid);
                    }
                }
                if (!newSteinerPoints.isEmpty()) {
                    tin.add(newSteinerPoints, null);
                    refine = true;
                    if (verbose) {
                        LOGGER.info("Refining Delaunay with {} points", newSteinerPoints.size());
                    }
                }
            }
        } while (refine);
        List<Vertex> verts = tin.getVertices();
        Map<Vertex, Integer> vertIndex = new HashMap<>();
        this.vertices = new ArrayList<>(verts.size());
        for (int i = 0; i < verts.size(); i++) {
            Vertex v = verts.get(i);
            vertIndex.put(v, i);
            this.vertices.add(toCoordinate(v));
        }
        // Collect triangles but with tin index of vertices
        Map<Integer, Integer> edgeIndexToTriangleIndex = new HashMap<>();
        for(SimpleTriangle t : tin.triangles()) {
            Triangle newTriangle = new Triangle(vertIndex.get(t.getVertexA()), vertIndex.get(t.getVertexB()), vertIndex.get(t.getVertexC()), 0);
            triangles.add(newTriangle);
            if(computeNeighbors) {
                edgeIndexToTriangleIndex.put(t.getEdgeA().getIndex(), triangles.size() - 1);
                edgeIndexToTriangleIndex.put(t.getEdgeB().getIndex(), triangles.size() - 1);
                edgeIndexToTriangleIndex.put(t.getEdgeC().getIndex(), triangles.size() - 1);
            }
        }
        if(computeNeighbors) {
            for(SimpleTriangle t : tin.triangles()) {
                Integer neighA = edgeIndexToTriangleIndex.get(t.getEdgeA().getDual().getIndex());
                Integer neighB = edgeIndexToTriangleIndex.get(t.getEdgeB().getDual().getIndex());
                Integer neighC = edgeIndexToTriangleIndex.get(t.getEdgeC().getDual().getIndex());
                neighbors.add(new Triangle(neighA != null ? neighA : -1,
                        neighB != null ? neighB : -1,
                        neighC != null ? neighC : -1));
            }
        }
        // Update triangle polygon association
        if(!polygonMap.isEmpty()) {
            GeometryFactory gf = new GeometryFactory();
            IntStream.range(0, triangles.size()).parallel().forEach(i -> {
                Triangle triangle = triangles.get(i);
                // Look for associated polygon area
                Coordinate inCenter = org.locationtech.jts.geom.Triangle.inCentre(toCoordinate(verts.get(triangle.getA())),
                        toCoordinate(verts.get(triangle.getB())), toCoordinate(verts.get(triangle.getC())));
                Point inCenterPoint = gf.createPoint(inCenter);
                int polygonIndex = findPolygonIndexByPoint(inCenterPoint);
                if (polygonIndex != -1) {
                    triangle.setAttribute(polygonIndex);
                }
            });
        }
    }


    /**
     * Append a polygon into the triangulation
     *
     * @param newPoly Polygon to append into the mesh, internal rings willb be inserted as holes.
     * @param buildingId Polygon attribute. {@link Triangle#getAttribute()}
     */
    @Override
    public void addPolygon(Polygon newPoly, int buildingId) throws LayerDelaunayError {
        final Coordinate[] coordinates = newPoly.getExteriorRing().getCoordinates();
        // Exterior ring must be CCW
        if(!Orientation.isCCW(coordinates)) {
            CoordinateArrays.reverse(coordinates);
        }
        if (coordinates.length >= 4) {
            List<Vertex> vertexList = new ArrayList<>();
            for(int vId = 0; vId < coordinates.length - 1 ; vId++) {
                vertexList.add(addCoordinate(coordinates[vId], buildingId));
            }
            PolygonConstraint polygonConstraint = new PolygonConstraint(vertexList);
            polygonConstraint.complete();
            if(polygonConstraint.isValid()) {
                constraints.add(polygonConstraint);
                polygonRtree.insert(newPoly.getEnvelopeInternal(), buildingId);
                polygonMap.put(buildingId, new PreparedPolygon(newPoly));
                // Append holes
                final int holeCount = newPoly.getNumInteriorRing();
                for (int holeIndex = 0; holeIndex < holeCount; holeIndex++) {
                    LineString holeLine = newPoly.getInteriorRingN(holeIndex);
                    final Coordinate[] hCoordinates = holeLine.getCoordinates();
                    // Holes must be CW
                    if(Orientation.isCCW(hCoordinates)) {
                        CoordinateArrays.reverse(hCoordinates);
                    }
                    vertexList = new ArrayList<>(hCoordinates.length);
                    for(int vId = 0; vId < hCoordinates.length - 1 ; vId++) {
                        vertexList.add(addCoordinate(hCoordinates[vId], buildingId));
                    }
                    polygonConstraint = new PolygonConstraint(vertexList);
                    polygonConstraint.complete();
                    if(polygonConstraint.isValid()) {
                        constraints.add(polygonConstraint);
                    }
                }
            }
        }
    }

    @Override
    public void setMinAngle(Double minAngle) {
        // TODO Auto-generated method stub

    }

    @Override
    public void hintInit(Envelope bBox, long polygonCount, long verticesCount) throws LayerDelaunayError {
    }

    @Override
    public List<Coordinate> getVertices() throws LayerDelaunayError {
        return this.vertices;
    }

    @Override
    public List<Triangle> getTriangles() throws LayerDelaunayError {
        return triangles;
    }

    private static Coordinate toCoordinate(Vertex v) {
        return new Coordinate(v.getX(), v.getY(), v.getZ());
    }

    /**
     * Append a vertex into the triangulation
     *
     * @param vertexCoordinate coordinate of the new vertex
     */
    @Override
    public void addVertex(Coordinate vertexCoordinate) throws LayerDelaunayError {
        addCoordinate(vertexCoordinate, 0);
    }

    /**
     * Set the maximum area in m²
     *
     * @param maxArea Maximum area in m²
     */
    @Override
    public void setMaxArea(Double maxArea) throws LayerDelaunayError {
        this.maxArea = Math.max(0, maxArea);
    }

    /**
     * Append a LineString into the triangulation
     * @param buildingID Associated ID building that will be available on points
     */
    //add buildingID to edge property and to points property
    public void addLineString(LineString lineToProcess, int buildingID) throws LayerDelaunayError {
        Coordinate[] coordinates = lineToProcess.getCoordinates();
        List<Vertex> vertexList = new ArrayList<>();
        for(Coordinate coordinate : coordinates) {
            vertexList.add(addCoordinate(coordinate, buildingID));
        }
        LinearConstraint linearConstraint = new LinearConstraint(vertexList);
        linearConstraint.complete();
        if(linearConstraint.isValid()) {
            constraints.add(linearConstraint);
        }
    }
    //add buildingID to edge property and to points property

    @Override
    public void reset() {
        // TODO Auto-generated method stub

    }

    /**
     * retrieve results Triangle link
     * @return list of triangles neighbor by their index.
     */
    @Override
    public List<Triangle> getNeighbors() throws LayerDelaunayError {
        if (computeNeighbors) {
            return neighbors;
        } else {
            throw new LayerDelaunayError("You must call setRetrieveNeighbors(True) before process delaunay triangulation");
        }
    }

    /**
     * Enable or Disable the collecting of triangles neighboring data.
     * @param retrieve
     */
    @Override
    public void setRetrieveNeighbors(boolean retrieve) {
        this.computeNeighbors = retrieve;

    }
}
