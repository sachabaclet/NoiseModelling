/**
 * NoiseModelling is a library capable of producing noise maps. It can be freely used either for research and education, as well as by experts in a professional use.
 * <p>
 * NoiseModelling is distributed under GPL 3 license. You can read a copy of this License in the file LICENCE provided with this software.
 * <p>
 * Official webpage : http://noise-planet.org/noisemodelling.html
 * Contact: contact@noise-planet.org
 */

package org.noise_planet.noisemodelling.pathfinder.path;

import org.locationtech.jts.algorithm.Intersection;
import org.locationtech.jts.algorithm.LineIntersector;
import org.locationtech.jts.algorithm.RobustLineIntersector;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineSegment;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.prep.PreparedGeometry;
import org.locationtech.jts.geom.prep.PreparedGeometryFactory;
import org.locationtech.jts.index.ItemVisitor;
import org.locationtech.jts.index.strtree.AbstractNode;
import org.locationtech.jts.index.strtree.Boundable;
import org.locationtech.jts.index.strtree.ItemBoundable;
import org.locationtech.jts.index.strtree.STRtree;
import org.locationtech.jts.io.WKTWriter;
import org.locationtech.jts.math.Vector2D;
import org.locationtech.jts.operation.buffer.BufferParameters;
import org.locationtech.jts.triangulate.quadedge.Vertex;
import org.noise_planet.noisemodelling.pathfinder.profilebuilder.ProfileBuilder;
import org.noise_planet.noisemodelling.pathfinder.profilebuilder.Wall;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;

public class MirrorReceiversCompute {
    private static final double DEFAULT_CIRCLE_POINT_ANGLE = Math.PI / 24;
    STRtree mirrorReceiverTree;
    public static final int DEFAULT_MIRROR_RECEIVER_CAPACITY = 50000;
    private int mirrorReceiverCapacity = DEFAULT_MIRROR_RECEIVER_CAPACITY;
    private final Coordinate receiverCoordinate;
    private final List<Wall> buildWalls;
    private final double maximumDistanceFromWall;
    private final double maximumPropagationDistance;
    int numberOfImageReceivers = 0;
    /** With the sources given to the constructor, the visitor that holds the images of each source */
    private Map<Coordinate, ReceiverImageVisitor> visitorBySource = null;

    public static Polygon createWallReflectionVisibilityCone(Coordinate receiverImage, LineSegment wall,
                                                             double maximumPropagationDistance,
                                                             double maximumDistanceFromWall) {
        double distanceMin = wall.distance(receiverImage);

        GeometryFactory factory = new GeometryFactory();
        if(distanceMin > maximumPropagationDistance) {
            return factory.createPolygon();
        }
        ArrayList<Coordinate> circleSegmentPoints = new ArrayList<>();

        Vector2D rP0 = new Vector2D(receiverImage, wall.p0).normalize();
        Vector2D rP1 = new Vector2D(receiverImage, wall.p1).normalize();
        double angleSign = rP0.angleTo(rP1) >= 0 ? 1 : -1;
        int numberOfStep = Math.max(1, (int)(Math.abs(rP0.angleTo(rP1)) / DEFAULT_CIRCLE_POINT_ANGLE));
        Coordinate lastWallIntersectionPoint = new Coordinate();
        for(int angleStep = 0 ; angleStep <= numberOfStep; angleStep++) {
            Vector2D newPointTranslationVector = rP0.rotate(DEFAULT_CIRCLE_POINT_ANGLE * angleSign * angleStep);
            if(angleStep == numberOfStep) {
                newPointTranslationVector = rP1;
            } else if(angleStep == 0) {
                newPointTranslationVector = rP0;
            }
            Coordinate newPoint = newPointTranslationVector.translate(receiverImage);
            Coordinate wallIntersectionPoint = Intersection.intersection(wall.p0, wall.p1, receiverImage, newPoint);
            if(wallIntersectionPoint != null) {
                double wallIntersectionPointDistance = wallIntersectionPoint.distance(receiverImage);
                if (wallIntersectionPointDistance < maximumPropagationDistance) {
                    double vectorLength = maximumPropagationDistance;
                    newPoint = newPointTranslationVector.multiply(vectorLength).translate(receiverImage);
                    if (circleSegmentPoints.isEmpty()) {
                        circleSegmentPoints.add(wallIntersectionPoint);
                    }
                    lastWallIntersectionPoint = wallIntersectionPoint;
                    circleSegmentPoints.add(newPoint);
                }
            }
        }
        if(!circleSegmentPoints.isEmpty()) {
            circleSegmentPoints.add(lastWallIntersectionPoint);
            circleSegmentPoints.add(circleSegmentPoints.get(0));
            Coordinate[] conePolygon = circleSegmentPoints.toArray(new Coordinate[0]);
            return factory.createPolygon(conePolygon);
        } else {
            return factory.createPolygon();
        }
    }
    /**
     * Generate all image receivers from the provided list of walls
     * @param buildWalls
     * @param receiverCoordinates
     * @param reflectionOrder
     */
    public MirrorReceiversCompute(List<Wall> buildWalls, Coordinate receiverCoordinates,
                                  int reflectionOrder, double maximumPropagationDistance,
                                  double maximumDistanceFromWall) {
        this(buildWalls, receiverCoordinates, reflectionOrder, maximumPropagationDistance, maximumDistanceFromWall,
                null, 0);
    }

    /**
     * Generate the image receivers from the provided list of walls
     * @param sources Positions of the sources that will be given to {@link #findCloseMirrorReceivers(Coordinate)},
     *                or null if unknown. When known, only the images that give a reflection path to one of these
     *                sources are created.
     * @param closeReceiverWallDistance If greater than 0, the walls closer than this distance to the receiver do not
     *                                  create first order images: the paths whose last reflection is on such a wall
     *                                  are ignored (see Scene#getCloseReceiverReflectionWallDistance)
     */
    public MirrorReceiversCompute(List<Wall> buildWalls, Coordinate receiverCoordinates,
                                  int reflectionOrder, double maximumPropagationDistance,
                                  double maximumDistanceFromWall, List<Coordinate> sources,
                                  double closeReceiverWallDistance) {
        GeometryFactory gf = new GeometryFactory();
        this.receiverCoordinate = receiverCoordinates;
        this.buildWalls = buildWalls;
        this.maximumDistanceFromWall = maximumDistanceFromWall;
        this.maximumPropagationDistance = maximumPropagationDistance;
        mirrorReceiverTree = new STRtree();
        // Create the geometry of each wall once, and index them, instead of creating a new
        // geometry for every parent image x wall combination in the loop below
        List<LineString> wallGeometries = new ArrayList<>(buildWalls.size());
        for (Wall wall : buildWalls) {
            wallGeometries.add(wall.getLineSegment().toGeometry(gf));
        }
        // With known sources, each source gets the visitor of findCloseMirrorReceivers, fed while the images are
        // created. The visitor rejects an image for a source if a wall of its chain is too far from the
        // source-receiver segment, so for each wall keep the sources that pass this test.
        ReceiverImageVisitor[] sourceVisitors = null;
        BitSet[] wallSources = null;
        if (sources != null) {
            sourceVisitors = new ReceiverImageVisitor[sources.size()];
            visitorBySource = new IdentityHashMap<>();
            for (int idSource = 0; idSource < sources.size(); idSource++) {
                sourceVisitors[idSource] = new ReceiverImageVisitor(buildWalls, sources.get(idSource),
                        receiverCoordinates, maximumDistanceFromWall, maximumPropagationDistance);
                visitorBySource.put(sources.get(idSource), sourceVisitors[idSource]);
            }
            wallSources = computeWallSources(buildWalls, receiverCoordinates, sourceVisitors, maximumDistanceFromWall);
        }
        STRtree wallsTree = null;
        if (reflectionOrder > 1) {
            wallsTree = new STRtree();
            for (int idWall = 0; idWall < buildWalls.size(); idWall++) {
                if (wallSources == null || !wallSources[idWall].isEmpty()) {
                    wallsTree.insert(wallGeometries.get(idWall).getEnvelopeInternal(), idWall);
                }
            }
            wallsTree.build();
        }
        LineIntersector lineIntersector = new RobustLineIntersector();
        ArrayList<MirrorReceiver> parentsToProcess = new ArrayList<>();
        // With known sources: the sources for which every wall of the parent image chain passes the distance test
        ArrayList<BitSet> parentsSources = new ArrayList<>();
        for(int currentDepth = 0; currentDepth < reflectionOrder; currentDepth++) {
            if(currentDepth == 0) {
                parentsToProcess.add(null);
                parentsSources.add(null);
            }
            final boolean lastDepth = currentDepth == reflectionOrder - 1;
            ArrayList<MirrorReceiver> nextParentsToProcess = new ArrayList<>();
            ArrayList<BitSet> nextParentsSources = new ArrayList<>();
            for (int idParent = 0; idParent < parentsToProcess.size(); idParent++) {
                MirrorReceiver parent = parentsToProcess.get(idParent);
                // For the first depth every wall can create an image. For the next depths only
                // the walls under the visibility cone of the parent image can, so ask the wall
                // index instead of testing every wall
                List<?> wallCandidates = null;
                PreparedGeometry parentCone = null;
                if (parent != null) {
                    List<Integer> parentWedgeWalls = new ArrayList<>();
                    queryParentWedge(wallsTree.getRoot(), parent.getImageReceiverVisibilityCone().getEnvelopeInternal(),
                            parent, parentWedgeWalls);
                    wallCandidates = parentWedgeWalls;
                    parentCone = PreparedGeometryFactory.prepare(parent.getImageReceiverVisibilityCone());
                }
                int candidateCount = parent == null ? buildWalls.size() : wallCandidates.size();
                for (int idCandidate = 0; idCandidate < candidateCount; idCandidate++) {
                    int wallIndex = parent == null ? idCandidate : (Integer) wallCandidates.get(idCandidate);
                    Wall wall = buildWalls.get(wallIndex);
                    // The tests below only reject images, from the cheapest to the most expensive one
                    if (parent != null && isOutsideParentWedge(parent, wall.getLineSegment().p0, wall.getLineSegment().p1)) {
                        continue; // cheap version of the visibility cone test below
                    }
                    BitSet imageSources = null;
                    if (wallSources != null) {
                        imageSources = (BitSet) wallSources[wallIndex].clone();
                        if (parent != null) {
                            imageSources.and(parentsSources.get(idParent));
                        }
                        if (imageSources.isEmpty()) {
                            continue; // no source can get a path from this image or from its children
                        }
                    }
                    Coordinate receiverImage;
                    if (parent != null) {
                        if(wall == parent.getWall()) {
                            continue;
                        } else {
                            receiverImage = parent.getReceiverPos();
                        }
                    } else {
                        if (wall.getLineSegment().distance(receiverCoordinates) < closeReceiverWallDistance) {
                            continue; // the paths whose last reflection is on this wall are ignored
                        }
                        receiverImage = receiverCoordinates;
                    }
                    //Calculate the coordinate of projection
                    Coordinate proj = wall.getLineSegment().project(receiverImage);
                    Coordinate rcvMirror = new Coordinate(2 * proj.x - receiverImage.x,
                            2 * proj.y - receiverImage.y, receiverImage.z);
                    if(wall.getLineSegment().distance(rcvMirror) > maximumPropagationDistance) {
                        // wall is too far from the receiver image, there is no receiver image
                        continue;
                    }
                    // Walls that belong to a building (polygon) does not create image receiver
                    // from the two sides of the wall
                    // Exterior polygons are CW we can check if the receiver is on the reflective side of the wall
                    // (on the exterior side of the wall)
                    if(wall.getType() == ProfileBuilder.IntersectionType.BUILDING &&
                            !wallPointTest(wall.getLineSegment(), receiverImage)) {
                        continue;
                    }
                    if (lastDepth && imageSources != null && !hasReflectionPoint(wall.getLineSegment(), rcvMirror,
                            sources, imageSources, lineIntersector)) {
                        continue; // every source visitor would reject this image, and it has no children
                    }
                    if(parent != null) {
                        // check if the wall is visible from the previous image receiver
                        if(!parentCone.intersects(wallGeometries.get(wallIndex))) {
                            continue; // this wall is out of the bound of the receiver visibility
                        }
                    }
                    // create the visibility cone of this receiver image
                    Polygon imageReceiverVisibilityCone = createWallReflectionVisibilityCone(rcvMirror,
                            wall.getLineSegment(), maximumPropagationDistance, maximumDistanceFromWall);
                    MirrorReceiver receiverResultNext = new MirrorReceiver(rcvMirror, parent, wall);
                    receiverResultNext.setImageReceiverVisibilityCone(imageReceiverVisibilityCone);
                    Envelope coneEnvelope = imageReceiverVisibilityCone.getEnvelopeInternal();
                    if (imageSources == null) {
                        mirrorReceiverTree.insert(coneEnvelope, receiverResultNext.copyWithoutCone());
                    } else {
                        // Same visits as the query of the tree with each source position
                        MirrorReceiver receiverResult = receiverResultNext.copyWithoutCone();
                        for (int idSource = imageSources.nextSetBit(0); idSource >= 0;
                             idSource = imageSources.nextSetBit(idSource + 1)) {
                            if (coneEnvelope.intersects(sources.get(idSource))) {
                                sourceVisitors[idSource].visitItem(receiverResult);
                            }
                        }
                    }
                    nextParentsToProcess.add(receiverResultNext);
                    nextParentsSources.add(imageSources);
                    numberOfImageReceivers++;
                    if(numberOfImageReceivers >= mirrorReceiverCapacity) {
                        return;
                    }
                }
            }
            parentsToProcess = nextParentsToProcess;
            parentsSources = nextParentsSources;
        }
        mirrorReceiverTree.build();
    }

    /**
     * @return For each wall, the sources for which the wall is not farther than maximumDistanceFromWall from the
     * source-receiver segment (the distance test of the source visitors)
     */
    private static BitSet[] computeWallSources(List<Wall> walls, Coordinate receiver,
                                               ReceiverImageVisitor[] sourceVisitors, double maximumDistanceFromWall) {
        // Cheap pre-test of the distance test, done in the frame of each receiver-source segment
        double[] directionX = new double[sourceVisitors.length];
        double[] directionY = new double[sourceVisitors.length];
        double[] length = new double[sourceVisitors.length];
        for (int idSource = 0; idSource < sourceVisitors.length; idSource++) {
            Coordinate source = sourceVisitors[idSource].source;
            length[idSource] = receiver.distance(source);
            if (length[idSource] > 0) {
                directionX[idSource] = (source.x - receiver.x) / length[idSource];
                directionY[idSource] = (source.y - receiver.y) / length[idSource];
            }
        }
        BitSet[] wallSources = new BitSet[walls.size()];
        for (int idWall = 0; idWall < walls.size(); idWall++) {
            wallSources[idWall] = new BitSet(sourceVisitors.length);
            LineSegment wallSegment = walls.get(idWall).getLineSegment();
            for (int idSource = 0; idSource < sourceVisitors.length; idSource++) {
                if (!isFarFromSegment(wallSegment, receiver, directionX[idSource], directionY[idSource],
                        length[idSource], maximumDistanceFromWall) &&
                        wallSegment.distance(sourceVisitors[idSource].sourceReceiverSegment) <= maximumDistanceFromWall) {
                    wallSources[idWall].set(idSource);
                }
            }
        }
        return wallSources;
    }

    /**
     * Cheap test of the distance between a wall and a segment
     * @param origin First end of the segment
     * @param directionX Unit direction of the segment, or 0 if its length is 0
     * @param length Length of the segment
     * @return True if both ends of the wall are on the same side of the line of the segment, or beyond the same end
     * of the segment, farther than the distance (plus 1e-6 m so that the rounding can not exclude a wall at this
     * distance)
     */
    private static boolean isFarFromSegment(LineSegment wall, Coordinate origin, double directionX, double directionY,
                                            double length, double distance) {
        double limit = distance + 1e-6;
        double along0 = directionX * (wall.p0.x - origin.x) + directionY * (wall.p0.y - origin.y);
        double along1 = directionX * (wall.p1.x - origin.x) + directionY * (wall.p1.y - origin.y);
        double across0 = directionX * (wall.p0.y - origin.y) - directionY * (wall.p0.x - origin.x);
        double across1 = directionX * (wall.p1.y - origin.y) - directionY * (wall.p1.x - origin.x);
        return (across0 > limit && across1 > limit) || (across0 < -limit && across1 < -limit) ||
                (along0 < -limit && along1 < -limit) || (along0 > length + limit && along1 > length + limit);
    }

    /**
     * Same walls, in the same order, as the query of the walls tree with the envelope, without the tree nodes that
     * are outside the wedge of the parent image
     */
    private static void queryParentWedge(Boundable node, Envelope envelope, MirrorReceiver parent,
                                         List<Integer> walls) {
        Envelope bounds = (Envelope) node.getBounds();
        if (bounds == null || !bounds.intersects(envelope)) {
            return;
        }
        if (node instanceof ItemBoundable) {
            walls.add((Integer) ((ItemBoundable) node).getItem());
        } else if (!isOutsideParentWedge(parent, new Coordinate(bounds.getMinX(), bounds.getMinY()),
                new Coordinate(bounds.getMinX(), bounds.getMaxY()), new Coordinate(bounds.getMaxX(), bounds.getMinY()),
                new Coordinate(bounds.getMaxX(), bounds.getMaxY()))) {
            for (Object child : ((AbstractNode) node).getChildBoundables()) {
                queryParentWedge((Boundable) child, envelope, parent, walls);
            }
        }
    }

    /**
     * The visibility cone of an image is inside the wedge from the image through its wall, beyond its wall.
     * @return True if all the points are outside this wedge of the parent image, by more than a margin that covers
     * the rounding of the cone vertices: a segment or a box with these points can not intersect the visibility cone
     */
    private static boolean isOutsideParentWedge(MirrorReceiver parent, Coordinate... points) {
        Coordinate image = parent.getReceiverPos();
        Coordinate p0 = parent.getWall().getLineSegment().p0;
        Coordinate p1 = parent.getWall().getLineSegment().p1;
        // the inside of each side of the wedge is the side of the other end of the wall, beyond the wall the inside
        // is the side opposite to the image
        return isOutside(image, p0, Math.signum(cross(image, p0, p1)), points) ||
                isOutside(image, p1, Math.signum(cross(image, p1, p0)), points) ||
                isOutside(p0, p1, -Math.signum(cross(p0, p1, image)), points);
    }

    /**
     * @param insideSide 1 if the inside is on the left of the line (a, b), -1 if on the right, 0 if unknown
     * @return True if all the points are outside the line (a, b) by more than 1e-6 m
     */
    private static boolean isOutside(Coordinate a, Coordinate b, double insideSide, Coordinate[] points) {
        double margin = 1e-6 * Math.sqrt((b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y));
        if (insideSide == 0) {
            return false;
        }
        for (Coordinate point : points) {
            if (cross(a, b, point) * insideSide >= -margin) {
                return false;
            }
        }
        return true;
    }

    /** @return The cross product of (b - a) and (p - a), positive if p is on the left of the line (a, b) */
    private static double cross(Coordinate a, Coordinate b, Coordinate p) {
        return (b.x - a.x) * (p.y - a.y) - (b.y - a.y) * (p.x - a.x);
    }

    /**
     * First test of the source visitor: the segment from the receiver image to the source must cross the wall
     * @return True if this test passes for at least one of the sources
     */
    private static boolean hasReflectionPoint(LineSegment wall, Coordinate receiverImage, List<Coordinate> sources,
                                              BitSet sourceIds, LineIntersector lineIntersector) {
        for (int idSource = sourceIds.nextSetBit(0); idSource >= 0; idSource = sourceIds.nextSetBit(idSource + 1)) {
            lineIntersector.computeIntersection(wall.p0, wall.p1, receiverImage, sources.get(idSource));
            if (lineIntersector.hasIntersection()) {
                return true;
            }
        }
        return false;
    }

    /**
     * Occlusion test between one wall and a viewer.
     * Simple Feature Access (ISO 19125-1) say that:
     * On polygon exterior ring are CCW, and interior rings are CW.
     * @param wall1 Wall segment
     * @param pt Observer
     * @return True if the wall is oriented to the point, false if the wall Occlusion Culling (transparent)
     */
    public static boolean wallPointTest(LineSegment wall1, Coordinate pt) {
        return org.locationtech.jts.algorithm.Orientation.isCCW(new Coordinate[]{wall1.getCoordinate(0),
                wall1.getCoordinate(1), pt, wall1.getCoordinate(0)});
    }

    public int getMirrorReceiverCapacity() {
        return mirrorReceiverCapacity;
    }

    public void setMirrorReceiverCapacity(int mirrorReceiverCapacity) {
        this.mirrorReceiverCapacity = mirrorReceiverCapacity;
    }

    public void exportVisibility(StringBuilder sb, double maxPropagationDistance,
                                 double maxPropagationDistanceFromWall, int t, List<MirrorReceiver> MirrorReceiverList, boolean includeHeader) {
        WKTWriter wktWriter = new WKTWriter();
        GeometryFactory factory = new GeometryFactory();
        if(includeHeader) {
            sb.append("the_geom,type,ref_index,ref_order,wall_id,t\n");
        }
        int refIndex = 0;
        for (MirrorReceiver res : MirrorReceiverList) {
            Polygon visibilityCone = createWallReflectionVisibilityCone(
                    res.getReceiverPos(), res.getWall().getLineSegment(),
                    maxPropagationDistance, maxPropagationDistanceFromWall);
            if(!visibilityCone.isEmpty()) {
                int refOrder=1;
                MirrorReceiver parent = res.getParentMirror();
                while (parent != null) {
                    refOrder++;
                    parent = parent.getParentMirror();
                }

                while(res != null) {
                    sb.append("\"");
                    sb.append(wktWriter.write(visibilityCone));
                    sb.append("\",0");
                    sb.append(",").append(refIndex);
                    sb.append(",").append(refOrder);
                    sb.append(",").append(res.getWall().getProcessedObstructionIndex());
                    sb.append(",").append(t).append("\n");
                    sb.append("\"");
                    sb.append(wktWriter.write(factory.createPoint(res.getReceiverPos()).buffer(0.1,
                            12, BufferParameters.CAP_ROUND)));
                    sb.append("\",4");
                    sb.append(",").append(refIndex);
                    sb.append(",").append(refOrder);
                    sb.append(",").append(res.getWall().getProcessedObstructionIndex());
                    sb.append(",").append(t).append("\n");
                    sb.append("\"");
                    sb.append(wktWriter.write(factory.createLineString(new Coordinate[]{res.getWall().line.p0, res.getWall().line.p1}).
                            buffer(0.05, 8, BufferParameters.CAP_SQUARE)));
                    sb.append("\",1");
                    sb.append(",").append(refIndex);
                    sb.append(",").append(refOrder);
                    sb.append(",").append(res.getWall().getProcessedObstructionIndex());
                    sb.append(",").append(t).append("\n");
                    res = res.getParentMirror();
                    if(res != null) {
                        visibilityCone = createWallReflectionVisibilityCone(
                                res.getReceiverPos(), res.getWall().getLineSegment(),
                                maxPropagationDistance, maxPropagationDistanceFromWall);
                    }
                    refOrder-=1;
                }
                refIndex+=1;
            }
        }
        sb.append("\"");
        sb.append(wktWriter.write(factory.createPoint(receiverCoordinate).buffer(0.1, 12, BufferParameters.CAP_ROUND)));
        sb.append("\",2");
        sb.append(",").append(t).append("\n");
    }

    public List<MirrorReceiver> findCloseMirrorReceivers(Coordinate sourcePosition) {
        if(Double.isNaN(sourcePosition.z)) {
            throw new IllegalArgumentException("Not supported NaN z value");
        }
        if (visitorBySource != null) {
            ReceiverImageVisitor sourceVisitor = visitorBySource.get(sourcePosition);
            if (sourceVisitor == null) {
                throw new IllegalArgumentException("The source position has not been given to the constructor");
            }
            return sourceVisitor.result;
        }
        Envelope env = new Envelope(sourcePosition);
        ReceiverImageVisitor receiverImageVisitor = new ReceiverImageVisitor(buildWalls, sourcePosition,
                receiverCoordinate, maximumDistanceFromWall, maximumPropagationDistance);
        mirrorReceiverTree.query(env, receiverImageVisitor);
        return receiverImageVisitor.result;
    }

    private static class ReceiverImageVisitor implements ItemVisitor {
        List<MirrorReceiver> result = new ArrayList<>();
        List<Wall> buildWalls;
        Coordinate source;
        Coordinate receiver;
        LineSegment sourceReceiverSegment;
        double maximumDistanceFromSegment;
        double maximumPropagationDistance;
        int visitedNode = 0;

        public ReceiverImageVisitor(List<Wall> buildWalls, Coordinate source, Coordinate receiver,
                                    double maximumDistanceFromSegment,
                                    double maximumPropagationDistance) {
            this.buildWalls = buildWalls;
            this.source = source;
            this.receiver = receiver;
            this.sourceReceiverSegment = new LineSegment(source, receiver);
            this.maximumDistanceFromSegment = maximumDistanceFromSegment;
            this.maximumPropagationDistance = maximumPropagationDistance;
        }

        @Override
        public void visitItem(Object item) {
            visitedNode++;
            // try to exclude walls without taking into account the topography and other factors
            // we intentionnaly do not check for wall height here as meteo conditions might virtually raise or lower the wall.

            MirrorReceiver receiverImage = (MirrorReceiver) item;
            // Check propagation distance
            if(receiverImage.getReceiverPos().distance3D(source) < maximumPropagationDistance) {
                // Check distance of walls
                MirrorReceiver currentReceiverImage = receiverImage;
                Coordinate reflectionPoint = source;
                while (currentReceiverImage != null) {
                    final Wall currentWall = currentReceiverImage.getWall();
                    final LineSegment currentWallLineSegment = currentWall.getLineSegment();
                    if (currentWallLineSegment.distance(sourceReceiverSegment) > maximumDistanceFromSegment) {
                        return;
                    }
                    // Check if reflection is placed on the wall segment
                    LineSegment srcMirrRcvLine = new LineSegment(currentReceiverImage.getReceiverPos(), reflectionPoint);
                    LineIntersector li = new RobustLineIntersector();
                    li.computeIntersection(currentWallLineSegment.p0, currentWallLineSegment.p1,
                            srcMirrRcvLine.p0, srcMirrRcvLine.p1);
                    if(!li.hasIntersection()) {
                        // No reflection on this wall
                        return;
                    } else {
                        // Set the height for the reflection point.
                        // intersect3D's height is actually sitting vertically between the two lines.
                        // But we want the reflection point to be on the line between the source and the receiver image
                        // So we need to recompute the Z value on the srcMirrRcvLine at intersection coordinates
                        Coordinate intersectionPoint = li.getIntersection(0);
                        double zIntersect = Vertex.interpolateZ(intersectionPoint, srcMirrRcvLine.p0, srcMirrRcvLine.p1);
                        reflectionPoint = new Coordinate(intersectionPoint.x, intersectionPoint.y, zIntersect);
                    }
                    currentReceiverImage = currentReceiverImage.getParentMirror();
                }
                // not rejected
                result.add(receiverImage);
            }
        }
    }
}
