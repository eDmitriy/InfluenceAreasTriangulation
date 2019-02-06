using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using EPPZ.Geometry.AddOns;
using Poly2Tri;
using Sirenix.OdinInspector;
using UnityEngine;
using Edge = EPPZ.Geometry.Model.Edge;


[ExecuteInEditMode]
public class InfluenceCircle : MonoBehaviour
{
    #region Vars
    
    public Transform circleSpawnTr;

    public Color lineColor = Color.yellow;
    public Material mat;

    //public int numberOfSides;
    public float polygonRadius;
    public float influence = 1;
    public int fraction = 1;

    public float area = 0;

    //public Vector2 polygonCenter;




    [HideInInspector]
    public List<BorderPoint> borderPoints = new List<BorderPoint>();


    public float offset = 0.0f;

/*    public enum UpdateMode { Awake, Update, LateUpdate };
    public UpdateMode update = UpdateMode.Awake;

    public enum Coordinates { World, Local }
    public Coordinates coordinates = Coordinates.World;*/

    EPPZ.Geometry.Model.Polygon _polygon;
    EPPZ.Geometry.Model.Polygon _offsetPolygon;
    EPPZ.Geometry.Model.Polygon polygon { get { return ( offset != 0.0f ) ? _offsetPolygon : _polygon; } }




    public static List<InfluenceCircle> allInfluenceCircles = new List<InfluenceCircle>();

    #endregion



    #region OnEnableDisable

    void OnEnable()
    {
        if( !allInfluenceCircles.Contains(this) ) allInfluenceCircles.Add(this);

        if (_polygon == null) _polygon = EPPZ.Geometry.Model.Polygon.PolygonWithSource(this);

    }

    void OnDisable()
    {
        if( allInfluenceCircles.Contains( this ) ) allInfluenceCircles.Remove( this );
    }

    #endregion



    #region GenerateCircle

    public void GenerateCircle( Vector3 centerPoint )
    {
        GenerateCircle( polygonRadius, InfluenceCirclesManager.NumSides(), centerPoint );
    }

    public List<Vector3> GenerateCircle(  float radius, int numSides, Vector3 centerPoint, bool raycast=true )
    {

        Vector2 startCorner = new Vector2(radius, 0) ;
        Vector2 previousCorner = startCorner;
        float lastRadius = radius;
        List<Vector3> points = new List<Vector3>();

        // For each corner after the starting corner...
        for( int i = 1; i < numSides; i++ )
        {
            // Calculate the angle of the corner in radians.
            float cornerAngle = 2f * Mathf.PI / (float)numSides * i;

            // Get the X and Y coordinates of the corner point.
            Vector2 currentCorner = new Vector2(Mathf.Cos(cornerAngle) * radius, Mathf.Sin(cornerAngle) * radius) ;

            Vector3 currCornerV3 = Vector2ToWorldPoint( currentCorner, centerPoint );
            if (raycast)
            {
/*                Vector3 origin = transform.position;
                Vector3 dir = (currCornerV3 - origin).normalized;
                float dist = Raycast(origin, dir, radius);
                //dist = Mathf.Lerp(dist, lastRadius, 0.5f);
                currCornerV3 = origin + dir * dist /*( radius - (radius - dist)/2)#1#;
                lastRadius = dist;*/
            }

            // Draw a side of the polygon by connecting the current corner to the previous one.
            /*Debug*/
            //Gizmos.DrawLine( currCornerV3, Vector2ToWorldPoint(previousCorner) );
            //Gizmos.DrawLine( origin, currCornerV3 );
            points.Add(currCornerV3);

            // Having used the current corner, it now becomes the previous corner.
            previousCorner = currentCorner;
        }

        // Draw the final side by connecting the last corner to the starting corner.
        /*Debug*/
        //Gizmos.DrawLine( Vector2ToWorldPoint(startCorner), Vector2ToWorldPoint(previousCorner) );

        return points;
    }

    #endregion



    #region Raycast

    private Vector3 Vector2ToWorldPoint( Vector2 currentCorner, Vector3 centerPoint )
    {
        return /*transform.position*/centerPoint + new Vector3( currentCorner.x, 0, currentCorner.y );
    }

/*    List<Vector3> GetCircleDirections()
    {
        var numSides = InfluenceCirclesManager.NumSides();/*numberOfSides;#1#
        var radius = polygonRadius;

        List<Vector3> directions = new List<Vector3>();

        for( int i = 1; i < numSides; i++ )
        {
            // Calculate the angle of the corner in radians.
            float cornerAngle = 2f * Mathf.PI / (float)numSides * i;

            // Get the X and Y coordinates of the corner point.
            Vector2 currentCorner = new Vector2(Mathf.Cos(cornerAngle) * radius, Mathf.Sin(cornerAngle) * radius) ;

            Vector3 currCornerV3 = Vector2ToWorldPoint(currentCorner);
            Vector3 origin = transform.position;
            Vector3 dir = (currCornerV3 - origin).normalized;

            directions.Add(dir);
        }

        return directions;
    }



    public float Raycast(Vector3 origin, Vector3 direction, float maxDistance)
    {
        var influenceCircles = allInfluenceCircles
            .Where(v=>v!=this &&  Vector3.Distance(v.transform.position, origin) < maxDistance+v.polygonRadius )
            .Where(v=> Vector3.Angle(direction, v.transform.position - origin)<89)
            //.OrderBy(v=> Vector3.Angle(direction, v.transform.position - origin))
            .ToList();

        List<float> intersectionDistancesList = new List<float>();

        foreach( InfluenceCircle influenceCircle in influenceCircles)
        {
            float f = RaycastLine( origin, direction, maxDistance, influenceCircle );
            intersectionDistancesList.Add(f);
        }

        if( intersectionDistancesList.Count>0 ) maxDistance = intersectionDistancesList.Min();
/*        InfluenceCircle circle = influenceCircles.FirstOrDefault();


        if (circle != null)
        {
            maxDistance = RaycastLine(origin, direction, maxDistance, circle);
        }#1#

        return maxDistance;
    }


    public float RaycastLine( Vector3 origin, Vector3 direction, float maxDistance, InfluenceCircle circle)
    {
        float maxDistTEmp = maxDistance;

        Vector3 circlePos = circle.transform.position;
        var dirList = circle.GetCircleDirections();
        var dirListSorted = dirList/*.OrderBy( v => OrderByAngle( direction, v, origin, circlePos ) ).Take(4).ToList()#1#;
        List<float> intersectionDistancesList = new List<float>();


        for (var i = 1; i < dirListSorted.Count; i++)
        {
            //Vector3 cDir = dirListSorted[i];
            Vector3 circleDir_1 = dirListSorted[i];//cDir;
            Vector3 circleDir_2 = dirListSorted[i-1];//cDir;

            Vector3 originLineEnd = origin + direction * maxDistance;
            Vector3 circleDirEnd_1 = circlePos + circleDir_1 * circle.polygonRadius;
            Vector3 circleDirEnd_2 = circlePos + circleDir_2 * circle.polygonRadius;

            Vector2 intersectionPoint = origin;


            if ( LineIntersection(
                new Vector2(origin.x, origin.z) * 10,
                new Vector2(originLineEnd.x, originLineEnd.z) * 10,
                //new Vector2(circlePos.x, circlePos.z) * 10, 
                new Vector2(circleDirEnd_1.x, circleDirEnd_1.z) * 10,
                new Vector2( circleDirEnd_2.x, circleDirEnd_2.z ) * 10,
                ref intersectionPoint )
            )
            {
                intersectionPoint /= 10;
                maxDistance = Vector3.Distance(new Vector3(intersectionPoint.x, 0, intersectionPoint.y), origin);
                intersectionDistancesList.Add(maxDistance);
            }
        }

        if ( intersectionDistancesList.Count>0 )
        {
            float min = intersectionDistancesList.Min();

/*            float diff = Vector3.Distance(origin, circlePos ) - min;
            min = maxDistance - diff;#1#

            return min;
        }


        return maxDistance;
    }

    float OrderByAngle(Vector3 direction, Vector3 circleDirection, Vector3 origin, Vector3 circlePos )
    {
        var dirAngleFromCircle =  Vector3.Angle( circleDirection, origin - circlePos/*, Vector3.up#1# );
        var dirAngleFromOrigin = Vector3.Angle(direction,  circlePos - origin /*, Vector3.up#1#);

        return Mathf.Abs(dirAngleFromCircle - dirAngleFromOrigin);
    }*/

    #endregion




    #region LineIntersection

    public static bool LineIntersection( Vector2 p1, Vector2 p2, Vector2 p3, Vector2 p4, ref Vector2 intersection )
    {
        float Ax,Bx,Cx,Ay,By,Cy,d,e,f,num,offset;
        float x1lo,x1hi,y1lo,y1hi;

        Ax = p2.x - p1.x;
        Bx = p3.x - p4.x;

        // X bound box test/
        if( Ax < 0 )
        {
            x1lo = p2.x; x1hi = p1.x;
        }
        else
        {
            x1hi = p2.x; x1lo = p1.x;
        }

        if( Bx > 0 )
        {
            if( x1hi < p4.x || p3.x < x1lo ) return false;
        }
        else
        {
            if( x1hi < p3.x || p4.x < x1lo ) return false;
        }

        Ay = p2.y - p1.y;
        By = p3.y - p4.y;

        // Y bound box test//
        if( Ay < 0 )
        {
            y1lo = p2.y; y1hi = p1.y;
        }
        else
        {
            y1hi = p2.y; y1lo = p1.y;
        }

        if( By > 0 )
        {
            if( y1hi < p4.y || p3.y < y1lo ) return false;
        }
        else
        {
            if( y1hi < p3.y || p4.y < y1lo ) return false;
        }

        Cx = p1.x - p3.x;
        Cy = p1.y - p3.y;
        d = By * Cx - Bx * Cy;  // alpha numerator//
        f = Ay * Bx - Ax * By;  // both denominator//

        // alpha tests//
        if( f > 0 )
        {
            if( d < 0 || d > f ) return false;
        }
        else
        {
            if( d > 0 || d < f ) return false;
        }

        e = Ax * Cy - Ay * Cx;  // beta numerator//

        // beta tests //
        if( f > 0 )
        {
            if( e < 0 || e > f ) return false;
        }
        else
        {
            if( e > 0 || e < f ) return false;
        }

        // check if they are parallel
        if( f == 0 ) return false;

        // compute intersection coordinates //
        num = d * Ax; // numerator //
        offset = same_sign( num, f ) ? f * 0.5f : -f * 0.5f;   // round direction //
        intersection.x = p1.x + ( num + offset ) / f;

        num = d * Ay;
        offset = same_sign( num, f ) ? f * 0.5f : -f * 0.5f;
        intersection.y = p1.y + ( num + offset ) / f;

        return true;
    }

    private static bool same_sign( float a, float b )
    {
        return ( ( a * b ) >= 0f );
    }


    public class LineCast_GrowResult
    {
        public float dist;
        public BorderPoint intersectionBP;
    }

    private LineCast_GrowResult LineCastForGrow( BorderPoint nearPoint, float towardShift, float shiftCorrection,
        Vector3 segmentNormal, List<BorderPoint> allSceneBorderPoints )
    {
        Vector3 desiredPos = nearPoint.pointA + segmentNormal * towardShift;
        Vector2 instersectPoint = Vector2.zero;
        float vectMult = 10f;

        BorderPoint intersectionBP = null;

        foreach( var bp in allSceneBorderPoints )
        {
            if (bp == nearPoint || Vector3.Angle(bp.normal, nearPoint.normal)<90)
            {
                continue;
            }

            if( LineIntersection(
                new Vector2( nearPoint.pointA.x, nearPoint.pointA.z ) * vectMult,
                new Vector2( desiredPos.x, desiredPos.z ) * vectMult,
                new Vector2( bp.pointA.x, bp.pointA.z ) * vectMult,
                new Vector2( bp.pointB.x, bp.pointB.z ) * vectMult,
                ref instersectPoint
            ) )
            {
                instersectPoint /= vectMult;
                float distance = Vector3.Distance(nearPoint.pointA, new Vector3(instersectPoint.x,0,instersectPoint.y)  );
                //if(distance<0.01f) continue;

                distance = Mathf.Clamp( distance - shiftCorrection/*0.2f*/, 0, float.MaxValue );
                towardShift = Mathf.Min( towardShift, distance );

                //Debug.DrawRay( new Vector3( instersectPoint.x, 0, instersectPoint.y ), Vector3.up, Color.magenta, 10 );
                intersectionBP = bp;
                break;
            }
        }

        return /*towardShift*/new LineCast_GrowResult()
        {
            dist = towardShift, intersectionBP = intersectionBP
        };
    }

    #endregion




    #region Triangulation

    public void GenerateTriangulation()
    {
        //Poly2Tri_Test.OrderBorderPoints(borderPoints);

        /*        float bpDebugShift = 0;
                foreach (BorderPoint bp in borderPoints )
                {
                    Debug.DrawLine(
                        bp.pointA + Vector3.up * bpDebugShift,
                        bp.pointB + Vector3.up * bpDebugShift,
                        Color.red, InfluenceCirclesManager._instance.loopWaitTime );

                    bpDebugShift += 0.2f;
                }*/
        UpdateModel();
        return;


        List<PolygonPoint> polygonPoints = new List<PolygonPoint>();
        PolygonPoint polygonPoint_Prev = null;

        foreach( var p in borderPoints )
        {
            Vector3 pPointA = p.pointA;
            pPointA = transform.InverseTransformPoint(pPointA);
            PolygonPoint polygonPoint = new PolygonPoint( pPointA.x, pPointA.z );
            polygonPoint.Previous = polygonPoint_Prev;
            if (polygonPoint_Prev != null) polygonPoint_Prev.Next = polygonPoint;
            polygonPoint_Prev = polygonPoint;

            polygonPoints.Add( polygonPoint );

            //Debug.DrawRay( pPointA , Vector3.up, Color.blue, 1);
        }

        if ( polygonPoints.Count>2 )
        {
            Polygon polygon = new Polygon(polygonPoints);
            polygon.Precision = TriangulationPoint.kVertexCodeDefaultPrecision;
/*            polygon.WindingOrder = Point2DList.WindingOrderType.CW;
            polygon.CalculateWindingOrder();*/
            //polygon.RemoveDuplicateNeighborPoints();
            //polygon.MergeParallelEdges(0.1f);
            //polygon.CalculateWindingOrder();

            GenerateTriangulation( polygon);
        }
    }

    private Polygon polygon_triangulated = null;

    public void GenerateTriangulation( /*List<Vector2> points2D*/Polygon polygon )
    {
        var triangulationResults = Poly2Tri_Test.Triangulate( polygon );
        polygon_triangulated = polygon;
        area = (float) polygon_triangulated.GetArea();


        #region BorderPoints

        borderPoints = new List<BorderPoint>();

        foreach (var bp in triangulationResults.borderPoints )
        {
            Vector3 pA_v3 = transform.TransformPoint( bp.pointA );
            Vector3 pB_v3 = transform.TransformPoint( bp.pointB );

            Vector3 segmentVector = (pA_v3-pB_v3).normalized;
            Vector3 segmentNormal = Quaternion.Euler(0, -90/*-90*/, 0) * segmentVector ;


            BorderPoint newBP = new BorderPoint()
            {
                pointA = pA_v3,
                pointB = pB_v3,
                normal = segmentNormal,
                circle = this,
            };
            borderPoints.Add( newBP );
        }

        #endregion

/*        float bpDebugShift = 0;
        foreach( BorderPoint bp in borderPoints )
        {
            Debug.DrawLine(
                bp.pointA + Vector3.up * bpDebugShift,
                bp.pointB + Vector3.up * bpDebugShift,
                Color.red, InfluenceCirclesManager._instance.loopWaitTime );

            bpDebugShift += 0.2f;
        }*/


        #region Setup Unity mesh

        // build a mesh from triangles in a Triangulation2D instance
        Mesh mesh = triangulationResults.mesh;// triangulation.Build();


        MeshFilter meshFilter = GetComponent<MeshFilter>();
        if (meshFilter == null) meshFilter = gameObject.AddComponent<MeshFilter>();
        if (meshFilter != null) meshFilter.sharedMesh = mesh;

        MeshRenderer meshRenderer = GetComponent<MeshRenderer>();
        if( meshRenderer == null ) meshRenderer = gameObject.AddComponent<MeshRenderer>();
        if (meshRenderer != null)
        {
            if (mat != null)
            {
                meshRenderer.sharedMaterial = new Material(mat);
                meshRenderer./*material*/sharedMaterial.SetColor( "_FillColor", lineColor);
            }
        }

/*        var meshColl = GetComponent<MeshCollider>();
        if( meshColl == null ) meshColl = gameObject.AddComponent<MeshCollider>();
        if( meshColl != null )
        {
            meshColl.sharedMesh = mesh;
        }*/

        #endregion

    }

    #endregion


    public void Grow(float growDist )
    {
        List<BorderPoint> vertices = new List<BorderPoint>( borderPoints);
        if( vertices.Count < 1 ) return;


        #region Calc grow points

        Vector3 centerOfMass = CenterOfMass_BorderPoints( );
        var allSceneBorderPoints = allInfluenceCircles
            .Where(v=>v!=this)
            .SelectMany( v => v.borderPoints ).ToList();

        var notBlocked_BP_Points = vertices
            .Where(v=>
                allSceneBorderPoints
                    .Where(c=>c.circle.fraction!=fraction)
                    .All(b=> Vector3.Distance(b.pointA, v.pointA)>0.5f)
                )
            .ToList();
        if(notBlocked_BP_Points.Count==0) return;

        BorderPoint closestPoint = notBlocked_BP_Points
            .OrderBy(v=>  Vector3.Distance(v.pointA, centerOfMass)).ToList()[0];
        Vector3 closestPoint_Start = closestPoint.pointA;

        var growPoints = borderPoints.Where( v => Vector3.Distance( v.pointA, closestPoint_Start )< growDist ).ToList();

        #endregion

        //subdivide
/*        if (CheckIfNeedSubdivide(growPoints))
        {
            growPoints = borderPoints.Where( v => Vector3.Distance( v.pointA, closestPoint_Start ) < growDist ).ToList();
        }*/


        List< BorderPoint > growedEdges = new List<BorderPoint>(borderPoints);


        foreach( BorderPoint nearPoint in growPoints )
        {
            Vector3 np = nearPoint.pointA;
            float distToStartPoint = Vector3.Distance(np, closestPoint_Start);
            float towardShift = InfluenceCirclesManager.SampleGrowCurve((growDist - distToStartPoint ) / growDist) *growDist;


            #region Normal

            Vector3 segmentNormal = nearPoint.normal;
/*            if (Vector3.Angle(segmentNormal, averageNormal) > 100)
            {
                segmentNormal *= -1;
            }
            averageNormal = Vector3.Lerp(averageNormal, segmentNormal, 0.3f);*/

            #endregion

            //raycast against others
            /*towardShift*/var lineCast_growResult = LineCastForGrow(
                nearPoint, towardShift, 0.2f, segmentNormal, allSceneBorderPoints 
                );
            towardShift = lineCast_growResult.dist;

            #region Connect similar fraction circles

            BorderPoint ibp = lineCast_growResult.intersectionBP;
            if ( ibp!=null && ibp.circle.fraction == fraction)
            {
                ConnectTwoCircles( growDist, nearPoint, ibp );
                break;
            }

            #endregion

            Vector3 newPos = nearPoint.pointA + segmentNormal * towardShift;


            #region raycast against self

            int indexOfnearPoint = borderPoints.IndexOf(nearPoint);

            var intersectionBP = LineCastForGrow( nearPoint, towardShift +0.2f, 0, segmentNormal, growedEdges ).intersectionBP;
            if ( intersectionBP!=null )
            {

                borderPoints.Remove( nearPoint );
/*                Debug.DrawRay( nearPoint.pointA, Vector3.up*10, 
                    Color.magenta, InfluenceCirclesManager._instance.loopWaitTime );*/

                continue;
            }

            var newPosLocal = transform.InverseTransformPoint(newPos);
            if ( polygon_triangulated != null 
                 && polygon_triangulated.IsPointInside( new TriangulationPoint( newPosLocal.x, newPosLocal.z ) ) )
            {
                borderPoints.Remove( nearPoint );

/*                Debug.DrawLine(nearPoint.pointA + Vector3.up*0.1f, newPos + Vector3.up * 0.1f,
                    Color.magenta, InfluenceCirclesManager._instance.loopWaitTime );
                Debug.DrawRay( nearPoint.pointA, Vector3.up * 10,
                    Color.magenta, InfluenceCirclesManager._instance.loopWaitTime );*/

                //Debug.LogError(this);
                continue;
            }

            #region AddGrowedEdges

            Vector3 segmentVector = (nearPoint.pointA-newPos).normalized;
            Vector3 newGrowedEdgeNormal = Quaternion.Euler(0, -90, 0) * segmentVector ;
            growedEdges.Add( new BorderPoint()
            {
                pointA = nearPoint.pointA,
                pointB = newPos,
                normal = newGrowedEdgeNormal
            } );
            segmentVector = (  newPos - nearPoint.pointA ).normalized;
            newGrowedEdgeNormal = Quaternion.Euler( 0, -90, 0 ) * segmentVector;
            growedEdges.Add( new BorderPoint()
            {
                pointA = newPos,
                pointB = nearPoint.pointA,
                normal = newGrowedEdgeNormal
            } );

            #endregion

            #endregion


            Debug.DrawRay( nearPoint.pointA + Vector3.up*0.25f, segmentNormal * towardShift, 
                Color.green, InfluenceCirclesManager._instance.loopWaitTime );
            
            nearPoint.pointA = newPos;
            //nearPoint.pointB += segmentNormal * towardShift;

/*            if (indexOfnearPoint > 1 && indexOfnearPoint + 1 < borderPoints.Count )
            {
                borderPoints[indexOfnearPoint - 1].pointB = newPos;
                //borderPoints[indexOf + 1].pointA = nearPoint.pointB;
            }*/

        }


/*        foreach (var borderPoint in borderPoints)
        {
            Debug.DrawRay(borderPoint.pointA, Vector3.up, Color.cyan, InfluenceCirclesManager._instance.loopWaitTime);
        }*/


        GenerateTriangulation();
    }


    #region Points order on connection

    private void MakeCorrectOrderOnConnect( InfluenceCircle ibpCircle, BorderPoint lastBP )
    {
        //BorderPoint lastBP = borderPoints.Last();
        Debug.DrawRay( lastBP.pointA + Vector3.up * 0.25f, Vector3.up * 20, Color.blue, 30 );

        while( CheckCorrectOrder( ibpCircle.borderPoints, lastBP ) /*ibpCircle.borderPoints[0] is not closest to last*/)
        {
            var intersectionFirstPoint = ibpCircle.borderPoints[0];
            ibpCircle.borderPoints.Insert( ibpCircle.borderPoints.Count /*- 1*/, intersectionFirstPoint );
            ibpCircle.borderPoints.RemoveAt(0);
        }
    }

    bool CheckCorrectOrder( List<BorderPoint> bpPoints, BorderPoint lastBP )
    {
        List<BorderPoint> bplist = new List<BorderPoint>(/*ibpCircle.borderPoints*/bpPoints);
        BorderPoint cp = bplist.OrderBy(v=>Vector3.Distance(v.pointA, lastBP.pointA)).First();
        return  bplist.IndexOf( cp ) != 0 ;
    }


    private void ConnectTwoCircles( float growDist, BorderPoint nearPoint, BorderPoint ibp )
    {
        InfluenceCircle ibpCircle = ibp.circle;

        List<BorderPoint> intersectPoints = ibpCircle.borderPoints
            .Where(v=>Vector3.Distance(v.pointA, nearPoint.pointA)<growDist/2).ToList();
        foreach( var intersectPoint in intersectPoints )
        {
            ibpCircle.borderPoints.Remove( intersectPoint );
        }
        int indexOf = borderPoints.IndexOf(nearPoint);

        //borderPoints.Remove(nearPoint);
        intersectPoints = borderPoints
            .Where( v => Vector3.Distance( v.pointA, nearPoint.pointA ) < growDist / 2 ).ToList();
        foreach( var intersectPoint in intersectPoints )
        {
            borderPoints.Remove( intersectPoint );
        }

        MakeCorrectOrderOnConnect( ibpCircle, nearPoint );
        borderPoints.InsertRange( indexOf, ibpCircle.borderPoints );

        /*                foreach( var bp in borderPoints )
                        {
                            Debug.DrawRay( bp.pointA, Vector3.up * 1, Color.red, 5 );
                        }*/

        Destroy( ibpCircle.gameObject );
    }


    #endregion



    private Vector3 CenterOfMass_BorderPoints(  )
    {
        Vector3 centerOfMass = Vector3.zero;
        foreach( var bp in borderPoints )
        {
            centerOfMass += bp.pointA;
        }
        centerOfMass /= borderPoints.Count;
        return centerOfMass;
    }


    #region Subdivision

    private bool CheckIfNeedSubdivide( List<BorderPoint> growPoints )
    {
        float segementSubdDist = InfluenceCirclesManager.SegementSubdDist();
        if( growPoints.Any( v => Vector3.Distance( v.pointA, v.pointB ) > segementSubdDist ) )
        {
            SubdSegments( segementSubdDist );
            return true;
        }

        return false;
    }

    private void SubdSegments( float segementSubdDist )
    {
        //float segementSubdDist = InfluenceCirclesManager.SegementSubdDist();

        for( var i = 1; i < borderPoints.Count-1; i++ )
        {
            //BorderPoint bp = borderPoints[i];
            var p1 = borderPoints[i].pointA;
            var p2 = borderPoints[i-1].pointA;

            if(Vector3.Distance(p1,p2) <  segementSubdDist/*0.3f*/) continue;

            var pLerp = Vector3.Lerp(p1,p2, 0.5f);

            int indexOf = i;//borderPoints.IndexOf( bp );
            BorderPoint newBorderPoint = new BorderPoint()
            {
                pointA = pLerp,
                pointB = p1,
                normal = borderPoints[i].normal,
                //fraction = borderPoints[i].fraction,
                circle = borderPoints[i].circle
            };
            //bp.pointB = pLerp;
            borderPoints.Insert( indexOf, newBorderPoint );
            i--;
            //i = indexOf+1;
        }
    }

        #endregion



    #region Buttons

    [Button(), ContextMenu( "GenerateInitial" )]
    public void GenerateInitial()
    {
        var points = GenerateCircle(polygonRadius, InfluenceCirclesManager.NumSides(), transform.position, false );


/*        List<Vector2> points2D = new List<Vector2>();

        foreach( Vector3 v3 in points )
        {
            var v3_local = transform.InverseTransformPoint(v3);
            var v2 = new Vector2(v3_local.x, v3_local.z);
            points2D.Add( v2 );
        }*/

/*        List<PolygonPoint> polygonPoints = new List<PolygonPoint>();
        foreach( Vector3 p in points )
        {
            //polygonPoints.Add( new PolygonPoint( p.x, p.z ) );

            Vector3 pPointA = p;
            pPointA = transform.InverseTransformPoint( pPointA );
            polygonPoints.Add( new PolygonPoint( pPointA.x, pPointA.z ) );
        }


        Polygon polygon = new Polygon(polygonPoints);

        GenerateTriangulation( /*points2D#1#polygon );*/
        borderPoints = new List<BorderPoint>();
        foreach( Vector3 p in points )
        {
            borderPoints.Add( new BorderPoint()
            {
                pointA = p
            });
        }
        GenerateTriangulation( );

    }

    [Button()]
    public void GenerateCircle()
    {
        var points = GenerateCircle( polygonRadius, InfluenceCirclesManager.NumSides(),transform.position, false );

        foreach (Vector3 p in points)
        {
            borderPoints.Add(new BorderPoint(){pointA = p});
        }

        GenerateTriangulation( );
    }

    [Button()]
    public void GenerateCircles_2()
    {
        var points = GenerateCircle( polygonRadius, InfluenceCirclesManager.NumSides(),transform.position, false );

        foreach( Vector3 p in points )
        {
            borderPoints.Add( new BorderPoint() { pointA = p } );
        }

        points = GenerateCircle( polygonRadius * 0.5f, InfluenceCirclesManager.NumSides(), transform.position, false );

        foreach( Vector3 p in points )
        {
            borderPoints.Add( new BorderPoint() { pointA = p } );
        }

        GenerateTriangulation();
    }


    [Button()]
    public void GenerateCircle_AtPoint() //broken
    {
        var points = GenerateCircle( polygonRadius, InfluenceCirclesManager.NumSides(), circleSpawnTr.position, false );
        
        foreach( Vector3 p in points )
        {
            borderPoints.Add( new BorderPoint() { pointA = p } );
        }
        
        //GenerateTriangulation();
    }


    [Button()]
    public void GenerateCircle_Poly2Tri()
    {
        var points = GenerateCircle( polygonRadius, InfluenceCirclesManager.NumSides(),transform.position, false );


        List<PolygonPoint> polygonPoints = new List<PolygonPoint>();
        foreach( Vector3 p in points )
        {
            polygonPoints.Add( new PolygonPoint( p.x, p.z ) );

            //Debug.DrawRay(p, Vector3.up, Color.white, 10);
        }


        Polygon polygon = new Polygon(polygonPoints);
        polygon.Precision = TriangulationPoint.kVertexCodeDefaultPrecision;
        Poly2Tri_Test.Triangulate( polygon );
    }

    #endregion


    void UpdateModel()
    {
        // Update polygon model with transforms, also update calculations.
        _polygon = EPPZ.Geometry.Model.Polygon.PolygonWithSource( this );
        _polygon.UpdatePointPositionsWithSource( this );

        if( offset != 0.0f )
        {
            //Model.Polygon.clipperArcTolerance = clipperArcTolerance;
            //Model.Polygon.clipperScale = clipperScale;

            _offsetPolygon = _polygon.SimplifiedNotRoundedOffsetPolygon(offset); 
            //.SimplifiedAndRoundedOffsetPolygon( offset );
        }



        #region Border points
        

        borderPoints = new List<BorderPoint>();
        polygon.EnumerateEdgesRecursive( ( Edge eachEdge ) =>
        {
            borderPoints.Add( new BorderPoint()
            {
                circle = this,
                pointA = transform.TransformPoint( new Vector3( eachEdge.a.x, 0, eachEdge.a.y ) ),
                pointB = transform.TransformPoint( new Vector3( eachEdge.b.x, 0, eachEdge.b.y ) ),
                normal = new Vector3( eachEdge.normal.x, 0, eachEdge.normal.y ).normalized

            } );

        } );


/*        if( InfluenceCirclesManager._instance != null )
        {
            float bpDebugShift = 0;
            float step = 3f / (float)borderPoints.Count;
            foreach( BorderPoint bp in borderPoints )
            {
                Debug.DrawLine(
                    bp.pointA + Vector3.up * bpDebugShift,
                    bp.pointB + Vector3.up * bpDebugShift,
                    Color.red, InfluenceCirclesManager._instance.loopWaitTime );

                bpDebugShift += step;
            }
        }*/

        #endregion

        

        #region Setup mesh

        //cast points to local space
/*        for( var i = 0; i < polygon.points.Length; i++ )
        {
            Vector3 localPoint = transform.InverseTransformPoint( new Vector3( polygon.points [ i ].x, 0, polygon.points [ i ].y) );
            polygon.points [ i ] = new Vector2( localPoint.x, localPoint.z );
        }*/


        MeshFilter meshFilter = GetComponent<MeshFilter>();
        meshFilter.sharedMesh = polygon.Mesh( lineColor, TriangulatorType.Dwyer );

        MeshRenderer meshRenderer = GetComponent<MeshRenderer>();
        if( meshRenderer == null ) meshRenderer = gameObject.AddComponent<MeshRenderer>();
        if( meshRenderer != null )
        {
            if( mat != null /*&& meshRenderer.sharedMaterial ==null*/)
            {
                meshRenderer.sharedMaterial = new Material( mat );
                //meshRenderer./*material*/sharedMaterial.SetColor( "_FillColor", lineColor );
                meshRenderer./*material*/sharedMaterial.SetColor( "_Color", lineColor );

            }
        }

        #endregion

    }

}



[Serializable]
public class BorderPoint
{
    public Vector3 pointA;
    public Vector3 pointB;
    public Vector3 normal;
    //public int fraction;
    public InfluenceCircle circle;

    public BorderPoint prevPoint;
    public BorderPoint nextPoint;
}

