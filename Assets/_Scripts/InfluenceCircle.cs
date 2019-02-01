using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using mattatz.Triangulation2DSystem;
using UnityEngine;

[ExecuteInEditMode]
public class InfluenceCircle : MonoBehaviour
{
    #region Vars

    public Color lineColor = Color.yellow;
    public Material mat;

    public int numberOfSides;
    public float polygonRadius;
    public float influence = 1;
    //public Vector2 polygonCenter;

    public static List<InfluenceCircle> allInfluenceCircles = new List<InfluenceCircle>();

    #endregion



    #region OnEnableDisable

    void OnEnable()
    {
        if( !allInfluenceCircles.Contains(this) ) allInfluenceCircles.Add(this);
    }

    void OnDisable()
    {
        if( allInfluenceCircles.Contains( this ) ) allInfluenceCircles.Remove( this );
    }

    #endregion


    Vector3 lastPos = Vector3.down*9000;
    public bool drawLines = true;
    void OnDrawGizmos()
    {
        Gizmos.color = lineColor/*Color.yellow*/;

/*        if ( IsPositionChanged() )
        {
            DebugDrawPolygon(polygonRadius, /*numberOfSides#1#InfluenceCirclesManager.NumSides());
            GenerateTriangulation();
        }*/

        if ( drawLines )
        {
            foreach( Vector3 point in points)
            {
                Gizmos.DrawLine( transform.position, point );
            }
        }
    }

    public bool IsPositionChanged()
    {
        if (Vector3.Distance(transform.position, lastPos) > 0.5f)
        {
            lastPos = transform.position;
            return true;
        }

        return false;
    }




    public List<Vector3> points = new List<Vector3>();

    public void DebugDrawPolygon()
    {
        DebugDrawPolygon( polygonRadius, InfluenceCirclesManager.NumSides() );
    }

    public void DebugDrawPolygon(  float radius, int numSides, bool raycast=true )
    {

        Vector2 startCorner = new Vector2(radius, 0) ;
        Vector2 previousCorner = startCorner;
        float lastRadius = radius;
        points = new List<Vector3>();

        // For each corner after the starting corner...
        for( int i = 1; i < numSides; i++ )
        {
            // Calculate the angle of the corner in radians.
            float cornerAngle = 2f * Mathf.PI / (float)numSides * i;

            // Get the X and Y coordinates of the corner point.
            Vector2 currentCorner = new Vector2(Mathf.Cos(cornerAngle) * radius, Mathf.Sin(cornerAngle) * radius) ;

            Vector3 currCornerV3 = Vector2ToWorldPoint(currentCorner);
            if (raycast)
            {
                Vector3 origin = transform.position;
                Vector3 dir = (currCornerV3 - origin).normalized;
                float dist = Raycast(origin, dir, radius);
                //dist = Mathf.Lerp(dist, lastRadius, 0.5f);
                currCornerV3 = origin + dir * dist /*( radius - (radius - dist)/2)*/;
                lastRadius = dist;
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
    }



    #region Raycast

    private Vector3 Vector2ToWorldPoint( Vector2 currentCorner )
    {
        return transform.position + new Vector3( currentCorner.x, 0, currentCorner.y );
    }

    List<Vector3> GetCircleDirections()
    {
        var numSides = numberOfSides;
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
        }*/

        return maxDistance;
    }


    public float RaycastLine( Vector3 origin, Vector3 direction, float maxDistance, InfluenceCircle circle)
    {
        float maxDistTEmp = maxDistance;

        Vector3 circlePos = circle.transform.position;
        var dirList = circle.GetCircleDirections();
        var dirListSorted = dirList/*.OrderBy( v => OrderByAngle( direction, v, origin, circlePos ) ).Take(4).ToList()*/;
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
            min = maxDistance - diff;*/

            return min;
        }


        return maxDistance;
    }

    float OrderByAngle(Vector3 direction, Vector3 circleDirection, Vector3 origin, Vector3 circlePos )
    {
        var dirAngleFromCircle =  Vector3.Angle( circleDirection, origin - circlePos/*, Vector3.up*/ );
        var dirAngleFromOrigin = Vector3.Angle(direction,  circlePos - origin /*, Vector3.up*/);

        return Mathf.Abs(dirAngleFromCircle - dirAngleFromOrigin);
    }

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

    #endregion

    [Serializable]
    public class BorderPoint
    {
        public Vector3 point;
    }
    public List<BorderPoint> borderPoints = new List<BorderPoint>();
    //[ContextMenu( "GenerateTriangulation" )]
    public void GenerateTriangulation( List<Vector2> points2D )
    {
        // construct Polygon2D 
        Polygon2D polygon = Polygon2D.Contour(points2D.ToArray());

        // construct Triangulation2D with Polygon2D and threshold angle (18f ~ 27f recommended)
        Triangulation2D triangulation = new Triangulation2D(polygon, 22.5f);

        #region BorderPoints

        borderPoints = new List<BorderPoint>();
        foreach (Segment2D segment in triangulation.Polygon.Segments )
        {
            var midPoint2d = segment.a.Coordinate;
            Vector3 midPoint3d = new Vector3(midPoint2d.x, 0, midPoint2d.y);
            midPoint3d = transform.TransformPoint(midPoint3d);
            //Debug.DrawRay( midPoint3d, Vector3.up * 20, Color.blue, 10 );

            borderPoints.Add( new BorderPoint(){point = midPoint3d } );
        }

        #endregion

        // build a mesh from triangles in a Triangulation2D instance
        Mesh mesh = triangulation.Build();

        //rotate mesh
        Vector3[] vertices = mesh.vertices;
        for (var i = 0; i < vertices.Length; i++)
        {
            var vertex = vertices[i];
            vertex = Quaternion.Euler(90, 0, 0) * vertex;
            vertices[i] = vertex;
        }
        mesh.vertices = vertices;
        mesh.MarkDynamic();


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
                meshRenderer./*material*/sharedMaterial.color = lineColor;
            }
        }
    }


    public void Grow(float growDist/*, float growPow*/ )
    {
        List<BorderPoint> vertices = new List<BorderPoint>( borderPoints);
        if( vertices.Count < 1 ) return;


        //find closest vertex 
        BorderPoint closestPoint = vertices.OrderBy(v=>  Vector3.Distance(v.point, transform.position)).ToList()[0];


        Vector3 closestPoint_Start = closestPoint.point;
/*        Vector3 dirFromCenter = (closestPoint.point - transform.position).normalized * 1000;
        closestPoint.point = Vector3.MoveTowards(closestPoint.point, closestPoint.point + dirFromCenter, growDist );*/
        //vertices[0] = closestPoint;
        //vertices.Add( closestPoint  + Quaternion.Euler(0,90,0)* dirFromCenter.normalized * growDist );


        var nearPoints = /*vertices*/borderPoints.Where( v => Vector3.Distance( v.point, closestPoint_Start )< growDist ).ToList();

        List<BorderPoint> subdPoints = new List<BorderPoint>();
        for( var i = 1; i < nearPoints.Count; i++ )
        {
            var p1 = nearPoints[i].point;
            var p2 = nearPoints[i-1].point;

            var pLerp = Vector3.Lerp(p1,p2, 0.5f);
            pLerp += ( pLerp - transform.position ).normalized * 0.05f;

            //Debug.DrawRay( pLerp, Vector3.up*5, Color.red, 10 );

            int indexOf = borderPoints.IndexOf( nearPoints [ i ] );
            BorderPoint newBorderPoint = new BorderPoint(){point = pLerp };
            borderPoints.Insert( indexOf, newBorderPoint );
            subdPoints.Add(newBorderPoint);
        }

        foreach (var subdPoint in subdPoints)
        {
            nearPoints.Add(subdPoint);
        }

        foreach (BorderPoint nearPoint in nearPoints)
        {
            Vector3 np = nearPoint.point;
            float distToStartPoint = Vector3.Distance(np, closestPoint_Start);
            float towardShift = InfluenceCirclesManager.SampleGrowCurve((growDist - distToStartPoint ) / growDist) *growDist;
            var dirFromCenter = ( np - transform.position ).normalized * towardShift;
            nearPoint.point += dirFromCenter;

        }





        List<Vector2> points2D = new List<Vector2>();
        foreach( var v3 in /*vertices*/borderPoints )
        {
            var v3_local = transform.InverseTransformPoint(v3.point);
            var v2 = new Vector2(v3_local.x, v3_local.z);
            points2D.Add( v2 );
        }


        GenerateTriangulation(points2D);

    }


    public void GenerateInitial()
    {
        DebugDrawPolygon(polygonRadius, InfluenceCirclesManager.NumSides(), false );


        List<Vector2> points2D = new List<Vector2>();

        foreach( Vector3 v3 in points )
        {
            var v3_local = transform.InverseTransformPoint(v3);
            var v2 = new Vector2(v3_local.x, v3_local.z);
            points2D.Add( v2 );
        }
        GenerateTriangulation( points2D );

    }

}
