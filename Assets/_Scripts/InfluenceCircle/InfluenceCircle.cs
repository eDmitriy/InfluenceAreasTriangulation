using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using EPPZ.Geometry;
using EPPZ.Geometry.AddOns;
using EPPZ.Geometry.Model;
using Sirenix.OdinInspector;
using UnityEngine;
using Edge = EPPZ.Geometry.Model.Edge;


[ExecuteInEditMode]
public class InfluenceCircle : MonoBehaviour
{
    #region Vars
    
    //public Transform circleSpawnTr;

    public Color lineColor = Color.yellow;
    public Material mat;

    public float polygonRadius;
    //public float influence = 1;
    public int fraction = 1;
    public float area = 0;

    [HideInInspector]
    public List< List<BorderPoint>> borderPoints_subPolygons = new List<List<BorderPoint>>();
    public List<BorderPoint> borderPoints = new List<BorderPoint>();
    //public List<Vector3> growPoints = new List<Vector3>();

    #region EPPZ_Polygon

    float offset = 0;
    
    EPPZ.Geometry.Model.Polygon _polygon;
    EPPZ.Geometry.Model.Polygon _offsetPolygon;
    EPPZ.Geometry.Model.Polygon polygon { get { return ( true/*offset != 0.0f*/ ) ? _offsetPolygon : _polygon; } }

    #endregion

    
    public static List<InfluenceCircle> allInfluenceCircles = new List<InfluenceCircle>();

    #endregion



    #region OnEnableDisable

    void OnEnable()
    {
        if( !allInfluenceCircles.Contains(this) ) allInfluenceCircles.Add(this);

        //if (_polygon == null) _polygon = EPPZ.Geometry.Model.Polygon.PolygonWithSource(this);

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

            Vector3 currCornerV3 = centerPoint + new Vector3(currentCorner.x, 0, currentCorner.y);
            //Vector2ToWorldPoint( currentCorner, centerPoint );
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



    #region LineIntersection

    public class LineCast_GrowResult
    {
        public float dist;
        public BorderPoint intersectionBP;
    }

    private LineCast_GrowResult LineCastForGrow( BorderPoint nearPoint, float towardShift, float shiftCorrection,
        Vector3 segmentNormal, List<BorderPoint> allSceneBorderPoints, bool checkNormals )
    {
        Vector3 desiredPos = nearPoint.pointA + segmentNormal * towardShift;
        //Vector2 instersectPoint = Vector2.zero;
        float vectMult = 100f;

        BorderPoint intersectionBP = null;

        foreach( var bp in allSceneBorderPoints )
        {
            if ( checkNormals && bp == nearPoint || Vector3.Angle(bp.normal, nearPoint.normal)<90)
            {
                continue;
            }


            Segment segment_1 = new Segment()
            {
                a = nearPoint.pointA.xz(),
                b = desiredPos.xz(),
                //normal = p1.normal.xz()
            };
            Segment segment_2 = new Segment()
            {
                a= bp.pointA.xz(),
                b = bp.pointB.xz(),
                normal = bp.normal.xz()
            };

            Vector2 intersectionPoint;
            if( segment_1.IntersectionWithSegmentWithAccuracy( segment_2, 0.1f, out intersectionPoint ) )
            {
                float distance = Vector3.Distance(nearPoint.pointA, new Vector3(intersectionPoint.x,0,intersectionPoint.y)  );

                distance = Mathf.Clamp( distance - shiftCorrection, 0, float.MaxValue );
                towardShift = Mathf.Min( towardShift, distance );

                intersectionBP = bp;
                break;
            }



/*            if( LineIntersection(
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

                distance = Mathf.Clamp( distance - shiftCorrection/*0.2f#1#, 0, float.MaxValue );
                towardShift = Mathf.Min( towardShift, distance );

                //Debug.DrawRay( new Vector3( instersectPoint.x, 0, instersectPoint.y ), Vector3.up, Color.magenta, 10 );
                intersectionBP = bp;
                break;
            }*/
        }

        return /*towardShift*/new LineCast_GrowResult()
        {
            dist = towardShift, intersectionBP = intersectionBP
        };
    }

    #endregion


    


    #region Grow

/*    public void GrowAllGrowPoints(float growDist)
    {
        for (var i = 0; i < growPoints.Count; i++)
        {
            Vector3 growPoint = growPoints[i];
            Grow(growDist, growPoint);
        }
    }*/


    public void Grow(float growDist, Vector3 centerOfMass )
    {
        #region Remove parallel border points

/*        for( var i = 0; i < borderPoints.Count; i++ )
        {
            var bp = borderPoints[i];
            BorderPoint firstOrDefault = borderPoints.FirstOrDefault(v =>v!= bp && BorderPointsNormalIntersected(bp, v, 20));
            if( firstOrDefault != null )
            {
                borderPoints.Remove( firstOrDefault );
                borderPoints.Remove( bp );

                i = -1;
            }
        }*/

        #endregion

        borderPoints_subPolygons = new List<List<BorderPoint>>(){borderPoints};

        List<BorderPoint> vertices = new List<BorderPoint>( borderPoints);
        if( vertices.Count < 1 ) return;


        #region Calc grow points

        var allSceneBorderPoints = allInfluenceCircles
            .Where(v=>v!=this)
            .SelectMany( v => v.borderPoints ).ToList();
        var allSceneBorderPoints_NotThisFraction = allSceneBorderPoints
            .Where(c => c.circle.fraction != fraction).ToList();


        var notBlocked_BP_Points = vertices
            .Where(v=>
                allSceneBorderPoints_NotThisFraction
                    .All(b=> Vector3.Distance(b.pointA, v.pointA)>0.5f/*PointNotBlocked(b, v)*/ )
            )
            .ToList();

/*        foreach (BorderPoint bp in notBlocked_BP_Points )
        {
            Debug.DrawRay(bp.pointA , Vector3.up,
                Color.red, InfluenceCirclesManager._instance.loopWaitTime );
        }*/

        if( notBlocked_BP_Points.Count==0) return;


        //Vector3 centerOfMass = CenterOfMass_BorderPoints( borderPoints);


        BorderPoint closestPoint = notBlocked_BP_Points
            .OrderBy(v=>  Vector3.Distance(v.pointA, centerOfMass)).ToList()[0];
        Vector3 closestPoint_Start = closestPoint.pointA;

        var growPoints = borderPoints.Where( v => Vector3.Distance( v.pointA, closestPoint_Start )< growDist ).ToList();

        #endregion

        //subdivide
        if (CheckIfNeedSubdivide(growPoints))
        {
            growPoints = borderPoints.Where( v => Vector3.Distance( v.pointA, closestPoint_Start ) < growDist ).ToList();
        }


        List< BorderPoint > growedEdges = new List<BorderPoint>(borderPoints);


        //GROW
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
                nearPoint, towardShift, 0.2f, segmentNormal, allSceneBorderPoints, false
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

            //int indexOfnearPoint = borderPoints.IndexOf(nearPoint);

            var intersectionBP = LineCastForGrow( 
                nearPoint, towardShift +0.2f, 0, segmentNormal, growedEdges, true ).intersectionBP;
            if ( intersectionBP!=null )
            {

                borderPoints.Remove( nearPoint );
/*                Debug.DrawRay( nearPoint.pointA, Vector3.up*10, 
                    Color.magenta, InfluenceCirclesManager._instance.loopWaitTime );*/

                continue;
            }

/*            var newPosLocal = transform.InverseTransformPoint(newPos);
            if ( polygon != null 
                 && polygon.ContainsPoint( new Vector2( newPosLocal.x, newPosLocal.z ) ) )
            {
                borderPoints.Remove( nearPoint );

/*                Debug.DrawLine(nearPoint.pointA + Vector3.up*0.1f, newPos + Vector3.up * 0.1f,
                    Color.magenta, InfluenceCirclesManager._instance.loopWaitTime );
                Debug.DrawRay( nearPoint.pointA, Vector3.up * 10,
                    Color.magenta, InfluenceCirclesManager._instance.loopWaitTime );#1#

                //Debug.LogError(this);
                continue;
            }*/

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
        }


        //GenerateTriangulation();
        UpdateModel();
    }


    public void GrowUniform( float growDist )
    {
        #region Calc grow points

        SubdSegments( InfluenceCirclesManager.SegementSubdDist() );


        var allSceneBorderPoints = allInfluenceCircles
            .Where(v=>v!=this)
            .SelectMany( v => v.borderPoints ).ToList();
        var allSceneBorderPoints_NotThisFraction = allSceneBorderPoints
            .Where(c => c.circle.fraction != fraction).ToList();


        var notBlocked_BP_Points = /*vertices*/borderPoints
            .Where(v=>
                allSceneBorderPoints_NotThisFraction
                    .All(b=> Vector3.Distance(b.pointA, v.pointA)>0.25f )
            )
            .ToList();


        if( notBlocked_BP_Points.Count == 0 ) return;


        var growPoints = new List<BorderPoint>(/*notBlocked_BP_Points*/notBlocked_BP_Points);

        List< BorderPoint > growedEdges = new List<BorderPoint>(borderPoints);

        #endregion

        borderPoints_subPolygons = new List<List<BorderPoint>>() { borderPoints };
        //List<BorderPoint> connectionBorderPoints = new List<BorderPoint>();

        //GROW
        foreach( BorderPoint nearPoint in growPoints )
        {
            Vector3 np = nearPoint.pointA;
            float towardShift = growDist;
            Vector3 segmentNormal = nearPoint.normal;


            //raycast against others

            var lineCast_growResult = LineCastForGrow(
                nearPoint, towardShift, 0.2f, segmentNormal, allSceneBorderPoints, false
            );
            towardShift = lineCast_growResult.dist;

            #region Connect similar fraction circles

            BorderPoint ibp = lineCast_growResult.intersectionBP;
            if( ibp != null && ibp.circle.fraction == fraction )
            {
                ConnectTwoCircles( growDist, nearPoint, ibp );
                break;
                //connectionBorderPoints.Add(ibp);
            }

            #endregion

            Vector3 newPos = nearPoint.pointA + segmentNormal * towardShift;


            #region raycast against self

            //int indexOfnearPoint = borderPoints.IndexOf(nearPoint);

            var intersectionResult = LineCastForGrow(
                nearPoint, towardShift +0.2f, 0, segmentNormal, growedEdges, true );
            if( intersectionResult.intersectionBP != null )
            {
                BorderPoint intBP = intersectionResult.intersectionBP;

/*                if (nearPoint.nextPoint != null) nearPoint.nextPoint.prevPoint = null;
                if (nearPoint.prevPoint != null) nearPoint.prevPoint.nextPoint = null;*/
                borderPoints.Remove( nearPoint );

/*                if( intBP.nextPoint != null ) intBP.nextPoint.prevPoint = null;
                if( intBP.prevPoint != null ) intBP.prevPoint.nextPoint = null;*/
                borderPoints.Remove( intBP );

/*                Debug.DrawLine(nearPoint.pointA + Vector3.up, nearPoint.pointB +Vector3.up, Color.red, 0.1f);
                Debug.DrawLine( intersectionResult.intersectionBP.pointA + Vector3.up,
                    intersectionResult.intersectionBP.pointB + Vector3.up, Color.red, 0.1f );*/

                continue;
                //newPos = nearPoint.pointA + segmentNormal * intersectionResult.dist*1.05f;
            }

            #region Is inside polygon

            /*            var newPosLocal = transform.InverseTransformPoint(newPos);
                        if ( polygon != null 
                             && polygon.ContainsPoint( new Vector2( newPosLocal.x, newPosLocal.z ) ) )
                        {
                            borderPoints.Remove( nearPoint );

            /*                Debug.DrawLine(nearPoint.pointA + Vector3.up*0.1f, newPos + Vector3.up * 0.1f,
                                Color.magenta, InfluenceCirclesManager._instance.loopWaitTime );
                            Debug.DrawRay( nearPoint.pointA, Vector3.up * 10,
                                Color.magenta, InfluenceCirclesManager._instance.loopWaitTime );#1#

                            //Debug.LogError(this);
                            continue;
                        }*/

            #endregion


            #region AddGrowedEdges

            Vector3 segmentVector = (nearPoint.pointA-newPos).normalized;
            Vector3 newGrowedEdgeNormal = Quaternion.Euler(0, -90, 0) * segmentVector ;
            growedEdges.Add( new BorderPoint()
            {
                pointA = nearPoint.pointA,
                pointB = newPos,
                normal = newGrowedEdgeNormal
            } );
            segmentVector = ( newPos - nearPoint.pointA ).normalized;
            newGrowedEdgeNormal = Quaternion.Euler( 0, -90, 0 ) * segmentVector;
            growedEdges.Add( new BorderPoint()
            {
                pointA = newPos,
                pointB = nearPoint.pointA,
                normal = newGrowedEdgeNormal
            } );

            #endregion

            #endregion


/*            Debug.DrawRay( nearPoint.pointA + Vector3.up * 0.25f, segmentNormal * towardShift,
                Color.green, InfluenceCirclesManager._instance.loopWaitTime );*/

            nearPoint.pointA = newPos;
        }

        
        //GenerateTriangulation();
        UpdateModel();
    }

    #endregion




    #region Points order on connection

    private void MakeCorrectOrderOnConnect( InfluenceCircle ibpCircle, BorderPoint lastBP )
    {
        //BorderPoint lastBP = borderPoints.Last();
        //Debug.DrawRay( lastBP.pointA + Vector3.up * 0.25f, Vector3.up * 20, Color.blue, 30 );

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


    private void ConnectTwoCircles( BorderPoint ibp)
    {
        InfluenceCircle ibpCircle = ibp.circle;
        borderPoints_subPolygons.Add( ibpCircle.borderPoints );
        Destroy( ibpCircle.gameObject );
    }


    private void ConnectTwoCircles( float growDist, BorderPoint nearPoint, BorderPoint ibp )
    {
        InfluenceCircle ibpCircle = ibp.circle;
        int indexOf = borderPoints.IndexOf(nearPoint);


        #region Remove points on intersection radius

/*        List<BorderPoint> intersectPoints = ibpCircle.borderPoints
            .Where(v=>Vector3.Distance(v.pointA, nearPoint.pointA)<growDist/2).ToList();
        foreach( var intersectPoint in intersectPoints )
        {
            ibpCircle.borderPoints.Remove( intersectPoint );
        }

        //borderPoints.Remove(nearPoint);
        intersectPoints = borderPoints
            .Where( v => Vector3.Distance( v.pointA, nearPoint.pointA ) < growDist / 2 ).ToList();
        foreach( var intersectPoint in intersectPoints )
        {
            borderPoints.Remove( intersectPoint );
        }*/

        #endregion


        MakeCorrectOrderOnConnect( ibpCircle, nearPoint );
        borderPoints.InsertRange( indexOf, ibpCircle.borderPoints );
        //growPoints.AddRange(ibpCircle.growPoints);


        //borderPoints_subPolygons.Add( ibpCircle.borderPoints );

        Destroy( ibpCircle.gameObject );
    }


    #endregion



    #region Utils

    private Vector3 CenterOfMass_BorderPoints(List<BorderPoint> borderPointsList  )
    {
        Vector3 centerOfMass = Vector3.zero;
        foreach( var bp in borderPointsList )
        {
            centerOfMass += bp.pointA;
        }
        centerOfMass /= borderPointsList.Count;
        return centerOfMass;
    }

    bool PointNotBlocked( BorderPoint p1, BorderPoint p2 )
    {
        //bool result = false;

        /*        float dist = 0.5f;
                bool result = Vector3.Distance(p1.pointA, p2.pointA) < dist
                       || Vector3.Distance(p1.pointA, p2.pointB) < dist
                       || Vector3.Distance(p1.pointB, p2.pointB) < dist
                       || Vector3.Distance(p1.pointB, p2.pointA) < dist;
                result = !result;

                if (result == true) return result;*/


        Vector3 startMidPoint = Vector3.Lerp(p1.pointA, p1.pointB, 0.5f);
        Vector3 desiredPos = p1.pointA/*startMidPoint*/ + p1.normal * 0.3f;
        //Vector2 instersectPoint = Vector2.zero;
        float mult = 1;

        Segment segment_1 = new Segment()
        {
            a = /*startMidPoint*/p1.pointA.xz(),
            b = desiredPos.xz(),
            //normal = p1.normal.xz()
        };
        Segment segment_2 = new Segment()
        {
            a= p2.pointA.xz(),
            b = p2.pointB.xz(),
            normal = p2.normal.xz()
        };

        Vector2 intersectionPoint;
        if( segment_1.IntersectionWithSegmentWithAccuracy( segment_2, 0.01f, out intersectionPoint ) )
        {
            return false;
        }


        /*        if (LineIntersection(
                    new Vector2( startMidPoint.x, startMidPoint.z) * mult,
                    new Vector2(desiredPos.x, desiredPos.z) * mult,
                    new Vector2(p2.pointA.x, p2.pointA.z) * mult,
                    new Vector2(p2.pointB.x, p2.pointB.z) * mult,
                    ref instersectPoint
                ))
                {
                    return false;
                }*/

        return true;
    }
    bool BorderPointsNormalIntersected( BorderPoint p1, BorderPoint p2, float angle )
    {
        if( Vector2.Angle( -p1.normal.xz(), p2.normal.xz() ) > angle ) return false;


        Vector3 startMidPoint = Vector3.Lerp(p1.pointA, p1.pointB, 0.5f);
        Vector3 desiredPos = startMidPoint + p1.normal * 0.3f;
        //Vector2 instersectPoint = Vector2.zero;
        float mult = 1;

        Segment segment_1 = new Segment()
        {
            a = startMidPoint.xz(),
            b = desiredPos.xz(),
            normal = (Quaternion.Euler(0,-90,0)* (desiredPos - startMidPoint).normalized).normalized.xz()
        };
        Segment segment_2 = new Segment()
        {
            a= p2.pointA.xz(),
            b = p2.pointB.xz(),
            normal = p2.normal.xz()
        };

        Vector2 intersectionPoint;
        if( segment_1.IntersectionWithSegmentWithAccuracy( segment_2, 0.1f, out intersectionPoint ) )
        {
            return true;
        }

        return false;
    }

    #endregion



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
                circle = borderPoints[i].circle,
/*                nextPoint = borderPoints[i],
                prevPoint = borderPoints[i-1]*/
            };
            //bp.pointB = pLerp;
            borderPoints.Insert( indexOf, newBorderPoint );
            //i--;
            i = 0;
        }
    }

        #endregion



    #region Buttons

    [Button(), ContextMenu( "GenerateInitial" )]
    public void GenerateInitial()
    {
        //growPoints.Add(transform.position);

        var points = GenerateCircle(polygonRadius, InfluenceCirclesManager.NumSides(), transform.position, false );


        borderPoints = new List<BorderPoint>();
        foreach( Vector3 p in points )
        {
            borderPoints.Add( new BorderPoint()
            {
                pointA = p
            });
        }

        borderPoints_subPolygons = new List<List<BorderPoint>>() { borderPoints };

        UpdateModel();
    }


    [ContextMenu( "GenerateInitial_ring" )]
    public void GenerateInitial_ring()
    {
        //growPoints.Add(transform.position);

        var points = GenerateCircle(polygonRadius, InfluenceCirclesManager.NumSides(), transform.position, false );
        var points_2 = GenerateCircle(polygonRadius*1.5f, InfluenceCirclesManager.NumSides(), transform.position, false );
        borderPoints_subPolygons = new List<List<BorderPoint>>();


        borderPoints = new List<BorderPoint>();
        foreach( Vector3 p in points )
        {
            borderPoints.Add( new BorderPoint()
            {
                pointA = p
            } );
        }

        borderPoints_subPolygons.Add( borderPoints );


        borderPoints = new List<BorderPoint>();
        foreach( Vector3 p in points_2 )
        {
            borderPoints.Add( new BorderPoint()
            {
                pointA = p
            } );
        }
        borderPoints.Reverse();

        borderPoints_subPolygons.Add( borderPoints);


        UpdateModel();
    }



    #endregion


    /*    public float clipperScale = 1000;
        public float clipperArcTolerance = 10e+3f; // 2 magnitude smaller*/

    void UpdateModel()
    {
        // Update polygon model with transforms, also update calculations.
        _polygon = Polygon.PolygonWithSource( this );
        //_polygon.UpdatePointPositionsWithSource( this );

        //if( offset != 0.0f )
        {
/*            EPPZ.Geometry.Model.Polygon.clipperArcTolerance = clipperArcTolerance;
            EPPZ.Geometry.Model.Polygon.clipperScale = clipperScale;*/
            _offsetPolygon = _polygon.SimplifiedNotRoundedOffsetPolygon(offset);
            //_offsetPolygon = _polygon.SimplifiedAndRoundedOffsetPolygon( offset );
            //.SimplifiedAndRoundedOffsetPolygon( offset );
        }



        #region Border points


        borderPoints = new List<BorderPoint>();
        //BorderPoint prevPoint = null;
        polygon.EnumerateEdgesRecursive( ( Edge eachEdge ) =>
        {
            BorderPoint newBorderPoint = new BorderPoint()
            {
                circle = this,
                pointA = transform.TransformPoint( new Vector3( eachEdge.a.x, 0, eachEdge.a.y ) ),
                pointB = transform.TransformPoint( new Vector3( eachEdge.b.x, 0, eachEdge.b.y ) ),
                normal = new Vector3( eachEdge.normal.x, 0, eachEdge.normal.y ).normalized,
                //prevPoint = prevPoint

            };
            borderPoints.Add( newBorderPoint );
/*            if (prevPoint != null) prevPoint.nextPoint = newBorderPoint;
            prevPoint = newBorderPoint;*/
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


        area = polygon.area;
    }

}



[Serializable]
public class BorderPoint
{
    public Vector3 pointA;
    public Vector3 pointB;
    public Vector3 normal;
    public InfluenceCircle circle;

/*    public BorderPoint prevPoint;
    public BorderPoint nextPoint;*/
}

