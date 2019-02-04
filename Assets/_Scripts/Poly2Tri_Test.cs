using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using UnityEngine;
using Poly2Tri;
using Sirenix.OdinInspector;

public class Poly2Tri_Test : MonoBehaviour {

	// Use this for initialization
/*	void Start ()
    {
        TestTriangulation();

    }*/

/*    [Button()]
    private void TestTriangulation()
    {
        List<TriangulationPoint> points = PointGenerator.UniformDistribution(200, 10);
        var pointSet = new /*Constrained#1#PointSet(points);
        //pointSet.CalculateWindingOrder();

        Triangulate( pointSet , transform);

    }*/



    public class TriangulationResult
    {
        public List<Vector3> vertices = new List<Vector3>();
        public List<int> triangles = new List<int>();

        public List<BorderPoint> borderPoints = new List<BorderPoint>();

        public Mesh mesh;
    }


    public static TriangulationResult Triangulate( Polygon polygon )
    {
        TriangulationResult triangulationResult = new TriangulationResult();

/*        polygon.WindingOrder = Point2DList.WindingOrderType.CW;
        polygon.CalculateWindingOrder();*/

        P2T.Triangulate( TriangulationAlgorithm.DTSweep, polygon );

        polygon.WindingOrder = Point2DList.WindingOrderType.CW;
        polygon.CalculateWindingOrder();
        polygon.RemoveDuplicateNeighborPoints();


        //float bpDebugShift = 0;
        foreach( var tri in polygon.Triangles )
        {
            Vector3[] pArr = new Vector3[3];

            pArr [ 0 ] = new Vector3( tri.Points [ 0 ].Xf, 0, tri.Points [ 0 ].Yf );
            pArr [ 1 ] = new Vector3( tri.Points [ 1 ].Xf, 0, tri.Points [ 1 ].Yf );
            pArr [ 2 ] = new Vector3( tri.Points [ 2 ].Xf, 0, tri.Points [ 2 ].Yf );

            /*            pArr [ 0 ] = rooTransform.TransformPoint( pArr [ 0 ] );
                        pArr [ 1 ] = rooTransform.TransformPoint( pArr [ 1 ] );
                        pArr [ 2 ] = rooTransform.TransformPoint( pArr [ 2 ] );*/




            for( int i = 0; i < 3; ++i )
            {
                Color color = Color.green;
                //Vector3 posShift = Vector3.zero;
                bool isBorder = false;
                //var curPen = penNormal;
                DTSweepConstraint edge = null;
                bool isConstrained = tri./*GetConstrainedEdgeCCW*/GetConstrainedEdgeCW(tri.Points[i]);
                bool hasConstrainedEdge = tri./*GetEdgeCCW*/GetEdgeCW(tri.Points[i], out edge);
                if( isConstrained || hasConstrainedEdge )
                {
                    if( isConstrained && hasConstrainedEdge )
                    {
                        //curPen = penConstrained;
                        color = Color.red;
                        isBorder = true;
                        //posShift = Vector3.up;
                    }
                    else if( isConstrained && !hasConstrainedEdge )
                    {
                        // this will happen when edges are split and is expected
                        //curPen = penErrorCase1;
                        //curPen = penConstrained;
                        color = Color.red;
                        //isBorder = true;
                        //posShift = Vector3.up;
                    }
                    else
                    {
                        //curPen = penErrorCase2;
                        //color = Color.magenta;
                        //posShift = Vector3.up;
                    }
                }

                //int nextIndex = (i < 2 ? i + 1 : 0);



                if( isBorder && edge != null)
                {
                    
                    triangulationResult.borderPoints.Add(
                        new BorderPoint()
                        {
                            pointA =  /*pArr [ i ]*/ new Vector3( edge.EdgeEnd.Xf, 0, edge.EdgeEnd.Yf ),
                            pointB = /* pArr [ nextIndex ] */ new Vector3( edge.EdgeStart.Xf, 0, edge.EdgeStart.Yf )
                        }
                    );
                }
                triangulationResult.vertices.Add( pArr [ i ] );
                triangulationResult.triangles.Add( triangulationResult.vertices.Count - 1 );

            }

        }

        triangulationResult.mesh = BuildMesh( triangulationResult.vertices, triangulationResult.triangles );

        

        OrderBorderPoints( triangulationResult.borderPoints );

        /*        bpDebugShift = 0;
                foreach (BorderPoint bp in triangulationResult.borderPoints )
                {
                    Debug.DrawLine(
                        rootTransform.TransformPoint( bp.pointA ) + Vector3.up * bpDebugShift,
                        rootTransform.TransformPoint( bp.pointB ) + Vector3.up * bpDebugShift,
                        Color.red, 10 );

                    bpDebugShift += 0.2f;
                }*/


        return triangulationResult;
    }

    public static void OrderBorderPoints( List<BorderPoint> borderPoints )
    {
        const float maxDist = 0.003f;
        //int endReachedCounter = 0;
        List<int> lastBreakedIndexes = new List<int>();

        for( int i = 1; i < borderPoints.Count; i++ )
        {
            //var bp_next = borderPoints[i];
            var bp_curr = borderPoints[i-1];

            if( bp_curr.nextPoint==null 
                //|| Vector3.Distance( bp_curr.pointB, bp_curr.nextPoint.pointA ) > maxDist 
                )
            {

                //if(bpListTemp.Count(v => v != bp_curr && v != bp_next)==0) break;

                BorderPoint newOrdered_Pb = borderPoints
                    .Where(v =>
                    {
                        bool b = false;
                        b = bp_curr != v;
                        if (bp_curr.nextPoint != null)
                        {
                            b &= bp_curr.nextPoint != v;
                        }

                        return b;
                    })
                    /*.OrderBy(v=> Vector3.Distance(v.pointA, bp_curr.pointB )  )
                    .First();*/
                    .FirstOrDefault(v => Vector3.Distance(v.pointA, bp_curr.pointB) < maxDist);
                //if( newOrdered_Pb ==null) continue;
                if( newOrdered_Pb == null )
                {
                    newOrdered_Pb = borderPoints
                        .Where( v => bp_curr != v )
                        .FirstOrDefault( v => Vector3.Distance( v.pointB, bp_curr.pointB ) < maxDist );
                    if( newOrdered_Pb ==null) continue;
                    else
                    {
                        Vector3 pointB_temp = newOrdered_Pb.pointB;
                        newOrdered_Pb.pointB = newOrdered_Pb.pointA;
                        newOrdered_Pb.pointA = pointB_temp;

                        //Debug.Log("Reversed BP");
                    }

                }

                //if( Vector3.Distance( bp_curr.pointB, newOrdered_Pb.pointA ) > 0.1f ) continue;

                //int indexOf = borderPoints.IndexOf(newOrdered_Pb);



                #region Swap

/*                List<BorderPoint> borderPoints_sorted = new List<BorderPoint>();
                //List<int> sortedPointsIndexes = new List<int>();
                BorderPoint newOrdered_Pb_SequenceItem = newOrdered_Pb;
                BorderPoint whileStartPoint = newOrdered_Pb.nextPoint;
                //sortedPointsIndexes.Add(indexOf);
                borderPoints_sorted.Add( newOrdered_Pb );
                borderPoints.Remove( newOrdered_Pb );

                while( ( newOrdered_Pb_SequenceItem = newOrdered_Pb_SequenceItem.nextPoint) != null )
                {
                    if( newOrdered_Pb_SequenceItem == whileStartPoint) break;
                    //sortedPointsIndexes.Add( borderPoints.IndexOf( newOrdered_Pb_SequenceItem ) );
                    borderPoints_sorted.Add( newOrdered_Pb_SequenceItem );
                    borderPoints.Remove(newOrdered_Pb_SequenceItem);

                }*/

                int indexOf = borderPoints.IndexOf(newOrdered_Pb);
                Swap( borderPoints, i, indexOf);
                /*                borderPoints.RemoveAt( indexOf );
                                borderPoints.Insert( i, newOrdered_Pb );*/
/*                int indexOf_curr = borderPoints.IndexOf(bp_curr) + 1;

                borderPoints.InsertRange( indexOf_curr, borderPoints_sorted );*/
/*                for( var index = 0; index < sortedPointsIndexes.Count; index++)
                {
                    int sortedPointsIndex = sortedPointsIndexes[index];
                    int insertIndex = i + index >= borderPoints.Count 
                        ? (i + index) - borderPoints.Count 
                        : i + index;
                    int sortedIndexToInsert = sortedPointsIndex >= borderPoints.Count 
                        ? sortedPointsIndex - borderPoints.Count 
                        : sortedPointsIndex;

                    Swap( 
                        borderPoints,
                        insertIndex,
                        sortedIndexToInsert
                        );
                }*/

                bp_curr.nextPoint = newOrdered_Pb;

                #endregion




/*                if( lastBreakedIndexes.Count==20 && lastBreakedIndexes.All(v=>v == lastBreakedIndexes[0] ) ) break;
                else
                {
                    if ( lastBreakedIndexes.Count>0 && lastBreakedIndexes.All(v=>v==indexOf))
                    {
                        lastBreakedIndexes.Add(indexOf);
                    }
                    else
                    {
                        lastBreakedIndexes.Clear();
                        lastBreakedIndexes.Add(indexOf);
                    }
                }*/

                i = 0;
            }
        }
    }

    public static void Swap<T>( IList<T> list, int indexA, int indexB )
    {
        if ( indexA> -1 && list.Count > indexA 
             && indexB >-1 && list.Count > indexB )
        {
            T tmp = list[indexA];
            list [ indexA ] = list [ indexB ];
            list [ indexB ] = tmp;
        }
    }


    class DistinctBorderPointComparer : IEqualityComparer<BorderPoint>
    {
        #region Implementation of IEqualityComparer<in BorderPoint>

        public bool Equals(BorderPoint bp_1, BorderPoint bp_2)
        {
            return bp_1 != null && bp_2!=null 
                                && Vector3.Distance(bp_1.pointB, bp_2.pointB) < 0.02f;
        }

        public int GetHashCode(BorderPoint obj)
        {
            return /*obj.GetHashCode() ^ *//*obj.pointA.GetHashCode() ^*/ obj.pointB.GetHashCode();
        }

        #endregion
    }



    static Mesh BuildMesh( List<Vector3> vertices, List<int> triangles )
    {
        Mesh mesh = new Mesh();

/*        for (var i = 0; i < vertices.Count; i++)
        {
            vertices [ i ] = Quaternion.Euler(180, 180,180) * /*rooTransform.InverseTransformPoint#1#( vertices [ i ] );
        }*/

        mesh.vertices = vertices.ToArray();
        mesh.triangles = triangles.ToArray();

        Vector3[] normals = new Vector3[vertices.Count];
        for (int i = 0; i < normals.Length; i++)
        {
            normals[i] = Vector3.up;
        }

        mesh.normals = normals;
        mesh.MarkDynamic();

        return mesh;
    }
}
