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


    public static TriangulationResult Triangulate( Polygon pointSet, Transform rootTransform )
    {
        TriangulationResult triangulationResult = new TriangulationResult();

        //pointSet.CalculateWindingOrder();
        //pointSet.DisplayFlipX = true;
        P2T.Triangulate( TriangulationAlgorithm.DTSweep, pointSet );

        pointSet.CalculateWindingOrder();

        float bpDebugShift = 0;
        foreach( var tri in pointSet.Triangles )
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
                bool isConstrained = tri.GetConstrainedEdgeCCW(tri.Points[i]);
                bool hasConstrainedEdge = tri.GetEdgeCCW(tri.Points[i], out edge);
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
                        isBorder = true;
                        //posShift = Vector3.up;
                    }
                    else
                    {
                        //curPen = penErrorCase2;
                        //color = Color.magenta;
                        //posShift = Vector3.up;
                    }
                }

                int nextIndex = (i < 2 ? i + 1 : 0);



                if( isBorder )
                {
                    /*                    Debug.DrawLine(
                                            rootTransform.TransformPoint( pArr [ i ] ) + Vector3.up * bpDebugShift,
                                            rootTransform.TransformPoint( pArr [ nextIndex ] ) + Vector3.up * bpDebugShift,
                                            color, 10 );

                                        bpDebugShift += 0.2f;*/

                    triangulationResult.borderPoints.Add(
                        new BorderPoint()
                        {
                            pointA = /*rootTransform.TransformPoint*/( pArr [ i ] ),
                            pointB = /*rootTransform.TransformPoint*/( pArr [ nextIndex ] )
                        }
                    );
                }
                triangulationResult.vertices.Add( pArr [ i ] );
                triangulationResult.triangles.Add( triangulationResult.vertices.Count - 1 );

            }

        }

        triangulationResult.mesh = BuildMesh( triangulationResult.vertices, triangulationResult.triangles, rootTransform );


        OrderBorderPoints( triangulationResult.borderPoints );
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
        //List<BorderPoint> borderPoints = triangulationResult.borderPoints;
        for( int i = 1; i < borderPoints.Count; i++ )
        {
            var bp = borderPoints[i];
            var bp_prev = borderPoints[i-1];

            if( Vector3.Distance( bp_prev.pointA, bp.pointB ) > 0.01f )
            {
                var nextPb = borderPoints
                    .OrderBy(v=>Vector3.Distance(v.pointB, bp_prev.pointA ) ).ToList()[0];
                int indexOf = borderPoints.IndexOf(nextPb);

                borderPoints.RemoveAt( indexOf );
                borderPoints.Insert( i, nextPb );
            }
        }
    }

    static Mesh BuildMesh( List<Vector3> vertices, List<int> triangles, Transform rooTransform )
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
