using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Poly2Tri;
using Sirenix.OdinInspector;

public class Poly2Tri_Test : MonoBehaviour {

	// Use this for initialization
	void Start ()
    {
        TestTriangulation();

    }

    [Button()]
    private void TestTriangulation()
    {
        List<TriangulationPoint> points = PointGenerator.UniformDistribution(200, 10);
        var pointSet = new /*Constrained*/PointSet(points);
        //pointSet.CalculateWindingOrder();

        Triangulate( pointSet );

    }

    public static void Triangulate( ITriangulatable pointSet )
    {
        P2T.Triangulate( TriangulationAlgorithm.DTSweep, pointSet );


        foreach( var tri in pointSet.Triangles )
        {
            Vector3[] pArr = new Vector3[3];

            pArr [ 0 ] = new Vector3( tri.Points [ 0 ].Xf, 0, tri.Points [ 0 ].Yf );
            pArr [ 1 ] = new Vector3( tri.Points [ 1 ].Xf, 0, tri.Points [ 1 ].Yf );
            pArr [ 2 ] = new Vector3( tri.Points [ 2 ].Xf, 0, tri.Points [ 2 ].Yf );

/*            pArr [ 0 ] = transform.TransformPoint( pArr [ 0 ] );
            pArr [ 1 ] = transform.TransformPoint( pArr [ 1 ] );
            pArr [ 2 ] = transform.TransformPoint( pArr [ 2 ] );*/




            for( int i = 0; i < 3; ++i )
            {
                Color color = Color.green;
                Vector3 posShift = Vector3.zero;

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
                        //posShift = Vector3.up;
                    }
                    else if( isConstrained && !hasConstrainedEdge )
                    {
                        // this will happen when edges are split and is expected
                        //curPen = penErrorCase1;
                        //curPen = penConstrained;
                        color = Color.red;
                        //posShift = Vector3.up;
                    }
                    else
                    {
                        //curPen = penErrorCase2;
                        color = Color.magenta;
                        //posShift = Vector3.up;
                    }
                }

                int nextIndex = (i < 2 ? i + 1 : 0);
                Debug.DrawLine( pArr [ i ] + posShift, pArr [ nextIndex ] + posShift, color, 10 );

            }

        }
    }
}
