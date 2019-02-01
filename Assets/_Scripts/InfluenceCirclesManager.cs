using System.Collections;
using System.Collections.Generic;
using System.Linq;
using Sirenix.OdinInspector;
using UnityEngine;


[ExecuteInEditMode]
public class InfluenceCirclesManager : MonoBehaviour
{
    public int numSides = 32;
    public float segmentSubdDist = 0.3f;

    public float loopWaitTime = 0.3f;

    public AnimationCurve growCurve = new AnimationCurve();
    public static InfluenceCirclesManager _instance;


    void OnEnable()
    {
        _instance = this;

        if (Application.isPlaying)
        {
            StartCoroutine(Loop());
        }
    }


    public int loopTicks = 0;
    IEnumerator Loop()
    {
        loopTicks = 0;
        while (true)
        {
            GrowCircles();
            loopTicks++;

            yield return new WaitForSeconds( loopWaitTime );
        }
    }



    public static int NumSides()
    {
        return _instance!=null ? _instance.numSides : 16;
    }
    public static float SegementSubdDist()
    {
        return _instance != null ?
            _instance.segmentSubdDist 
            : 0.3f;
    }


    public static float SampleGrowCurve(float t)
    {
        AnimationCurve curve = _instance.growCurve;
        float m = curve.Evaluate(t);
        return m;
    }



    /*    void OnDrawGizmos()
        {
            if (InfluenceCircle.allInfluenceCircles.Any(v => v.IsPositionChanged()))
            {
                foreach (var circle in InfluenceCircle.allInfluenceCircles )
                {
                    circle.DebugDrawPolygon();
                }


                foreach( var circle in InfluenceCircle.allInfluenceCircles )
                {
                    List<Vector2> points2D = new List<Vector2>();

                    foreach( Vector3 v3 in circle.points )
                    {
                        var v3_local = transform.InverseTransformPoint(v3);
                        var v2 = new Vector2(v3_local.x, v3_local.z);
                        points2D.Add( v2 );
                    }

                    circle.GenerateTriangulation(points2D);
                }
            }
        }*/


    [Header("Grow")]

    public float growDist = 1;
    //public float growPow = 2;

    [ContextMenu( "GrowCircles" )]
    [Button()]
    void GrowCircles()
    {
        foreach( var circle in InfluenceCircle.allInfluenceCircles )
        {
            circle.Grow( growDist/*, growPow*/ );
        }
    }

    public int growCount = 10;
    [Button()]
    void GrowCircles_NTimes()
    {
        for (int i = 0; i < growCount; i++)
        {
            GrowCircles();
        }
    }



    [ContextMenu( "GenerateInitial" )]
    [Button()]
    void GenerateInitial()
    {
        foreach( var circle in InfluenceCircle.allInfluenceCircles )
        {
            circle.GenerateInitial();
        }
    }
}
