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
    public int newCircleIterationsCount = 100;

    public AnimationCurve growCurve = new AnimationCurve();
    public Transform circleTransf;

    public static InfluenceCirclesManager _instance;



    void OnEnable()
    {
        _instance = this;

        if (Application.isPlaying)
        {

            StartCoroutine(Loop());
        }
    }

    void Update()
    {
        if (Input.GetMouseButtonDown(0))
        {
            Plane plane = new Plane(Vector3.up, Vector3.zero);
            Ray ray = Camera.main.ScreenPointToRay(Input.mousePosition);
            float raycastDist = 0;

            if (plane.Raycast(ray, out raycastDist))
            {
                Vector3 point = ray.GetPoint(raycastDist);

                if (circleTransf != null)
                {
                    Transform instPrefab = Instantiate(circleTransf, point, Quaternion.identity, transform);
                    InfluenceCircle circle = instPrefab.GetComponent<InfluenceCircle>();

                    if (circle != null)
                    {
                        circle.GenerateInitial();

                        for (int i = 0; i < newCircleIterationsCount; i++)
                        {
                            if(circle==null) break;
                            //circle.GrowAllGrowPoints(growDist);
                            circle.GrowUniform( growDist );

                        }
                    }
                }
            }
        }
    }


    public int loopTicks = 0;
    IEnumerator Loop()
    {
        loopTicks = 0;

        yield return 10;

        GenerateInitial();

        while (true)
        {
            yield return new WaitForSeconds( loopWaitTime );


            GrowCircles();
            loopTicks++;

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
        for (var i = 0; i < InfluenceCircle.allInfluenceCircles.Count; i++)
        {
            var circle = InfluenceCircle.allInfluenceCircles[i];
            //circle.GrowAllGrowPoints(growDist);
            circle.GrowUniform( growDist );

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
