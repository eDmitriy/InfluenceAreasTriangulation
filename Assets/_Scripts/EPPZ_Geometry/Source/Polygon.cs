//
// Copyright (c) 2017 Geri Borbás http://www.twitter.com/_eppz
//
//  Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
//  The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
using UnityEngine;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using UnityEngine.Serialization;



namespace EPPZ.Geometry.Source
{


	public class Polygon : MonoBehaviour
	{


		[UnityEngine.Serialization.FormerlySerializedAs("pointTransforms")]
		//public Transform[] points;
        public List<Transform> points = new List<Transform>();
	    public Transform [ ] GetPoints
	    {
	        get
	        {
/*	            points.Clear();

                foreach( Transform child in transform )
	            {
	                if (child != transform)
	                {
	                    points.Add( child );
	                }
	            }*/

	            return points.ToArray();
	        }
            set { points = value.ToList(); }
	    }



        public float offset = 0.0f;

		public enum UpdateMode { Awake, Update, LateUpdate };
		public UpdateMode update = UpdateMode.Awake;	

		public enum Coordinates { World, Local }
		public Coordinates coordinates = Coordinates.World;	

		Model.Polygon _polygon;
		Model.Polygon _offsetPolygon;		
		public Model.Polygon polygon { get { return (offset != 0.0f) ? _offsetPolygon : _polygon; } }




	    void Awake()
		{
			// Construct a polygon model from transforms (if not created by a root polygon already).
			//if (_polygon == null) _polygon = Model.Polygon.PolygonWithSource(this);
		    if (_polygon == null)
		    {
		        _polygon = Model.Polygon.PolygonWithSource( this );
		    }

            if( offset != 0.0f) _offsetPolygon = _polygon.SimplifiedAndRoundedOffsetPolygon( offset )/*OffsetPolygon(offset)*/;
		}

		void Update()
		{
			if (update == UpdateMode.Update)
			{ UpdateModel(); }
		}

		void LateUpdate()
		{
			if (update == UpdateMode.LateUpdate)
			{ UpdateModel(); }
		}


        //public float clipperArcTolerance = 10e+3f; // 2 magnitude smaller
        //public float clipperScale = 10e+5f;

        [ContextMenu( "UpdateModel" )]
        void UpdateModel()
		{
			// Update polygon model with transforms, also update calculations.
            _polygon.UpdatePointPositionsWithSource(this);
		    //_polygon = _polygon.SimplifiedAndRoundedOffsetPolygon( 0.05f )/*UnionPolygon() */;

		    if (offset != 0.0f)
		    {
                //Model.Polygon.clipperArcTolerance = clipperArcTolerance;
                //Model.Polygon.clipperScale = clipperScale;

                _offsetPolygon = _polygon./*OffsetPolygon(offset)*/SimplifiedAndRoundedOffsetPolygon( offset );
		    }
		}
	}
}
