#ifndef RAYTRACINGRENDERER_H
#define RAYTRACINGRENDERER_H

#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <math.h>


//Viewplane window size
#define XMAX 1.0
#define XMIN -1.0
#define YMAX 1.0
#define YMIN -1.0

//Number of spheres and lights
#define NSPHERE 3
#define NLIGHT 2



// Structure declaration_______________________________________________________

//3D data structure to hold 3D points
typedef struct {
	float x;
	float y;
	float z;
}
Point;

//Colour data structure
typedef struct {
	float r;
	float g;
	float b;
}
Colour;

typedef struct _colourGrid_t{
    Colour **data;
    int width;
    int height;
} *ColourGrid;

//light data structure
typedef struct {
	float r;
	float g;
	float b;
	float x;
	float y;
	float z;
}
Light;

//Sphere data structure
typedef struct {
	float x, y, z;//location
	float r, g, b;//colour
	float radius; // radius
}
Sphere;

//Sphere3D data structure with K values for materials
typedef struct {
	float x, y, z;//location
	float radius; // radius
	Colour ka, kd, ks; //the material coeficients
	float n; // the specular reflexion for this sphere
	float ni ; // the refraction index of material

}
Sphere3D;

//create the datatypes for a type Vector and for a type Matrix of Vectors
typedef float Vector[3] ; //means a Vector is a pointer to float
typedef float ** Matrix; //a Matrix is a pointer to a pointer to float



class raytracingRenderer
{
    public:
        raytracingRenderer();
        virtual ~raytracingRenderer();
        void render();
        static ColourGrid getPixels();
    protected:
    private:
        //normalize fonction, to normalize any vector
        float * normalize( Vector );

        //normal fonction, to get the normal vector to a surface
        float *  normal( Point, Sphere3D);

        //viewVector fonction, to get the "Va" vector
        float *  viewVector( Point ) ;

        float * invert(Vector v);

        //lightVector fonction, to get the "La" vector
        float *  lightVector( Point , Light);

        //reflectVector function, to get the "Ra" vector
        float *  reflectVector( Vector, Vector);

        //refractVector function
        float * refractVector(Point, Vector, Vector, float, float );

        //dot fonction to compute scalar product.
        float  dot( Vector , Vector);

        // the reflection intensity function
        Colour reflectionIntensity( Vector, Point, Colour, Sphere3D, Vector);

        //recursive RayTracing function
        void rayTrace( Point , Vector, int , int  );

        int intersectionTest (Point, Vector, Sphere3D );


};

#endif // RAYTRACINGRENDERER_H
