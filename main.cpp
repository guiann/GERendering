//Raycasting - involves following the path of rays from the COP - centre of camera projection through the centre of
//scene pixels.  Follow each of these rays to find that if it intersects an object to represent the colour of the
//corresponding scene pixel on the viewing plane.  At the end of this process scene pixles will form an image on the
//screen

#include <GL/glut.h>
#include <stdio.h>
#include <stdlib.h>
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

// Global variables____________________________________________________________

//array to hold three elements type sphere
Sphere3D sphere [3];
Light lights[2];

//Centre of camera projection co-ordinates
Point COP = {0, 0, 3 };

//define an ambiant light for the scene
Colour Ia = {0.3, 0.3, 0.3};

//Pixel array for our viewing plane
Colour pixels [512][512];

//count the number of recursive call for rayTrace
int count;

// Functions prototypes________________________________________________________

//Display function
void display();

//normalize fonction, to normalize any vector
float * normalize( Vector );

//normal fonction, to get the normal vector to a surface
float *  normal( Point, Sphere3D);

//viewVector fonction, to get the "Va" vector
float *  viewVector( Point ) ;

//lightVector fonction, to get the "La" vector
float *  lightVector( Point , Light);

//reflectVector function, to get the "Ra" vector
float *  reflectVector( Vector, Vector) ;

//refractVector function
float * refractVector(Point, Vector, Vector, float, float );

//dot fonction to compute scalar product.
float  dot( Vector , Vector);

// the reflection intensity function
Colour reflectionIntensity( Vector, Point, Colour, Sphere3D, Vector);

//recursive RayTracing function
void rayTrace( Point , Vector, int , int  );

int intersectionTest (Point, Vector, Sphere3D );



// Functions bodies____________________________________________________________

int main(int argc, char **argv){

	//Creating graphics window
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);//GLUT DOUBLE prevents flickering of screen
	glutInitWindowSize(512, 512);
	glutInitWindowPosition(0, 0);
	glutCreateWindow("RayTracing");
	glutDisplayFunc(display);
	glMatrixMode(GL_PROJECTION);

	glLoadIdentity();
	glOrtho(0.0, 512.0, 0.0, 512.0, -100.0, 100.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	//Lights information

	//light 1 coordinates and colour
	lights[0].x = -3;
	lights[0].y = 2;
	lights[0].z = 0;
	lights[0].r = 1;
	lights[0].g = 1;
	lights[0].b = 1;

	//light 2 coordinates and colour
	lights[1].x = 5;
	lights[1].y = -5;
	lights[1].z = 3;
	lights[1].r = 1;
	lights[1].g = 1;
	lights[1].b = 1;

	//Sphere information.  Co-ordinates, radius and materials specifications_______

	//Sphere 1 co-ordinates, radius and colour
	sphere[0].x = 0.0;
	sphere[0].y = -0.5;
	sphere[0].z = -2.0;

	sphere[0].radius = 0.6;

	//Sphere 1 materials
	sphere[0].ka.r = 0.6;
	sphere[0].ka.g = 0.0;
	sphere[0].ka.b = 0.0;

	sphere[0].kd.r = 0.5;
	sphere[0].kd.g = 0.0;
	sphere[0].kd.b = 0.0;

	sphere[0].ks.r = 1;
	sphere[0].ks.g = 1;
	sphere[0].ks.b = 1;

	sphere[0].n = 120;
	sphere[0].ni = 1.309;

	//Sphere 2 co-ordinates, radius and materials specifications_______
	sphere[1].x = 1.0;
	sphere[1].y = 0.0;
	sphere[1].z = -3.5;

	sphere[1].radius = 1.0;

	//Sphere 2 materials
	sphere[1].ka.r = 0.5;
	sphere[1].ka.g = 0.5;
	sphere[1].ka.b = 0.0;

	sphere[1].kd.r = 0.5;
	sphere[1].kd.g = 0.5;
	sphere[1].kd.b = 0.0;

	sphere[1].ks.r = 0.9;
	sphere[1].ks.g = 0.9;
	sphere[1].ks.b = 0.9;

	sphere[1].n =120;
	sphere[1].ni = -1;

	//Sphere 3 co-ordinates, radius and colour
	sphere[2].x = -1.0;
	sphere[2].y = 1.0;
	sphere[2].z = -3.0;

	sphere[2].radius = 1.0;

	//Sphere 3 materials
	sphere[2].ka.r = 0.0;
	sphere[2].ka.g = 0.0;
	sphere[2].ka.b = 0.5;

	sphere[2].kd.r = 0.0;
	sphere[2].kd.g = 0.0;
	sphere[2].kd.b = 0.5;

	sphere[2].ks.r = 0.9;
	sphere[2].ks.g = 0.9;
	sphere[2].ks.b = 0.9;

	sphere[2].n = 120;
	sphere[2].ni = -1;


	int i,j;
	float dx,dy,dz;
	Vector direction;


	//Determing the width and height of each pixel
	//width and height of each scene pixel
	float pixelwidth, pixelheight;

	//Centre of pixels type point
	Point centreofpix;

	pixelwidth = (XMAX - XMIN) / 512;
	pixelheight = (YMAX - YMIN) / 512;



	//Incrementing through our scene
	for ( i = 0; i < 512; i++){
		for ( j = 0; j < 512; j++){

			pixels[i][j].r= 0.0;
			pixels[i][j].g= 0.0;
			pixels[i][j].b= 0.0;

			//Determining the centre of our pixels
			centreofpix.x = XMIN + pixelwidth * (i + 0.5);
			centreofpix.y = YMIN + pixelheight * (j + 0.5);
			centreofpix.z = 0;

			//dx, dy, and dz are the values of the directional vector of the first ray
			dx = centreofpix.x - COP.x;
			dy = centreofpix.y - COP.y;
			dz = centreofpix.z - COP.z;

			direction[0] = dx/ sqrt(dx*dx + dy*dy + dz*dz);
			direction[1] = dy/ sqrt(dx*dx + dy*dy + dz*dz);
			direction[2] = dz/ sqrt(dx*dx + dy*dy + dz*dz);

			count=0;
			//calling rayTrace function
			rayTrace(COP, direction, i , j );

			pixels[i][j].r = pixels[i][j].r/count;
			pixels[i][j].g = pixels[i][j].g/count;
			pixels[i][j].b = pixels[i][j].b/count;

			if(pixels[i][j].r < 0) pixels[i][j].r=0;
			if(pixels[i][j].g < 0) pixels[i][j].g=0;
			if(pixels[i][j].b < 0) pixels[i][j].b=0;

			if(pixels[i][j].r > 1) pixels[i][j].r=1;
			if(pixels[i][j].g > 1) pixels[i][j].g=1;
			if(pixels[i][j].b > 1) pixels[i][j].b=1;
		}
	}
	glutMainLoop();

  return 0;
}

//Display function creating a graphics window
void display()
{

	int r, c;

	// clear the offscreen buffer
	glClearColor(0.0, 0.0, 0.0, 0.0);/*Background colour*/
	glClear(GL_COLOR_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	//Displaying scene pixels and forming an image
	glBegin(GL_POINTS);

	for (r = 0; r < 512; r++)
		for (c = 0; c < 512; c++)
		{
			glColor3f(pixels[r][c].r, pixels[r][c].g, pixels[r][c].b);
			glVertex2f(r, c);
		}

	glEnd();

	glFlush();

	// swap the visible on the offscreen buffer - required to prevent flickering.
	//Function used in conjunction with GLUT DOUBLE
	glutSwapBuffers();

	// ask for glut to recall display() immediately
	glutPostRedisplay();

	return;
}


/***************************************************************************/
/*we'll need to normalize several vectors*/
float * normalize( Vector v){
	float norm= sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
	int i;
	for(	i=0 ; i<=2 ; i++){
		v[i] = v[i]/norm;
	}
	return v;
}

/*Compute and return the normal vector to the surface
given:
Point P - the Point on the sphere
Sphere O -
*/
float *  normal( Point P, Sphere3D O) {
	float * result = (float*) malloc (sizeof(float)*3);
	result[0]= P.x - O.x ;
	result[1]= P.y - O.y ;
	result[2]= P.z - O.z ;

	return normalize(result);
}

/*Compute and return the "Va" vector from the point to the camera ( COP )
using the value of the global variable COP.
given -
Point P
*/
float *  viewVector( Point P) {
	float * result = (float*) malloc (sizeof(float)*3);
	result[0]= COP.x - P.x  ;
	result[1]= COP.y - P.y  ;
	result[2]= COP.z - P.z  ;

	return normalize(result);
}


/*Compute and return the "La" vector from the point P to the light
given -
Point P
Point light
*/
float *  lightVector( Point P, Light light) {
	float * result = (float *) malloc (sizeof(float)*3);
	result[0]= light.x - P.x  ;
	result[1]= light.y - P.y  ;
	result[2]= light.z - P.z  ;

	return normalize(result);
}


/*Compute and return the "Ra" vector of the refraction
given -
Vector Na
Vector La
*/
float *  reflectVector( Vector Na, Vector La) {
	float * result = (float *)malloc (sizeof(float)*3);
	float scalarProduct = Na[0]*La[0]+Na[1]*La[1]+Na[2]*La[2];

	result[0]= 2*scalarProduct*Na[0]-La[0] ;
	result[1]= 2*scalarProduct*Na[1]-La[1] ;
	result[2]= 2*scalarProduct*Na[2]-La[2] ;

	return normalize(result);
}


/*compute the refraction vector
given :
Point p - the point oin the sphere
Vector N - the normal vector to the surface
Vector L - the light vector (from p to the light source)
float ni - the refraction index of the first material
float nj - the refraction index of the second material
*/
float* refractVector(Point p, Vector N, Vector L, float ni, float nj ){
	float *result = (float*) malloc(sizeof(float) * 3 );
	float n = ni/ nj;
	float c =  - dot(N,L);

	int i;
	for(i=0; i<3 ; i++){
		result[i] = n*L[i] + N[i]*( n*c- sqrt(1+ n*n * ((c*c) -1)));
	}

	return result;
}

/****************************************************************************************/

//scalar product
float  dot( Vector Na, Vector La) {
	return Na[0]*La[0]+Na[1]*La[1]+Na[2]*La[2];
}

/* return -v
given, Vector v
*/
float * invert(Vector v){
	float * result = (float *)malloc (sizeof(float)*3);
	result[0] = - v[0];
	result[1] = - v[1];
	result[2] = - v[2];

	return result;
}

/*****************************************************************************************/
/*Compute the reflection intensity on an object from a light source and add it to the existing colour
given
Vector Na - the normal vector to the surface
Point p - the point on the sphere
Colour Ia - ambiant intensity
Colour Ka - rgb values for the Ka coeficient
Colour Kd -	rgb values for the Kd coeficient
Colour Ks -	rgb values for the Ks coeficient
float n - the specular intensity
 Colour init - the existing colour of the pixel
*/
Colour reflectionIntensity( Vector Na, Point p, Colour Ia, Sphere3D ball,Vector Va){

	Colour result; //result will be the result of our calculation
	float *La, *Ra; //the computed value of the light and reflection vectors

	result.r = Ia.r*ball.ka.r;
	result.g = Ia.g*ball.ka.g;
	result.b = Ia.b*ball.ka.b;



	int k;
	for(k=0 ; k < NLIGHT; k++){
		La = lightVector( p, lights[k]);
		if( !intersectionTest(p,La , ball) && dot(Na,La) >=0 ){
			Ra = reflectVector(Na, La);
			result.r = result.r + lights[k].r*((ball.kd.r*(dot(Na,La))));
			result.g = result.g + lights[k].g*((ball.kd.g*(dot(Na,La))));
			result.b = result.b + lights[k].b*((ball.kd.b*(dot(Na,La))));

			if(dot(Va,Ra)>0){
				result.r = result.r +lights[k].r* (ball.ks.r*(pow(dot(Va,Ra),ball.n)));
				result.g = result.g +lights[k].g* (ball.ks.g*(pow(dot(Va,Ra),ball.n)));
				result.b = result.b +lights[k].b* (ball.ks.b*(pow(dot(Va,Ra),ball.n)));
			}
		}
	}

	return  result;
}

/*****************************************************************************************/


/*RayTracing recursive algorithm - involves following the path of rays from the COP - centre of camera projection through the centre of
scene pixels.  Follow each of these rays to find that if it intersects an object to represent the colour of the
corresponding scene pixel on the viewing plane.  At the end of this process scene pixles will form an image on the
screen
given :
Point p
Vector direction
*/
void rayTrace(Point toto , Vector direction, int i , int j){
	count++;//increment the number of recursive call for rayTrace

	//dx, dy, and dz are the values of the directional vector of the ray, (x0, y0, z0)
	float dx, dy,dz;

	//A, B, and C of the quadratic equation
	float A, B, C;

	//Point on ray, or point of intersection
	float t;

	//Used for the co-ordinates of the intersection point
	float x, y, z;

	//distance between intersection point and toto
	float distance;

	int k;
	double depth = 1000;//To compare point of first intersection with this large value

	//Intersetion test carried out according to the number of spheres or objects
	for ( k = 0; k < NSPHERE; k++){
		//A, B, and C of the quadratic equation
		A = direction[0] * direction[0] + direction[1] * direction[1] + direction[2] * direction[2];

		B = (toto.x - sphere[k].x) * direction[0] + (toto.y - sphere[k].y) * direction[1] + (toto.z - sphere[k].z) * direction[2];

		C = toto.x * toto.x + toto.y * toto.y + toto.z * toto.z- 2 * (toto.x * sphere[k].x + toto.y * sphere[k].y + toto.z * sphere[k].z) +  sphere[k].x * sphere[k].x + sphere[k].y * sphere[k].y + sphere[k].z * sphere[k].z - sphere[k].radius * sphere[k].radius;

		//Ray tangent to sphere and intersecting in two places to form image of spheres on viweing plane
		if (B * B - A * C >= 0){

			//Determining the point of intersection
			t = - B - sqrt(B * B - A * C) / A;

			//Determining the co-ordinates of the intersection point
			x = toto.x + t * direction[0];
			y = toto.y + t * direction[1];
			z = toto.z + t * direction[2];
			Point p = {x,y,z};

			//Using pythagoras to work out distance between intersection point and toto
			distance = sqrt((toto.x - x) * (toto.x - x) + (toto.y - y) * (toto.y - y) + (toto.z - z )*( toto.z - z));

			//Comparing distance with depth enabling pixels nearest to toto to form nearest image of sphere
			if (distance < depth){

				float * Na = normal( p , sphere[k] );
				float * Va = viewVector ( p);



				pixels[i][j].r += reflectionIntensity( Na, p, Ia, sphere[k], Va).r;
				pixels[i][j].g += reflectionIntensity( Na, p, Ia, sphere[k], Va).g;
				pixels[i][j].b += reflectionIntensity( Na, p, Ia, sphere[k], Va).b;

				float * next = reflectVector(Na,invert(direction));
				if( intersectionTest(p, next, sphere[k]) ){

					rayTrace(p, next, i,j);
				}
				depth = distance ;
			}
		}
	}
}


/***********************************************************************************************/

/* Test the intersection between a segment (from a point and a light source) and each sphere in the scene and return the parameter t
(when t<0 the sphere isn't between the points, so there isn't any shadow)
given -
Point p1 - the point on the sphere "ball"
Vector La - the Light vector
*/
int intersectionTest (Point p, Vector v, Sphere3D ball){

   	int shadow = 0;

	// A, B, and C of the quadratic equation
	float A, B, C;

	// the parameter to find the intersection point
	float t;

	int k;
	for (k = 0; k < NSPHERE; k++) {
	   if (sphere[k].x != ball.x && sphere[k].y != ball.y && sphere[k].z != ball.z) {

			A = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];

			B = (p.x - sphere[k].x) * v[0] + (p.y - sphere[k].y) * v[1] + (p.z - sphere[k].z) * v[2];

			C = p.x * p.x + p.y * p.y + p.z * p.z - 2 * (p.x * sphere[k].x + p.y * sphere[k].y + p.z * sphere[k].z) + sphere[k].x * sphere[k].x + sphere[k].y * sphere[k].y + sphere[k].z * sphere[k].z - sphere[k].radius * sphere[k].radius;


			if (B * B - A * C >= 0) {
				shadow = 1;

				float t = - B - sqrt(B * B - A * C) / A; // find the intersection point

				if (t < 0) {
					shadow = 0;
 				}

			}
		}
	}


	return shadow;
}
