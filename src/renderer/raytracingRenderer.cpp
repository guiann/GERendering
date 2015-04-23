#include "raytracingRenderer.h"

//Raycasting - involves following the path of rays from the COP - center of camera projection through the center of
//scene pixelsImgMap.  Follow each of these rays to find that if it intersects an object to represent the color of the
//corresponding scene pixel on the viewing plane.  At the end of this process scene pixelsImgMap will form an image on the
//screen



// Global variables____________________________________________________________

//Pixel array for our viewing plane
ColourGrid pixelsImgMap;

//array to hold three elements type sphere
Sphere3D sphere [3];
Light lights[2];

//Centre of camera projection co-ordinates
Point COP = {0, 0, 3 };

//define an ambiant light for the scene
Colour Ia = {0.3, 0.3, 0.3};



//count the number of recursive call for rayTrace
int count;


raytracingRenderer::raytracingRenderer()
{

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
}

raytracingRenderer::~raytracingRenderer()
{
    //dtor
}

ColourGrid raytracingRenderer::getPixels(){

     ColourGrid grid = (ColourGrid) malloc(sizeof(struct _colourGrid_t));

    Colour **data = (Colour **) malloc(sizeof(Colour *) * 512);

    for (int i = 0; i < 512; i++) {
        data[i] = (Colour*) malloc(sizeof(Colour) * 512);

        memset(data[i], 0, sizeof(int) * 512);
    }

    grid->data = data;
    grid->width = 512;
    grid->height = 512;

    return grid;
}
void raytracingRenderer::render()
{
	int i,j;
	float dx,dy,dz;
	Vector direction;


	//Determing the width and height of each pixel
	//width and height of each scene pixel
	float pixelwidth, pixelheight;

	//Centre of type point
	Point centreofpix;

	pixelwidth = (XMAX - XMIN) / 512;
	pixelheight = (YMAX - YMIN) / 512;



	//Incrementing through our scene
	for ( i = 0; i < 512; i++){
		for ( j = 0; j < 512; j++){

			pixelsImgMap->data[i][j].r= 0.0;
			pixelsImgMap->data[i][j].g= 0.0;
			pixelsImgMap->data[i][j].b= 0.0;

			//Determining the centre of our pixelsImgMap
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

			pixelsImgMap->data[i][j].r = pixelsImgMap->data[i][j].r/count;
			pixelsImgMap->data[i][j].g = pixelsImgMap->data[i][j].g/count;
			pixelsImgMap->data[i][j].b = pixelsImgMap->data[i][j].b/count;

			if(pixelsImgMap->data[i][j].r < 0) pixelsImgMap->data[i][j].r=0;
			if(pixelsImgMap->data[i][j].g < 0) pixelsImgMap->data[i][j].g=0;
			if(pixelsImgMap->data[i][j].b < 0) pixelsImgMap->data[i][j].b=0;

			if(pixelsImgMap->data[i][j].r > 1) pixelsImgMap->data[i][j].r=1;
			if(pixelsImgMap->data[i][j].g > 1) pixelsImgMap->data[i][j].g=1;
			if(pixelsImgMap->data[i][j].b > 1) pixelsImgMap->data[i][j].b=1;
		}
	}
}


/*we'll need to normalize several vectors*/
float * raytracingRenderer::normalize( Vector v){
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
float *  raytracingRenderer::normal( Point P, Sphere3D O) {
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
float *  raytracingRenderer::viewVector( Point P) {
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
float *  raytracingRenderer::lightVector( Point P, Light light) {
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
float *  raytracingRenderer::reflectVector( Vector Na, Vector La) {
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
float* raytracingRenderer::refractVector(Point p, Vector N, Vector L, float ni, float nj ){
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
float  raytracingRenderer::dot( Vector Na, Vector La) {
	return Na[0]*La[0]+Na[1]*La[1]+Na[2]*La[2];
}

/* return -v
given, Vector v
*/
float * raytracingRenderer::invert(Vector v){
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
Colour raytracingRenderer::reflectionIntensity( Vector Na, Point p, Colour Ia, Sphere3D ball,Vector Va){

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
scene pixelsImgMap.  Follow each of these rays to find that if it intersects an object to represent the colour of the
corresponding scene pixel on the viewing plane.  At the end of this process scene pixles will form an image on the
screen
given :
Point p
Vector direction
*/
void raytracingRenderer::rayTrace(Point toto , Vector direction, int i , int j){
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

			//Comparing distance with depth enabling pixelsImgMap nearest to toto to form nearest image of sphere
			if (distance < depth){

				float * Na = normal( p , sphere[k] );
				float * Va = viewVector ( p);



				pixelsImgMap->data[i][j].r += reflectionIntensity( Na, p, Ia, sphere[k], Va).r;
				pixelsImgMap->data[i][j].g += reflectionIntensity( Na, p, Ia, sphere[k], Va).g;
				pixelsImgMap->data[i][j].b += reflectionIntensity( Na, p, Ia, sphere[k], Va).b;

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
int raytracingRenderer::intersectionTest (Point p, Vector v, Sphere3D ball){

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
