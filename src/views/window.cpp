#include "window.h"

window::window(int argc, char **argv)
{
    initGlutWindow(argc, argv);
}

window::~window()
{
    //dtor
}

void window::initGlutWindow(int argc, char **argv)
{

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
}

void window::glutMainLoop()
{
    glutMainLoop();
}

//Display function creating a graphics window
void window::display()
{

	int r, c;

	// clear the offscreen buffer
	glClearColor(0.0, 0.0, 0.0, 0.0);/*Background colour*/
	glClear(GL_COLOR_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	//Displaying scene pixels and forming an image
	glBegin(GL_POINTS);

    ColourGrid pixels = raytracingRenderer::getPixels();

	for (r = 0; r < 512; r++)
		for (c = 0; c < 512; c++)
		{
			glColor3f(pixels->data[r][c].r, pixels->data[r][c].g, pixels->data[r][c].b);
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
