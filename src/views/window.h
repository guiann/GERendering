

#ifndef WINDOW_H
#define WINDOW_H

#include <GL/glut.h>
#include <stdio.h>
#include <stdlib.h>

#include "../renderer/raytracingRenderer.h"

class window
{
    public:
        window(int argc, char **argv);
        virtual ~window();
        void glutMainLoop();
    protected:
    private:


        //GLUT init function
        void initGlutWindow(int argc, char **argv);

        //GLUT Display function
        static void display();
};

#endif // WINDOW_H
