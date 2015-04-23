#include "views/window.h"
#include "renderer/raytracingRenderer.h"



int main(int argc, char **argv){


    window* w = new window(argc, argv);

    raytracingRenderer* renderer = new raytracingRenderer();
    renderer->render();

    w->glutMainLoop();

    return 0;
}



