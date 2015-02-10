/////////////////////////////////////////////////////////////
// FILE:      mygl.cpp
// CONTAINS:  your implementations of various GL functions
////////////////////////////////////////////////////////////.

#define _USE_MATH_DEFINES

#include <cstdlib>
#include <cfloat>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

#ifdef __APPLE__
#  include <GLUT/glut.h>
#elif _WIN32
#  include "GL/glut.h"
#else
#  include <GL/glut.h>
#endif

#include "image.hpp"
#include "mygl.hpp"

// Flags that are defined/set in a3.cpp.
extern bool perspectiveCorrectTextures;
extern bool gridIsVisible;
extern bool drawPolysAsPoints;

// The dimensions of the virtual window we are drawing into.
int virtualWidth;
int virtualHeight;

// The current projection and model-view matrices.
Matrix projectionMatrix;
Matrix modelViewMatrix;

// Current color, texture, and texture coordinates.
Vector currentColor;
Image *currentTexture;
double currentTextureCoord[2];
Vector currentNormal;
int currentPolygonMode;
int currentVertexMode;

struct Vertex {
	Vector position, texture_coord, normal, color;
	Vertex(Vector position, Vector texture_coord, Vector normal, Vector color) : position(position), texture_coord(texture_coord), normal(normal), color(color) {}
};
std::vector<Vertex> currentVertices;

// A class to simplify lookup of two-dimensional zBuffer data from an array
// of doubles. This zBuffer MUST have reshape(...) called before use.
class ZBuffer {
private:

    int height_, size_, allocated_;
    double *data_;

public:

    ZBuffer() : allocated_(0), data_(NULL) {}

    // Reshape the zBuffer because the window was reshaped.
    void reshape(int w, int h) {
        height_ = h;
        size_ = w * h;
        // Only need to bother reallocating the array when we need a larger
        // array.
        if (size_ > allocated_) {
            delete [] data_;
            allocated_ = 4 * size_;
            data_ = new double[allocated_];
        }
    }

    // Clear out all of the depth values.
    void clear() {
        for (int i = 0; i < size_; i++) {
            data_[i] = DBL_MAX;
        }
    }

    // Get a pointer to the given row of depth data, which can be indexed
    // again to get/set depth values.
    // E.g.: `zBuffer[x][y] = newDepthValue`.
    double* operator[](int index) {
        return data_ + index * height_;
    }

} zBuffer;


// A function to set a pixel value on the screen. This is the entry point that
// you MUST use to draw to the screen.
void setPixel(int x, int y, double r, double g, double b)
{
    if (x < 0 || x >= virtualWidth || y < 0 || y >= virtualHeight) {
        std::cerr << "attempting to set a pixel that is off-screen;" << x << ", " << y << std::endl;
        return;
    }
    glColor3d(r, g, b);
    glBegin(GL_POINTS);
    glVertex2f(x + 0.5, y + 0.5);
    glEnd();
}


// Draws the virtual pixel grid.
void drawPixelGrid() {

    // Dark gray.
    glColor4d(0.15, 0.15, 0.15, 1.0);

    glBegin(GL_LINES);
	  
    // Draw vertical grid lines.
    for (int x = 0; x <= virtualWidth; x++) {
        glVertex3f(x, 0.0, 1.0);
        glVertex3f(x, virtualHeight, 1.0);
    }

    // Draw horizontal grid lines.
    for (int y = 0; y <= virtualHeight; y++) {
        glVertex3f(0.0, y, 1.0);
        glVertex3f(virtualWidth, y, 1.0);
    }
    
    glEnd();
}


void myViewport(int w, int h) {
    virtualWidth = w;
    virtualHeight = h;
    zBuffer.reshape(w, h);
}


void myClear() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    if (gridIsVisible) {
        drawPixelGrid();
    }
    zBuffer.clear();
}


void myLoadIdentity() {
    modelViewMatrix = Matrix::identity();
    projectionMatrix = Matrix::identity();
}

void myOrtho() {
    projectionMatrix *= Matrix(
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0,-1, 0, // switches Z axis
        0, 0, 0, 1
    );	
}

void myBindTexture(char const *name) {
    if (name) {
        currentTexture = &Image::fromFile(name);
    } else {
        currentTexture = NULL;
    }
}

void myPolygonMode(int mode) {
    currentPolygonMode = mode;
}


void myBegin(int mode) {
	currentVertexMode = mode;
	currentVertices.clear();
}


void myColor(double r, double g, double b) {
    currentColor = Vector(r, g, b, 1.0);
}


void myVertex(double x, double y, double z) {
	Vector position = Vector(x, y, z, 1.0);
	Vector texture_coord = Vector(currentTextureCoord[0], currentTextureCoord[1], 0.0, 0.0);
	currentVertices.push_back(Vertex(position, texture_coord, currentNormal, currentColor));
}


void drawPoint(Vector position, Vector color) {
	if (abs(position[0]) <= 1.0 && abs(position[1]) <= 1.0) {
		int x = position[0] * virtualWidth/2.0 + virtualWidth/2.0,
			y = position[1] * virtualHeight/2.0 + virtualHeight/2.0;
		if (position[2] < zBuffer[x][y]) {
			setPixel(x, y, color[0], color[1], color[2]);
		}
	}
}

void drawTriangle(Vector position0, Vector color0, Vector position1, Vector color1, Vector position2, Vector color2) {
	
	int x0 = position0[0] * virtualWidth/2.0 + virtualWidth/2.0,
		y0 = position0[1] * virtualHeight/2.0 + virtualHeight/2.0;
	int x1 = position1[0] * virtualWidth/2.0 + virtualWidth/2.0,
		y1 = position1[1] * virtualHeight/2.0 + virtualHeight/2.0;
	int x2 = position2[0] * virtualWidth/2.0 + virtualWidth/2.0,
		y2 = position2[1] * virtualHeight/2.0 + virtualHeight/2.0;

	int max_x = (x0 > x1) ? ((x0 > x2) ? x0 : x2) : ((x1 > x2) ? x1 : x2);
		max_x = (max_x > virtualWidth) ? virtualWidth : max_x;
	int min_x = (x0 < x1) ? ((x0 < x2) ? x0 : x2) : ((x1 < x2) ? x1 : x2);
		min_x = (min_x < 0) ? 0 : min_x;
	int max_y = (y0 > y1) ? ((y0 > y2) ? y0 : y2) : ((y1 > y2) ? y1 : y2);
		max_y = (max_y > virtualHeight) ? virtualHeight : max_y;
	int min_y = (y0 < y1) ? ((y0 < y2) ? y0 : y2) : ((y1 < y2) ? y1 : y2);
		min_y = (min_y < 0) ? 0 : min_y;

	int dx0 = (x1 - x0),
		dx1 = (x2 - x1),
		dx2 = (x0 - x2);
	int dy0 = (y1 - y0),
		dy1 = (y2 - y1),
		dy2 = (y0 - y2);

//	double A = Vector((x1 - x0), (y1 - x0), 0.0, 1.0).cross(Vector((x2 - x0), (y2 - y0), 0.0, 1.0)).length();
//	std::cout << "A: " << A << "\n";

	for (int i = min_x; i <= max_x; i++) {
		int l0 = dy0 * (i - x0) - dx0 * (min_y - y0),
			l1 = dy1 * (i - x1) - dx1 * (min_y - y1),
			l2 = dy2 * (i - x2) - dx2 * (min_y - y2);
		for (int j = min_y; j <= max_y; j++) {
			if ((l0 < 0 && l1 < 0 && l2 < 0) || (l0 >= 0 && l1 >= 0 && l2 >= 0)) {
				setPixel(i, j, color0[0], color0[1], color0[2]);
//				double a0 = Vector((x1 - i), (y1 - j), 0.0, 1.0).cross(Vector((x2 - i), (y2 - j), 0.0, 1.0)).length() / A;
//				double a1 = Vector((x2 - i), (y2 - j), 0.0, 1.0).cross(Vector((x0 - i), (y0 - j), 0.0, 1.0)).length() / A;
//				double a2 = Vector((x0 - i), (y0 - j), 0.0, 1.0).cross(Vector((x1 - i), (y1 - j), 0.0, 1.0)).length() / A;
////				std::cout << "As: " << a0 << "," << a1 << "," << a2;
//
//				double z = a0 * position0[2] + a1 * position1[2] + a2 * position2[2];
////				std::cout << "\tz: " << z << "\n";
//				if (z < zBuffer[i][j]) {
//					setPixel(i, j,	a0 * color0[0] + a1 * color1[0] + a2 * color2[0],
//									a0 * color0[1] + a1 * color1[1] + a2 * color2[1],
//									a0 * color0[2] + a1 * color1[2] + a2 * color2[2]);
//					zBuffer[i][j] = z;
//				}
			}
			l0 -= dx0;
			l1 -= dx1;
			l2 -= dx2;
		}
	}
}

void myEnd() {
	switch (currentVertexMode) {
	case GL_POINT:
		for (int i = 0; i < currentVertices.size(); i++) {
			Vertex currentVertex = currentVertices[i];
			Vector position = projectionMatrix * modelViewMatrix * currentVertex.position;
			position = position/position[3];
			drawPoint(position, currentVertex.color);
		}
		break;

	case GL_LINE:
		for (int i = 0; i < currentVertices.size(); i += 2) {

			Vertex currentVertex0 = currentVertices[i],
				currentVertex1 = currentVertices[i+1];
			Vector position0 = projectionMatrix * modelViewMatrix * currentVertex0.position,
				position1 = projectionMatrix * modelViewMatrix * currentVertex1.position;
			position0 = position0/position0[3],
				position1 = position1/position1[3];

			switch (currentPolygonMode) {
			case GL_POINT:
				drawPoint(position0, currentVertex0.color);
				drawPoint(position1, currentVertex1.color);
				break;

			case GL_LINE:
			case GL_FILL:
				;
				//drawLine(p0, p1);
			}
		}
		break;

	case GL_TRIANGLES:
		for (int i = 0; i < currentVertices.size(); i += 3) {

			Vertex currentVertex0 = currentVertices[i],
				currentVertex1 = currentVertices[i+1],
				currentVertex2 = currentVertices[i+2];
			Vector position0 = projectionMatrix * modelViewMatrix * currentVertex0.position,
				position1 = projectionMatrix * modelViewMatrix * currentVertex1.position,
				position2 = projectionMatrix * modelViewMatrix * currentVertex2.position;
			position0 = position0/position0[3],
				position1 = position1/position1[3],
				position2 = position2/position2[3];

			switch (currentPolygonMode) {
			case GL_POINT:
				drawPoint(position0, currentVertex0.color);
				drawPoint(position1, currentVertex1.color);
				drawPoint(position2, currentVertex2.color);
				break;

			case GL_LINE:
				//drawLine(p0, p1);
				//drawLine(p1, p2);
				//drawLine(p2, p0);
				//break;

			case GL_FILL:
				drawTriangle(position0, currentVertex0.color, position1, currentVertex1.color, position2, currentVertex2.color);
			}
		}
		break;

	case GL_POLYGON:

	default:
		std::cerr << "unsupported vertex mode";
	}
}


void myTranslate(double tx, double ty, double tz) {

	modelViewMatrix *= Matrix(	1.0, 0.0, 0.0, tx,
								0.0, 1.0, 0.0, ty,
								0.0, 0.0, 1.0, tz,
								0.0, 0.0, 0.0, 1.0	);
}


void myRotate(double angle, double axisX, double axisY, double axisZ) {

	modelViewMatrix *= Matrix::rotation(angle * M_PI/180.0, Vector(axisX, axisY, axisZ));
}


void myScale(double sx, double sy, double sz) {

	modelViewMatrix *= Matrix(	sx, 0.0, 0.0, 0.0,
								0.0, sy, 0.0, 0.0,
								0.0, 0.0, sz, 0.0,
								0.0, 0.0, 0.0, 1.0	);
}


void myFrustum(double left, double right, double bottom, double top, double mynear, double myfar) {

	projectionMatrix = Matrix(	2.0*mynear/(right-left), 0.0, (right+left)/(right-left), 0.0,
								0.0, 2.0*mynear/(top-bottom), (top+bottom)/(top-bottom), 0.0,
								0.0, 0.0, -(myfar+mynear)/(myfar-mynear), -2.0*myfar*mynear/(myfar-mynear),
								0.0, 0.0, -1.0, 0.0	);
}


void myLookAt(double eyeX, double eyeY, double eyeZ,
              double cenX, double cenY, double cenZ,
              double  upX, double  upY, double  upZ) {

	Vector eye = Vector(eyeX, eyeY, eyeZ);
	Vector cen = Vector(cenX, cenY, cenZ);
	Vector up = Vector(upX, upY, upZ);

	Vector w = (eye - cen).normalized();
	Vector u = up.cross(w).normalized();
	Vector v = w.cross(u);

	modelViewMatrix *= Matrix(	u[0], u[1], u[2], (u * -1.0).dot(eye),
								v[0], v[1], v[2], (v * -1.0).dot(eye),
								w[0], w[1], w[2], (w * -1.0).dot(eye),
								0.0, 0.0, 0.0, 1.0	);
}


void myTexCoord(double s, double t) {
    // @@@@ YOUR CODE HERE
}


void myNormal(double x, double y, double z) {
    // @@@@ YOUR CODE HERE
}


void myEnableLighting() {
    // @@@@ YOUR CODE HERE
}


void myDisableLighting() {
    // @@@@ YOUR CODE HERE
}
