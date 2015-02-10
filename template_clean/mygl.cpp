/////////////////////////////////////////////////////////////
// FILE:      mygl.cpp
// CONTAINS:  your implementations of various GL functions
////////////////////////////////////////////////////////////.

#define _USE_MATH_DEFINES
#define MAX_LIGHTS 4

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

bool enableLighting = false;
Vector currentEye;

struct Vertex {
	Vector position, texture_coord, normal, color;
	Vertex(Vector position, Vector texture_coord, Vector normal, Vector color) : position(position), texture_coord(texture_coord), normal(normal), color(color) {}
};
std::vector<Vertex> currentVertices;

struct Light {
	bool positionEnabled, ambientEnabled, diffuseEnabled, specularEnabled;
	Vector position, ambient, diffuse, specular;
	
	Light() {
		positionEnabled = ambientEnabled = diffuseEnabled = specularEnabled = false;
		position = Vector(0.0, 0.0, 0.0, 0.0);
	}
};
Light currentLights[MAX_LIGHTS];
bool lightEnabled[MAX_LIGHTS] = {false};

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


Vector computeLighting(Vector normal) {
	double r_intensity = 0.0,
		g_intensity = 0.0,
		b_intensity = 0.0;
	for (int i = 0; i < MAX_LIGHTS; i++) {
		if (lightEnabled[i]) {
			Light light = currentLights[i];
			if (light.ambientEnabled) {
				r_intensity += light.ambient[0];
				g_intensity += light.ambient[1];
				b_intensity += light.ambient[2];
			}
			if (light.diffuseEnabled) {
				double diffuse_intensity = normal.dot(light.position);
				r_intensity += light.diffuse[0] * diffuse_intensity;
				g_intensity += light.diffuse[1] * diffuse_intensity;
				b_intensity += light.diffuse[2] * diffuse_intensity;
			}
			if (light.specularEnabled) {
				Vector reflection = normal * 2 * (light.position.dot(normal));
				double specular_intensity = currentEye.dot(reflection.normalized());
				if (specular_intensity > 0) {
					specular_intensity = pow(specular_intensity, 4);
					r_intensity += light.specular[0] * specular_intensity;
					g_intensity += light.specular[1] * specular_intensity;
					b_intensity += light.specular[2] * specular_intensity;
				}
			}
		}
	}
	return Vector(r_intensity, g_intensity, b_intensity, 1.0);
}

void drawPoint(Vertex vertex) {
	Vector position = projectionMatrix * modelViewMatrix * vertex.position;
	position = position/position[3];
	Vector color = vertex.color;
	Vector normal = vertex.normal;

	if (abs(position[0]) <= 1.0 && abs(position[1]) <= 1.0) {
		int x = position[0] * virtualWidth/2.0 + virtualWidth/2.0,
			y = position[1] * virtualHeight/2.0 + virtualHeight/2.0;
		if (position[2] < zBuffer[x][y]) {
			if (enableLighting) {
				Vector lighting = computeLighting(normal);
				setPixel(x, y, color[0] * lighting[0], color[1] * lighting[0], color[2] * lighting[0]);
			}
			else
				setPixel(x, y, color[0], color[1], color[2]);
		}
	}
}

void drawLine(Vertex vertex0, Vertex vertex1) {

	Vector position0 = projectionMatrix * modelViewMatrix * vertex0.position,
		position1 = projectionMatrix * modelViewMatrix * vertex1.position;
	position0 = position0/position0[3],
		position1 = position1/position1[3];

	int clip0 = (position0[0] >= 0.0) | ((position0[1] >= 0.0) << 1) | ((position0[0] < 1.0) << 2) | ((position0[1] < 1.0) << 3);
	int clip1 = (position1[0] >= 0.0) | ((position1[1] >= 0.0) << 1) | ((position1[0] < 1.0) << 2) | ((position1[1] < 1.0) << 3);
	if (clip0 & clip1)
		return;

	int x0 = position0[0] * virtualWidth/2.0 + virtualWidth/2.0,
		y0 = position0[1] * virtualHeight/2.0 + virtualHeight/2.0;
	int x1 = position1[0] * virtualWidth/2.0 + virtualWidth/2.0,
		y1 = position1[1] * virtualHeight/2.0 + virtualHeight/2.0;
	int x = x0,
		y = y0;
	int	dx = x1 - x0,
		dy = y1 - y0;

	// TODO: implement clipping
	// account for the points being out of order

	Vector color0 = vertex0.color,
		color1 = vertex1.color;
	Vector normal0 = vertex0.normal,
		normal1 = vertex1.normal;
	double z0 = position0[2],
		z1 = position1[2];
	Vector color = color0;
	Vector normal = normal0;
	double z = z0;

	if (enableLighting) {
		if (z < zBuffer[x][y]) {
			Vector lighting = computeLighting(normal);
			setPixel(x, y, color[0] * lighting[0], color[1] * lighting[0], color[2] * lighting[0]);
			zBuffer[x][y] = z;
		}

		if (abs(dx) >= abs(dy)) {
			int d = 2 * dy - dx;
			int d_e = 2 * dy,
				d_ne = 2 * (dy - dx);
			Vector dcolor = (color1 - color0) / dx;
			Vector dnormal = (normal1 - normal0) / dx;
			double dz = (z1 - z0) / dx;

			if (dy >= 0) {
				while (x++ <= x1) {
					if (d < 0)
						d += d_e;
					else {
						d += d_ne;
						y++;
					}
					color = color + dcolor;
					normal = normal + dnormal;
					z += dz;

					if (z < zBuffer[x][y]) {
						Vector lighting = computeLighting(normal);
						setPixel(x, y, color[0] * lighting[0], color[1] * lighting[0], color[2] * lighting[0]);
						zBuffer[x][y] = z;
					}
				}
			} else {	// dy < 0
				while (x++ <= x1) {
					if (d < 0)
						d += d_e;
					else {
						d += d_ne;
						y--;
					}
					color = color + dcolor;
					normal = normal + dnormal;
					z += dz;

					if (z < zBuffer[x][y]) {
						Vector lighting = computeLighting(normal);
						setPixel(x, y, color[0] * lighting[0], color[1] * lighting[0], color[2] * lighting[0]);
						zBuffer[x][y] = z;
					}
				}
			}
		} else {	// dy > dx
			int d = 2 * dx - dy;
			int d_e = 2 * dx,
				d_ne = 2 * (dx - dy);
			Vector dcolor = (color1 - color0) / dy;
			Vector dnormal = (normal1 - normal0) / dy;
			double dz = (z1 - z0) / dy;

			if (dx >= 0) {
				while (y++ <= y1) {
					if (d < 0)
						d += d_e;
					else {
						d += d_ne;
						x++;
					}
					color = color + dcolor;
					normal = normal + dnormal;
					z += dz;
					
					if (z < zBuffer[x][y]) {
						Vector lighting = computeLighting(normal);
						setPixel(x, y, color[0] * lighting[0], color[1] * lighting[0], color[2] * lighting[0]);
						zBuffer[x][y] = z;
					}
				}
			} else {	// dx < 0
				while (y++ <= y1) {
					if (d < 0)
						d += d_e;
					else {
						d += d_ne;
						x--;
					}
					color = color + dcolor;
					normal = normal + dnormal;
					z += dz;
					
					if (z < zBuffer[x][y]) {
						Vector lighting = computeLighting(normal);
						setPixel(x, y, color[0] * lighting[0], color[1] * lighting[0], color[2] * lighting[0]);
						zBuffer[x][y] = z;
					}
				}
			}
		}
	} else {	// no lighting
		if (z < zBuffer[x][y]) {
			setPixel(x, y, color[0], color[1], color[2]);
			zBuffer[x][y] = z;
		}
	
		if (dx >= dy) {
			int d = 2 * dy - dx;
			int d_e = 2 * dy,
				d_ne = 2 * (dy - dx);
			Vector dcolor = (color1 - color0) / dx;
			double dz = (z1 - z0) / dx;

			if (dy >= 0) {
				while (x++ <= x1) {
					if (d < 0)
						d += d_e;
					else {
						d += d_ne;
						y++;
					}
					color = color + dcolor;
					z += dz;

					if (z < zBuffer[x][y]) {
						setPixel(x, y, color[0], color[1], color[2]);
						zBuffer[x][y] = z;
					}
				}
			} else {	// dy < 0
				while (x++ <= x1) {
					if (d < 0)
						d += d_e;
					else {
						d += d_ne;
						y--;
					}
					color = color + dcolor;
					z += dz;

					if (z < zBuffer[x][y]) {
						setPixel(x, y, color[0], color[1], color[2]);
						zBuffer[x][y] = z;
					}
				}
			}
		} else {	// dy > dx
			int d = 2 * dx - dy;
			int d_e = 2 * dx,
				d_ne = 2 * (dx - dy);
			Vector dcolor = (color1 - color0) / dy;
			double dz = (z1 - z0) / dy;

			if (dx >= 0) {
				while (y++ <= y1) {
					if (d < 0)
						d += d_e;
					else {
						d += d_ne;
						x++;
					}
					color = color + dcolor;
					z += dz;

					if (z < zBuffer[x][y]) {
						setPixel(x, y, color[0], color[1], color[2]);
						zBuffer[x][y] = z;
					}
				}
			} else {	// dx < 0
				while (y++ <= y1) {
					if (d < 0)
						d += d_e;
					else {
						d += d_ne;
						x--;
					}
					color = color + dcolor;
					z += dz;

					if (z < zBuffer[x][y]) {
						setPixel(x, y, color[0], color[1], color[2]);
						zBuffer[x][y] = z;
					}
				}
			}
		}
	}
}

void drawTriangle(Vertex vertex0, Vertex vertex1, Vertex vertex2) {

	Vector homogenous0 = modelViewMatrix * vertex0.position,
		homogenous1 = modelViewMatrix * vertex1.position,
		homogenous2 = modelViewMatrix * vertex2.position;	
	Vector position0 = projectionMatrix * homogenous0,
		position1 = projectionMatrix * homogenous1,
		position2 = projectionMatrix * homogenous2;
	position0 = position0/position0[3],
		position1 = position1/position1[3],
		position2 = position2/position2[3];

	int x0 = position0[0] * virtualWidth/2.0 + virtualWidth/2.0,
		y0 = position0[1] * virtualHeight/2.0 + virtualHeight/2.0;
	int x1 = position1[0] * virtualWidth/2.0 + virtualWidth/2.0,
		y1 = position1[1] * virtualHeight/2.0 + virtualHeight/2.0;
	int x2 = position2[0] * virtualWidth/2.0 + virtualWidth/2.0,
		y2 = position2[1] * virtualHeight/2.0 + virtualHeight/2.0;

	int max_x = (x0 > x1) ? ((x0 > x2) ? x0 : x2) : ((x1 > x2) ? x1 : x2);
		max_x = (max_x >= virtualWidth) ? (virtualWidth - 1) : max_x;
	int min_x = (x0 < x1) ? ((x0 < x2) ? x0 : x2) : ((x1 < x2) ? x1 : x2);
		min_x = (min_x < 0) ? 0 : min_x;
	int max_y = (y0 > y1) ? ((y0 > y2) ? y0 : y2) : ((y1 > y2) ? y1 : y2);
		max_y = (max_y >= virtualHeight) ? (virtualHeight - 1) : max_y;
	int min_y = (y0 < y1) ? ((y0 < y2) ? y0 : y2) : ((y1 < y2) ? y1 : y2);
		min_y = (min_y < 0) ? 0 : min_y;

	if ((min_x == max_x) || (min_y == max_y))
		return;

	Matrix A = Matrix(	x0, y0, 1.0, 1.0,
						x1, y1, 1.0, 1.0,
						x2, y2, 1.0, 1.0,
						0.0, 0.0, 0.0, 1.0	);
	Matrix A_inv = A.inverse();

	//Matrix A = Matrix(	Vector(homogenous0[0], homogenous1[0], homogenous2[0], 1.0),
	//					Vector(homogenous0[1], homogenous1[1], homogenous2[1], 1.0),
	//					Vector(homogenous0[2], homogenous1[2], homogenous2[2], 1.0),
	//					Vector(0.0, 0.0, 0.0, 1.0)	);
	//Matrix A_inv = A.inverse();
	//Vector a_w = A_inv * Vector(1.0, 1.0, 1.0, 1.0);
	//Vector a_z = A_inv * Vector(position0[2] * homogenous0[3], position1[2] * homogenous1[3], position2[2] * homogenous2[3], 1.0);

	//homogenous0.print();
	//std::cout << "\n";
	//homogenous1.print();
	//std::cout << "\n";
	//homogenous2.print();
	//std::cout << "\n";
	//modelViewMatrix.print(true);
	//std::cout << "<- modelView\n";
	//projectionMatrix.print(true);
	//std::cout << "<- projection\n";
	//A.print(true);
	//std::cout << "<- A\n";
	//A_inv.print(true);
	//std::cout << "<- A_inv\n";
	//a_w.print();
	//std::cout << "\na_z: ";
	//a_z.print();
	//std::cout << "\n: ";

	Vector a_z = A_inv * Vector(position0[2], position1[2], position2[2], 1.0);
	Vector a_r = A_inv * Vector(vertex0.color[0], vertex1.color[0], vertex2.color[0], 1.0);
	Vector a_g = A_inv * Vector(vertex0.color[1], vertex1.color[1], vertex2.color[1], 1.0);
	Vector a_b = A_inv * Vector(vertex0.color[2], vertex1.color[2], vertex2.color[2], 1.0);

	int dx0 = (x1 - x0),
		dx1 = (x2 - x1),
		dx2 = (x0 - x2);
	int dy0 = (y1 - y0),
		dy1 = (y2 - y1),
		dy2 = (y0 - y2);

	if (!enableLighting) {

		for (int i = min_x; i <= max_x; i++) {

			int l0 = dy0 * (i - x0) - dx0 * (min_y - y0),
				l1 = dy1 * (i - x1) - dx1 * (min_y - y1),
				l2 = dy2 * (i - x2) - dx2 * (min_y - y2);
			double z = a_z[0] * i + a_z[1] * min_y + a_z[2],
				r = a_r[0] * i + a_r[1] * min_y + a_r[2],
				g = a_g[0] * i + a_g[1] * min_y + a_g[2],
				b = a_b[0] * i + a_b[1] * min_y + a_b[2];
			//double w_inv = a_w[0] * i + a_w[1] * min_y + a_w[2];

			for (int j = min_y; j <= max_y; j++) {
				if ((l0 <= 0 && l1 <= 0 && l2 <= 0) || (l0 >= 0 && l1 >= 0 && l2 >= 0)) {
					if (z < zBuffer[i][j]) {
						setPixel(i, j, r, g, b);
						zBuffer[i][j] = z;

						//double w = 1.0 / w_inv;
						//setPixel(i, j, r*w, g*w, b*w);
						//zBuffer[i][j] = z;
						//std::cout << "w: " << w <<  ", x: " << i << ", y: " << j << ", z: " << z << "\n";
						//std::cout << "r: " << r << ", g: " << g << ", b: " << b << "\n";
					}
				}
				l0 -= dx0;
				l1 -= dx1;
				l2 -= dx2;
				z += a_z[1];
				r += a_r[1];
				g += a_g[1];
				b += a_b[1];

				//w_inv += a_w[1];
			}
		}
	}
	else {

		// lighting enabled; need to perform lighting computations
		Vector a_nx = A_inv * Vector(vertex0.normal[0], vertex1.normal[0], vertex2.normal[0], 1.0);
		Vector a_ny = A_inv * Vector(vertex0.normal[1], vertex1.normal[1], vertex2.normal[1], 1.0);
		Vector a_nz = A_inv * Vector(vertex0.normal[2], vertex1.normal[2], vertex2.normal[2], 1.0);

		for (int i = min_x; i <= max_x; i++) {

			int l0 = dy0 * (i - x0) - dx0 * (min_y - y0),
				l1 = dy1 * (i - x1) - dx1 * (min_y - y1),
				l2 = dy2 * (i - x2) - dx2 * (min_y - y2);
			double z = a_z[0] * i + a_z[1] * min_y + a_z[2],
				r = a_r[0] * i + a_r[1] * min_y + a_r[2],
				g = a_g[0] * i + a_g[1] * min_y + a_g[2],
				b = a_b[0] * i + a_b[1] * min_y + a_b[2],
				nx = a_nx[0] * i + a_nx[1] * min_y + a_nx[2],
				ny = a_ny[0] * i + a_ny[1] * min_y + a_ny[2],
				nz = a_nz[0] * i + a_nz[1] * min_y + a_nz[2];

			for (int j = min_y; j <= max_y; j++) {
				if ((l0 <= 0 && l1 <= 0 && l2 <= 0) || (l0 > 0 && l1 > 0 && l2 > 0)) {
					if (z < zBuffer[i][j]) {
						Vector lighting = computeLighting(Vector(nx, ny, nz, 1.0));
						setPixel(i, j, r * lighting[0], g * lighting[1], b * lighting[2]);
						zBuffer[i][j] = z;

//						std::cout << "normal: ";
//						Vector normal = Vector(nx, ny, nz, 1.0);
//						normal.print();
//						std::cout << "\n";
					}
				}
				l0 -= dx0;
				l1 -= dx1;
				l2 -= dx2;
				r += a_r[1];
				g += a_g[1];
				b += a_b[1];
				nx += a_nx[1];
				ny += a_ny[1];
				nz += a_ny[1];
				z += a_z[1];
			}
		}
	}
}

void myEnd() {
	switch (currentVertexMode) {
	case GL_POINT:
		for (int i = 0; i < currentVertices.size(); i++) {
			drawPoint(currentVertices[i]);
		}
		break;

	case GL_LINE:
		for (int i = 0; i < currentVertices.size(); i += 2) {
			Vertex currentVertex0 = currentVertices[i],
				currentVertex1 = currentVertices[i+1];

			switch (currentPolygonMode) {
			case GL_POINT:
				drawPoint(currentVertex0);
				drawPoint(currentVertex1);
				break;

			case GL_LINE:
			case GL_FILL:
				drawLine(currentVertex0, currentVertex1);
				break;
				
			default:
				std::cerr << "invalid polygon mode\n";
			}
		}
		break;

	case GL_TRIANGLES:
		for (int i = 0; i < currentVertices.size(); i += 3) {
			Vertex currentVertex0 = currentVertices[i],
				currentVertex1 = currentVertices[i+1],
				currentVertex2 = currentVertices[i+2];

			switch (currentPolygonMode) {
			case GL_POINT:
				drawPoint(currentVertex0);
				drawPoint(currentVertex1);
				drawPoint(currentVertex2);
				break;

			case GL_LINE:
				drawLine(currentVertex0, currentVertex1);
				drawLine(currentVertex1, currentVertex2);
				drawLine(currentVertex2, currentVertex0);
				break;

			case GL_FILL:
				drawTriangle(currentVertex0, currentVertex1, currentVertex2);
				break;

			default:
				std::cerr << "invalid polygon mode\n";
			}
		}
		break;

	case GL_POLYGON:
		switch (currentPolygonMode) {
		case GL_POINT:
			for (int i = 0; i < currentVertices.size(); i++) {
				drawPoint(currentVertices[i]);
			}
			break;

		case GL_LINE: {
			Vertex currentVertex0 = currentVertices[0];

			for (int i = 1; i < currentVertices.size(); i++) {
				Vertex currentVertex1 = currentVertices[i];
				drawLine(currentVertex0, currentVertex1);
				currentVertex0 = currentVertex1;
			}
			drawLine(currentVertex0, currentVertices[0]);
			break;
		}

		case GL_FILL: {
			Vertex currentVertex0 = currentVertices[0],
				currentVertex1 = currentVertices[1];

			for (int i = 2; i < currentVertices.size(); i++) {
				Vertex currentVertex2 = currentVertices[i];
				drawTriangle(currentVertex0, currentVertex1, currentVertex2);
				currentVertex1 = currentVertex2;
			}
			break;
		}

		default:
			std::cerr << "invalid polygon mode: " << currentPolygonMode << "\n";
		}
		break;

	default:
		std::cerr << "unsupported vertex mode: " << currentVertexMode << "\n";
		std::cerr << "GL_POLYGON: " << GL_POLYGON;
		std::cerr << "\tGL_TRIANGLES: " << GL_TRIANGLES;
		std::cerr << "\tGL_LINE: " << GL_LINE;
		std::cerr << "\tGL_POINT: " << GL_POINT << "\n";
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

	currentEye = eye.normalized();
	modelViewMatrix *= Matrix(	u[0], u[1], u[2], (u * -1.0).dot(eye),
								v[0], v[1], v[2], (v * -1.0).dot(eye),
								w[0], w[1], w[2], (w * -1.0).dot(eye),
								0.0, 0.0, 0.0, 1.0	);
}


void myTexCoord(double s, double t) {
	currentTextureCoord[0] = s;
	currentTextureCoord[1] = t;
}


void myNormal(double x, double y, double z) {
	currentNormal[0] = x;
	currentNormal[1] = y;
	currentNormal[2] = z;
	currentNormal[3] = 1.0;
	currentNormal.normalize();
}


void myEnableLight(int lightNo) {
	if (lightNo < MAX_LIGHTS)
		lightEnabled[lightNo] = true;
	else
		std::cerr << "myEnableLight: invalid lightNo: " << lightNo << "\n";
}

void myLight(int lightNo, int field, Vector value) {

	if (!(lightNo < MAX_LIGHTS)) {
		std::cerr << "myEnableLight: invalid lightNo: " << lightNo << "\n";
		return;
	}

	Light *light = &(currentLights[lightNo]);
	switch (field) {
	case GL_POSITION:
		light->positionEnabled = true;
		light->position = value;
		break;
	case GL_AMBIENT:
		light->ambientEnabled = true;
		light->ambient = value;
		break;
	case GL_DIFFUSE:
		if (light->positionEnabled) {
			light->diffuseEnabled = true;
			light->diffuse = value;
		}
		else
			std::cerr << "mylighting: attempt to set GL_DIFFUSE with position not set for lightNo: " << lightNo << "\n";
		break;
	case GL_SPECULAR:
		if (light->positionEnabled) {
			light->specularEnabled = true;
			light->specular = value;
		}
		else
			std::cerr << "mylighting: attempt to set GL_SPECULAR with position not set for lightNo: " << lightNo << "\n";
		break;
	default:
		std::cerr << "unsupported GL property name: " << field << "\n";
	}
}


void myEnableLighting() {
	enableLighting = true;
	myLight(0, GL_POSITION, Vector(-3.0, 4.0, 0.0, 1.0));
	myLight(0, GL_AMBIENT, Vector(0.4, 0.2, 0.2, 1.0));
	myLight(0, GL_DIFFUSE, Vector(0.2, 0.4, 0.2, 1.0));
//	myLight(0, GL_SPECULAR, Vector(0.2, 0.2, 0.4, 1.0));
	myEnableLight(0);
}


void myDisableLighting() {
	enableLighting = false;
	for (int i = 0; i < MAX_LIGHTS; i++)
		lightEnabled[i] = false;
}
