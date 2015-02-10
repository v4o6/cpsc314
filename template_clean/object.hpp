#ifndef OBJECT_H
#define OBJECT_H

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>

#include "linalg.hpp"
#include "mygl.hpp"

class Object
{
public:

    Object() {}

    // Parses an Object from a file, or retrieves it from the cache if we
    // have seen it before.
    static Object& fromFile(std::string const &filename);

    // Is there data here?
    bool good() const { return polygons_.size(); }

    // Draw the object.
    void draw(bool monochrome = false) const;

protected:
    // A class to represent a single vertex of a polygon. The ints stored within
    // are indices into the positions/texCoords/normals/colors vectors of the
    // Object that it belongs to.
    struct Vertex {
        // Indices into positions, texCoods, normals, and colors vectors.
        int pi, ti, ni, ci;
        Vertex(int pi, int ti, int ni, int ci) : pi(pi), ti(ti), ni(ni), ci(ci) {}
    };

    // Draw this single vertex using myGL commands.
    inline void drawVertex_(const Vertex & v) const;
	mutable bool drawMonochrome_;

    // Storage for positions/texCoords/normals/colors. Looked up by index.
    std::vector<Vector> positions_;
    std::vector<Vector> texCoords_;
    std::vector<Vector> normals_;
    std::vector<Vector> colors_;

    // Polygons are a set of vertices.
    std::vector<std::vector<Vertex> > polygons_;

    // A mapping of filenames to already loaded Objects.
    static std::map<std::string,Object> cache_;

    // Read OBJ data from a given file.
    bool readOBJ(std::string const &filename);

    // Calculates missing normals.
    void fillNormals();

    // Fills in missing colors with normals.
    void fillColors();

    // Centers and resizes the object to be nearly unit size.
    void normalize();

};


#endif

