CPSC 314 Assignment 3
Myron Yao 80226046 v4o6

Compiling and Running on *nix
=============================

The provided Makefile will compile this assignment on a number of Linux/Unix
base systems. All you need to do is run:

    make

To run the executable:

    ./a3

To run the executable and display a different OBJ file (which you will need to
find yourself), pass it as the first argument:

    ./a3 stanford_bunny.obj

Compiling and Running on Windows
=============================

Use the provided Visual Studio 2010 project file


Interface
=========

You can manipulate the camera by dragging with the mouse. There are several
variations:

- normal drag: orbit the camera around the origin.
- alt-drag: move the camera towards/away from the camera's focus. (only effective if perspective is on)
- shift-drag: pan the camera parallel to the viewport.

There are a several key bindings to control the running executable:

- q OR esc: exit.
- a through l: display a different scene.
- < OR >: decrease or increase the virtual pixel size.
- /: disable the grid and set pixel size to 1.
- .: toggle visibility of the grid.
- t: toggle perspective correct textures.
- o: toggle orthographic/perspective projection


Scenarios
=========

The pre-existing scenarios are:

A: Three points, one each of red, green, and blue.

B: A white, transformed triangle. If implemented in assignment order it will
   be displayed as three points.

C: The front faces of a unit cube centered at the origin, with the faces
   coloured red, green, and blue. If implemented in assignment order it will
   be displayed as 7 points (there are no backfaces so we are missing the eighth).

D: A red triangle, and a blue triangle (which is partially off screen to test
   proper clipping; you will be warned if you try to draw off screen).

E: A single triangle (made from the points in scenario A). Colours should
   interpolate.

F: An intersecting square and triangle.

G: A loaded OBJ, centered at the origin.

H: A solid with 2 pentagons and 5 quadrangles

I: A loaded OBJ, enabling lighting

J: A loaded OBJ, calling myPolygonMode(GL_LINE) (i.e., wire-frame mode)

K: The front faces of a texture mapped unit cube centered at the origin.

L: Empty scene.
