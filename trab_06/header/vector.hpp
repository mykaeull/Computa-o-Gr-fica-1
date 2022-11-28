#ifndef VECTOR_HPP
#define VECTOR_HPP

typedef struct Vector
{
    double x, y, z, w;

    Vector(double xx, double yy, double zz, double ww)
    {
        x = xx, y = yy, z = zz, w = ww;
    }
    Vector() {}
} Vector;

#endif