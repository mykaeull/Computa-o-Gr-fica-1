#ifndef VIEWPORT_HPP
#define VIEWPORT_HPP

typedef struct Viewport
{
    double w, h, d;

    Viewport(double vw, double vh, double vd)
    {
        w = vw, h = vh, d = vd;
    }
    Viewport() {}
} Viewport;

#endif