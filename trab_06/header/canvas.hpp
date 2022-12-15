#ifndef CANVAS_HPP
#define CANVAS_HPP
#include "./viewport.hpp"
#include "./vector.hpp"
#include "./color.hpp"

typedef struct Canvas
{
    double w, h, dx, dy;
    Viewport &vp;
    Color bg;

    Canvas(double wc, double hc, Viewport &vpc, Color bgc) : vp(vpc)
    {
        w = wc, h = hc;
        dx = vpc.w / wc, dy = vpc.h / hc;
        vp = vpc;
        bg = bgc;
    }

    Vector canvas_to_viewport(double x, double y)
    {
        return Vector(-vp.w / 2. + dx / 2. + dx * x, vp.h / 2. - dy / 2. - dy * y, -vp.d, 0);
    }

    void update_data()
    {
        dx = vp.w / w, dy = vp.h / h;
    }
}

Canvas;

#endif