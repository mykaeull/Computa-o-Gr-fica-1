#ifndef COLOR_HPP
#define COLOR_HPP

typedef struct Color
{
    double r, g, b;

    Color(double rr, double gg, double bb)
    {
        r = rr, g = gg, b = bb;
    }
    Color() {}
} Color;

#endif