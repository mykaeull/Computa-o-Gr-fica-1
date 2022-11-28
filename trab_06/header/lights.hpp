#ifndef LIGHTS_HPP
#define LIGHTS_HPP
#include "./vector.hpp"
#include <bits/stdc++.h>

using namespace std;

typedef struct Light
{
    Vector intensity;
    Vector position;
    Vector direction;
    double grau;
    string type;

    Light(Vector intensity_l, Vector pisition_l, string type_l) // point or directional
    {
        intensity = intensity_l;
        if (type_l == "point")
            position = pisition_l;
        if (type_l == "directional")
            direction = pisition_l;
        type = type_l;
    }

    Light(Vector intensity_l, Vector pisition_l, Vector direction_l, double grau_l, string type_l) // spot
    {
        intensity = intensity_l;
        position = pisition_l;
        direction = direction_l;
        grau = grau_l;
        type = type_l;
    }

    Light() {}
} Light;

#endif