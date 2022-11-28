#ifndef OBJECTS_HPP
#define OBJECTS_HPP
#include "./vector.hpp"
#include "./mesh_struct.hpp"
#include <bits/stdc++.h>

using namespace std;

typedef struct Object
{
    string type;
    double radius;
    Vector center;
    Vector p_pi;
    Vector normal;
    Vector base;
    double edge;
    vector<Vector> LV;
    vector<VertexIndex> LA;
    vector<EdgeIndex> LF;
    Vector u; // cylinder
    double h;
    double specular;
    Vector k_d;
    Vector k_e;
    Vector k_a;
    string position;
    Vector vc;

    Object(string object_type, double r, Vector c, double s, Vector K_d, Vector K_e, Vector K_a) // Sphere
    {
        type = object_type;
        radius = r;
        center = c;
        specular = s;
        k_d = K_d;
        k_e = K_e;
        k_a = K_a;
    }

    Object(string object_type, Vector P_pi, Vector n, double s, Vector K_d, Vector K_e, Vector K_a, string object_position) // Plane
    {
        type = object_type;
        p_pi = P_pi;
        normal = n;
        specular = s;
        k_d = K_d;
        k_e = K_e;
        k_a = K_a;
        position = object_position;
    }

    Object(string object_type, double r, Vector b, double height, Vector uu, double s, Vector K_d, Vector K_e, Vector K_a) // cylinder or cone
    {
        type = object_type;
        radius = r;
        base = b; // centro da base
        h = height;
        u = uu;
        specular = s;
        k_d = K_d;
        k_e = K_e;
        k_a = K_a;
        vc = Vector(b.x + (height * uu.x), b.y + (height * uu.y), b.z + (height * uu.z), 1); // case cone
    }

    Object(string object_type, double edge_obj, double s, Vector b, Vector K_d, Vector K_e, Vector K_a) // cube
    {
        type = object_type;
        base = b; // centro da base
        edge = edge_obj;
        specular = s;
        k_d = K_d;
        k_e = K_e;
        k_a = K_a;
    }

    Object(string object_type, double edge_obj, double s, Vector b, Vector K_d, Vector K_e, Vector K_a, vector<Vector> LV_cube, vector<VertexIndex> LA_cube, vector<EdgeIndex> LF_cube) // cube
    {
        type = object_type;
        base = b; // centro da base
        edge = edge_obj;
        specular = s;
        k_d = K_d;
        k_e = K_e;
        k_a = K_a;
        LV = LV_cube;
        LA = LA_cube;
        LF = LF_cube;
    }

    Object(string object_type, double s, Vector N, Vector K_d, Vector K_e, Vector K_a)
    { // face of cube
        type = object_type;
        normal = N;
        specular = s;
        k_d = K_d;
        k_e = K_e;
        k_a = K_a;
    }

    Object()
    {
        radius = -1;
    }

} Object;

typedef struct Closest_Object
{
    double closest_t;
    Object object;

    Closest_Object(double t, Object obj)
    {
        closest_t = t;
        object = obj;
    }
    Closest_Object() {}

} Closest_Object;

#endif