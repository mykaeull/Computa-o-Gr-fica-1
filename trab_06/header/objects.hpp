#ifndef OBJECTS_HPP
#define OBJECTS_HPP
#include "./vector.hpp"
#include "./mesh_struct.hpp"
#include <bits/stdc++.h>
#include "./color.hpp"

#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "../stb_image.h"
#include "../stb_image_write.h"

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
    bool texture;
    Color **matrix_img;
    int height;
    int width;

    tuple<vector<Vector>, vector<VertexIndex>, vector<EdgeIndex>> cube_maping(Vector b, double edge_obj)
    {
        Vector centro = Vector(b.x, b.y, b.z, 1);
        Vector A = Vector(centro.x - (edge_obj / 2.), centro.y - (edge_obj / 2.), centro.z + (edge_obj / 2.), 1);
        Vector B = Vector(centro.x - (edge_obj / 2.), centro.y - (edge_obj / 2.), centro.z - (edge_obj / 2.), 1);
        Vector C = Vector(centro.x + (edge_obj / 2.), centro.y - (edge_obj / 2.), centro.z - (edge_obj / 2.), 1);
        Vector D = Vector(centro.x + (edge_obj / 2.), centro.y - (edge_obj / 2.), centro.z + (edge_obj / 2.), 1);
        Vector E = Vector(centro.x - (edge_obj / 2.), centro.y + (edge_obj / 2.), centro.z + (edge_obj / 2.), 1);
        Vector F = Vector(centro.x - (edge_obj / 2.), centro.y + (edge_obj / 2.), centro.z - (edge_obj / 2.), 1);
        Vector G = Vector(centro.x + (edge_obj / 2.), centro.y + (edge_obj / 2.), centro.z - (edge_obj / 2.), 1);
        Vector H = Vector(centro.x + (edge_obj / 2.), centro.y + (edge_obj / 2.), centro.z + (edge_obj / 2.), 1);
        vector<Vector> LV;
        LV.push_back(A);
        LV.push_back(B);
        LV.push_back(C);
        LV.push_back(D);
        LV.push_back(E);
        LV.push_back(F);
        LV.push_back(G);
        LV.push_back(H);
        vector<VertexIndex> LA;
        LA.push_back(VertexIndex(0, 1));
        LA.push_back(VertexIndex(1, 2));
        LA.push_back(VertexIndex(2, 3));
        LA.push_back(VertexIndex(3, 0));
        LA.push_back(VertexIndex(4, 5));
        LA.push_back(VertexIndex(5, 6));
        LA.push_back(VertexIndex(6, 7));
        LA.push_back(VertexIndex(7, 4));
        LA.push_back(VertexIndex(0, 4));
        LA.push_back(VertexIndex(1, 5));
        LA.push_back(VertexIndex(2, 6));
        LA.push_back(VertexIndex(3, 7));
        LA.push_back(VertexIndex(2, 7));
        LA.push_back(VertexIndex(5, 7));
        LA.push_back(VertexIndex(5, 2));
        LA.push_back(VertexIndex(1, 4));
        LA.push_back(VertexIndex(1, 3));
        LA.push_back(VertexIndex(3, 4));
        vector<EdgeIndex> LF;
        LF.push_back(EdgeIndex(6, 10, 12));
        LF.push_back(EdgeIndex(12, 2, 11));
        LF.push_back(EdgeIndex(7, 4, 13));
        LF.push_back(EdgeIndex(13, 5, 6));
        LF.push_back(EdgeIndex(5, 14, 10));
        LF.push_back(EdgeIndex(9, 1, 14));
        LF.push_back(EdgeIndex(4, 15, 9));
        LF.push_back(EdgeIndex(8, 0, 15));
        LF.push_back(EdgeIndex(1, 16, 2));
        LF.push_back(EdgeIndex(3, 16, 10));
        LF.push_back(EdgeIndex(11, 17, 7));
        LF.push_back(EdgeIndex(3, 8, 17));
        auto t = make_tuple(LV, LA, LF);
        return t;
    }

    Vector convert_rgb(Vector rgb_color)
    {
        // cout << rgb_color.x / 255 << " " << rgb_color.y / 255 << " " << rgb_color.z / 255 << "\n";
        return Vector(rgb_color.x / 255., rgb_color.y / 255., rgb_color.z / 255., 0.);
    }

    Object(string object_type, double r, Vector c, double s, Vector K_d, Vector K_e, Vector K_a) // Sphere
    {
        type = object_type;
        radius = r;
        center = c;
        specular = s;
        k_d = convert_rgb(K_d);
        k_e = convert_rgb(K_e);
        k_a = convert_rgb(K_a);
    }

    Object(string object_type, Vector P_pi, Vector n, double s, Vector K_d, Vector K_e, Vector K_a, bool have_texture, const char *path_image) // Plane
    {
        type = object_type;
        p_pi = P_pi;
        normal = n;
        specular = s;
        k_d = convert_rgb(K_d);
        k_e = convert_rgb(K_e);
        k_a = convert_rgb(K_a);
        texture = have_texture;
        if (have_texture)
        {
            int w, h, chan;
            Color **matrix_image;

            unsigned char *image = stbi_load(path_image, &w, &h, &chan, 3);
            matrix_image = new Color *[h];
            for (int i = 0; i < h; i++)
            {
                matrix_image[i] = new Color[w];
            }
            for (int i = 0; i < h; i++)
            {
                for (int j = 0; j < w; j++)
                {
                    matrix_image[i][j].r = (double)image[j * 3 + w * i * 3];
                    matrix_image[i][j].g = (double)image[j * 3 + w * i * 3 + 1];
                    matrix_image[i][j].b = (double)image[j * 3 + w * i * 3 + 2];
                }
            }

            height = h;
            width = w;

            stbi_image_free(image);

            matrix_img = matrix_image;
        }
    }

    Object(string object_type, double r, Vector b, double height, Vector uu, double s, Vector K_d, Vector K_e, Vector K_a) // cylinder or cone
    {
        type = object_type;
        radius = r;
        base = b; // centro da base
        h = height;
        u = uu;
        specular = s;
        k_d = convert_rgb(K_d);
        k_e = convert_rgb(K_e);
        k_a = convert_rgb(K_a);
        vc = Vector(b.x + (height * uu.x), b.y + (height * uu.y), b.z + (height * uu.z), 1); // case cone
    }

    Object(string object_type, double edge_obj, double s, Vector b, Vector K_d, Vector K_e, Vector K_a) // cube
    {
        type = object_type;
        base = b; // centro da base
        edge = edge_obj;
        specular = s;
        k_d = convert_rgb(K_d);
        k_e = convert_rgb(K_e);
        k_a = convert_rgb(K_a);
        auto t = cube_maping(b, edge_obj);
        LV = get<0>(t);
        LA = get<1>(t);
        LF = get<2>(t);
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