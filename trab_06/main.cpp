#include <bits/stdc++.h>
#include "./header/vector.hpp"
#include "./header/viewport.hpp"
#include "./header/color.hpp"
#include "./header/canvas.hpp"
#include "./header/mesh_struct.hpp"
#include "./header/objects.hpp"
#include "./header/lights.hpp"
#include "./header/scene.hpp"

#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "stb_image.h"
#include "stb_image_write.h"

typedef unsigned char int8;
#define CHANNEL_NUM 3

using namespace std;

double dot(Vector a, Vector b)
{
    return (a.x * b.x) + (a.y * b.y) + (a.z * b.z) + (a.w * b.w);
}

Vector vector_mult(Vector a, Vector b, bool p)
{
    double row1 = (a.y * b.z) - (a.z * b.y);
    double row2 = (a.z * b.x) - (a.x * b.z);
    double row3 = (a.x * b.y) - (a.y * b.x);
    double x = p == 0 ? 0 : 1; // 0: vetor; 1: ponto
    return Vector(row1, row2, row3, x);
}

Vector sum_vector(Vector a, Vector b, bool p)
{
    double x = p == 0 ? 0 : 1; // 0: vetor; 1: ponto
    return Vector(a.x + b.x, a.y + b.y, a.z + b.z, x);
}

Vector sub_vector(Vector a, Vector b, bool p)
{
    double x = p == 0 ? 0 : 1; // 0: vetor; 1: ponto
    return Vector(a.x - b.x, a.y - b.y, a.z - b.z, x);
}

double length(Vector v)
{
    return sqrt((v.x * v.x) + (v.y * v.y) + (v.z * v.z));
}

Vector norm_vector(Vector v)
{
    double length_v = length(v);
    return Vector(v.x / length_v, v.y / length_v, v.z / length_v, 0);
}

vector<vector<double>> matrix_I()
{
    vector<vector<double>> I;

    I.push_back({1., 0., 0., 0.});
    I.push_back({0., 1., 0., 0.});
    I.push_back({0., 0., 1., 0.});
    I.push_back({0., 0., 0., 1.});

    return I;
}

vector<vector<double>> mult_matrix(vector<vector<double>> M, vector<vector<double>> I)
{
    vector<vector<double>> X;

    double a = (M[0][0] * I[0][0]) + (M[0][1] * I[1][0]) + (M[0][2] * I[2][0]) + (M[0][3] * I[3][0]);
    double b = (M[0][0] * I[0][1]) + (M[0][1] * I[1][1]) + (M[0][2] * I[2][1]) + (M[0][3] * I[3][1]);
    double c = (M[0][0] * I[0][2]) + (M[0][1] * I[1][2]) + (M[0][2] * I[2][2]) + (M[0][3] * I[3][2]);
    double d = (M[0][0] * I[0][3]) + (M[0][1] * I[1][3]) + (M[0][2] * I[2][3]) + (M[0][3] * I[3][3]);
    X.push_back({a, b, c, d});
    a = (M[1][0] * I[0][0]) + (M[1][1] * I[1][0]) + (M[1][2] * I[2][0]) + (M[1][3] * I[3][0]);
    b = (M[1][0] * I[0][1]) + (M[1][1] * I[1][1]) + (M[1][2] * I[2][1]) + (M[1][3] * I[3][1]);
    c = (M[1][0] * I[0][2]) + (M[1][1] * I[1][2]) + (M[1][2] * I[2][2]) + (M[1][3] * I[3][2]);
    d = (M[1][0] * I[0][3]) + (M[1][1] * I[1][3]) + (M[1][2] * I[2][3]) + (M[1][3] * I[3][3]);
    X.push_back({a, b, c, d});
    a = (M[2][0] * I[0][0]) + (M[2][1] * I[1][0]) + (M[2][2] * I[2][0]) + (M[2][3] * I[3][0]);
    b = (M[2][0] * I[0][1]) + (M[2][1] * I[1][1]) + (M[2][2] * I[2][1]) + (M[2][3] * I[3][1]);
    c = (M[2][0] * I[0][2]) + (M[2][1] * I[1][2]) + (M[2][2] * I[2][2]) + (M[2][3] * I[3][2]);
    d = (M[2][0] * I[0][3]) + (M[2][1] * I[1][3]) + (M[2][2] * I[2][3]) + (M[2][3] * I[3][3]);
    X.push_back({a, b, c, d});
    a = (M[3][0] * I[0][0]) + (M[3][1] * I[1][0]) + (M[3][2] * I[2][0]) + (M[3][3] * I[3][0]);
    b = (M[3][0] * I[0][1]) + (M[3][1] * I[1][1]) + (M[3][2] * I[2][1]) + (M[3][3] * I[3][1]);
    c = (M[3][0] * I[0][2]) + (M[3][1] * I[1][2]) + (M[3][2] * I[2][2]) + (M[3][3] * I[3][2]);
    d = (M[3][0] * I[0][3]) + (M[3][1] * I[1][3]) + (M[3][2] * I[2][3]) + (M[3][3] * I[3][3]);
    X.push_back({a, b, c, d});

    return X;
}

Vector transform_W_to_C(vector<vector<double>> M, Vector P, bool p)
{
    Vector a;
    Vector b;
    Vector c;
    Vector d;
    Vector new_P;
    double x = p == 0 ? 0 : 1; // 0: vetor; 1: ponto

    a = Vector(P.x * M[0][0], P.x * M[1][0], P.x * M[2][0], P.x * M[3][0]);
    b = Vector(P.y * M[0][1], P.y * M[1][1], P.y * M[2][1], P.y * M[3][1]);
    c = Vector(P.z * M[0][2], P.z * M[1][2], P.z * M[2][2], P.z * M[3][2]);
    d = Vector(P.w * M[0][3], P.w * M[1][3], P.w * M[2][3], P.w * x);
    new_P = sum_vector(sum_vector(a, b, 1), sum_vector(c, d, 1), 1);
    // Vector new_P_norm = Vector(new_P.x / length(new_P), new_P.y / length(new_P), new_P.z / length(new_P), new_P.w / length(new_P));
    return new_P;
}

vector<vector<double>> matrix_C_to_W(Vector O, Vector look_at, Vector view_up)
{
    Vector O_sub_lookAt = sub_vector(O, look_at, 1);

    Vector kc = Vector(O_sub_lookAt.x / length(O_sub_lookAt), O_sub_lookAt.y / length(O_sub_lookAt), O_sub_lookAt.z / length(O_sub_lookAt), 0);

    Vector view_up_mult_kc = vector_mult(view_up, kc, 0);

    Vector ic = Vector(view_up_mult_kc.x / length(view_up_mult_kc), view_up_mult_kc.y / length(view_up_mult_kc), view_up_mult_kc.z / length(view_up_mult_kc), 0);

    Vector kc_mult_ic = vector_mult(kc, ic, 0);

    Vector jc = Vector(kc_mult_ic.x / length(kc_mult_ic), kc_mult_ic.y / length(kc_mult_ic), kc_mult_ic.z / length(kc_mult_ic), 0);

    vector<vector<double>> M_C_to_W;

    M_C_to_W.push_back({ic.x, ic.y, ic.z, -dot(ic, O)});
    M_C_to_W.push_back({jc.x, jc.y, jc.z, -dot(jc, O)});
    M_C_to_W.push_back({kc.x, kc.y, kc.z, -dot(kc, O)});
    M_C_to_W.push_back({0., 0., 0., 1.});

    return mult_matrix(M_C_to_W, matrix_I());
}

vector<vector<double>> matrix_translatef(double x, double y, double z, Vector P)
{
    vector<vector<double>> M_Traslatef;

    M_Traslatef.push_back({1., 0., 0., x - P.x});
    M_Traslatef.push_back({0., 1., 0., y - P.y});
    M_Traslatef.push_back({0., 0., 1., z - P.z});
    M_Traslatef.push_back({0., 0., 0., 1.});

    return mult_matrix(M_Traslatef, matrix_I());
}

vector<vector<double>> matrix_rotation(string axis, double angle)
{
    vector<vector<double>> M_Rotation;
    if (axis == "x")
    {
        M_Rotation.push_back({1., 0., 0., 0.});
        M_Rotation.push_back({0., cos(angle), -sin(angle), 0.});
        M_Rotation.push_back({0., sin(angle), cos(angle), 0.});
        M_Rotation.push_back({0., 0., 0., 1.});

        return mult_matrix(M_Rotation, matrix_I());
    }
    else if (axis == "y")
    {
        M_Rotation.push_back({cos(angle), 0., sin(angle), 0.});
        M_Rotation.push_back({0., 1, 0., 0.});
        M_Rotation.push_back({-sin(angle), 0., cos(angle), 0.});
        M_Rotation.push_back({0., 0., 0., 1.});

        return mult_matrix(M_Rotation, matrix_I());
    }

    M_Rotation.push_back({cos(angle), -sin(angle), 0., 0.});
    M_Rotation.push_back({sin(angle), cos(angle), 0., 0.});
    M_Rotation.push_back({0., 0., 1., 0.});
    M_Rotation.push_back({0., 0., 0., 1.});

    return mult_matrix(M_Rotation, matrix_I());
}

vector<vector<double>> matrix_scale(double x, double y, double z)
{
    vector<vector<double>> M_Scale;

    M_Scale.push_back({x, 0., 0., 0.});
    M_Scale.push_back({0., y, 0., 0.});
    M_Scale.push_back({0., 0., z, 0.});
    M_Scale.push_back({0., 0., 0., 1.});

    return mult_matrix(M_Scale, matrix_I());
}

// vector<vector<double>> matrices_transform(vector<vector<vector<double>>> matrix_list)
// {
//     vector<vector<double>> M;

//     M.push_back({1., 0., 0., 0.});
//     M.push_back({0., 1., 0., 0.});
//     M.push_back({0., 0., 1., 0.});
//     M.push_back({0., 0., 0., 1.});

//     for (int i = 0; i < matrix_list.size(); i++)
//     {
//         M = mult_matrix(matrix_list[i], M);
//     }

//     return M;
// }

int main()
{
    Viewport vp(60., 60., 30.);

    Canvas canva(500., 500., vp, Color(100, 100, 100));

    vector<Object> objects;

    vector<Light> lights;

    vector<vector<double>> M_C_to_W;
    vector<vector<double>> M_translatef;
    vector<vector<double>> M_rotation;
    vector<vector<double>> M_scale;

    Vector O = Vector(0., 0., 0., 1);

    M_C_to_W = matrix_C_to_W(O, Vector(0, -60, -200, 0), Vector(0., 100., 0., 1));
    M_translatef = matrix_translatef(60., -150., -165., Vector(0, -150., -165., 1.));
    M_rotation = matrix_rotation("z", 3.141592 / 4.);
    // matrix_list.push_back(matrix_scale(2., 2., 1.));

    // vector<vector<double>> M_I;

    // M_I = mult_matrix(matrix_C_to_W(O, Vector(0, -60, -200, 0), Vector(0., 100., 0., 1)), I);

    Object plane1("plane", transform_W_to_C(M_C_to_W, Vector(0, -150, 0, 1), 1), norm_vector(transform_W_to_C(M_C_to_W, Vector(0., 1., 0., 0), 0)), 1., Vector(0.933, 0.933, 0.933, 0), Vector(0.933, 0.933, 0.933, 0), Vector(0.933, 0.933, 0.933, 0), "floor");
    Object plane2("plane", transform_W_to_C(M_C_to_W, Vector(200, -150, 0, 1), 1), norm_vector(transform_W_to_C(M_C_to_W, Vector(-1., 0., 0., 0), 0)), 1., Vector(0.686, 0.933, 0.933, 0), Vector(0.686, 0.933, 0.933, 0), Vector(0.686, 0.933, 0.933, 0), "right");
    Object plane3("plane", transform_W_to_C(M_C_to_W, Vector(200, -150, -400, 1), 1), norm_vector(transform_W_to_C(M_C_to_W, Vector(0., 0., 1., 0), 0)), 1., Vector(0.686, 0.933, 0.933, 0), Vector(0.686, 0.933, 0.933, 0), Vector(0.686, 0.933, 0.933, 0), "front");
    Object plane4("plane", transform_W_to_C(M_C_to_W, Vector(-200, -150, 0, 1), 1), norm_vector(transform_W_to_C(M_C_to_W, Vector(1., 0., 0., 0), 0)), 1., Vector(0.686, 0.933, 0.933, 0), Vector(0.686, 0.933, 0.933, 0), Vector(0.686, 0.933, 0.933, 0), "left");
    Object plane5("plane", transform_W_to_C(M_C_to_W, Vector(0, 150, 0, 1), 1), norm_vector(transform_W_to_C(M_C_to_W, Vector(0., -1., 0., 0), 0)), 1., Vector(0.933, 0.933, 0.933, 0), Vector(0.933, 0.933, 0.933, 0), Vector(0.933, 0.933, 0.933, 0), "ceil");

    objects.push_back(plane1);
    objects.push_back(plane2);
    objects.push_back(plane3);
    objects.push_back(plane4);
    objects.push_back(plane5);

    Object cylinder1("cylinder", 5., transform_W_to_C(M_C_to_W, Vector(0, -150, -200., 1.), 1), 90., norm_vector(transform_W_to_C(M_C_to_W, Vector(0, 1., 0, 0), 0)), 10., Vector(0.824, 0.706, 0.549, 0), Vector(0.824, 0.706, 0.549, 0), Vector(0.824, 0.706, 0.549, 0));
    objects.push_back(cylinder1);

    Object cone1("cone", 90, Vector(0, -60, -200, 1.), 150, Vector(0, 1., 0, 0), 10., Vector(0., 1., 0.498, 0), Vector(0., 1., 0.498, 0), Vector(0., 1., 0.498, 0));
    cone1.base = transform_W_to_C(M_C_to_W, cone1.base, 1);
    cone1.u = norm_vector(transform_W_to_C(M_C_to_W, cone1.u, 0));
    cone1.vc = transform_W_to_C(M_C_to_W, cone1.vc, 1);
    objects.push_back(cone1);

    Object sphere1("sphere", 5., transform_W_to_C(M_C_to_W, Vector(0, 95, -200., 1.), 1), 10., Vector(0.854, 0.647, 0.125, 0), Vector(0.854, 0.647, 0.125, 0), Vector(0.854, 0.647, 0.125, 0));
    objects.push_back(sphere1);

    Object cube("cube", 40., 10., Vector(0, -150., -165., 1.), Vector(1., 0.078, 0.576, 0), Vector(1., 0.078, 0.576, 0), Vector(1., 0.078, 0.576, 0));
    cube.base = transform_W_to_C(M_C_to_W, transform_W_to_C(M_rotation, cube.base, 1), 1);
    cube.LV[0] = transform_W_to_C(M_C_to_W, transform_W_to_C(M_rotation, cube.LV[0], 1), 1);
    cube.LV[1] = transform_W_to_C(M_C_to_W, transform_W_to_C(M_rotation, cube.LV[1], 1), 1);
    cube.LV[2] = transform_W_to_C(M_C_to_W, transform_W_to_C(M_rotation, cube.LV[2], 1), 1);
    cube.LV[3] = transform_W_to_C(M_C_to_W, transform_W_to_C(M_rotation, cube.LV[3], 1), 1);
    cube.LV[4] = transform_W_to_C(M_C_to_W, transform_W_to_C(M_rotation, cube.LV[4], 1), 1);
    cube.LV[5] = transform_W_to_C(M_C_to_W, transform_W_to_C(M_rotation, cube.LV[5], 1), 1);
    cube.LV[6] = transform_W_to_C(M_C_to_W, transform_W_to_C(M_rotation, cube.LV[6], 1), 1);
    cube.LV[7] = transform_W_to_C(M_C_to_W, transform_W_to_C(M_rotation, cube.LV[7], 1), 1);
    objects.push_back(cube);

    Light point_light(Vector(0.7, 0.7, 0.7, 0), transform_W_to_C(M_C_to_W, Vector(-100, 140., -20., 1.), 1), "point");
    Light ambient_light(Vector(0.3, 0.3, 0.3, 0), Vector(0, 0, 0, 0), "ambient");
    Light directional_light(Vector(0.5, 0.5, 0.5, 0), norm_vector(transform_W_to_C(M_C_to_W, Vector(0, -150., -165., 0), 0)), "directional");
    Light spot_light(Vector(0.7, 0.7, 0.7, 0), transform_W_to_C(M_C_to_W, Vector(0, 140, 0, 1.), 1), norm_vector(transform_W_to_C(M_C_to_W, Vector(0, -150., -165., 0), 0)), 0.5, "spot");

    // Vector test = transform_W_to_C(M_rotation, Vector(20., 20., 20., 1.), 1);

    // cout << test.x << " " << test.y << " " << test.z << " " << test.w << "\n";

    lights.push_back(point_light);
    lights.push_back(ambient_light);
    // lights.push_back(directional_light);
    // lights.push_back(spot_light);

    Scene scene(objects, canva, lights);

    int w, h, chan;
    Color **matrix_image;

    unsigned char *image = stbi_load("madeira.png", &w, &h, &chan, 3);
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

    stbi_image_free(image);

    ofstream out("out.ppm");

    out << "P3";
    out << endl;
    out << canva.w << endl;
    out << canva.h << endl;
    out << 255 << endl;

    int8 *rgb_image = new int8[(int)canva.w * (int)canva.h * CHANNEL_NUM];
    for (int y = 0, c = 0; y < canva.h; y++)
    {
        for (int x = 0; x < canva.w; x++)
        {
            Vector D = canva.canvas_to_viewport(x, y);

            Color pixel_image = matrix_image[y % h][x % w];

            Color color = scene.trace_ray(O, D, 0.0, INFINITY, pixel_image);

            rgb_image[c++] = min((int)color.r, 255);
            rgb_image[c++] = min((int)color.g, 255);
            rgb_image[c++] = min((int)color.b, 255);

            out << min((int)color.r, 255) << endl;
            out << min((int)color.g, 255) << endl;
            out << min((int)color.b, 255) << endl;
        }
    }

    out.close();

    stbi_write_png("out.png", canva.w, canva.h, CHANNEL_NUM, rgb_image, canva.w * CHANNEL_NUM);

    stbi_image_free(rgb_image);

    return 0;
}