#include <bits/stdc++.h>
#include "./header/vector.hpp"
#include "./header/viewport.hpp"
#include "./header/color.hpp"
#include "./header/canvas.hpp"
#include "./header/mesh_struct.hpp"
#include "./header/objects.hpp"
#include "./header/lights.hpp"
#include "./header/scene.hpp"
#include <SFML/Graphics.hpp>
#include <thread>

typedef unsigned char int8;
#define CHANNEL_NUM 4

using namespace std;

Object *old_obj = nullptr;
Object *current_object = nullptr;

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

vector<vector<double>> matrix_W_to_C(Vector O, Vector look_at, Vector view_up)
{
    Vector K = sub_vector(O, look_at, 0);

    double Lk = length(K);

    Vector kc = Vector(K.x / Lk, K.y / Lk, K.z / Lk, 0);

    Vector O_to_up = sub_vector(view_up, O, 0);

    Vector I = vector_mult(O_to_up, kc, 0);

    double LI = length(I);

    Vector ic = Vector(I.x / LI, I.y / LI, I.z / LI, 0);

    Vector jc = vector_mult(kc, ic, 0);

    vector<vector<double>> M_W_to_C;

    M_W_to_C.push_back({ic.x, ic.y, ic.z, -dot(ic, O)});
    M_W_to_C.push_back({jc.x, jc.y, jc.z, -dot(jc, O)});
    M_W_to_C.push_back({kc.x, kc.y, kc.z, -dot(kc, O)});
    M_W_to_C.push_back({0., 0., 0., 1.});

    return mult_matrix(M_W_to_C, matrix_I());
}

vector<vector<double>> matrix_C_to_W(Vector E, Vector look_at, Vector view_up)
{
    Vector K = sub_vector(E, look_at, 0);

    double Lk = length(K);

    Vector kc = Vector(K.x / Lk, K.y / Lk, K.z / Lk, 0);

    Vector O_to_up = sub_vector(view_up, E, 0);

    Vector I = vector_mult(O_to_up, kc, 0);

    double LI = length(I);

    Vector ic = Vector(I.x / LI, I.y / LI, I.z / LI, 0);

    Vector jc = vector_mult(kc, ic, 0);

    vector<vector<double>> M_C_to_W;

    M_C_to_W.push_back({ic.x, jc.x, kc.x, E.x});
    M_C_to_W.push_back({ic.y, jc.y, kc.y, E.y});
    M_C_to_W.push_back({ic.z, jc.z, kc.z, E.z});
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

vector<vector<double>> matrix_shear(string axis, double angle)
{
    vector<vector<double>> M_shear;
    if (axis == "x")
    {
        M_shear.push_back({1., tan(angle), 0., 0.});
        M_shear.push_back({0., 1., 0., 0.});
        M_shear.push_back({0., 0., 1., 0.});
        M_shear.push_back({0., 0., 0., 1.});

        return mult_matrix(M_shear, matrix_I());
    }

    M_shear.push_back({1., 0., 0., 0.});
    M_shear.push_back({tan(angle), 1., 0., 0.});
    M_shear.push_back({0., 0., 1., 0.});
    M_shear.push_back({0., 0., 0., 1.});

    return mult_matrix(M_shear, matrix_I());
}

void transform_plane(Object *plane, vector<vector<double>> M)
{
    plane->center = transform_W_to_C(M, plane->center, 1);
    plane->normal = norm_vector(transform_W_to_C(M, plane->normal, 0));
}

void transform_cylinder(Object *cylinder, vector<vector<double>> M)
{
    cylinder->center = transform_W_to_C(M, cylinder->center, 1);
    cylinder->u = norm_vector(transform_W_to_C(M, cylinder->u, 0));
}

void transform_cone(Object *cone, vector<vector<double>> M)
{
    cone->center = transform_W_to_C(M, cone->center, 1);
    cone->u = norm_vector(transform_W_to_C(M, cone->u, 0));
    cone->vc = transform_W_to_C(M, cone->vc, 1);
}

void transform_sphere(Object *sphere, vector<vector<double>> M)
{
    sphere->center = transform_W_to_C(M, sphere->center, 1);
}

void transform_cube(Object *cube, vector<vector<double>> M)
{
    cube->center = transform_W_to_C(M, cube->center, 1);
    for (int i = 0; i < cube->LV.size(); i++)
    {
        cube->LV[i] = transform_W_to_C(M, cube->LV[i], 1);
    }
}

auto switch_transform(string type)
{
    if (type == "plane")
    {
        return transform_plane;
    }
    else if (type == "cylinder")
    {
        return transform_cylinder;
    }
    else if (type == "cone")
    {
        return transform_cone;
    }
    else if (type == "sphere")
    {
        return transform_sphere;
    }
    return transform_cube;
}

Vector E1 = Vector(0., 500., 800., 1);

Vector A1 = Vector(0., 115., -550., 1);

Vector V1 = Vector(0., 501., 800., 1);

double canva_width = 500.;
double canva_height = 500.;

int main()
{
    Viewport vp(60., 60., 30.); // 60., 60., 30. - 1100., 1000., -900.

    Canvas canva(canva_width, canva_height, vp, Color(100, 100, 100));

    vector<Object *> objects;

    vector<Light *> lights;

    vector<vector<double>> M_W_to_C;

    // visão boa de cima

    M_W_to_C = matrix_W_to_C(E1, A1, V1);

    Scene scene(canva);

    // --------------------------------------

    // Vector O = Vector(0., 0., 0., 1);

    // M_W_to_C = matrix_W_to_C(Vector(0., 250., 0., 1), Vector(40., 115., -550., 0), Vector(0., 251., 0., 1));

    // -----------------------------------

    // Vector O = Vector(-300., -170., 800., 1);

    // M_W_to_C = matrix_W_to_C(Vector(0., 0., 0., 1), Vector(40., 115., -550., 0), Vector(0., 1., 0., 1));

    Object *down_plane = new Object("plane", Vector(0, -150, 0, 1), Vector(0., 1., 0., 0), 1., Vector(236., 219., 182., 0), Vector(236., 219., 182., 0), Vector(236., 219., 182., 0), true, "madeira5.png");
    transform_plane(down_plane, M_W_to_C);
    Object *right_plane = new Object("plane", Vector(600, -150, 0, 1), Vector(-1., 0., 0., 0), 1., Vector(236., 219., 182., 0), Vector(236., 219., 182., 0), Vector(236., 219., 182., 0), false, "");
    transform_plane(right_plane, M_W_to_C);
    Object *back_plane = new Object("plane", Vector(600, -150, -600, 1), Vector(0., 0., 1., 0), 1., Vector(236., 219., 182., 0), Vector(236., 219., 182., 0), Vector(236., 219., 182., 0), false, "");
    transform_plane(back_plane, M_W_to_C);
    Object *left_plane = new Object("plane", Vector(-600, -150, 0, 1), Vector(1., 0., 0., 0), 1., Vector(236., 219., 182., 0), Vector(236., 219., 182., 0), Vector(236., 219., 182., 0), false, "");
    transform_plane(left_plane, M_W_to_C);

    scene.objects.push_back(down_plane);
    scene.objects.push_back(right_plane);
    scene.objects.push_back(back_plane);
    scene.objects.push_back(left_plane);

    // Object *tree_trunk = new Object("cylinder", 5., Vector(0, -150, 350., 1.), 90., Vector(0, 1., 0, 0), 10., Vector(127., 80., 47., 0), Vector(127., 80., 47., 0), Vector(127., 80., 47., 0));
    // transform_cylinder(tree_trunk, M_W_to_C);
    // scene.objects.push_back(tree_trunk);

    // Object *tree_base = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(73., 17., 13., 0), Vector(73., 17., 13., 0), Vector(73., 17., 13., 0));
    // transform_cube(tree_base, matrix_scale(100. / 10., 30. / 10., 1.));
    // transform_cube(tree_base, matrix_translatef(0, -135, 350., tree_base->center));
    // transform_cube(tree_base, M_W_to_C);
    // scene.objects.push_back(tree_base);

    // Object *tree_top = new Object("cone", 75, Vector(0, -60, 350, 1.), 150, Vector(0, 1., 0, 0), 10., Vector(37., 83., 19., 0), Vector(37., 83., 19., 0), Vector(37., 83., 19., 0));
    // transform_cone(tree_top, M_W_to_C);
    // scene.objects.push_back(tree_top);

    // Object *tree_ball = new Object("sphere", 5., Vector(0, 95, 350., 1.), 10., Vector(247., 212., 36., 0), Vector(247., 212., 36., 0), Vector(247., 212., 36., 0));
    // transform_sphere(tree_ball, M_W_to_C);
    // scene.objects.push_back(tree_ball);

    // Object *door = new Object("cube", 20., 10., Vector(-400, 100., -600., 1.), Vector(128., 69., 22., 0), Vector(128., 69., 22., 0), Vector(128., 69., 22., 0));
    // transform_cube(door, matrix_scale(200. / 20., 425. / 20., 1.));
    // transform_cube(door, matrix_translatef(-500., 65., -590., door->center));
    // transform_cube(door, M_W_to_C);
    // scene.objects.push_back(door);

    // Object *table_leg1 = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(218., 161., 17., 0), Vector(218., 161., 17., 0), Vector(218., 161., 17., 0));
    // transform_cube(table_leg1, matrix_scale(15. / 10., 190. / 10., 15. / 10.));
    // transform_cube(table_leg1, matrix_translatef(-165., -55., -580., table_leg1->center));
    // transform_cube(table_leg1, M_W_to_C);
    // scene.objects.push_back(table_leg1);

    // Object *table_leg2 = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(218., 161., 17., 0), Vector(218., 161., 17., 0), Vector(218., 161., 17., 0));
    // transform_cube(table_leg2, matrix_scale(15. / 10., 190. / 10., 15. / 10.));
    // transform_cube(table_leg2, matrix_translatef(165., -55., -580., table_leg2->center));
    // transform_cube(table_leg2, M_W_to_C);
    // scene.objects.push_back(table_leg2);

    // Object *table_leg3 = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(218., 161., 17., 0), Vector(218., 161., 17., 0), Vector(218., 161., 17., 0));
    // transform_cube(table_leg3, matrix_scale(15. / 10., 190. / 10., 15. / 10.));
    // transform_cube(table_leg3, matrix_translatef(-165., -55., -480., table_leg3->center));
    // transform_cube(table_leg3, M_W_to_C);
    // scene.objects.push_back(table_leg3);

    // Object *table_leg4 = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(218., 161., 17., 0), Vector(218., 161., 17., 0), Vector(218., 161., 17., 0));
    // transform_cube(table_leg4, matrix_scale(15. / 10., 190. / 10., 15. / 10.));
    // transform_cube(table_leg4, matrix_translatef(165., -55., -480., table_leg4->center));
    // transform_cube(table_leg4, M_W_to_C);
    // scene.objects.push_back(table_leg4);

    // Object *table_cover = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(218., 161., 17., 0), Vector(218., 161., 17., 0), Vector(218., 161., 17., 0));
    // transform_cube(table_cover, matrix_scale(350. / 10., 15. / 10., 120. / 10.));
    // transform_cube(table_cover, matrix_translatef(0., 50., -530., table_cover->center));
    // transform_cube(table_cover, M_W_to_C);
    // scene.objects.push_back(table_cover);

    // Object *pc_base = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(179., 183., 186., 0), Vector(179., 183., 186., 0), Vector(179., 183., 186., 0));
    // transform_cube(pc_base, matrix_scale(50. / 10., 5. / 10., 30. / 10.));
    // transform_cube(pc_base, matrix_translatef(40., 60., -550., pc_base->center));
    // transform_cube(pc_base, M_W_to_C);
    // scene.objects.push_back(pc_base);

    // Object *pc_suport_bar = new Object("cylinder", 5., Vector(40., 60., -550., 1.), 25., Vector(0, 1., 0, 0), 10., Vector(179., 183., 186., 0), Vector(179., 183., 186., 0), Vector(179., 183., 186., 0));
    // transform_cylinder(pc_suport_bar, M_W_to_C);
    // scene.objects.push_back(pc_suport_bar);

    // Object *pc_screen = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(38., 38., 38., 0), Vector(38., 38., 38., 0), Vector(38., 38., 38., 0));
    // transform_cube(pc_screen, matrix_scale(150. / 10., 80. / 10., 15. / 10.));
    // transform_cube(pc_screen, matrix_translatef(40., 115., -550., pc_screen->center));
    // transform_cube(pc_screen, M_W_to_C);
    // scene.objects.push_back(pc_screen);

    // Object *pc_keyboard = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(98., 102., 111., 0), Vector(98., 102., 111., 0), Vector(98., 102., 111., 0));
    // transform_cube(pc_keyboard, matrix_scale(100. / 10., 5. / 10., 30. / 10.));
    // transform_cube(pc_keyboard, matrix_translatef(40., 60., -510., pc_keyboard->center));
    // transform_cube(pc_keyboard, M_W_to_C);
    // scene.objects.push_back(pc_keyboard);

    // Object *pc_mouse = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(98., 102., 111., 0), Vector(98., 102., 111., 0), Vector(98., 102., 111., 0));
    // transform_cube(pc_mouse, matrix_scale(10. / 10., 5. / 10., 10. / 10.));
    // transform_cube(pc_mouse, matrix_translatef(125., 60., -510., pc_mouse->center));
    // transform_cube(pc_mouse, M_W_to_C);
    // scene.objects.push_back(pc_mouse);

    // Object *pc_cpu = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(62., 70., 73., 0), Vector(62., 70., 73., 0), Vector(62., 70., 73., 0));
    // transform_cube(pc_cpu, matrix_scale(50. / 10., 80. / 10., 150. / 10.));
    // transform_cube(pc_cpu, matrix_translatef(-90., 95., -570., pc_cpu->center));
    // transform_cube(pc_cpu, M_W_to_C);
    // scene.objects.push_back(pc_cpu);

    // Object *chair_leg1 = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(218., 161., 17., 0), Vector(218., 161., 17., 0), Vector(218., 161., 17., 0));
    // transform_cube(chair_leg1, matrix_scale(15. / 10., 120. / 10., 15. / 10.));
    // transform_cube(chair_leg1, matrix_translatef(-60., -90., -430., chair_leg1->center));
    // transform_cube(chair_leg1, M_W_to_C);
    // scene.objects.push_back(chair_leg1);

    // Object *chair_leg2 = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(218., 161., 17., 0), Vector(218., 161., 17., 0), Vector(218., 161., 17., 0));
    // transform_cube(chair_leg2, matrix_scale(15. / 10., 120. / 10., 15. / 10.));
    // transform_cube(chair_leg2, matrix_translatef(60., -90., -430., chair_leg2->center));
    // transform_cube(chair_leg2, M_W_to_C);
    // scene.objects.push_back(chair_leg2);

    // Object *chair_leg3 = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(218., 161., 17., 0), Vector(218., 161., 17., 0), Vector(218., 161., 17., 0));
    // transform_cube(chair_leg3, matrix_scale(15. / 10., 120. / 10., 15. / 10.));
    // transform_cube(chair_leg3, matrix_translatef(-60., -90., -330., chair_leg3->center));
    // transform_cube(chair_leg3, M_W_to_C);
    // scene.objects.push_back(chair_leg3);

    // Object *chair_leg4 = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(218., 161., 17., 0), Vector(218., 161., 17., 0), Vector(218., 161., 17., 0));
    // transform_cube(chair_leg4, matrix_scale(15. / 10., 120. / 10., 15. / 10.));
    // transform_cube(chair_leg4, matrix_translatef(60., -90., -330., chair_leg4->center));
    // transform_cube(chair_leg4, M_W_to_C);
    // scene.objects.push_back(chair_leg4);

    // Object *chair_cover = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(180., 81., 226., 0), Vector(180., 81., 226., 0), Vector(180., 81., 226., 0));
    // transform_cube(chair_cover, matrix_scale(135. / 10., 15. / 10., 120. / 10.));
    // transform_cube(chair_cover, matrix_translatef(0., -25., -380., chair_cover->center));
    // transform_cube(chair_cover, M_W_to_C);
    // scene.objects.push_back(chair_cover);

    // Object *chair_back = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(185., 105., 60., 0), Vector(185., 105., 60., 0), Vector(185., 105., 60., 0));
    // transform_cube(chair_back, matrix_scale(135. / 10., 140. / 10., 15. / 10.));
    // transform_cube(chair_back, matrix_translatef(0., 50., -330., chair_back->center));
    // transform_cube(chair_back, M_W_to_C);
    // scene.objects.push_back(chair_back);

    // Object *comfortable_leg1 = new Object("cylinder", 10., Vector(550, -150., -580., 1.), 10., Vector(0, 1., 0, 0), 10., Vector(185., 105., 60., 0), Vector(185., 105., 60., 0), Vector(185., 105., 60., 0));
    // transform_cylinder(comfortable_leg1, M_W_to_C);
    // scene.objects.push_back(comfortable_leg1);

    // Object *comfortable_leg2 = new Object("cylinder", 10., Vector(410., -150., -580., 1.), 10., Vector(0, 1., 0, 0), 10., Vector(185., 105., 60., 0), Vector(185., 105., 60., 0), Vector(185., 105., 60., 0));
    // transform_cylinder(comfortable_leg2, M_W_to_C);
    // scene.objects.push_back(comfortable_leg2);

    // Object *comfortable_leg3 = new Object("cylinder", 10., Vector(550., -150., -480., 1.), 10., Vector(0, 1., 0, 0), 10., Vector(185., 105., 60., 0), Vector(185., 105., 60., 0), Vector(185., 105., 60., 0));
    // transform_cylinder(comfortable_leg3, M_W_to_C);
    // scene.objects.push_back(comfortable_leg3);

    // Object *comfortable_leg4 = new Object("cylinder", 10., Vector(410., -150., -480., 1.), 10., Vector(0, 1., 0, 0), 10., Vector(185., 105., 60., 0), Vector(185., 105., 60., 0), Vector(185., 105., 60., 0));
    // transform_cylinder(comfortable_leg4, M_W_to_C);
    // scene.objects.push_back(comfortable_leg4);

    // Object *comfortable = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(216., 207., 198., 0), Vector(216., 207., 198., 0), Vector(216., 207., 198., 0));
    // transform_cube(comfortable, matrix_scale(160. / 10., 150. / 10., 160. / 10.));
    // transform_cube(comfortable, matrix_translatef(480., -65., -530., comfortable->center));
    // transform_cube(comfortable, M_W_to_C);
    // scene.objects.push_back(comfortable);

    // Object *comfortable_drawer_top = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(185., 105., 60., 0), Vector(185., 105., 60., 0), Vector(185., 105., 60., 0));
    // transform_cube(comfortable_drawer_top, matrix_scale(130. / 10., 50. / 10., 10. / 10.));
    // transform_cube(comfortable_drawer_top, matrix_translatef(480., -30., -445., comfortable_drawer_top->center));
    // transform_cube(comfortable_drawer_top, M_W_to_C);
    // scene.objects.push_back(comfortable_drawer_top);

    // Object *comfortable_drawer_bottom = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(185., 105., 60., 0), Vector(185., 105., 60., 0), Vector(185., 105., 60., 0));
    // transform_cube(comfortable_drawer_bottom, matrix_scale(130. / 10., 50. / 10., 10. / 10.));
    // transform_cube(comfortable_drawer_bottom, matrix_translatef(480., -100., -445., comfortable_drawer_bottom->center));
    // transform_cube(comfortable_drawer_bottom, M_W_to_C);
    // scene.objects.push_back(comfortable_drawer_bottom);

    // Object *luminary_suport = new Object("cylinder", 40., Vector(480, 10., -530., 1.), 10., Vector(0, 1., 0, 0), 10., Vector(183., 99., 32., 0), Vector(183., 99., 32., 0), Vector(183., 99., 32., 0));
    // transform_cylinder(luminary_suport, M_W_to_C);
    // scene.objects.push_back(luminary_suport);

    // Object *luminary = new Object("sphere", 50., Vector(480, 50, -530., 1.), 10., Vector(180., 81., 226., 0), Vector(180., 81., 226., 0), Vector(180., 81., 226., 0));
    // transform_sphere(luminary, M_W_to_C);
    // scene.objects.push_back(luminary);

    Object *bed_leg1 = new Object("cylinder", 15., Vector(550, -150., -100., 1.), 25., Vector(0, 1., 0, 0), 10., Vector(26., 27., 25., 0), Vector(26., 27., 25., 0), Vector(26., 27., 25., 0));
    transform_cylinder(bed_leg1, M_W_to_C);
    scene.objects.push_back(bed_leg1);

    Object *bed_leg2 = new Object("cylinder", 15., Vector(350., -150., -100., 1.), 25., Vector(0, 1., 0, 0), 10., Vector(26., 27., 25., 0), Vector(26., 27., 25., 0), Vector(26., 27., 25., 0));
    transform_cylinder(bed_leg2, M_W_to_C);
    scene.objects.push_back(bed_leg2);

    Object *bed_leg3 = new Object("cylinder", 15., Vector(550., -150., 300., 1.), 25., Vector(0, 1., 0, 0), 10., Vector(26., 27., 25., 0), Vector(26., 27., 25., 0), Vector(26., 27., 25., 0));
    transform_cylinder(bed_leg3, M_W_to_C);
    scene.objects.push_back(bed_leg3);

    Object *bed_leg4 = new Object("cylinder", 15., Vector(350., -150., 300., 1.), 25., Vector(0, 1., 0, 0), 10., Vector(26., 27., 25., 0), Vector(26., 27., 25., 0), Vector(26., 27., 25., 0));
    transform_cylinder(bed_leg4, M_W_to_C);
    scene.objects.push_back(bed_leg4);

    Object *bed_bottom = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(26., 27., 25., 0), Vector(26., 27., 25., 0), Vector(26., 27., 25., 0));
    transform_cube(bed_bottom, matrix_scale(250. / 10., 60. / 10., 450. / 10.));
    transform_cube(bed_bottom, matrix_translatef(450., -100., 100., bed_bottom->center));
    transform_cube(bed_bottom, M_W_to_C);

    scene.objects.push_back(bed_bottom);

    Object *bed_top = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(189., 189., 189., 0), Vector(189., 189., 189., 0), Vector(189., 189., 189., 0));
    transform_cube(bed_top, matrix_scale(250. / 10., 20. / 10., 450. / 10.));
    transform_cube(bed_top, matrix_translatef(450., -60., 100., bed_top->center));
    transform_cube(bed_top, M_W_to_C);
    scene.objects.push_back(bed_top);

    Object *bed_pillow = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(250., 249., 247., 0), Vector(250., 249., 247., 0), Vector(250., 249., 247., 0));
    transform_cube(bed_pillow, matrix_scale(120. / 10., 15. / 10., 50. / 10.));
    transform_cube(bed_pillow, matrix_translatef(450., -45., 265., bed_pillow->center));
    transform_cube(bed_pillow, M_W_to_C);
    scene.objects.push_back(bed_pillow);

    // Object *closet_leg1 = new Object("cylinder", 15., Vector(-550, -150., 0., 1.), 25., Vector(0, 1., 0, 0), 10., Vector(26., 27., 25., 0), Vector(26., 27., 25., 0), Vector(26., 27., 25., 0));
    // transform_cylinder(closet_leg1, M_W_to_C);
    // scene.objects.push_back(closet_leg1);

    // Object *closet_leg2 = new Object("cylinder", 15., Vector(-400., -150., 0., 1.), 25., Vector(0, 1., 0, 0), 10., Vector(26., 27., 25., 0), Vector(26., 27., 25., 0), Vector(26., 27., 25., 0));
    // transform_cylinder(closet_leg2, M_W_to_C);
    // scene.objects.push_back(closet_leg2);

    // Object *closet_leg3 = new Object("cylinder", 15., Vector(-550., -150., 300., 1.), 25., Vector(0, 1., 0, 0), 10., Vector(26., 27., 25., 0), Vector(26., 27., 25., 0), Vector(26., 27., 25., 0));
    // transform_cylinder(closet_leg3, M_W_to_C);
    // scene.objects.push_back(closet_leg3);

    // Object *closet_leg4 = new Object("cylinder", 15., Vector(-400., -150., 300., 1.), 25., Vector(0, 1., 0, 0), 10., Vector(26., 27., 25., 0), Vector(26., 27., 25., 0), Vector(26., 27., 25., 0));
    // transform_cylinder(closet_leg4, M_W_to_C);
    // scene.objects.push_back(closet_leg4);

    // Object *closet = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(56., 56., 56., 0), Vector(56., 56., 56., 0), Vector(56., 56., 56., 0));
    // transform_cube(closet, matrix_scale(200. / 10., 400. / 10., 350. / 10.));
    // transform_cube(closet, matrix_translatef(-475., 75., 150., closet->center));
    // transform_cube(closet, M_W_to_C);
    // scene.objects.push_back(closet);

    // Object *closet_door_1 = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(192., 132., 71., 0), Vector(192., 132., 71., 0), Vector(192., 132., 71., 0));
    // transform_cube(closet_door_1, matrix_scale(10. / 10., 375. / 10., 140. / 10.));
    // transform_cube(closet_door_1, matrix_translatef(-370., 75., 225., closet_door_1->center));
    // transform_cube(closet_door_1, M_W_to_C);
    // scene.objects.push_back(closet_door_1);

    // Object *closet_door_2 = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(192., 132., 71., 0), Vector(192., 132., 71., 0), Vector(192., 132., 71., 0));
    // transform_cube(closet_door_2, matrix_scale(10. / 10., 375. / 10., 140. / 10.));
    // transform_cube(closet_door_2, matrix_translatef(-370., 75., 65., closet_door_2->center));
    // transform_cube(closet_door_2, M_W_to_C);
    // scene.objects.push_back(closet_door_2);

    Light *point_light = new Light(Vector(0.5, 0.5, 0.5, 0), transform_W_to_C(M_W_to_C, Vector(0, 1200, 0, 1.), 1), "point");
    Light *ambient_light = new Light(Vector(0.3, 0.3, 0.3, 0), Vector(0, 0, 0, 0), "ambient");
    Light *directional_light = new Light(Vector(0.0, 0.0, 0.0, 0), norm_vector(transform_W_to_C(M_W_to_C, Vector(0, 1., 0., 0), 0)), "directional");
    Light *spot_light = new Light(Vector(0.0, 0.0, 0.0, 0), transform_W_to_C(M_W_to_C, Vector(0, 800, 0, 1.), 1), norm_vector(transform_W_to_C(M_W_to_C, Vector(0, -1., 0., 0), 0)), 0.5, "spot");

    scene.lights.push_back(directional_light);
    scene.lights.push_back(point_light);
    scene.lights.push_back(ambient_light);
    // scene.lights.push_back(spot_light);

    // ofstream out("out.ppm");

    // out << "P3";
    // out << endl;
    // out << canva.w << endl;
    // out << canva.h << endl;
    // out << 255 << endl;

    unsigned char *rgba_image = new int8[(int)canva.w * (int)canva.h * CHANNEL_NUM + 1];

    auto scene_generator = [&canva, &scene, &rgba_image]()
    {
        for (int y = 0, c = 0; y < canva.w; y++) // 500 * 500 * 4
        {
            for (int x = 0; x < canva.h; x++)
            {
                Vector D = canva.canvas_to_viewport(x, y);

                pair<Color, Object *> color_and_obj;

                if (scene.projection == PERSPECT)
                {
                    color_and_obj = scene.trace_ray(Vector(0., 0., 0., 1), D, 0.0, INFINITY, y, x); // Vector(0., 0., 0., 1), D - D, Vector(0., 0., -1., 1)
                }
                else
                {
                    color_and_obj = scene.trace_ray(D, Vector(0., 0., -1., 1), 0.0, INFINITY, y, x);
                }

                rgba_image[c++] = min((int)color_and_obj.first.r, 255);
                rgba_image[c++] = min((int)color_and_obj.first.g, 255);
                rgba_image[c++] = min((int)color_and_obj.first.b, 255);
                rgba_image[c++] = 255;

                scene.point_to_obj[x][y] = color_and_obj.second;

                // mico <3
                // out << min((int)color.r, 255) << endl;
                // out << min((int)color.g, 255) << endl;
                // out << min((int)color.b, 255) << endl;
            }
        }
    };

    sf::RenderWindow window(sf::VideoMode(canva.w, canva.h), "SFML works!");

    sf::Mouse mouse;

    sf::Image img;
    sf::Texture tx;
    sf::Sprite spt;

    auto render = [scene_generator, &img, &tx, &spt, canva, rgba_image]()
    {
        scene_generator();
        img.create(canva.w, canva.h, rgba_image);
        tx.loadFromImage(img);
        spt.setTexture(tx);
        old_obj = nullptr;
    };

    render();

    auto check_picking = [&scene](int x, int y) -> Object *
    {
        if (x >= scene.canva.w || x < 0 || y >= scene.canva.h || y < 0)
        {
            return nullptr;
        }

        auto obj = scene.point_to_obj[x][y];

        if (obj != nullptr)
        {
            return obj;
        }

        return nullptr;
    };

    auto menu_object = [render, &scene](Object *obj)
    {
        system("cls");
        auto transform = switch_transform(obj->type);
        cout << "Objeto selecionado: " << obj->type << "\n\n";
        cout << "Coordenadas do centro: " << obj->center.x << ", " << obj->center.y << ", " << obj->center.z << "\n";
        cout << "K_a: " << obj->k_a.x << ", " << obj->k_a.y << ", " << obj->k_a.z << "\n";
        cout << "K_e: " << obj->k_e.x << ", " << obj->k_e.y << ", " << obj->k_e.z << "\n";
        cout << "K_d: " << obj->k_d.x << ", " << obj->k_d.y << ", " << obj->k_d.z << "\n\n";

        cout << "[1] - Mudar coordenadas do centro\n";
        cout << "[2] - Rotacionar objeto\n";
        cout << "[3] - Escalar objeto\n";
        cout << "[4] - Cisalhar objeto\n";
        cout << "[5] - Mudar propriedade ambiente\n";
        cout << "[6] - Mudar propriedade especular\n";
        cout << "[7] - Mudar propriedade difusa\n";
        cout << "[0] - Sair\n\n";
        double x, y, z;
        double xs = 0, ys = 0, zs = 0;
        string axis;
        double angle;
        Vector p;
        int op = -1;
        while (op == -1)
        {
            cout << "Digite a opcao desejada: ";
            cin >> op;
            if (op < 0 || op > 7)
                op == -1;
        }
        switch (op)
        {
        case 1: // coordenadas centro
            cout << "Digite o x: ";
            cin >> x;
            cout << "Digite o y: ";
            cin >> y;
            cout << "Digite o z: ";
            cin >> z;
            transform(obj, matrix_C_to_W(E1, A1, V1));
            transform(obj, matrix_translatef(x, y, z, obj->center));
            transform(obj, matrix_W_to_C(E1, A1, V1));
            system("cls");
            render();
            break;
        case 2: // rotacionar
            cout << "Digite o eixo no qual deseja rotacionar: ";
            cin >> axis;
            cout << "Digite o angulo (em radianos): ";
            cin >> angle;
            transform(obj, matrix_C_to_W(E1, A1, V1));
            p = obj->center;
            transform(obj, matrix_rotation(axis, angle));
            transform(obj, matrix_translatef(p.x, p.y, p.z, obj->center));
            transform(obj, matrix_W_to_C(E1, A1, V1));
            system("cls");
            render();
            break;
        case 3: // escala
            if (obj->type == "sphere")
            {
                cout << "Digite o quanto deseja multiplicar o raio da esfera: ";
                cin >> xs;
                obj->radius = obj->radius * xs;
                system("cls");
                render();
                break;
            }
            else if (obj->type == "cylinder")
            {
                cout << "Digite o quanto deseja multiplicar o raio do cilindro: ";
                cin >> xs;
                cout << "Digite o quanto deseja multiplicar na altura do cilindro: ";
                cin >> ys;
                obj->radius = obj->radius * xs;
                obj->h = obj->h * ys;
                system("cls");
                render();
                break;
            }
            else if (obj->type == "cone")
            {
                cout << "Digite o quanto deseja multiplicar o raio do cone: ";
                cin >> xs;
                cout << "Digite o quanto deseja multiplicar na altura do cone: ";
                cin >> ys;
                obj->radius = obj->radius * xs;
                obj->h = obj->h * ys;
                system("cls");
                render();
                break;
            }
            else if (obj->type == "cube")
            {
                cout << "Digite o tamanho da coordenada x: ";
                cin >> xs;
                cout << "Digite o tamanho da coordenada y: ";
                cin >> ys;
                cout << "Digite o tamanho da coordenada z: ";
                cin >> zs;
                transform(obj, matrix_C_to_W(E1, A1, V1));
                p = obj->center;
                transform(obj, matrix_scale(xs, ys, zs));
                transform(obj, matrix_translatef(p.x, p.y, p.z, obj->center));
                transform(obj, matrix_W_to_C(E1, A1, V1));
                system("cls");
                render();
                break;
            }

        case 4: // cisalhamento
            if (obj->type != "cube")
            {
                cout << "Opcao invalida para esse objeto\n";
                break;
            }
            cout << "Digite o eixo no qual deseja cisalhar: ";
            cin >> axis;
            cout << "Digite o angulo (em radianos): ";
            cin >> angle;
            transform(obj, matrix_C_to_W(E1, A1, V1));
            p = obj->center;
            transform(obj, matrix_shear(axis, angle));
            transform(obj, matrix_translatef(p.x, p.y, p.z, obj->center));
            transform(obj, matrix_W_to_C(E1, A1, V1));
            system("cls");
            render();
            break;
        case 5:
            cout << "Digite o r: ";
            cin >> x;
            cout << "Digite o g: ";
            cin >> y;
            cout << "Digite o b: ";
            cin >> z;
            obj->k_a.x = x;
            obj->k_a.y = y;
            obj->k_a.z = z;
            system("cls");
            render();
            break;
        case 6: // especular
            cout << "Digite o r: ";
            cin >> x;
            cout << "Digite o g: ";
            cin >> y;
            cout << "Digite o b: ";
            cin >> z;
            obj->k_e.x = x;
            obj->k_e.y = y;
            obj->k_e.z = z;
            system("cls");
            render();
            break;
        case 7: // difusa
            cout << "Digite o r: ";
            cin >> x;
            cout << "Digite o g: ";
            cin >> y;
            cout << "Digite o b: ";
            cin >> z;
            obj->k_d.x = x;
            obj->k_d.y = y;
            obj->k_d.z = z;
            system("cls");
            render();
            break;
        case 0:
            system("cls");
            // render();
            break;
        }
    };

    auto menu = [render, &scene, &canva, &vp, &window]()
    {
        system("cls");

        // auto transform = switch_transform(obj->type);
        // cout << "Objeto selecionado: " << obj->type << "\n\n";
        // cout << "Coordenadas do centro: " << obj->center.x << ", " << obj->center.y << ", " << obj->center.z << "\n";
        // cout << "K_a: " << obj->k_a.x << ", " << obj->k_a.y << ", " << obj->k_a.z << "\n";
        // cout << "K_e: " << obj->k_e.x << ", " << obj->k_e.y << ", " << obj->k_e.z << "\n";
        // cout << "K_d: " << obj->k_d.x << ", " << obj->k_d.y << ", " << obj->k_d.z << "\n\n";

        cout << "[1] - Mudar posicao camera\n";
        cout << "[2] - Mudar tamanho Canvas\n";
        cout << "[3] - Mudar distancia focal\n";
        cout << "[4] - Acender e apagar luzes\n";
        cout << "[5] - Mudar luzes (posicao, intensidade, direcao, angulo de abertura)\n";
        cout << "[6] - Mudar projecao perspectiva/ortogonal\n";
        cout << "[0] - Sair\n\n";
        Vector E2 = E1;
        Vector A2 = A1;
        Vector V2 = V1;
        double x, y, z;
        double xe, ye, ze;
        double xa, ya, za;
        double xv, yv, zv;
        double xj, yj;
        double dist;
        double grau;
        Vector closed_light(0., 0., 0., 0.);
        int op = -1;
        int op_perspective = -1;
        int op_light = -1;
        int op_lights1 = -1;
        int op_lights2 = -1;

        while (op == -1)
        {
            cout << "Digite a opcao desejada: ";
            cin >> op;
            if (op < 0 || op > 6)
                op == -1;
        }
        switch (op)
        {
        case 1:
            for (int i = 0; i < scene.objects.size(); i++)
            {

                auto transform = switch_transform(scene.objects[i]->type);
                transform(scene.objects[i], matrix_C_to_W(E2, A2, V2));

                cout << "Coordenadas do Olho do Observador";
                cout << "Digite o x: ";
                cin >> xe;
                cout << "Digite o y: ";
                cin >> ye;
                cout << "Digite o z: ";
                cin >> ze;

                cout << "Coordenadas do LookAt";
                cout << "Digite o x: ";
                cin >> xa;
                cout << "Digite o y: ";
                cin >> ya;
                cout << "Digite o z: ";
                cin >> za;

                cout << "Coordenadas do ViewUp";
                cout << "Digite o x: ";
                cin >> xv;
                cout << "Digite o y: ";
                cin >> yv;
                cout << "Digite o z: ";
                cin >> zv;

                E1 = Vector(xe, ye, ze, 1.);
                A1 = Vector(xa, ya, za, 1.);
                V1 = Vector(xv, yv, zv, 1.);

                transform(scene.objects[i], matrix_W_to_C(E1, A1, V1));
            }
            system("cls");
            render();
            break;
        case 2: // tamanho janela
            cout << "Digite x da janela: ";
            cin >> xj;
            cout << "Digite y da janela: ";
            cin >> yj;
            vp.w = xj;
            vp.h = yj;
            system("cls");
            render();
            break;

        case 3: // distância focal
            cout << "Digite nova distancia focal: ";
            cin >> dist;
            vp.d = dist;
            system("cls");
            render();
            break;

        case 4: // apagar luzes
            system("cls");
            cout << "Luzes acesas: \n\n";
            for (int i = 0; i < scene.lights.size(); i++)
            {
                if (scene.lights[i]->type == "point" && (scene.lights[i]->intensity.x != 0. && scene.lights[i]->intensity.y != 0. && scene.lights[i]->intensity.z != 0.))
                {
                    cout << "Pontual\n";
                }
                if (scene.lights[i]->type == "directional" && (scene.lights[i]->intensity.x != 0. && scene.lights[i]->intensity.y != 0. && scene.lights[i]->intensity.z != 0.))
                {
                    cout << "Direcional\n";
                }
                if (scene.lights[i]->type == "spot" && (scene.lights[i]->intensity.x != 0. && scene.lights[i]->intensity.y != 0. && scene.lights[i]->intensity.z != 0.))
                {
                    cout << "Spot\n";
                }
            }
            cout << "[1] - Luz pontual\n";
            cout << "[2] - Luz direcional\n\n";
            while (op_light == -1)
            {
                cout << "Escolher uma luz: ";
                cin >> op_light;
                if (op_light < 0 || op_light > 2)
                    op_light == -1;
            }
            cout << "intensidade x: ";
            cin >> x;
            cout << "intensidade y: ";
            cin >> y;
            cout << "intensidade z: ";
            cin >> z;
            // Vector l(x, y, z, 0);
            if (op_light == 1)
            {
                for (int i = 0; i < scene.lights.size(); i++)
                {
                    if (scene.lights[i]->type == "point")
                    {
                        scene.lights[i]->intensity = Vector(x, y, z, 0);
                    }
                }
            }
            else if (op_light == 2)
            {
                for (int i = 0; i < scene.lights.size(); i++)
                {
                    if (scene.lights[i]->type == "directional")
                    {
                        scene.lights[i]->intensity = Vector(x, y, z, 0);
                    }
                }
            }
            render();
            break;
        case 5: // mudar luzes
            system("cls");
            cout << "[1] - Luz pontual\n";
            cout << "[2] - Luz direcional\n";
            cout << "[3] - Luz spot\n\n";
            while (op_lights1 == -1)
            {
                cout << "Escolher uma luz: ";
                cin >> op_lights1;
                if (op_lights1 < 0 || op_lights1 > 3)
                    op_lights1 == -1;
            }

            system("cls");
            cout << "[1] - Mudar intensidade\n";
            cout << "[2] - Mudar posição\n";
            cout << "[3] - Mudar direção\n";
            cout << "[4] - Mudar ângulo de abertura\n";
            while (op_lights2 == -1)
            {
                cout << "Escolher uma opcao: ";
                cin >> op_lights2;
                if (op_lights2 < 0 || op_lights2 > 4)
                    op_lights2 == -1;
            }

            if (op_lights2 == 1)
            {
                cout << "Digite intensidade x: ";
                cin >> x;
                cout << "Digite intensidade y: ";
                cin >> y;
                cout << "Digite intensidade z: ";
                cin >> z;
                Vector l(x, y, z, 0);

                if (op_lights1 == 1)
                {
                    for (int i = 0; i < scene.lights.size(); i++)
                    {
                        if (scene.lights[i]->type == "point")
                        {
                            scene.lights[i]->intensity = l;
                        }
                    }
                }

                else if (op_lights1 == 2)
                {
                    for (int i = 0; i < scene.lights.size(); i++)
                    {
                        if (scene.lights[i]->type == "directional")
                        {
                            scene.lights[i]->intensity = l;
                        }
                    }
                }

                else if (op_lights1 == 3)
                {
                    for (int i = 0; i < scene.lights.size(); i++)
                    {
                        if (scene.lights[i]->type == "spot")
                        {
                            scene.lights[i]->intensity = l;
                        }
                    }
                }
            }

            else if (op_lights2 == 2)
            {
                cout << "Digite posicao x: ";
                cin >> x;
                cout << "Digite posicao y: ";
                cin >> y;
                cout << "Digite posicao z: ";
                cin >> z;
                Vector l(x, y, z, 0);

                if (op_lights1 == 1)
                {
                    for (int i = 0; i < scene.lights.size(); i++)
                    {
                        if (scene.lights[i]->type == "point")
                        {
                            scene.lights[i]->position = l;
                        }
                    }
                }

                else if (op_lights1 == 2)
                {
                    for (int i = 0; i < scene.lights.size(); i++)
                    {
                        if (scene.lights[i]->type == "directional")
                        {
                            scene.lights[i]->position = l;
                        }
                    }
                }

                else if (op_lights1 == 3)
                {
                    for (int i = 0; i < scene.lights.size(); i++)
                    {
                        if (scene.lights[i]->type == "spot")
                        {
                            scene.lights[i]->position = l;
                        }
                    }
                }
            }

            else if (op_lights2 == 3)
            {
                if (op_lights1 == 1)
                {
                    cout << "Opcao invalida para essa luz";
                }

                else if (op_lights1 == 2)
                {
                    cout << "Digite direcao x: ";
                    cin >> x;
                    cout << "Digite direcao y: ";
                    cin >> y;
                    cout << "Digite direcao z: ";
                    cin >> z;
                    Vector l(x, y, z, 0);
                    for (int i = 0; i < scene.lights.size(); i++)
                    {
                        if (scene.lights[i]->type == "directional")
                        {
                            scene.lights[i]->direction = l;
                        }
                    }
                }

                else if (op_lights1 == 3)
                {
                    cout << "Digite direcao x: ";
                    cin >> x;
                    cout << "Digite direcao y: ";
                    cin >> y;
                    cout << "Digite direcao z: ";
                    cin >> z;
                    Vector l(x, y, z, 0);
                    for (int i = 0; i < scene.lights.size(); i++)
                    {
                        if (scene.lights[i]->type == "spot")
                        {
                            scene.lights[i]->direction = l;
                        }
                    }
                }
            }

            else if (op_lights2 == 4)
            {
                if (op_lights1 == 1)
                {
                    cout << "Opcao invalida para essa luz";
                }

                else if (op_lights1 == 2)
                {
                    cout << "Opcao invalida para essa luz";
                }

                else if (op_lights1 == 3)
                {
                    cout << "Digite angulo de abertura: ";
                    cin >> grau;
                    for (int i = 0; i < scene.lights.size(); i++)
                    {
                        if (scene.lights[i]->type == "spot")
                        {
                            scene.lights[i]->grau = grau;
                        }
                    }
                }
            }

            render();
            break;
        case 6: // projeções
            system("cls");
            cout << "[1] - Projecao ortografica\n";
            cout << "[2] - Projecao perspectiva\n\n";
            while (op_perspective == -1)
            {
                cout << "Digite a opcao desejada: ";
                cin >> op_perspective;
                if (op_perspective < 0 || op_perspective > 2)
                    op_perspective == -1;
            }
            if (op_perspective == 1)
            {
                cout << "x do canva: ";
                cin >> x;
                cout << "y do canva: ";
                cin >> y;
                cout << "d do canva: ";
                cin >> z;
                vp.w = x;
                vp.h = y;
                vp.d = z;
                canva.update_data();
                scene.projection = ORTO;
            }
            else if (op_perspective == 2)
            {
                cout << "x do canva: ";
                cin >> x;
                cout << "y do canva: ";
                cin >> y;
                cout << "d do canva: ";
                cin >> z;
                vp.w = x;
                vp.h = y;
                vp.d = z;
                canva.update_data();
                scene.projection = PERSPECT;
            }
            system("cls");
            render();
            break;

        case 0:
            system("cls");
            // render();
            break;
        }
    };

    while (window.isOpen())
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        // CUBO

        if (mouse.sf::Mouse::isButtonPressed(sf::Mouse::Left))
        {
            auto pos = mouse.getPosition(window);
            current_object = check_picking(pos.x, pos.y);

            if (current_object != old_obj && current_object != nullptr)
            {
                old_obj = current_object;
                std::thread t(menu_object, ref(current_object));
                t.join();
            }
        }

        if (mouse.sf::Mouse::isButtonPressed(sf::Mouse::Right))
        {
            std::thread t(menu);
            t.join();
        }

        // window.clear();
        window.draw(spt);
        window.display();
    }

    // out.close();

    stbi_write_png("out.png", canva.w, canva.h, CHANNEL_NUM, rgba_image, canva.w * CHANNEL_NUM);

    stbi_image_free(rgba_image);

    return 0;
}