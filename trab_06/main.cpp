#include <bits/stdc++.h>
#include "./header/vector.hpp"
#include "./header/viewport.hpp"
#include "./header/color.hpp"
#include "./header/canvas.hpp"
#include "./header/mesh_struct.hpp"
#include "./header/objects.hpp"
#include "./header/lights.hpp"
#include "./header/scene.hpp"

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
    plane->p_pi = transform_W_to_C(M, plane->p_pi, 1);
    plane->normal = norm_vector(transform_W_to_C(M, plane->normal, 0));
}

void transform_cylinder(Object *cylinder, vector<vector<double>> M)
{
    cylinder->base = transform_W_to_C(M, cylinder->base, 1);
    cylinder->u = norm_vector(transform_W_to_C(M, cylinder->u, 0));
}

void transform_cone(Object *cone, vector<vector<double>> M)
{
    cone->base = transform_W_to_C(M, cone->base, 1);
    cone->u = norm_vector(transform_W_to_C(M, cone->u, 0));
    cone->vc = transform_W_to_C(M, cone->vc, 1);
}

void transform_sphere(Object *sphere, vector<vector<double>> M)
{
    sphere->center = transform_W_to_C(M, sphere->center, 1);
}

void transform_cube(Object *cube, vector<vector<double>> M)
{
    cube->base = transform_W_to_C(M, cube->base, 1);
    for (int i = 0; i < 8; i++)
    {
        cube->LV[i] = transform_W_to_C(M, cube->LV[i], 1);
    }
}

int main()
{
    Viewport vp(60., 60., 30.);

    Canvas canva(500., 500., vp, Color(100, 100, 100));

    vector<Object> objects;

    vector<Light> lights;

    vector<vector<double>> M_C_to_W;

    // visão boa de cima
    Vector O = Vector(0., 0., -150., 1);

    M_C_to_W = matrix_C_to_W(Vector(0., 150., 0., 1), Vector(40., 115., -550., 0), Vector(0., 1., 0., 1));

    // -----------------------------------

    // Vector O = Vector(0., 0., 0., 1);

    // M_C_to_W = matrix_C_to_W(Vector(0., 250., 0., 1), Vector(40., 115., -550., 0), Vector(0., 251., 0., 1));

    // -----------------------------------

    // Vector O = Vector(-300., -170., 800., 1);

    // M_C_to_W = matrix_C_to_W(Vector(0., 0., 0., 1), Vector(40., 115., -550., 0), Vector(0., 1., 0., 1));

    Object *down_plane = new Object("plane", Vector(0, -150, 0, 1), Vector(0., 1., 0., 0), 1., Vector(236., 219., 182., 0), Vector(236., 219., 182., 0), Vector(236., 219., 182., 0), true, "madeira5.png");
    transform_plane(down_plane, M_C_to_W);
    Object *right_plane = new Object("plane", Vector(600, -150, 0, 1), Vector(-1., 0., 0., 0), 1., Vector(236., 219., 182., 0), Vector(236., 219., 182., 0), Vector(236., 219., 182., 0), false, "");
    transform_plane(right_plane, M_C_to_W);
    Object *back_plane = new Object("plane", Vector(600, -150, -600, 1), Vector(0., 0., 1., 0), 1., Vector(236., 219., 182., 0), Vector(236., 219., 182., 0), Vector(236., 219., 182., 0), false, "");
    transform_plane(back_plane, M_C_to_W);
    Object *left_plane = new Object("plane", Vector(-600, -150, 0, 1), Vector(1., 0., 0., 0), 1., Vector(236., 219., 182., 0), Vector(236., 219., 182., 0), Vector(236., 219., 182., 0), false, "");
    transform_plane(left_plane, M_C_to_W);

    objects.push_back(*down_plane);
    objects.push_back(*right_plane);
    objects.push_back(*back_plane);
    objects.push_back(*left_plane);

    // Object *tree_trunk = new Object("cylinder", 5., Vector(0, -150, 300., 1.), 90., Vector(0, 1., 0, 0), 10., Vector(127., 80., 47., 0), Vector(127., 80., 47., 0), Vector(127., 80., 47., 0));
    // transform_cylinder(tree_trunk, M_C_to_W);
    // objects.push_back(*tree_trunk);

    // Object *tree_base = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(73., 17., 13., 0), Vector(73., 17., 13., 0), Vector(73., 17., 13., 0));
    // transform_cube(tree_base, matrix_scale(100. / 10., 30. / 10., 1.));
    // transform_cube(tree_base, matrix_translatef(0, -135, 300., tree_base->base));
    // transform_cube(tree_base, M_C_to_W);
    // objects.push_back(*tree_base);

    // Object *tree_top = new Object("cone", 75, Vector(0, -60, 300, 1.), 150, Vector(0, 1., 0, 0), 10., Vector(37., 83., 19., 0), Vector(37., 83., 19., 0), Vector(37., 83., 19., 0));
    // transform_cone(tree_top, M_C_to_W);
    // objects.push_back(*tree_top);

    // Object *tree_ball = new Object("sphere", 5., Vector(0, 95, 300., 1.), 10., Vector(247., 212., 36., 0), Vector(247., 212., 36., 0), Vector(247., 212., 36., 0));
    // transform_sphere(tree_ball, M_C_to_W);
    // objects.push_back(*tree_ball);

    // Object *door = new Object("cube", 20., 10., Vector(-400, 100., -600., 1.), Vector(128., 69., 22., 0), Vector(128., 69., 22., 0), Vector(128., 69., 22., 0));
    // transform_cube(door, matrix_scale(200. / 20., 425. / 20., 1.));
    // transform_cube(door, matrix_translatef(-500., 65., -590., door->base));
    // transform_cube(door, M_C_to_W);
    // objects.push_back(*door);

    Object *table_leg1 = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(218., 161., 17., 0), Vector(218., 161., 17., 0), Vector(218., 161., 17., 0));
    transform_cube(table_leg1, matrix_scale(15. / 10., 385. / 10., 15. / 10.));
    transform_cube(table_leg1, matrix_translatef(-165., -150., -580., table_leg1->base));
    transform_cube(table_leg1, M_C_to_W);
    objects.push_back(*table_leg1);

    Object *table_leg2 = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(218., 161., 17., 0), Vector(218., 161., 17., 0), Vector(218., 161., 17., 0));
    transform_cube(table_leg2, matrix_scale(15. / 10., 385. / 10., 15. / 10.));
    transform_cube(table_leg2, matrix_translatef(165., -150., -580., table_leg2->base));
    transform_cube(table_leg2, M_C_to_W);
    objects.push_back(*table_leg2);

    Object *table_leg3 = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(218., 161., 17., 0), Vector(218., 161., 17., 0), Vector(218., 161., 17., 0));
    transform_cube(table_leg3, matrix_scale(15. / 10., 385. / 10., 15. / 10.));
    transform_cube(table_leg3, matrix_translatef(-165., -150., -480., table_leg3->base));
    transform_cube(table_leg3, M_C_to_W);
    objects.push_back(*table_leg3);

    Object *table_leg4 = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(218., 161., 17., 0), Vector(218., 161., 17., 0), Vector(218., 161., 17., 0));
    transform_cube(table_leg4, matrix_scale(15. / 10., 385. / 10., 15. / 10.));
    transform_cube(table_leg4, matrix_translatef(165., -150., -480., table_leg4->base));
    transform_cube(table_leg4, M_C_to_W);
    objects.push_back(*table_leg4);

    Object *table_cover = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(218., 161., 17., 0), Vector(218., 161., 17., 0), Vector(218., 161., 17., 0));
    transform_cube(table_cover, matrix_scale(350. / 10., 15. / 10., 120. / 10.));
    transform_cube(table_cover, matrix_translatef(0., 50., -530., table_cover->base));
    transform_cube(table_cover, M_C_to_W);
    objects.push_back(*table_cover);

    Object *pc_base = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(179., 183., 186., 0), Vector(179., 183., 186., 0), Vector(179., 183., 186., 0));
    transform_cube(pc_base, matrix_scale(50. / 10., 5. / 10., 30. / 10.));
    transform_cube(pc_base, matrix_translatef(40., 60., -550., pc_base->base));
    transform_cube(pc_base, M_C_to_W);
    objects.push_back(*pc_base);

    Object *pc_suport_bar = new Object("cylinder", 5., Vector(40., 60., -550., 1.), 25., Vector(0, 1., 0, 0), 10., Vector(179., 183., 186., 0), Vector(179., 183., 186., 0), Vector(179., 183., 186., 0));
    transform_cylinder(pc_suport_bar, M_C_to_W);
    objects.push_back(*pc_suport_bar);

    Object *pc_screen = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(38., 38., 38., 0), Vector(38., 38., 38., 0), Vector(38., 38., 38., 0));
    transform_cube(pc_screen, matrix_scale(150. / 10., 80. / 10., 15. / 10.));
    transform_cube(pc_screen, matrix_translatef(40., 115., -550., pc_screen->base));
    transform_cube(pc_screen, M_C_to_W);
    objects.push_back(*pc_screen);

    Object *pc_keyboard = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(98., 102., 111., 0), Vector(98., 102., 111., 0), Vector(98., 102., 111., 0));
    transform_cube(pc_keyboard, matrix_scale(100. / 10., 5. / 10., 30. / 10.));
    transform_cube(pc_keyboard, matrix_translatef(40., 60., -510., pc_keyboard->base));
    transform_cube(pc_keyboard, M_C_to_W);
    objects.push_back(*pc_keyboard);

    Object *pc_mouse = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(98., 102., 111., 0), Vector(98., 102., 111., 0), Vector(98., 102., 111., 0));
    transform_cube(pc_mouse, matrix_scale(10. / 10., 5. / 10., 10. / 10.));
    transform_cube(pc_mouse, matrix_translatef(125., 60., -510., pc_mouse->base));
    transform_cube(pc_mouse, M_C_to_W);
    objects.push_back(*pc_mouse);

    Object *pc_cpu = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(62., 70., 73., 0), Vector(62., 70., 73., 0), Vector(62., 70., 73., 0));
    transform_cube(pc_cpu, matrix_scale(50. / 10., 80. / 10., 150. / 10.));
    transform_cube(pc_cpu, matrix_translatef(-90., 95., -570., pc_cpu->base));
    transform_cube(pc_cpu, M_C_to_W);
    objects.push_back(*pc_cpu);

    // Object *comfortable_leg1 = new Object("cylinder", 10., Vector(550, -150., -580., 1.), 10., Vector(0, 1., 0, 0), 10., Vector(185., 105., 60., 0), Vector(185., 105., 60., 0), Vector(185., 105., 60., 0));
    // transform_cylinder(comfortable_leg1, M_C_to_W);
    // objects.push_back(*comfortable_leg1);

    // Object *comfortable_leg2 = new Object("cylinder", 10., Vector(410., -150., -580., 1.), 10., Vector(0, 1., 0, 0), 10., Vector(185., 105., 60., 0), Vector(185., 105., 60., 0), Vector(185., 105., 60., 0));
    // transform_cylinder(comfortable_leg2, M_C_to_W);
    // objects.push_back(*comfortable_leg2);

    // Object *comfortable_leg3 = new Object("cylinder", 10., Vector(550., -150., -480., 1.), 10., Vector(0, 1., 0, 0), 10., Vector(185., 105., 60., 0), Vector(185., 105., 60., 0), Vector(185., 105., 60., 0));
    // transform_cylinder(comfortable_leg3, M_C_to_W);
    // objects.push_back(*comfortable_leg3);

    // Object *comfortable_leg4 = new Object("cylinder", 10., Vector(410., -150., -480., 1.), 10., Vector(0, 1., 0, 0), 10., Vector(185., 105., 60., 0), Vector(185., 105., 60., 0), Vector(185., 105., 60., 0));
    // transform_cylinder(comfortable_leg4, M_C_to_W);
    // objects.push_back(*comfortable_leg4);

    // Object *comfortable = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(216., 207., 198., 0), Vector(216., 207., 198., 0), Vector(216., 207., 198., 0));
    // transform_cube(comfortable, matrix_scale(160. / 10., 150. / 10., 160. / 10.));
    // transform_cube(comfortable, matrix_translatef(480., -65., -530., comfortable->base));
    // transform_cube(comfortable, M_C_to_W);
    // objects.push_back(*comfortable);

    // Object *comfortable_drawer_top = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(185., 105., 60., 0), Vector(185., 105., 60., 0), Vector(185., 105., 60., 0));
    // transform_cube(comfortable_drawer_top, matrix_scale(130. / 10., 50. / 10., 10. / 10.));
    // transform_cube(comfortable_drawer_top, matrix_translatef(480., -30., -445., comfortable_drawer_top->base));
    // transform_cube(comfortable_drawer_top, M_C_to_W);
    // objects.push_back(*comfortable_drawer_top);

    // Object *comfortable_drawer_bottom = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(185., 105., 60., 0), Vector(185., 105., 60., 0), Vector(185., 105., 60., 0));
    // transform_cube(comfortable_drawer_bottom, matrix_scale(130. / 10., 50. / 10., 10. / 10.));
    // transform_cube(comfortable_drawer_bottom, matrix_translatef(480., -100., -445., comfortable_drawer_bottom->base));
    // transform_cube(comfortable_drawer_bottom, M_C_to_W);
    // objects.push_back(*comfortable_drawer_bottom);

    // Object *luminary_suport = new Object("cylinder", 40., Vector(480, 10., -530., 1.), 10., Vector(0, 1., 0, 0), 10., Vector(183., 99., 32., 0), Vector(183., 99., 32., 0), Vector(183., 99., 32., 0));
    // transform_cylinder(luminary_suport, M_C_to_W);
    // objects.push_back(*luminary_suport);

    // Object *luminary = new Object("sphere", 50., Vector(480, 50, -530., 1.), 10., Vector(180., 81., 226., 0), Vector(180., 81., 226., 0), Vector(180., 81., 226., 0));
    // transform_sphere(luminary, M_C_to_W);
    // objects.push_back(*luminary);

    // Object *bed_leg1 = new Object("cylinder", 15., Vector(550, -150., -100., 1.), 25., Vector(0, 1., 0, 0), 10., Vector(26., 27., 25., 0), Vector(26., 27., 25., 0), Vector(26., 27., 25., 0));
    // transform_cylinder(bed_leg1, M_C_to_W);
    // objects.push_back(*bed_leg1);

    // Object *bed_leg2 = new Object("cylinder", 15., Vector(350., -150., -100., 1.), 25., Vector(0, 1., 0, 0), 10., Vector(26., 27., 25., 0), Vector(26., 27., 25., 0), Vector(26., 27., 25., 0));
    // transform_cylinder(bed_leg2, M_C_to_W);
    // objects.push_back(*bed_leg2);

    // Object *bed_leg3 = new Object("cylinder", 15., Vector(550., -150., 300., 1.), 25., Vector(0, 1., 0, 0), 10., Vector(26., 27., 25., 0), Vector(26., 27., 25., 0), Vector(26., 27., 25., 0));
    // transform_cylinder(bed_leg3, M_C_to_W);
    // objects.push_back(*bed_leg3);

    // Object *bed_leg4 = new Object("cylinder", 15., Vector(350., -150., 300., 1.), 25., Vector(0, 1., 0, 0), 10., Vector(26., 27., 25., 0), Vector(26., 27., 25., 0), Vector(26., 27., 25., 0));
    // transform_cylinder(bed_leg4, M_C_to_W);
    // objects.push_back(*bed_leg4);

    // Object *bed_bottom = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(26., 27., 25., 0), Vector(26., 27., 25., 0), Vector(26., 27., 25., 0));
    // transform_cube(bed_bottom, matrix_scale(250. / 10., 60. / 10., 450. / 10.));
    // transform_cube(bed_bottom, matrix_translatef(450., -100., 100., bed_bottom->base));
    // transform_cube(bed_bottom, M_C_to_W);
    // objects.push_back(*bed_bottom);

    // Object *bed_top = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(189., 189., 189., 0), Vector(189., 189., 189., 0), Vector(189., 189., 189., 0));
    // transform_cube(bed_top, matrix_scale(250. / 10., 20. / 10., 450. / 10.));
    // transform_cube(bed_top, matrix_translatef(450., -60., 100., bed_top->base));
    // transform_cube(bed_top, M_C_to_W);
    // objects.push_back(*bed_top);

    // Object *bed_pillow = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(250., 249., 247., 0), Vector(250., 249., 247., 0), Vector(250., 249., 247., 0));
    // transform_cube(bed_pillow, matrix_scale(100. / 10., 15. / 10., 80. / 10.));
    // transform_cube(bed_pillow, matrix_translatef(450., -45., 250., bed_pillow->base));
    // transform_cube(bed_pillow, M_C_to_W);
    // objects.push_back(*bed_pillow);

    // Object *closet_leg1 = new Object("cylinder", 15., Vector(-550, -150., 0., 1.), 25., Vector(0, 1., 0, 0), 10., Vector(26., 27., 25., 0), Vector(26., 27., 25., 0), Vector(26., 27., 25., 0));
    // transform_cylinder(closet_leg1, M_C_to_W);
    // objects.push_back(*closet_leg1);

    // Object *closet_leg2 = new Object("cylinder", 15., Vector(-400., -150., 0., 1.), 25., Vector(0, 1., 0, 0), 10., Vector(26., 27., 25., 0), Vector(26., 27., 25., 0), Vector(26., 27., 25., 0));
    // transform_cylinder(closet_leg2, M_C_to_W);
    // objects.push_back(*closet_leg2);

    // Object *closet_leg3 = new Object("cylinder", 15., Vector(-550., -150., 300., 1.), 25., Vector(0, 1., 0, 0), 10., Vector(26., 27., 25., 0), Vector(26., 27., 25., 0), Vector(26., 27., 25., 0));
    // transform_cylinder(closet_leg3, M_C_to_W);
    // objects.push_back(*closet_leg3);

    // Object *closet_leg4 = new Object("cylinder", 15., Vector(-400., -150., 300., 1.), 25., Vector(0, 1., 0, 0), 10., Vector(26., 27., 25., 0), Vector(26., 27., 25., 0), Vector(26., 27., 25., 0));
    // transform_cylinder(closet_leg4, M_C_to_W);
    // objects.push_back(*closet_leg4);

    // Object *closet = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(56., 56., 56., 0), Vector(56., 56., 56., 0), Vector(56., 56., 56., 0));
    // transform_cube(closet, matrix_scale(200. / 10., 400. / 10., 350. / 10.));
    // transform_cube(closet, matrix_translatef(-475., 75., 150., closet->base));
    // transform_cube(closet, M_C_to_W);
    // objects.push_back(*closet);

    // Object *closet_door_1 = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(192., 132., 71., 0), Vector(192., 132., 71., 0), Vector(192., 132., 71., 0));
    // transform_cube(closet_door_1, matrix_scale(10. / 10., 375. / 10., 140. / 10.));
    // transform_cube(closet_door_1, matrix_translatef(-370., 75., 225., closet_door_1->base));
    // transform_cube(closet_door_1, M_C_to_W);
    // objects.push_back(*closet_door_1);

    // Object *closet_door_2 = new Object("cube", 10., 10., Vector(0., 0., 0., 1.), Vector(192., 132., 71., 0), Vector(192., 132., 71., 0), Vector(192., 132., 71., 0));
    // transform_cube(closet_door_2, matrix_scale(10. / 10., 375. / 10., 140. / 10.));
    // transform_cube(closet_door_2, matrix_translatef(-370., 75., 65., closet_door_2->base));
    // transform_cube(closet_door_2, M_C_to_W);
    // objects.push_back(*closet_door_2);

    Light *point_light = new Light(Vector(0.5, 0.5, 0.5, 0), transform_W_to_C(M_C_to_W, Vector(0, 1200, 0, 1.), 1), "point");
    Light *ambient_light = new Light(Vector(0.3, 0.3, 0.3, 0), Vector(0, 0, 0, 0), "ambient");
    Light *directional_light = new Light(Vector(0.5, 0.5, 0.5, 0), norm_vector(transform_W_to_C(M_C_to_W, Vector(0, -150., -165., 0), 0)), "directional");
    Light *spot_light = new Light(Vector(0.7, 0.7, 0.7, 0), transform_W_to_C(M_C_to_W, Vector(0, 140, 0, 1.), 1), norm_vector(transform_W_to_C(M_C_to_W, Vector(0, -150., -165., 0), 0)), 0.5, "spot");

    lights.push_back(*point_light);
    lights.push_back(*ambient_light);
    // lights.push_back(*directional_light);
    // lights.push_back(*spot_light);

    Scene scene(objects, canva, lights);

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

            Color color = scene.trace_ray(O, D, 0.0, INFINITY, y, x);

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