#include <bits/stdc++.h>

#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "stb_image.h"
#include "stb_image_write.h"

typedef unsigned char int8;
#define CHANNEL_NUM 3

using namespace std;

int teste = 0;

typedef struct Vector
{
    double x, y, z;

    Vector(double xx, double yy, double zz)
    {
        x = xx, y = yy, z = zz;
    }
    Vector() {}
} Vector;

typedef struct Viewport
{
    double w, h, d;

    Viewport(double vw, double vh, double vd)
    {
        w = vw, h = vh, d = vd;
    }
    Viewport() {}
} Viewport;

typedef struct Color
{
    double r, g, b;

    Color(double rr, double gg, double bb)
    {
        r = rr, g = gg, b = bb;
    }
    Color() {}
} Color;

typedef struct Canvas
{
    double w, h, dx, dy;
    Viewport vp;
    Color bg;

    Canvas(double wc, double hc, Viewport vpc, Color bgc)
    {
        w = wc, h = hc;
        dx = vpc.w / wc, dy = vpc.h / hc;
        vp = vpc;
        bg = bgc;
    }
    Canvas() {}

    Vector canvas_to_viewport(double x, double y)
    {
        return Vector(-vp.w / 2. + dx / 2. + dx * x, vp.h / 2. - dy / 2. - dy * y, -vp.d);
    }
} Canvas;

typedef struct Object
{
    string type;
    double radius;
    Vector center;
    Vector p_pi;
    Vector normal;
    Vector base;
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
        vc = Vector(b.x + (height * uu.x), b.y + (height * uu.y), b.z + (height * uu.z)); // case cone
    }

    Object() { radius = -1; }

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

typedef struct Light
{
    Vector intensity;
    Vector position;
    string type;

    Light(Vector intensity_l, Vector pisition_l, string type_l)
    {
        intensity = intensity_l;
        position = pisition_l;
        type = type_l;
    }

    Light() {}
} Light;

typedef struct Scene
{
    vector<Object> objects;
    Canvas canva;
    vector<Light> lights;

    Scene(vector<Object> objects_scene, Canvas c, vector<Light> lights_scene)
    {
        objects = objects_scene;
        canva = c;
        lights = lights_scene;
    }
    Scene() {}

    double dot(Vector a, Vector b)
    {
        return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
    }

    Vector vector_mult(Vector a, Vector b)
    {
        double row1 = (a.y * b.z) - (a.z * b.y);
        double row2 = (a.z * b.x) - (a.x * b.z);
        double row3 = (a.x * b.y) - (a.y * b.x);
        return Vector(row1, row2, row3);
    }

    double calc_insensity(double i, double x, double length_a, double length_b, double specular, double k, int flag)
    {
        if (flag == 0)
        {
            return (i * k) * pow((x / (length_a * length_b)), specular);
        }
        return (i * k) * (x / length_a * length_b);
    }

    double length(Vector v)
    {
        return sqrt((v.x * v.x) + (v.y * v.y) + (v.z * v.z));
    }

    Vector compute_lighting(Vector pi, Vector N, Vector V, Object object)
    {
        Vector i(0.0, 0.0, 0.0);

        for (int j = 0; j < lights.size(); j++)
        {
            if (lights[j].type == "ambient")
            {
                i.x += lights[j].intensity.x * object.k_a.x;
                i.y += lights[j].intensity.y * object.k_a.y;
                i.z += lights[j].intensity.z * object.k_a.z;
            }
            else
            {
                Vector L = Vector(lights[j].position.x - pi.x, lights[j].position.y - pi.y, lights[j].position.z - pi.z);
                double Pf_Pi = length(L);
                L = Vector(L.x / Pf_Pi, L.y / Pf_Pi, L.z / Pf_Pi);
                double n_dot_l = dot(N, L);

                // DIFUSA
                if (n_dot_l > 0)
                {
                    i.x += calc_insensity(lights[j].intensity.x, n_dot_l, length(N), length(L), object.specular, object.k_d.x, 1);
                    i.y += calc_insensity(lights[j].intensity.y, n_dot_l, length(N), length(L), object.specular, object.k_d.y, 1);
                    i.z += calc_insensity(lights[j].intensity.z, n_dot_l, length(N), length(L), object.specular, object.k_d.z, 1);
                }

                // ESPECULAR
                if (object.specular != -1)
                {
                    Vector R = Vector(((2 * dot(L, N)) * N.x) - L.x, ((2 * dot(L, N)) * N.y) - L.y, ((2 * dot(L, N)) * N.z) - L.z);
                    double r_dot_v = dot(R, V);
                    if (r_dot_v > 0)
                    {
                        i.x += calc_insensity(lights[j].intensity.x, r_dot_v, length(R), length(V), object.specular, object.k_e.x, 0);
                        i.y += calc_insensity(lights[j].intensity.y, r_dot_v, length(R), length(V), object.specular, object.k_e.y, 0);
                        i.z += calc_insensity(lights[j].intensity.z, r_dot_v, length(R), length(V), object.specular, object.k_e.z, 0);
                    }
                }
            }
        }

        return i;
    }

    pair<double, double> intersect_ray_sphere(Vector p0, Vector D, Object sphere)
    {
        double r = sphere.radius;

        Vector w = Vector(p0.x - sphere.center.x, p0.y - sphere.center.y, p0.z - sphere.center.z);

        double a = (D.x * D.x) + (D.y * D.y) + (D.z * D.z);
        double b = 2 * ((w.x * D.x) + (w.y * D.y) + (w.z * D.z));
        double c = ((w.x * w.x) + (w.y * w.y) + (w.z * w.z)) - (r * r);

        double delta = b * b - 4 * a * c;
        if (delta < 0)
        {
            return {INFINITY, INFINITY};
        }

        return {(-b + sqrt(delta)) / (2. * a), (-b - sqrt(delta)) / (2. * a)};
    }

    double intersect_ray_plane(Vector p0, Vector D, Object plane)
    {
        Vector w = Vector(p0.x - plane.p_pi.x, p0.y - plane.p_pi.y, p0.z - plane.p_pi.z);
        Vector normal = plane.normal;

        double num = ((w.x * normal.x) + (w.y * normal.y) + (w.z * normal.z));
        double den = ((D.x * normal.x) + (D.y * normal.y) + (D.z * normal.z));

        double ti = -num / den;

        if (ti < 0)
        {
            return INFINITY;
        }

        return ti;
    }

    bool is_in_shell(Vector Pi, Object object)
    {
        Vector Pi_sub_base = Vector(Pi.x - object.base.x, Pi.y - object.base.y, Pi.z - object.base.z);
        return (dot(Pi_sub_base, object.u) < 0 || dot(Pi_sub_base, object.u) > object.h);
    }

    bool is_in_base(Vector P, Object plane, Object object)
    {
        Vector CP = Vector(P.x - plane.p_pi.x, P.y - plane.p_pi.y, P.z - plane.p_pi.z);
        double CP_length = length(CP);
        double cp_x_n = dot(CP, object.u);
        bool is_zero = (cp_x_n > 0.0) && (cp_x_n < 0.0);
        bool is_radius = CP_length <= object.radius;
        return (is_zero && is_radius);
    }

    double intersect_ray_base(Vector p0, Vector D, Object plane, Object object) // cylinder or cone
    {
        double t;
        t = intersect_ray_plane(p0, D, plane);
        Vector P = Vector(p0.x + (t * D.x), p0.y + (t * D.y), p0.z + (t * D.z));
        return is_in_base(P, plane, object) ? t : INFINITY;
    }

    pair<double, double> intersect_ray_cylinder(Vector p0, Vector D, Object cylinder)
    {
        double r = cylinder.radius;
        double t1, t2;

        Vector p0_sub_base = Vector(p0.x - cylinder.base.x, p0.y - cylinder.base.y, p0.z - cylinder.base.z);

        Vector v = Vector(p0_sub_base.x - (dot(p0_sub_base, cylinder.u)) * cylinder.u.x, p0_sub_base.y - (dot(p0_sub_base, cylinder.u)) * cylinder.u.y, p0_sub_base.z - (dot(p0_sub_base, cylinder.u)) * cylinder.u.z);
        Vector w = Vector(D.x - ((dot(D, cylinder.u)) * cylinder.u.x), D.y - ((dot(D, cylinder.u)) * cylinder.u.y), D.z - ((dot(D, cylinder.u)) * cylinder.u.z));

        double a = dot(w, w);
        double b = 2 * dot(v, w);
        double c = dot(v, v) - (r * r);

        double delta = b * b - 4 * a * c;
        if (a == 0 || delta < 0)
        {
            return {INFINITY, INFINITY};
        }

        t1 = (-b + sqrt(delta)) / (2. * a);
        t2 = (-b - sqrt(delta)) / (2. * a);

        Vector Pi1 = Vector(p0.x + (D.x * t1), p0.y + (D.y * t1), p0.z + (D.z * t1));
        t1 = is_in_shell(Pi1, cylinder) ? INFINITY : t1;
        Vector Pi2 = Vector(p0.x + (D.x * t2), p0.y + (D.y * t2), p0.z + (D.z * t2));
        t2 = is_in_shell(Pi2, cylinder) ? INFINITY : t2;

        return {t1, t2};
    }

    pair<double, double> intersect_ray_cone(Vector p0, Vector D, Object cone)
    {
        double r = cone.radius;
        double h = cone.h;
        double cos_sqr_teta = (h * h) / ((r * r) + (h * h));
        double t1, t2;

        Vector v = Vector(cone.vc.x - p0.x, cone.vc.y - p0.y, cone.vc.z - p0.z);

        double a = (dot(D, cone.u) * dot(D, cone.u)) - dot(D, D) * cos_sqr_teta;
        double b = (dot(v, D) * cos_sqr_teta - dot(v, cone.u) * dot(D, cone.u)) * 2;
        double c = (dot(v, cone.u) * dot(v, cone.u)) - dot(v, v) * cos_sqr_teta;

        double delta = b * b - 4 * a * c;

        if (delta < 0)
        {
            return {INFINITY, INFINITY};
        }

        t1 = (-b + sqrt(delta)) / (2 * a);
        t2 = (-b - sqrt(delta)) / (2 * a);

        Vector Pi1 = Vector(p0.x + (D.x * t1), p0.y + (D.y * t1), p0.z + (D.z * t1));
        t1 = is_in_shell(Pi1, cone) ? INFINITY : t1;
        Vector Pi2 = Vector(p0.x + (D.x * t2), p0.y + (D.y * t2), p0.z + (D.z * t2));
        t2 = is_in_shell(Pi2, cone) ? INFINITY : t2;

        return {t1, t2};
    }

    bool has_shadow(Vector pi, Vector L, double length_Pf_Pi)
    {
        bool shadow = false;
        double s1, s2;

        for (int i = 0; i < objects.size(); i++)
        {
            if (objects[i].type == "sphere")
            {
                tie(s1, s2) = intersect_ray_sphere(pi, L, objects[i]);
                if (s1 > 0 && s1 < length_Pf_Pi)
                {
                    shadow = true;
                }
                if (s2 > 0 && s2 < length_Pf_Pi)
                {
                    shadow = true;
                }
            }
            if (objects[i].type == "plane")
            {
                s1 = intersect_ray_plane(pi, L, objects[i]);
                if (s1 > 0 && s1 < length_Pf_Pi)
                {
                    shadow = true;
                }
            }
            if (objects[i].type == "cylinder")
            {
                tie(s1, s2) = intersect_ray_cylinder(pi, L, objects[i]);
                if (s1 > 0 && s1 < length_Pf_Pi)
                {
                    shadow = true;
                }
                if (s2 > 0 && s2 < length_Pf_Pi)
                {
                    shadow = true;
                }
            }
            if (objects[i].type == "cone")
            {
                tie(s1, s2) = intersect_ray_cone(pi, L, objects[i]);
                if (s1 > 0 && s1 < length_Pf_Pi)
                {
                    shadow = true;
                }
                if (s2 > 0 && s2 < length_Pf_Pi)
                {
                    shadow = true;
                }
            }
        }

        return shadow;
    }

    Color define_color(double closest_t, Object object, Vector p0, Vector D)
    {
        Vector pi = Vector(p0.x + (D.x * closest_t), p0.y + (D.y * closest_t), p0.z + (D.z * closest_t));
        Vector N;
        if (object.type == "sphere")
        {
            N = Vector(pi.x - object.center.x, pi.y - object.center.y, pi.z - object.center.z);
            N = Vector(N.x / object.radius, N.y / object.radius, N.z / object.radius);
        }
        if (object.type == "plane")
        {
            N = object.normal;
        }
        if (object.type == "cylinder")
        {
            Vector CP = Vector(pi.x - object.base.x, pi.y - object.base.y, pi.z - object.base.z);
            double CP_dot_u = dot(CP, object.u);
            Vector AP = Vector(CP.x - (object.u.x * (CP_dot_u)), CP.y - (object.u.y * (CP_dot_u)), CP.z - (object.u.z * (CP_dot_u)));
            double AP_length = length(AP);
            N = Vector(AP.x / AP_length, AP.y / AP_length, AP.z / AP_length);
        }
        if (object.type == "cone")
        {
            Vector w = Vector(object.vc.x - pi.x, object.vc.y - pi.y, object.vc.z - pi.z);
            Vector n_barra = vector_mult(w, object.u);
            N = vector_mult(n_barra, w);
            double N_length = length(N);
            N = Vector(N.x / N_length, N.y / N_length, N.z / N_length);
        }

        Vector L = Vector(lights[0].position.x - pi.x, lights[0].position.y - pi.y, lights[0].position.z - pi.z);
        double length_Pf_Pi = length(L);
        L = Vector(L.x / length_Pf_Pi, L.y / length_Pf_Pi, L.z / length_Pf_Pi);

        if (has_shadow(pi, L, length_Pf_Pi))
        {
            return Color(255 * (lights[1].intensity.x * object.k_a.x), 255 * (lights[1].intensity.y * object.k_a.y), 255 * (lights[1].intensity.z * object.k_a.z));
        }

        Vector i = compute_lighting(pi, N, Vector(-D.x, -D.y, -D.z), object);

        return Color(255 * i.x, 255 * i.y, 255 * i.z);
    }

    Color trace_ray(Vector p0, Vector D, double t_min, double t_max)
    {
        double length_D = length(D);
        D = Vector(D.x / length_D, D.y / length_D, D.z / length_D);
        double closest_t_sphere = INFINITY;
        double closest_t_plane = INFINITY;
        double closest_t_cylinder = INFINITY;
        double closest_t_cone = INFINITY;
        Object closest_sphere;
        Object closest_plane;
        Object closest_cylinder;
        Object closest_cone;
        double EPS = 0.01;
        double t1, t2;
        for (int i = 0; i < objects.size(); i++)
        {
            if (objects[i].type == "sphere")
            {
                tie(t1, t2) = intersect_ray_sphere(p0, D, objects[i]);
                if ((t1 > t_min && t1 < t_max) && t1 < closest_t_sphere)
                {
                    closest_t_sphere = t1;
                    closest_sphere = objects[i];
                }
                if ((t2 > t_min && t2 < t_max) && t2 < closest_t_sphere)
                {
                    closest_t_sphere = t2;
                    closest_sphere = objects[i];
                }
            }
            if (objects[i].type == "plane")
            {
                t1 = intersect_ray_plane(p0, D, objects[i]);
                if ((t1 > t_min && t1 < t_max) && t1 < closest_t_plane)
                {
                    closest_t_plane = t1;
                    closest_plane = objects[i];
                }
            }
            if (objects[i].type == "cylinder")
            {
                tie(t1, t2) = intersect_ray_cylinder(p0, D, objects[i]);
                if ((t1 > t_min && t1 < t_max) && t1 < closest_t_cylinder)
                {
                    closest_t_cylinder = t1;
                    closest_cylinder = objects[i];
                }
                if ((t2 > t_min && t2 < t_max) && t2 < closest_t_cylinder)
                {
                    closest_t_cylinder = t2;
                    closest_cylinder = objects[i];
                }
                Vector N_base = Vector(-closest_cylinder.u.x, -closest_cylinder.u.y, -closest_cylinder.u.z);
                Vector P_pi_base = closest_cylinder.base;
                Object plane_base("plane", P_pi_base, N_base, objects[i].specular, objects[i].k_d, objects[i].k_e, objects[i].k_a, "base");
                Vector N_cover = closest_cylinder.u;
                Vector P_pi_cover = Vector(closest_cylinder.base.x + (closest_cylinder.u.x * closest_cylinder.h), closest_cylinder.base.y + (closest_cylinder.u.y * closest_cylinder.h), closest_cylinder.base.z + (closest_cylinder.u.z * closest_cylinder.h));
                Object plane_top("plane", P_pi_cover, N_cover, objects[i].specular, objects[i].k_d, objects[i].k_e, objects[i].k_a, "cover");
                t1 = intersect_ray_base(p0, D, plane_base, closest_cylinder);
                t2 = intersect_ray_base(p0, D, plane_top, closest_cylinder);

                if ((t1 > t_min && t1 < t_max) && t1 < closest_t_cylinder)
                {
                    closest_t_cylinder = t1;
                    closest_cylinder = plane_base;
                }
                if ((t2 > t_min && t2 < t_max) && t2 < closest_t_cylinder)
                {
                    closest_t_cylinder = t2;
                    closest_cylinder = plane_top;
                }
            }
            if (objects[i].type == "cone")
            {
                tie(t1, t2) = intersect_ray_cone(p0, D, objects[i]);
                if ((t1 > t_min && t1 < t_max) && t1 < closest_t_cone)
                {
                    closest_t_cone = t1;
                    closest_cone = objects[i];
                }
                if ((t2 > t_min && t2 < t_max) && t2 < closest_t_cone)
                {
                    closest_t_cone = t2;
                    closest_cone = objects[i];
                }
                Vector N_base = Vector(closest_cone.u.x, closest_cone.u.y, closest_cone.u.z);
                Vector P_pi_base = closest_cone.base;
                Object plane_base("plane", P_pi_base, N_base, objects[i].specular, objects[i].k_d, objects[i].k_e, objects[i].k_a, "base");
                t1 = intersect_ray_base(p0, D, plane_base, closest_cone);
                if ((t1 > t_min && t1 < t_max) && t1 < closest_t_cone)
                {
                    closest_t_cone = t1;
                    closest_cone = plane_base;
                }
            }
        }

        if (closest_t_sphere == INFINITY && closest_t_plane == INFINITY && closest_t_cylinder == INFINITY && closest_t_cone == INFINITY)
        {
            return canva.bg;
        }

        vector<Closest_Object> closest_objects;
        closest_objects.push_back(Closest_Object(closest_t_sphere, closest_sphere));
        closest_objects.push_back(Closest_Object(closest_t_plane, closest_plane));
        closest_objects.push_back(Closest_Object(closest_t_cylinder, closest_cylinder));
        closest_objects.push_back(Closest_Object(closest_t_cone, closest_cone));

        double smaller = INFINITY;
        double closest_t;
        Object closest_object;
        for (int i = 0; i < closest_objects.size(); i++) // pega o objeto com menor t (objeto mais próximo)
        {
            if (closest_objects[i].closest_t < smaller)
            {
                smaller = closest_objects[i].closest_t;
                closest_t = closest_objects[i].closest_t;
                closest_object = closest_objects[i].object;
            }
        }

        return define_color(closest_t - EPS, closest_object, p0, D);
    }
} Scene;

int main()
{
    Viewport vp(60., 60., 30.);

    Canvas canva(500., 500., vp, Color(100, 100, 100));

    vector<Object> objects;
    vector<Light> lights;

    double raio = 40.;

    Object plane1("plane", Vector(0, -150, 0), Vector(0., 1., 0.), 1., Vector(0.933, 0.933, 0.933), Vector(0.933, 0.933, 0.933), Vector(0.933, 0.933, 0.933), "floor");
    Object plane2("plane", Vector(200, -150, 0), Vector(-1., 0., 0.), 1., Vector(0.686, 0.933, 0.933), Vector(0.686, 0.933, 0.933), Vector(0.686, 0.933, 0.933), "background");
    Object plane3("plane", Vector(200, -150, -400), Vector(0., 0., 1.), 1., Vector(0.686, 0.933, 0.933), Vector(0.686, 0.933, 0.933), Vector(0.686, 0.933, 0.933), "floor");
    Object plane4("plane", Vector(-200, -150, 0), Vector(1., 0., 0.), 1., Vector(0.686, 0.933, 0.933), Vector(0.686, 0.933, 0.933), Vector(0.686, 0.933, 0.933), "floor");
    Object plane5("plane", Vector(0, 150, 0), Vector(0., -1., 0.), 1., Vector(0.933, 0.933, 0.933), Vector(0.933, 0.933, 0.933), Vector(0.933, 0.933, 0.933), "floor");

    objects.push_back(plane1);
    objects.push_back(plane2);
    objects.push_back(plane3);
    objects.push_back(plane4);
    objects.push_back(plane5);

    Object cylinder1("cylinder", 5., Vector(0, -150, -200.), 90., Vector(0, 1., 0), 10., Vector(0.824, 0.706, 0.549), Vector(0.824, 0.706, 0.549), Vector(0.824, 0.706, 0.549));
    objects.push_back(cylinder1);

    Object cone1("cone", 90, Vector(0, -60, -200), 150, Vector(0, 1., 0), 10., Vector(0., 1., 0.498), Vector(0., 1., 0.498), Vector(0., 1., 0.498));
    objects.push_back(cone1);

    Object sphere1("sphere", 5., Vector(0, 95, -200.), 10., Vector(0.854, 0.647, 0.125), Vector(0.854, 0.647, 0.125), Vector(0.854, 0.647, 0.125));
    objects.push_back(sphere1);

    Light point_light(Vector(0.7, 0.7, 0.7), Vector(-100, 140., -20.), "point");
    Light ambient_light(Vector(0.3, 0.3, 0.3), Vector(0, 0, 0), "ambient");

    lights.push_back(point_light);
    lights.push_back(ambient_light);

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

            Color color = scene.trace_ray(Vector(0., 0., 0.), D, 0.0, INFINITY);

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