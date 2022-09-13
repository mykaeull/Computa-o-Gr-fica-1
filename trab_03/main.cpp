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

typedef struct Sphere
{
    double radius;
    Vector center;
    double specular;
    Vector k_d;
    Vector k_e;
    Vector k_a;

    Sphere(double r, Vector c, double s, Vector K_d, Vector K_e, Vector K_a)
    {
        radius = r;
        center = c;
        specular = s;
        k_d = K_d;
        k_e = K_e;
        k_a = K_a;
    }
    Sphere() { radius = -1; }

} Sphere;

typedef struct Plane
{
    Vector p_pi;
    Vector normal;
    double specular;
    Vector k_d;
    Vector k_e;
    Vector k_a;

    Plane(Vector P_pi, Vector n, double s, Vector K_d, Vector K_e, Vector K_a)
    {
        p_pi = P_pi;
        normal = n;
        specular = s;
        k_d = K_d;
        k_e = K_e;
        k_a = K_a;
    }
    Plane() {}

} Plane;

typedef struct Object
{
    string type;
    double radius;
    Vector center;
    Vector p_pi;
    Vector normal;
    double specular;
    Vector k_d;
    Vector k_e;
    Vector k_a;

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

    Object(string object_type, Vector P_pi, Vector n, double s, Vector K_d, Vector K_e, Vector K_a) // Plane
    {
        type = object_type;
        p_pi = P_pi;
        normal = n;
        specular = s;
        k_d = K_d;
        k_e = K_e;
        k_a = K_a;
    }

    Object() { radius = -1; }

} Object;

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
        }

        return shadow;
    }

    Color trace_ray(Vector p0, Vector D, double t_min, double t_max, Color pixel_image)
    {
        double length_D = length(D);
        D = Vector(D.x / length_D, D.y / length_D, D.z / length_D);
        double closest_t_sphere = INFINITY;
        double closest_t_plane = INFINITY;
        Object closest_sphere;
        Object closest_plane;
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
        }

        if (closest_t_sphere == INFINITY && closest_t_plane == INFINITY)
        {
            return canva.bg;
        }

        if (closest_t_sphere < closest_t_plane)
        {
            double closest_t_sphere_adjusted = closest_t_sphere - EPS;
            Vector pi = Vector(p0.x + (D.x * closest_t_sphere_adjusted), p0.y + (D.y * closest_t_sphere_adjusted), p0.z + (D.z * closest_t_sphere_adjusted));
            Vector N = Vector(pi.x - closest_sphere.center.x, pi.y - closest_sphere.center.y, pi.z - closest_sphere.center.z);
            N = Vector(N.x / closest_sphere.radius, N.y / closest_sphere.radius, N.z / closest_sphere.radius);

            Vector L = Vector(lights[0].position.x - pi.x, lights[0].position.y - pi.y, lights[0].position.z - pi.z);
            double length_Pf_Pi = length(L);
            L = Vector(L.x / length_Pf_Pi, L.y / length_Pf_Pi, L.z / length_Pf_Pi);

            if (has_shadow(pi, L, length_Pf_Pi))
            {
                return Color(255 * (lights[1].intensity.x * closest_sphere.k_a.x), 255 * (lights[1].intensity.y * closest_sphere.k_a.y), 255 * (lights[1].intensity.z * closest_sphere.k_a.z));
            }

            Vector i = compute_lighting(pi, N, Vector(-D.x, -D.y, -D.z), closest_sphere);

            return Color(255 * i.x, 255 * i.y, 255 * i.z);
        }

        double closest_t_plane_adjusted = closest_t_plane - EPS;

        Vector pi = Vector(p0.x + (D.x * closest_t_plane_adjusted), p0.y + (D.y * closest_t_plane_adjusted), p0.z + (D.z * closest_t_plane_adjusted));

        Vector L = Vector(lights[0].position.x - pi.x, lights[0].position.y - pi.y, lights[0].position.z - pi.z);
        double length_Pf_Pi = length(L);
        L = Vector(L.x / length_Pf_Pi, L.y / length_Pf_Pi, L.z / length_Pf_Pi);

        if (has_shadow(pi, L, length_Pf_Pi))
        {
            return Color(255 * (lights[1].intensity.x * closest_plane.k_a.x), 255 * (lights[1].intensity.y * closest_plane.k_a.y), 255 * (lights[1].intensity.z * closest_plane.k_a.z));
        }

        Vector i = compute_lighting(pi, closest_plane.normal, Vector(-D.x, -D.y, -D.z), closest_plane);

        // return pixel_image;
        return Color(255 * i.x, 255 * i.y, 255 * i.z);
    }
} Scene;

int main()
{
    Viewport vp(60., 60., 30.);

    Canvas canva(500., 500., vp, Color(100, 100, 100));

    vector<Object> objects;
    vector<Light> lights;

    Object sphere1("sphere", 40., Vector(0, 0, -100.), 10., Vector(0.7, 0.2, 0.2), Vector(0.7, 0.2, 0.2), Vector(0.7, 0.2, 0.2)); // blue
    // Sphere sphere2("sphere", 0.4, Vector(1., 0, -1.), Color(255, 0, 0), 10., Vector(0.7, 0.2, 0.2), Vector(0.7, 0.2, 0.2), Vector(0.7, 0.2, 0.2));  // red
    // Sphere sphere3("sphere", 0.4, Vector(-1., 0, -1.), Color(0, 255, 0), 10., Vector(0.7, 0.2, 0.2), Vector(0.7, 0.2, 0.2), Vector(0.7, 0.2, 0.2)); // green

    objects.push_back(sphere1);
    // spheres.push_back(sphere2);
    // spheres.push_back(sphere3);

    Object plane1("plane", Vector(0, -40., 0), Vector(0., 1., 0.), 1., Vector(0.2, 0.7, 0.2), Vector(0, 0, 0), Vector(0.2, 0.7, 0.2));
    Object plane2("plane", Vector(0, 0, -200.), Vector(0., 0., 1.), 1., Vector(0.3, 0.3, 0.7), Vector(0, 0, 0), Vector(0.3, 0.3, 0.7));

    objects.push_back(plane1);
    objects.push_back(plane2);

    Light point_light(Vector(0.7, 0.7, 0.7), Vector(0, 60., -30.), "point");
    Light ambient_light(Vector(0.3, 0.3, 0.3), Vector(0, 0, 0), "ambient");

    lights.push_back(point_light);
    lights.push_back(ambient_light);

    Scene scene(objects, canva, lights);

    int w, h, chan;
    Color **matrix_image;

    unsigned char *rgb_image = stbi_load("lol.png", &w, &h, &chan, 3);
    matrix_image = new Color *[h];
    for (int i = 0; i < h; i++)
    {
        matrix_image[i] = new Color[w];
    }
    for (int i = 0; i < h; i++)
    {
        for (int j = 0; j < w; j++)
        {
            matrix_image[i][j].r = (double)rgb_image[j * 3 + w * i * 3];
            matrix_image[i][j].g = (double)rgb_image[j * 3 + w * i * 3 + 1];
            matrix_image[i][j].b = (double)rgb_image[j * 3 + w * i * 3 + 2];
        }
    }
    stbi_image_free(rgb_image);

    ofstream out("out.ppm");

    out << "P3";
    out << endl;
    out << canva.w << endl;
    out << canva.h << endl;
    out << 255 << endl;

    int8 *sphere_image = new int8[(int)canva.w * (int)canva.h * CHANNEL_NUM];
    for (int y = 0, c = 0; y < canva.h; y++)
    {
        for (int x = 0; x < canva.w; x++)
        {
            Vector D = canva.canvas_to_viewport(x, y);

            Color pixel_image = matrix_image[y][x];

            Color color = scene.trace_ray(Vector(0., 0., 0.), D, 0.0, INFINITY, pixel_image);

            sphere_image[c++] = min((int)color.r, 255);
            sphere_image[c++] = min((int)color.g, 255);
            sphere_image[c++] = min((int)color.b, 255);

            out << min((int)color.r, 255) << endl;
            out << min((int)color.g, 255) << endl;
            out << min((int)color.b, 255) << endl;
        }
    }

    out.close();

    stbi_write_png("out.png", canva.w, canva.h, CHANNEL_NUM, sphere_image, canva.w * CHANNEL_NUM);

    stbi_image_free(sphere_image);

    return 0;
}