#include <bits/stdc++.h>

#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "stb_image.h"
#include "stb_image_write.h"

typedef unsigned char int8;
#define CHANNEL_NUM 3

using namespace std;

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
    Color bg;
    double specular;

    Sphere(double r, Vector c, Color color, double s)
    {
        radius = r;
        center = c;
        bg = color;
        specular = s;
    }
    Sphere() { radius = -1; }

} Sphere;

typedef struct Light
{
    Vector intensity;
    Vector position;

    Light(Vector intensity_l, Vector pisition_l)
    {
        intensity = intensity_l;
        position = pisition_l;
    }

    Light() {}
} Light;

typedef struct Scene
{
    vector<Sphere> spheres;
    Canvas canva;
    Light light;

    Scene(vector<Sphere> sspheres, Canvas c, Light light_scene)
    {
        spheres = sspheres;
        canva = c;
        light = light_scene;
    }
    Scene() {}

    double dot(Vector a, Vector b)
    {
        return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
    }

    double calc_insensity(Vector v, double x, double length_a, double length_b, Sphere sphere, int flag)
    {
        if (flag == 0)
        {
            // cout << v.x * pow(x / length_a * length_b, sphere.specular) << "\n";
            return v.x * pow((x / (length_a * length_b)), sphere.specular);
        }
        return v.x * (x / length_a * length_b);
    }

    double length(Vector v)
    {
        return sqrt((v.x * v.x) + (v.y * v.y) + (v.z * v.z));
    }

    double compute_lighting(Vector pi, Vector N, Vector V, Sphere sphere)
    {
        double i = 0.0; // intensidade
        Vector L = Vector(light.position.x - pi.x, light.position.y - pi.y, light.position.z - pi.z);
        L = Vector(L.x / length(L), L.y / length(L), L.z / length(L));
        double n_dot_l = dot(N, L);

        if (n_dot_l > 0)
        {
            i += calc_insensity(light.intensity, n_dot_l, length(N), length(L), sphere, 1);
        }

        if (sphere.specular != -1)
        {
            Vector R = Vector(((2 * dot(L, N)) * N.x) - L.x, ((2 * dot(L, N)) * N.y) - L.y, ((2 * dot(L, N)) * N.z) - L.z);
            double r_dot_v = dot(R, V);
            if (r_dot_v > 0)
            {
                i += calc_insensity(light.intensity, r_dot_v, length(R), length(V), sphere, 0);
            }
        }

        return i;
    }

    pair<double, double> intersect_ray_sphere(Vector p0, Vector D, Sphere sphere)
    {
        double r = sphere.radius;

        Vector co = Vector(p0.x - sphere.center.x, p0.y - sphere.center.y, p0.z - sphere.center.z);

        double a = (D.x * D.x) + (D.y * D.y) + (D.z * D.z);
        double b = 2 * ((co.x * D.x) + (co.y * D.y) + (co.z * D.z));
        double c = ((co.x * co.x) + (co.y * co.y) + (co.z * co.z)) - (r * r);

        double delta = b * b - 4 * a * c;
        if (delta < 0)
        {
            return {INFINITY, INFINITY};
        }

        return {(-b + sqrt(delta)) / (2. * a), (-b - sqrt(delta)) / (2. * a)};
    }

    Color trace_ray(Vector p0, Vector D, double t_min, double t_max)
    {
        double closest_t = INFINITY;
        Sphere closest_sphere;
        double t1, t2;
        for (int i = 0; i < spheres.size(); i++)
        {
            tie(t1, t2) = intersect_ray_sphere(p0, D, spheres[i]);
            if ((t1 >= t_min && t1 <= t_max) && t1 < closest_t)
            {
                closest_t = t1;
                closest_sphere = spheres[i];
            }
            if ((t2 >= t_min && t2 <= t_max) && t2 < closest_t)
            {
                closest_t = t2;
                closest_sphere = spheres[i];
            }
        }

        if (closest_sphere.radius == -1)
            return canva.bg;

        Vector pi = Vector(p0.x + D.x * closest_t, p0.y + D.y * closest_t, p0.z + D.z * closest_t);
        Vector N = Vector(pi.x - closest_sphere.center.x, pi.y - closest_sphere.center.y, pi.z - closest_sphere.center.z);
        N = Vector(N.x / closest_sphere.radius, N.y / closest_sphere.radius, N.z / closest_sphere.radius);

        double i = compute_lighting(pi, N, Vector(-D.x, -D.y, -D.z), closest_sphere); // intensidade

        // cout << i << "\n";

        return Color(closest_sphere.bg.r * i, closest_sphere.bg.g * i, closest_sphere.bg.b * i);

        // return closest_sphere.bg;
    }
} Scene;

int main()
{
    Viewport vp(2., 2., 2.);

    Canvas canva(500., 500., vp, Color(100, 100, 100));

    vector<Sphere> spheres;

    // Sphere sphere1(1, Vector(0, 0, -(vp.d + 1)), Color(255, 0, 0));    // red
    Sphere sphere2(1, Vector(0, 0, -(vp.d + 1)), Color(0, 0, 255), 1000.); // blue
    // Sphere sphere3(1, Vector(-2.2, 0, -(vp.d + 1)), Color(0, 255, 0)); // green

    // spheres.push_back(sphere1);
    spheres.push_back(sphere2);
    // spheres.push_back(sphere3);

    Light point_light(Vector(0.7, 0.7, 0.7), Vector(0, 5., 0));

    Scene scene(spheres, canva, point_light);

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

            Color color = scene.trace_ray(Vector(0., 0., 0.), D, 1.0, INFINITY);

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