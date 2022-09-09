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
    Sphere(double r, Vector c, Color color)
    {
        radius = r;
        center = c;
        bg = color;
    }
    Sphere() { radius = -1; }

} Sphere;

typedef struct Scene
{
    vector<Sphere> spheres;
    Canvas canva;

    Scene(vector<Sphere> sspheres, Canvas c)
    {
        spheres = sspheres;
        canva = c;
    }
    Scene() {}

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

        return {(-b + sqrt(delta)) / 2. * a, (-b - sqrt(delta)) / 2. * a};
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

        return closest_sphere.bg;
    }
} Scene;

int main()
{
    Viewport vp(2., 2., 2.);

    Canvas canva(500., 500., vp, Color(100, 100, 100));

    vector<Sphere> spheres;

    // Sphere sphere1(1, Vector(0, 0, -(vp.d + 1)), Color(255, 0, 0));    // red
    Sphere sphere2(1, Vector(0, 0, -(vp.d + 1)), Color(0, 0, 255)); // blue
    // Sphere sphere3(1, Vector(-2.2, 0, -(vp.d + 1)), Color(0, 255, 0)); // green

    // spheres.push_back(sphere1);
    spheres.push_back(sphere2);
    // spheres.push_back(sphere3);

    Scene scene(spheres, canva);

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

            out << color.r << endl;
            out << color.g << endl;
            out << color.b << endl;
        }
    }

    out.close();

    stbi_write_png("out.png", canva.w, canva.h, CHANNEL_NUM, sphere_image, canva.w * CHANNEL_NUM);

    stbi_image_free(sphere_image);

    return 0;
}