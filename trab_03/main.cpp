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
    Vector k_d;
    Vector k_e;
    Vector k_a;

    Sphere(double r, Vector c, Color color, double s, Vector K_d, Vector K_e, Vector K_a)
    {
        radius = r;
        center = c;
        bg = color;
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
    Color bg;
    double specular;
    Vector k_d;
    Vector k_e;
    Vector k_a;

    Plane(Vector P_pi, Vector n, Color color, double s, Vector K_d, Vector K_e, Vector K_a)
    {
        p_pi = P_pi;
        normal = n;
        bg = color;
        specular = s;
        k_d = K_d;
        k_e = K_e;
        k_a = K_a;
    }
    Plane() {}

} Plane;

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
    vector<Sphere> spheres;
    vector<Plane> planes;
    Canvas canva;
    vector<Light> lights;

    Scene(vector<Sphere> spheres_scene, vector<Plane> planes_scene, Canvas c, vector<Light> lights_scene)
    {
        spheres = spheres_scene;
        planes = planes_scene;
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

    Vector compute_lighting(Vector pi, Vector N, Vector V, auto object) // NAO SEI SE FUNCIONA (SE FUNCIONAR, TENTAR vector<auto> objetos)
    {
        // double i = 0.0; // intensidade
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
                L = Vector(L.x / length(L), L.y / length(L), L.z / length(L));
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

    pair<double, double> intersect_ray_sphere(Vector p0, Vector D, Sphere sphere)
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

    double intersect_ray_plane(Vector p0, Vector D, Plane plane)
    {
        Vector w = Vector(p0.x - plane.p_pi.x, p0.y - plane.p_pi.y, p0.z - plane.p_pi.z);
        Vector normal = plane.normal;

        double ti = -((w.x * normal.x) + (w.y * normal.y) + (w.z * normal.z)) / ((D.x * normal.x) + (D.y * normal.y) + (D.z * normal.z));

        if (ti < 0)
        {
            return INFINITY;
        }

        return ti;
    }

    Color trace_ray(Vector p0, Vector D, double t_min, double t_max)
    {
        double closest_t_sphere = INFINITY;
        double closest_t_plane = INFINITY;
        Sphere closest_sphere;
        Plane closest_plane;
        double t1, t2;
        for (int i = 0; i < spheres.size(); i++)
        {
            tie(t1, t2) = intersect_ray_sphere(p0, D, spheres[i]);
            if ((t1 >= t_min && t1 <= t_max) && t1 < closest_t_sphere)
            {
                closest_t_sphere = t1;
                closest_sphere = spheres[i];
            }
            if ((t2 >= t_min && t2 <= t_max) && t2 < closest_t_sphere)
            {
                closest_t_sphere = t2;
                closest_sphere = spheres[i];
            }
        }

        for (int i = 0; i < planes.size(); i++)
        {
            t1 = intersect_ray_plane(p0, D, planes[i]);
            if ((t1 >= t_min && t1 <= t_max) && t1 < closest_t_plane)
            {
                closest_t_plane = t1;
                closest_plane = planes[i];
            }
        }

        if (closest_t_sphere == INFINITY && closest_t_plane == INFINITY)
        {
            return canva.bg;
        }

        if (closest_t_sphere < closest_t_plane)
        {
            Vector pi = Vector(p0.x + (D.x * closest_t_sphere), p0.y + (D.y * closest_t_sphere), p0.z + (D.z * closest_t_sphere));
            Vector N = Vector(pi.x - closest_sphere.center.x, pi.y - closest_sphere.center.y, pi.z - closest_sphere.center.z);
            N = Vector(N.x / closest_sphere.radius, N.y / closest_sphere.radius, N.z / closest_sphere.radius);

            // Vector(-D.x, -D.y, -D.z) -> É PRA NORMALIZAR? v = -dr/||dr||
            Vector i = compute_lighting(pi, N, Vector(-D.x, -D.y, -D.z), closest_sphere); // intensidade

            return Color(closest_sphere.bg.r * i.x, closest_sphere.bg.g * i.y, closest_sphere.bg.b * i.z);
        }

        Vector pi = Vector(p0.x + (D.x * closest_t_plane), p0.y + (D.y * closest_t_plane), p0.z + (D.z * closest_t_plane));
        Vector i = compute_lighting(pi, closest_plane.normal, Vector(-D.x, -D.y, -D.z), closest_plane);

        return Color(closest_plane.bg.r * i.x, closest_plane.bg.g * i.y, closest_plane.bg.b * i.z);
    }
} Scene;

int main()
{
    Viewport vp(0.6, 0.6, 0.3);

    Canvas canva(500., 500., vp, Color(100, 100, 100));

    vector<Sphere> spheres;
    vector<Plane> planes;
    vector<Light> lights;

    Sphere sphere1(0.4, Vector(0, 0, -1.), Color(0, 0, 255), 10., Vector(0.7, 0.2, 0.2), Vector(0.7, 0.2, 0.2), Vector(0.7, 0.2, 0.2)); // blue
    // Sphere sphere2(0.4, Vector(1., 0, -1.), Color(255, 0, 0), 10., Vector(0.7, 0.2, 0.2), Vector(0.7, 0.2, 0.2), Vector(0.7, 0.2, 0.2));  // red
    // Sphere sphere3(0.4, Vector(-1., 0, -1.), Color(0, 255, 0), 10., Vector(0.7, 0.2, 0.2), Vector(0.7, 0.2, 0.2), Vector(0.7, 0.2, 0.2)); // green

    spheres.push_back(sphere1);
    // spheres.push_back(sphere2);
    // spheres.push_back(sphere3);

    Plane plane1(Vector(0, -0.4, 0), Vector(0, 1., 0), Color(30., 40., 50.), 1., Vector(0.2, 0.7, 0.2), Vector(0, 0, 0), Vector(0.2, 0.7, 0.2));
    Plane plane2(Vector(0, 0, -2.), Vector(0, 0, 1.), Color(80., 90., 100.), 1., Vector(0.3, 0.3, 0.7), Vector(0, 0, 0), Vector(0.3, 0.3, 0.7));

    planes.push_back(plane1);
    planes.push_back(plane2);

    Light point_light(Vector(0.7, 0.7, 0.7), Vector(0, 0.6, -0.3), "point");
    Light ambient_light(Vector(0.3, 0.3, 0.3), Vector(0, 0, 0), "ambient");

    lights.push_back(point_light);
    lights.push_back(ambient_light);

    Scene scene(spheres, planes, canva, lights);

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