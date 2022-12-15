#ifndef SCENE_HPP
#define SCENE_HPP
#include <bits/stdc++.h>
#include "./canvas.hpp"
#include "./objects.hpp"
#include "./lights.hpp"
#include "./vector.hpp"
#include "./color.hpp"
#include "./canvas.hpp"
#include "./mesh_struct.hpp"

using namespace std;

#define ORTO 2
#define PERSPECT 3

typedef struct Scene
{
    int projection = PERSPECT;
    vector<Object *> objects;
    Canvas &canva;
    vector<Light *> lights;
    vector<vector<Object *>> point_to_obj;

    Scene(Canvas &c) : canva(c)
    {
        point_to_obj = vector<vector<Object *>>(c.h, vector<Object *>(c.w, nullptr));
    }

    double dot(Vector a, Vector b)
    {
        return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
    }

    Vector vector_mult(Vector a, Vector b)
    {
        double row1 = (a.y * b.z) - (a.z * b.y);
        double row2 = (a.z * b.x) - (a.x * b.z);
        double row3 = (a.x * b.y) - (a.y * b.x);
        return Vector(row1, row2, row3, 0);
    }

    Vector sum_vector(Vector a, Vector b)
    {
        return Vector(a.x + b.x, a.y + b.y, a.z + b.z, 0);
    }

    Vector sub_vector(Vector a, Vector b)
    {
        return Vector(a.x - b.x, a.y - b.y, a.z - b.z, 0);
    }

    double calc_intensity(double i, double x, double length_a, double length_b, double specular, double k, int flag)
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

    pair<double, double> intersect_ray_sphere(Vector p0, Vector D, Object *sphere)
    {
        double r = sphere->radius;

        Vector w = Vector(p0.x - sphere->center.x, p0.y - sphere->center.y, p0.z - sphere->center.z, 0);

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

    double intersect_ray_plane(Vector p0, Vector D, Object *plane)
    {
        Vector w = Vector(p0.x - plane->center.x, p0.y - plane->center.y, p0.z - plane->center.z, 0);
        Vector normal = plane->normal;

        double num = ((w.x * normal.x) + (w.y * normal.y) + (w.z * normal.z));
        double den = ((D.x * normal.x) + (D.y * normal.y) + (D.z * normal.z));

        double ti = -num / den;

        if (ti < 0)
        {
            return INFINITY;
        }

        return ti;
    }

    bool is_in_shell(Vector Pi, Object *object)
    {
        Vector Pi_sub_base = Vector(Pi.x - object->center.x, Pi.y - object->center.y, Pi.z - object->center.z, 0);
        return (dot(Pi_sub_base, object->u) < 0 || dot(Pi_sub_base, object->u) > object->h);
    }

    bool is_in_base(Vector P, Object *plane, Object *object)
    {
        Vector CP = Vector(P.x - plane->center.x, P.y - plane->center.y, P.z - plane->center.z, 0);
        double CP_length = length(CP);
        double cp_x_n = dot(CP, object->u);
        bool is_zero = (cp_x_n > 0.0) && (cp_x_n < 0.0);
        bool is_radius = CP_length <= object->radius;
        return (is_zero && is_radius);
    }

    double intersect_ray_base(Vector p0, Vector D, Object *plane, Object *object) // cylinder or cone
    {
        double t;
        t = intersect_ray_plane(p0, D, plane);
        Vector P = Vector(p0.x + (t * D.x), p0.y + (t * D.y), p0.z + (t * D.z), 0);
        return is_in_base(P, plane, object) ? t : INFINITY;
    }

    pair<double, double> intersect_ray_cylinder(Vector p0, Vector D, Object *cylinder)
    {
        double r = cylinder->radius;
        double t1, t2;

        Vector p0_sub_base = Vector(p0.x - cylinder->center.x, p0.y - cylinder->center.y, p0.z - cylinder->center.z, 0);

        Vector v = Vector(p0_sub_base.x - (dot(p0_sub_base, cylinder->u)) * cylinder->u.x, p0_sub_base.y - (dot(p0_sub_base, cylinder->u)) * cylinder->u.y, p0_sub_base.z - (dot(p0_sub_base, cylinder->u)) * cylinder->u.z, 0);
        Vector w = Vector(D.x - ((dot(D, cylinder->u)) * cylinder->u.x), D.y - ((dot(D, cylinder->u)) * cylinder->u.y), D.z - ((dot(D, cylinder->u)) * cylinder->u.z), 0);

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

        Vector Pi1 = Vector(p0.x + (D.x * t1), p0.y + (D.y * t1), p0.z + (D.z * t1), 0);
        t1 = is_in_shell(Pi1, cylinder) ? INFINITY : t1;
        Vector Pi2 = Vector(p0.x + (D.x * t2), p0.y + (D.y * t2), p0.z + (D.z * t2), 0);
        t2 = is_in_shell(Pi2, cylinder) ? INFINITY : t2;

        return {t1, t2};
    }

    pair<double, double> intersect_ray_cone(Vector p0, Vector D, Object *cone)
    {
        double r = cone->radius;
        double h = cone->h;
        double cos_sqr_teta = (h * h) / ((r * r) + (h * h));
        double t1, t2;

        Vector v = Vector(cone->vc.x - p0.x, cone->vc.y - p0.y, cone->vc.z - p0.z, 0);

        double a = (dot(D, cone->u) * dot(D, cone->u)) - dot(D, D) * cos_sqr_teta;
        double b = (dot(v, D) * cos_sqr_teta - dot(v, cone->u) * dot(D, cone->u)) * 2;
        double c = (dot(v, cone->u) * dot(v, cone->u)) - dot(v, v) * cos_sqr_teta;

        double delta = b * b - 4 * a * c;

        if (delta < 0)
        {
            return {INFINITY, INFINITY};
        }

        t1 = (-b + sqrt(delta)) / (2 * a);
        t2 = (-b - sqrt(delta)) / (2 * a);

        Vector Pi1 = Vector(p0.x + (D.x * t1), p0.y + (D.y * t1), p0.z + (D.z * t1), 0);
        t1 = is_in_shell(Pi1, cone) ? INFINITY : t1;
        Vector Pi2 = Vector(p0.x + (D.x * t2), p0.y + (D.y * t2), p0.z + (D.z * t2), 0);
        t2 = is_in_shell(Pi2, cone) ? INFINITY : t2;

        return {t1, t2};
    }

    pair<double, Vector> intersect_ray_cube(Vector p0, Vector D, Object *cube)
    {
        // Object mesh_cube = mesh(cube);
        int v1, v2, v3;
        Vector P1, P2, P3;
        Vector r1, r2;
        Vector n_face, N_face;
        Vector Pi;
        double C1, C2, C3;
        double t = INFINITY;
        double closest_t = INFINITY;
        Vector normal_closest_face;
        double EPS = 0.0000001;
        for (int i = 0; i < cube->LF.size(); i++)
        {
            int idA1 = cube->LF[i].a1;
            int idA2 = cube->LF[i].a2;
            int idA3 = cube->LF[i].a3;

            int idV11 = cube->LA[idA1].v1 + 1;
            int idV12 = cube->LA[idA1].v2 + 1;
            int idV21 = cube->LA[idA2].v1 + 1;
            int idV22 = cube->LA[idA2].v2 + 1;

            int n1 = idV11 * idV12;
            int n = n1 / idV21;
            if (n == idV11 || n == idV12)
            {
                v1 = idV21;
                v2 = idV22;
                v3 = n;
            }
            else
            {
                v1 = idV22;
                v2 = idV21;
                v3 = n1 / v1;
            }

            v1 = v1 - 1;
            v2 = v2 - 1;
            v3 = v3 - 1;

            P1.x = cube->LV[v1].x;
            P1.y = cube->LV[v1].y;
            P1.z = cube->LV[v1].z;

            P2.x = cube->LV[v2].x;
            P2.y = cube->LV[v2].y;
            P2.z = cube->LV[v2].z;

            P3.x = cube->LV[v3].x;
            P3.y = cube->LV[v3].y;
            P3.z = cube->LV[v3].z;

            r1 = sub_vector(P2, P1);
            r2 = sub_vector(P3, P1);

            N_face = vector_mult(r1, r2);
            double N_length = length(N_face);
            n_face = Vector(-N_face.x / N_length, -N_face.y / N_length, -N_face.z / N_length, 0);

            t = -(dot(sub_vector(p0, P1), n_face)) / dot(D, n_face);
            Pi = Vector(p0.x + (t * D.x), p0.y + (t * D.y), p0.z + (t * D.z), 0);
            C1 = (dot(vector_mult(sub_vector(P3, Pi), sub_vector(P1, Pi)), n_face)) / dot(N_face, n_face);
            C2 = (dot(vector_mult(sub_vector(P1, Pi), sub_vector(P2, Pi)), n_face)) / dot(N_face, n_face);
            C3 = (dot(vector_mult(sub_vector(P2, Pi), sub_vector(P3, Pi)), n_face)) / dot(N_face, n_face);

            if (C1 >= 0. - EPS && C2 >= 0. - EPS && C3 >= 0. - EPS)
            {
                if (closest_t > t)
                {
                    closest_t = t;
                    normal_closest_face = n_face;
                }
            }
        }

        return {closest_t, normal_closest_face};
    }

    bool has_shadow(Vector pi, Vector L, double length_Pf_Pi)
    {
        bool shadow = false;
        double s1, s2;
        Vector normal_face;

        for (int i = 0; i < objects.size(); i++)
        {
            if (objects[i]->type == "sphere")
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
            if (objects[i]->type == "plane")
            {
                s1 = intersect_ray_plane(pi, L, objects[i]);
                if (s1 > 0 && s1 < length_Pf_Pi)
                {
                    shadow = true;
                }
            }
            if (objects[i]->type == "cylinder")
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
            if (objects[i]->type == "cone")
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
            if (objects[i]->type == "cube")
            {
                tie(s1, normal_face) = intersect_ray_cube(pi, L, objects[i]);
                if (s1 > 0 && s1 < length_Pf_Pi)
                {
                    shadow = true;
                }
            }
        }

        return shadow;
    }

    pair<Vector, bool> compute_lighting(Vector pi, Vector N, Vector V, Object *object, int y, int x)
    {
        Vector i(0.0, 0.0, 0.0, 0);
        Vector L;
        Vector ds;             // case spot
        Vector intensity_spot; // case spot
        Light *ambient;

        for (int k = 0; k < lights.size(); k++)
        {
            if (lights[k]->type == "ambient")
            {
                ambient = lights[k];
            }
        }

        for (int j = 0; j < lights.size(); j++)
        {
            if (lights[j]->type == "ambient")
            {
                i.x += lights[j]->intensity.x * object->k_a.x;
                i.y += lights[j]->intensity.y * object->k_a.y;
                i.z += lights[j]->intensity.z * object->k_a.z;
            }
            else
            {
                if (lights[j]->type == "point")
                {
                    L = Vector(lights[j]->position.x - pi.x, lights[j]->position.y - pi.y, lights[j]->position.z - pi.z, 0);
                }
                else if (lights[j]->type == "directional")
                {
                    L = lights[j]->direction;
                }
                double Pf_Pi = length(L);
                L = Vector(L.x / Pf_Pi, L.y / Pf_Pi, L.z / Pf_Pi, 0);

                if (lights[j]->type == "spot")
                {
                    double direction_length = length(lights[j]->direction);
                    ds = Vector(-(lights[j]->direction.x / direction_length), -(lights[j]->direction.y / direction_length), -(lights[j]->direction.z / direction_length), 0);
                    L = Vector(lights[j]->position.x - pi.x, lights[j]->position.y - pi.y, lights[j]->position.z - pi.z, 0);
                    Pf_Pi = length(L);
                    L = Vector(L.x / Pf_Pi, L.y / Pf_Pi, L.z / Pf_Pi, 0);
                    double l_dot_ds = dot(L, ds);
                    if (l_dot_ds < cos(lights[j]->grau))
                    {
                        return {Vector(255 * (ambient->intensity.x * object->k_a.x), 255 * (ambient->intensity.y * object->k_a.y), 255 * (ambient->intensity.z * object->k_a.z), 0), false};
                    }
                    intensity_spot = Vector(lights[j]->intensity.x * l_dot_ds, lights[j]->intensity.y * l_dot_ds, lights[j]->intensity.z * l_dot_ds, 0);
                }

                if (has_shadow(pi, L, Pf_Pi))
                {
                    if (object->type == "plane")
                    {
                        if (object->texture)
                        {
                            int h = object->height;
                            int w = object->width;
                            return {Vector(object->matrix_img[y % h][x % w].r * (ambient->intensity.x * object->k_a.x), object->matrix_img[y % h][x % w].g * (ambient->intensity.y * object->k_a.y), object->matrix_img[y % h][x % w].b * (ambient->intensity.z * object->k_a.z), 0), true};
                        }
                    }
                    return {Vector(255 * (ambient->intensity.x * object->k_a.x), 255 * (ambient->intensity.y * object->k_a.y), 255 * (ambient->intensity.z * object->k_a.z), 0), true};
                }

                double n_dot_l = dot(N, L);

                // DIFUSA
                if (n_dot_l > 0)
                {
                    if (lights[j]->type == "spot")
                    {
                        i.x += calc_intensity(intensity_spot.x, n_dot_l, length(N), length(L), object->specular, object->k_d.x, 1);
                        i.y += calc_intensity(intensity_spot.y, n_dot_l, length(N), length(L), object->specular, object->k_d.y, 1);
                        i.z += calc_intensity(intensity_spot.z, n_dot_l, length(N), length(L), object->specular, object->k_d.z, 1);
                    }
                    else
                    {
                        i.x += calc_intensity(lights[j]->intensity.x, n_dot_l, length(N), length(L), object->specular, object->k_d.x, 1);
                        i.y += calc_intensity(lights[j]->intensity.y, n_dot_l, length(N), length(L), object->specular, object->k_d.y, 1);
                        i.z += calc_intensity(lights[j]->intensity.z, n_dot_l, length(N), length(L), object->specular, object->k_d.z, 1);
                    }
                }

                // ESPECULAR
                if (object->specular != -1)
                {
                    Vector R = Vector(((2 * dot(L, N)) * N.x) - L.x, ((2 * dot(L, N)) * N.y) - L.y, ((2 * dot(L, N)) * N.z) - L.z, 0);
                    double r_dot_v = dot(R, V);
                    if (r_dot_v > 0)
                    {
                        if (lights[j]->type == "spot")
                        {
                            i.x += calc_intensity(intensity_spot.x, r_dot_v, length(R), length(V), object->specular, object->k_e.x, 0);
                            i.y += calc_intensity(intensity_spot.y, r_dot_v, length(R), length(V), object->specular, object->k_e.y, 0);
                            i.z += calc_intensity(intensity_spot.z, r_dot_v, length(R), length(V), object->specular, object->k_e.z, 0);
                        }
                        else
                        {
                            i.x += calc_intensity(lights[j]->intensity.x, r_dot_v, length(R), length(V), object->specular, object->k_e.x, 0);
                            i.y += calc_intensity(lights[j]->intensity.y, r_dot_v, length(R), length(V), object->specular, object->k_e.y, 0);
                            i.z += calc_intensity(lights[j]->intensity.z, r_dot_v, length(R), length(V), object->specular, object->k_e.z, 0);
                        }
                    }
                }
            }
        }

        return {Vector(255 * i.x, 255 * i.y, 255 * i.z, 0), false};
    }

    Color define_color(double closest_t, Object *object, Vector p0, Vector D, int y, int x)
    {
        Vector pi = Vector(p0.x + (D.x * closest_t), p0.y + (D.y * closest_t), p0.z + (D.z * closest_t), 0);
        Vector N;
        if (object->type == "sphere")
        {
            N = Vector(pi.x - object->center.x, pi.y - object->center.y, pi.z - object->center.z, 0);
            N = Vector(N.x / object->radius, N.y / object->radius, N.z / object->radius, 0);
        }
        if (object->type == "plane")
        {
            N = object->normal;
        }
        if (object->type == "cylinder")
        {
            Vector CP = Vector(pi.x - object->center.x, pi.y - object->center.y, pi.z - object->center.z, 0);
            double CP_dot_u = dot(CP, object->u);
            Vector AP = Vector(CP.x - (object->u.x * (CP_dot_u)), CP.y - (object->u.y * (CP_dot_u)), CP.z - (object->u.z * (CP_dot_u)), 0);
            double AP_length = length(AP);
            N = Vector(AP.x / AP_length, AP.y / AP_length, AP.z / AP_length, 0);
        }
        if (object->type == "cone")
        {
            Vector w = Vector(object->vc.x - pi.x, object->vc.y - pi.y, object->vc.z - pi.z, 0);
            Vector n_barra = vector_mult(w, object->u);
            N = vector_mult(n_barra, w);
            double N_length = length(N);
            N = Vector(N.x / N_length, N.y / N_length, N.z / N_length, 0);
        }
        if (object->type == "cube")
        {
            N = object->normal;
        }

        // Vector L = Vector(lights[0].position.x - pi.x, lights[0].position.y - pi.y, lights[0].position.z - pi.z, 0);
        // double length_Pf_Pi = length(L);
        // L = Vector(L.x / length_Pf_Pi, L.y / length_Pf_Pi, L.z / length_Pf_Pi, 0);

        // if (has_shadow(pi, L, length_Pf_Pi))
        // {
        //     if (object.type == "plane")
        //     {
        //         if (object.texture)
        //         {
        //             int h = object.height;
        //             int w = object.width;
        //             return Color(object.matrix_img[y % h][x % w].r * (lights[1].intensity.x * object.k_a.x), object.matrix_img[y % h][x % w].g * (lights[1].intensity.y * object.k_a.y), object.matrix_img[y % h][x % w].b * (lights[1].intensity.z * object.k_a.z));
        //         }
        //     }
        //     return Color(255 * (lights[1].intensity.x * object.k_a.x), 255 * (lights[1].intensity.y * object.k_a.y), 255 * (lights[1].intensity.z * object.k_a.z));
        // }

        Vector i;
        bool shadow;
        tie(i, shadow) = compute_lighting(pi, N, Vector(-D.x, -D.y, -D.z, 0), object, y, x);

        if (object->type == "plane")
        {
            if (object->texture)
            {
                int h = object->height;
                int w = object->width;

                if (shadow)
                {
                    return Color(i.x, i.y, i.z);
                }

                return Color(object->matrix_img[y % h][x % w].r * (i.x / 255.), object->matrix_img[y % h][x % w].g * (i.y / 255.), object->matrix_img[y % h][x % w].b * (i.z / 255.));
            }
        }

        return Color(i.x, i.y, i.z);
    }

    struct info_intercept
    {
        Color color;
        Object *object = nullptr;
    };

    pair<Color, Object *> trace_ray(Vector p0, Vector D, double t_min, double t_max, int y, int x)
    {
        double length_D = length(D);
        D = Vector(D.x / length_D, D.y / length_D, D.z / length_D, 0);
        double closest_t_sphere = INFINITY;
        double closest_t_plane = INFINITY;
        double closest_t_cylinder = INFINITY;
        double closest_t_cone = INFINITY;
        double closest_t_cube = INFINITY;
        Object *closest_sphere = nullptr;
        Object *closest_plane = nullptr;
        Object *closest_cylinder = nullptr;
        Object *closest_cone = nullptr;
        Object *closest_cube = nullptr;
        double EPS = 0.01;
        double t1, t2;
        Vector normal_face;
        for (int i = 0; i < objects.size(); i++)
        {
            if (objects[i]->type == "sphere")
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
            if (objects[i]->type == "plane")
            {
                t1 = intersect_ray_plane(p0, D, objects[i]);
                if ((t1 > t_min && t1 < t_max) && t1 < closest_t_plane)
                {
                    closest_t_plane = t1;
                    closest_plane = objects[i];
                }
            }
            if (objects[i]->type == "cylinder")
            {
                tie(t1, t2) = intersect_ray_cylinder(p0, D, objects[i]);
                if ((t1 > t_min && t1 < t_max) && t1 < closest_t_cylinder)
                {
                    closest_t_cylinder = t1;
                    closest_cylinder = objects[i];
                }
                else if ((t2 > t_min && t2 < t_max) && t2 < closest_t_cylinder)
                {
                    closest_t_cylinder = t2;
                    closest_cylinder = objects[i];
                }
                else
                {
                    continue;
                }
                Vector N_base = Vector(-closest_cylinder->u.x, -closest_cylinder->u.y, -closest_cylinder->u.z, 0);
                Vector P_pi_base = closest_cylinder->center;
                Object *plane_base = new Object("plane", P_pi_base, N_base, objects[i]->specular, objects[i]->k_d, objects[i]->k_e, objects[i]->k_a, false, "");
                Vector N_cover = closest_cylinder->u;
                Vector P_pi_cover = Vector(closest_cylinder->center.x + (closest_cylinder->u.x * closest_cylinder->h), closest_cylinder->center.y + (closest_cylinder->u.y * closest_cylinder->h), closest_cylinder->center.z + (closest_cylinder->u.z * closest_cylinder->h), 0);
                Object *plane_top = new Object("plane", P_pi_cover, N_cover, objects[i]->specular, objects[i]->k_d, objects[i]->k_e, objects[i]->k_a, false, "");
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
            if (objects[i]->type == "cone")
            {
                tie(t1, t2) = intersect_ray_cone(p0, D, objects[i]);
                if ((t1 > t_min && t1 < t_max) && t1 < closest_t_cone)
                {
                    closest_t_cone = t1;
                    closest_cone = objects[i];
                }
                else if ((t2 > t_min && t2 < t_max) && t2 < closest_t_cone)
                {
                    closest_t_cone = t2;
                    closest_cone = objects[i];
                }
                else
                {
                    continue;
                }
                Vector N_base = Vector(closest_cone->u.x, closest_cone->u.y, closest_cone->u.z, 0);
                Vector P_pi_base = closest_cone->center;
                Object *plane_base = new Object("plane", P_pi_base, N_base, objects[i]->specular, objects[i]->k_d, objects[i]->k_e, objects[i]->k_a, false, "");
                t1 = intersect_ray_base(p0, D, plane_base, closest_cone);
                if ((t1 > t_min && t1 < t_max) && t1 < closest_t_cone)
                {
                    closest_t_cone = t1;
                    closest_cone = plane_base;
                }
            }
            if (objects[i]->type == "cube")
            {
                tie(t1, normal_face) = intersect_ray_cube(p0, D, objects[i]);
                if ((t1 > t_min && t1 < t_max) && t1 < closest_t_cone)
                {
                    closest_t_cube = t1;
                    Object *cube_face = new Object("cube", objects[i]->specular, normal_face, objects[i]->k_d, objects[i]->k_e, objects[i]->k_a);
                    objects[i]->normal = normal_face;
                    closest_cube = objects[i];
                }
            }
        }

        if (closest_t_sphere == INFINITY && closest_t_plane == INFINITY && closest_t_cylinder == INFINITY && closest_t_cone == INFINITY && closest_t_cube == INFINITY)
        {
            return {canva.bg, nullptr};
        }

        vector<Closest_Object *> closest_objects;

        closest_objects.push_back(new Closest_Object(closest_t_sphere, closest_sphere));
        closest_objects.push_back(new Closest_Object(closest_t_plane, closest_plane));
        closest_objects.push_back(new Closest_Object(closest_t_cylinder, closest_cylinder));
        closest_objects.push_back(new Closest_Object(closest_t_cone, closest_cone));
        closest_objects.push_back(new Closest_Object(closest_t_cube, closest_cube));

        double smaller = INFINITY;
        double closest_t;
        Object *closest_object = nullptr;
        for (int i = 0; i < closest_objects.size(); i++) // pega o objeto com menor t (objeto mais próximo)
        {
            if (closest_objects[i]->closest_t < smaller)
            {
                smaller = closest_objects[i]->closest_t;
                closest_t = closest_objects[i]->closest_t;
                closest_object = closest_objects[i]->object;
            }
        }

        return {define_color(closest_t - EPS, closest_object, p0, D, y, x), closest_object};
    }
} Scene;

#endif