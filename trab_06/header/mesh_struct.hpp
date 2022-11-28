#ifndef MESH_STRUCT_HPP
#define MESH_STRUCT_HPP

typedef struct VertexIndex
{
    int v1, v2;
    VertexIndex(int vv1, int vv2)
    {
        v1 = vv1;
        v2 = vv2;
    }
    VertexIndex() {}
} VertexIndex;

typedef struct EdgeIndex
{
    int a1, a2, a3;
    EdgeIndex(int aa1, int aa2, int aa3)
    {
        a1 = aa1;
        a2 = aa2;
        a3 = aa3;
    }
    EdgeIndex() {}
} EdgeIndex;

#endif