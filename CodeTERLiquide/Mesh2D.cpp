#ifndef _MESH_2D_CPP

#include "Mesh2D.h"
#include<iostream>
#include<fstream>
#include <assert.h>

using namespace std;
using namespace Eigen;

Vertice::Vertice()
{
  _vert[0] = -10000; _vert[1] = -10000; _ref = -1;
}

Vertice::Vertice(double x, double y, int ref) : _ref(ref)
{
  _vert[0] = x; _vert[1] = y;
}

void Vertice::Print() const
{
  cout << "[x, y] = [" << _vert[0] << " " << _vert[1] << "];" << endl;
  cout << "ref = " << _ref << endl;
}

Edge::Edge()
{
  _edg[0] = -1; _edg[1] = -1; _ref = -1;
}

Edge::Edge(int edg1, int edg2, int ref) : _ref(ref)
{
  // sort
  if (edg1 > edg2)
  {
    _edg[0] = edg2;
    _edg[1] = edg1;
  }
  else
  {
    _edg[0] = edg1;
    _edg[1] = edg2;
  }
  _t1 = -1;
  _t2 = -1;
}

void Edge::Print() const
{
  cout << "[pt1, pt2] = [" << _edg[0] << " " << _edg[1] << "];" << endl;
  cout << "[t1, t2] = [" << _t1 << " " << _t2 << "];" << endl;
  cout << "ref = " << _ref << endl;
}

Triangle::Triangle()
{
  _tri[0] = -1; _tri[1] = -1; _tri[2] = -1; _ref = -1;
}

Triangle::Triangle(int tri1, int tri2, int tri3, int ref) : _ref(ref)
{
  _tri[0] = tri1; _tri[1] = tri2; _tri[2] = tri3;
}

void Triangle::Print() const
{
  cout << "[pt1, pt2, pt3] = [" << _tri[0] << " " << _tri[1] << " " << _tri[2] << "];" << endl;
  cout << "ref = " << _ref << endl;
}

Mesh2D::Mesh2D()
{
}

void Mesh2D::BuildTrianglesCenterAndArea()
{

  _tri_center.resize(_triangles.size(),2);
  _tri_area.resize(_triangles.size());
  _p.resize(_triangles.size(),3);
  _x.resize(_triangles.size(),3);
  _y.resize(_triangles.size(),3);

  double x1(0.),y1(0.),x2(0.),y2(0.),x3(0.),y3(0.);

  for (int i = 0; i < _triangles.size(); i++)
  {
    _p(i,0) = _triangles[i].GetTriangle()(0);
    _x(i,0) = _vertices[_p(i,0)].GetVertice()(0);
    _y(i,0) = _vertices[_p(i,0)].GetVertice()(1);
    _p(i,1) = _triangles[i].GetTriangle()(1);
    _x(i,1) = _vertices[_p(i,1)].GetVertice()(0);
    _y(i,1) = _vertices[_p(i,1)].GetVertice()(1);
    _p(i,2) = _triangles[i].GetTriangle()(2);
    _x(i,2) = _vertices[_p(i,2)].GetVertice()(0);
    _y(i,2) = _vertices[_p(i,2)].GetVertice()(1);

    //Construction de la matrice contenant les coordonnées des centres de chaque triangle
    _tri_center.row(i) << (_x(i,1)+_x(i,0)+_x(i,2))/3.,(_y(i,1)+_y(i,0)+_y(i,2))/3.;
    //Construction du vecteur contenant les aires de chaque triangle
    _tri_area(i) = abs((_x(i,1)-_x(i,0))*(_y(i,2)-_y(i,0))-(_y(i,1)-_y(i,0))*(_x(i,2)-_x(i,0)))/2.;
  }
}


void Mesh2D::ReadMesh(string name_mesh)
{
  ifstream mesh_file(name_mesh.data());
  if (!mesh_file.is_open())
  {
    cout << "Unable to open file " << name_mesh << endl;
    abort();
  }
  else
  {
    cout << "-------------------------------------------------" << endl;
    cout << "Reading mesh: " << name_mesh << endl;
  }

  string file_line;
  vector<Edge> edges_boundary;
  int dim = 3;

  while (!mesh_file.eof())
  {
    getline(mesh_file, file_line);
    if (file_line.find("Dimension") != std::string::npos)
    {
      mesh_file >> dim;
    }
    else if (file_line.find("Vertices") != std::string::npos)
    {
      int nb_vertices(0);
      mesh_file >> nb_vertices;
      cout << "Number of vertices  (" << nb_vertices << ")" << endl;
      _vertices.resize(nb_vertices);
      for (int i = 0 ; i < nb_vertices ; ++i)
      {
        double x,y,z; int ref;
        mesh_file >> x >> y >> z >> ref;
        _vertices[i] = Vertice(x, y, ref);
      }
    }
    else if (file_line.find("Edges") != std::string::npos)
    {
      int nb_edges(0);
      mesh_file >> nb_edges;
      cout << "Number of edges (" << nb_edges << ")" << endl;
      edges_boundary.resize(nb_edges);
      int n1, n2, ref;
      for (int i = 0 ; i < nb_edges ; ++i)
      {
        mesh_file >> n1 >> n2 >> ref;
        n1--; n2--;
        edges_boundary[i] = Edge(n1, n2, ref);
      }
    }
    else if (file_line.find("Triangles") != std::string::npos)
    {
      int nb_triangles(0);
      mesh_file >> nb_triangles;
      cout << "Number of triangles (" << nb_triangles << ")" << endl;
      _triangles.resize(nb_triangles);
      for (int i = 0 ; i < nb_triangles ; ++i)
      {
        int tri1, tri2, tri3, ref;
        mesh_file >> tri1 >> tri2 >> tri3 >> ref;
        tri1--; tri2--; tri3--;
        _triangles[i] = Triangle(tri1, tri2, tri3, ref);
      }
    }
  }


  cout << "-----------Triangles center and area-------------" << endl;
  BuildTrianglesCenterAndArea();


  cout << "-------------------------------------------------" << endl;
}

#define _MESH_2D_CPP
#endif
