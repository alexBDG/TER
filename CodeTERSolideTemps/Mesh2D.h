#ifndef _MESH_2D_H

#include<vector>
#include<string>
#include "Dense"
#include "Sparse"

class Vertice
{
  private:
    Eigen::Vector2d _vert;
    int _ref;
  public:
    Vertice();
    Vertice(double x, double y, int ref);
    void Print() const;
    const Eigen::Vector2d GetVertice() const {return _vert;};
};

class Edge
{
  private:
    Eigen::Vector2i _edg;
    int _ref;
    int _t1, _t2;
  public:
    Edge();
    Edge(int edg1, int edg2, int ref);
    void Print() const;
    void AddTriangle(int t)
    {
      if (_t1 == -1)
        _t1 = t;
      else
        _t2 = t;
    }
    const Eigen::Vector2i& GetEdge() const { return _edg;}
    int GetT1() const { return _t1; };
    int GetT2() const { return _t2; };
    int GetReference() const { return _ref;};
};

class Triangle
{
  private:
    Eigen::Vector3i _tri;
    int _ref;
  public:
    Triangle();
    Triangle(int t1, int t2, int t3, int ref);
    void Print() const;
    const Eigen::Vector3i& GetTriangle() const { return _tri; }
};

class Mesh2D
{
private:
  // liste de tous les sommets
  std::vector<Vertice> _vertices;
  // liste de tous les triangles
  std::vector<Triangle> _triangles;
  // centre de tous les triangles
  Eigen::Matrix<double, Eigen::Dynamic, 2> _tri_center;
  // aire de tous les triangles
  Eigen::VectorXd _tri_area;
  //liste des 3 somments asociés à chaque triangles
  Eigen::Matrix<double, Eigen::Dynamic, 3> _p;
  //liste des coordonnées des sommets associés a chaque triangle
  Eigen::Matrix<double, Eigen::Dynamic, 3> _x;
  Eigen::Matrix<double, Eigen::Dynamic, 3> _y;
  // liste de toutes les arêtes
  std::vector<Edge> _edges;
  // liste de toutes les normales unitaires !!!
  Eigen::Matrix<double, Eigen::Dynamic, 2> _edg_normal;
  // liste de toutes les longueurs d'arêtes
  Eigen::VectorXd _edg_length;
  // centre des aretes
  Eigen::Matrix<double, Eigen::Dynamic, 2> _edg_center;

public:
  Mesh2D();
  void ReadMesh(std::string name_mesh);
  void BuildTrianglesCenterAndArea();
  void BuildPhi(const double x, const double y, int i);

  const std::vector<Vertice> & GetVertices() const {return _vertices;};

  const std::vector<Triangle> & GetTriangles() const {return _triangles;};
  const Eigen::Matrix<double, Eigen::Dynamic, 2> & GetTrianglesCenter() const {return _tri_center;};
  const Eigen::VectorXd & GetTrianglesArea() const  {return _tri_area;};

  Eigen::Matrix<double, Eigen::Dynamic, 3> GetP() {return _p;};
  Eigen::Matrix<double, Eigen::Dynamic, 3> GetX() {return _x;};
  Eigen::Matrix<double, Eigen::Dynamic, 3> GetY() {return _y;};

  const std::vector<Edge> & GetEdges() const {return _edges;};
  const Eigen::VectorXd & GetEdgesLength() const {return _edg_length;};
  const Eigen::Matrix<double, Eigen::Dynamic, 2> & GetEdgesNormal() const {return _edg_normal;};
  const Eigen::Matrix<double, Eigen::Dynamic, 2> & GetEdgesCenter() const {return _edg_center;};

protected:
  void AddSingleEdge(const Edge& edge, int ne, std::vector<int>& head_minv,
		     std::vector<int>& next_edge, int& nb_edges);
};

#define _MESH_2D_H
#endif
