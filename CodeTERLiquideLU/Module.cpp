#ifndef _MODULE_CPP

#include "Module.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <chrono>

using namespace Eigen;
using namespace std;

Calcul::Calcul(Mesh2D* maillage) :
_vertices(maillage->GetVertices()),
_tri_area(maillage->GetTrianglesArea()),
_triangles(maillage->GetTriangles()),
_tri_center(maillage->GetTrianglesCenter()),
_edges(maillage->GetEdges())
{
  _c = 100.;
  _v = 148000;
  _maillage = maillage;
  nb_pts = _vertices.size();
  _Sommets = maillage->GetP();
  _X = maillage->GetX();
  _Y = maillage->GetY();
  system("mkdir -p ./Resultats");
  ofstream sensor;
  sensor.open("Resultats/sensor.dat", ios::out);
  sensor.close();
}

void Calcul::BuildMK(double dt)
{
  _M.resize(_vertices.size(),_vertices.size());
  _K.resize(_vertices.size(),_vertices.size());

  cout << "Création de M et K : " << endl;
  for (int u = 0; u < _triangles.size(); u++)
  {
    cout.flush();
    cout << "Progression : " << (double)u/((double)_triangles.size())*100 << "% \r";

    //Dérivés dues au changment de variables
    double dxpx, dxpy, dypx, dypy;
    dxpx = (_Y(u,2)-_Y(u,0))/((_Y(u,2)-_Y(u,0))*(_X(u,1)-_X(u,0))+(_X(u,0)-_X(u,2))*(_Y(u,1)-_Y(u,0)));
    dxpy = (_X(u,0)-_X(u,2))/((_Y(u,2)-_Y(u,0))*(_X(u,1)-_X(u,0))+(_X(u,0)-_X(u,2))*(_Y(u,1)-_Y(u,0)));
    dypx = (_Y(u,1)-_Y(u,0))/((_X(u,2)-_X(u,0))*(_Y(u,1)-_Y(u,0))+(_X(u,0)-_X(u,1))*(_Y(u,2)-_Y(u,0)));
    dypy = (_X(u,0)-_X(u,1))/((_X(u,2)-_X(u,0))*(_Y(u,1)-_Y(u,0))+(_X(u,0)-_X(u,1))*(_Y(u,2)-_Y(u,0)));


    for (int i = 0 ; i < 3 ; i++)
    {
      for (int j = 0 ; j < 3 ; j++)
      {

	if(i==j){
	  _M.coeffRef(_Sommets(u,i),_Sommets(u,j)) += _tri_area(u)/6.;
	  if(i==0){
	    _K.coeffRef(_Sommets(u,i),_Sommets(u,j)) += _tri_area(u)*(pow(dxpx+dypx,2)+pow(dxpy+dypy,2));
	  }
	  else if(i==1){
	    _K.coeffRef(_Sommets(u,i),_Sommets(u,j)) += _tri_area(u)*(dxpx*dxpx+dxpy*dxpy);
	  }
	  else if(i==2){
	    _K.coeffRef(_Sommets(u,i),_Sommets(u,j)) += _tri_area(u)*(dypx*dypx+dypy*dypy);
	  }
	}

	else{
	  _M.coeffRef(_Sommets(u,i),_Sommets(u,j)) += _tri_area(u)/12.;
	  if(((i==0)and(j==1))or((j==0)and(i==1))){
	    _K.coeffRef(_Sommets(u,i),_Sommets(u,j)) += _tri_area(u)*((-dxpx-dypx)*dxpx+(-dxpy-dypy)*dxpy);
	  }
	  if(((i==0)and(j==2))or((j==0)and(i==2))){
	    _K.coeffRef(_Sommets(u,i),_Sommets(u,j)) += _tri_area(u)*((-dxpx-dypx)*dypx+(-dxpy-dypy)*dypy);
	  }
	  if(((i==1)and(j==2))or((j==1)and(i==2))){
	    _K.coeffRef(_Sommets(u,i),_Sommets(u,j)) += _tri_area(u)*(dxpx*dypx+dxpy*dypy);
	  }
	}

      }
    }
  }
  Eigen::SparseMatrix<double, Eigen::ColMajor> A(_M+_v*_v*dt*dt*_K);
  cout << "-------------------------------------------------" << endl;
  cout << "Modification des limites de A : " << endl;
  BoundaryConditionA(A);
  _solver.analyzePattern(A);
  _solver.factorize(A);
}

void Calcul::BoundaryConditionA(SparseMatrix<double, ColMajor> A)
{
  SparseVector<double> V(nb_pts);
  V.setZero();
  A = A.transpose();

  for (int i = 0; i < _edges.size(); i++)
  {
    cout.flush();
    cout << "Progression : " << (double)i/(double)_edges.size()*100 << "%       \r";

    if ((_edges[i].GetReference()==4)or(_edges[i].GetReference()==2))//or(_edges[i].GetReference()==1)or(_edges[i].GetReference()==3))
    {
      A.col(_edges[i].GetEdge()[0]) = V;
      A.coeffRef(_edges[i].GetEdge()[0],_edges[i].GetEdge()[0]) = 1.;
    }
  }
}

void Calcul::BoundaryConditionB(SparseVector<double>& B, const double t)
{
  for (int i = 0; i < _edges.size(); i++)
  {
//////////////////////////////////////Côté gauche/////////////////////////
    if (_edges[i].GetReference()==4)
    {
      B.coeffRef(_edges[i].GetEdge()[0]) = f(((_vertices[_edges[i].GetEdge()[0]]).GetVertice())[0],((_vertices[_edges[i].GetEdge()[0]]).GetVertice())[1],t);
      if (_edges[i].GetEdge()[0] == 0)
      {
        B.coeffRef(_edges[i].GetEdge()[1]) = f(((_vertices[_edges[i].GetEdge()[1]]).GetVertice())[0],((_vertices[_edges[i].GetEdge()[1]]).GetVertice())[1],t);
      }
    }
//////////////////////////////////////Côté droit//////////////////////////
    else if (_edges[i].GetReference()==2)
    {
      B.coeffRef(_edges[i].GetEdge()[0]) = 0.;
      if (_edges[i].GetEdge()[0] == 2)
      {
        B.coeffRef(_edges[i].GetEdge()[1]) = 0.;
      }
    }
//////////////////////////////////////Côté haut///////////////////////////
    else if (_edges[i].GetReference()==3)
    {
      B.coeffRef(_edges[i].GetEdge()[0]) = 0.;
      if (_edges[i].GetEdge()[0] == 3)
      {
        B.coeffRef(_edges[i].GetEdge()[1]) = 0.;
      }
    }
//////////////////////////////////////Côté bas////////////////////////////
    else if (_edges[i].GetReference()==1)
    {
      B.coeffRef(_edges[i].GetEdge()[0]) = 0.;
      if (_edges[i].GetEdge()[0] == 1)
      {
        B.coeffRef(_edges[i].GetEdge()[1]) = 0.;
      }
    }
  }
}


void Calcul::BuildF(double t,double dt,SparseVector<double> & sol,SparseVector<double> & sol1)
{
  VectorXd dX0(_vertices.size());
/*
  SparseVector<double> B(_vertices.size());
  B = SecondMember(t);
  B.setZero();
*/
  _f.resize(_vertices.size());
  _f.setZero();

  SparseVector<double> C(_vertices.size());

  C = _M*(2*sol - sol1);// + dt*dt*B;

  BoundaryConditionB(C, t);

  VectorXd dC(_vertices.size());
  dC = MatrixXd(C);

  dX0 = _solver.solve(dC);

  SparseVector<double> spX0;
  spX0 = dX0.sparseView();

  _f = spX0;
}


SparseVector<double> Calcul::InitialCondition(const double t)
{
	SparseVector<double> sol0(_vertices.size());
  sol0.setZero();
  for (int i = 0; i < _vertices.size(); i++)
  {
  	sol0.coeffRef(i) = FunctionInitial(((_vertices[i]).GetVertice())[0],((_vertices[i]).GetVertice())[1],t);
  }
	return sol0;
}

double Calcul::FunctionInitial(const double x, const double y, const double t) const
{
  double _x0 = 0;
  double _y0 = 0;/*
  if (((x-_x0)*(x-_x0)+(y-_y0)*(y-_y0))<1)
  {
    return 1.;
  }
  else
  {
    return 0.;
  }*/
  return 0.;
}


SparseVector<double> Calcul::SecondMember(const double t)
{
  SparseVector<double> B(_vertices.size());
  for (int u = 0; u < _vertices.size(); u++)
  {
    for (int i = 0 ; i < 3 ; i++)
    {
      B.coeffRef(u/*_Sommets(u,i)*/) = f(((_vertices[u]).GetVertice())[0],((_vertices[u]).GetVertice())[1],t);//_tri_center(u,0),_tri_center(u,1),t)*Phi(_tri_center(u,0),_tri_center(u,1),u,i)*_tri_area(u);
    }
  }
  return B;
}

double Calcul::f(const double x, const double y, const double t) const
{
  double pi = 3.1415;
  if (t<0.000005)
  {
    return sin(t*2.*pi/0.000005)*exp(-(t-0.0000025)*(t-0.0000025));
  }
  else
  {
    return 0.;
  }
}

// Sauvegarde la solution
void Calcul::SaveSol(SparseVector<double> sol, int n, const std::string nom)
{
	string name_file = "Resultats/" + nom + std::to_string(n) + ".vtk";

  int nb_vert = _vertices.size();

  assert((sol.size() == _vertices.size()) && "The size of the solution vector is not the same than the number of _triangles !");

	ofstream solution;
	solution.open(name_file, ios::out);
	solution.precision(7);

  solution << "# vtk DataFile Version 3.0 " << endl;
  solution << "2D Unstructured Grid" << endl;
  solution << "ASCII" << endl;
  solution << "DATASET UNSTRUCTURED_GRID" << endl;

  solution << "POINTS " << nb_vert << " float " << endl;
  for (int i = 0 ; i < nb_vert ; ++i)
  {
    solution << ((_vertices[i]).GetVertice())[0] << " " << ((_vertices[i]).GetVertice())[1] << " 0." << endl;
  }
  solution << endl;

  solution << "CELLS " << _triangles.size() << " " << _triangles.size()*4 << endl;
  for (int i = 0 ; i < _triangles.size() ; ++i)
  {
    solution << 3 << " " << ((_triangles[i]).GetTriangle())[0] << " " << ((_triangles[i]).GetTriangle())[1]
    << " " << ((_triangles[i]).GetTriangle())[2] << endl;
  }
  solution << endl;

  solution << "CELL_TYPES " << _triangles.size() << endl;
  for (int i = 0 ; i < _triangles.size() ; ++i)
  {
    solution << 5 << endl;
  }
  solution << endl;

  solution << "POINT_DATA " << _vertices.size() << endl;
  solution << "SCALARS z float" << endl;
  solution << "LOOKUP_TABLE default" << endl;
  for (int i = 0 ; i < _vertices.size() ; ++i)
  {
    solution << float(sol.coeffRef(i)) << endl;
  }
  solution << endl;

	solution.close();
}

void Calcul::SaveSensor(SparseVector<double> sol, double t)
{

  	ofstream sensor;
  	sensor.open("Resultats/sensor.dat", ios::app);
  	sensor.precision(7);

    double val1_x(0.),val2_x(0.), x, y, iter1(0.), iter2(0.);
    double xg(-4.9),yg(0.),xd(4.9),yd(0.),r(0.05);
    for (int i = 0 ; i < nb_pts ; ++i)
    {
      x = ((_vertices[i]).GetVertice())[0];
      y = ((_vertices[i]).GetVertice())[1];
      if ((x-xg)*(x-xg)+(y-yg)*(y-yg) < r)
      {
        val1_x += float(sol.coeffRef(i));
        iter1++;
      }
      if ((x-xd)*(x-xd)+(y-yd)*(y-yd) < r)
      {
        val2_x += float(sol.coeffRef(i));
        iter2++;
      }
    }
    sensor << t << " " << val1_x/iter1 << " " << val2_x/iter2 << endl;

  	sensor.close();
}


#define _MODULE_CPP
#endif
