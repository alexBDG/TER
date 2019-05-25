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
_tri_center(maillage->GetTrianglesCenter())
{
  _c = 100.;
  _v = 1.5;
  _maillage = maillage;
 
  _Sommets = maillage->GetP();
  _X = maillage->GetX();
  _Y = maillage->GetY();
  system("mkdir -p ./Resultats");
  ofstream norme;
  norme.open("Resultats/norme.dat", ios::out);
  norme.close();
}


void Calcul::BuildMK(double dt)
{
  
  _M.resize(_vertices.size(),_vertices.size());
  _K.resize(_vertices.size(),_vertices.size());
  ofstream flux("matrice");
  
  cout << "Création de M et K : " << endl;
  for (int u = 0; u < _triangles.size(); u++)
  {
    cout.flush();
    cout << "Progression : " << (double)u/((double)_triangles.size()-1.)*100 << "% \r";

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

        if(i==j)
        {
          _M.coeffRef(_Sommets(u,i),_Sommets(u,j)) += _tri_area(u)/6.;
          if(i==0)
          {
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)) += _tri_area(u)*(pow(dxpx+dypx,2)+pow(dxpy+dypy,2));
          }
          else if(i==1)
          {
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)) += _tri_area(u)*(dxpx*dxpx+dxpy*dxpy);
          }
          else if(i==2)
          {
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)) += _tri_area(u)*(dypx*dypx+dypy*dypy);
          }
        }

        else
        {
          _M.coeffRef(_Sommets(u,i),_Sommets(u,j)) += _tri_area(u)/12.;
          if(((i==0)and(j==1))or((j==0)and(i==1)))
          {
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)) += _tri_area(u)*((-dxpx-dypx)*dxpx+(-dxpy-dypy)*dxpy);
          }
          if(((i==0)and(j==2))or((j==0)and(i==2)))
          {
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)) += _tri_area(u)*((-dxpx-dypx)*dypx+(-dxpy-dypy)*dypy);
          }
          if(((i==1)and(j==2))or((j==1)and(i==2)))
          {
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)) += _tri_area(u)*(dxpx*dypx+dxpy*dypy);
          }


	}
	flux << _Sommets(u,i) << " " << 5533-_Sommets(u,j) << endl;
      }

    }
    
    

    
  }

  auto start_mk = chrono::high_resolution_clock::now();
  _solver.compute(_M+_v*_v*dt*dt*_K);
  auto finish_mk = chrono::high_resolution_clock::now();
  double t_mk = chrono::duration_cast<chrono::milliseconds>(finish_mk-start_mk).count();
  cout << "\n decomposition de choleky dure " << t_mk/1000. << endl;
}


void Calcul::BuildF(double t,double dt,SparseVector<double> & sol,SparseVector<double> & sol1)
{
  VectorXd dX0(_vertices.size());

  SparseVector<double> B(_vertices.size());
  B = SecondMember(t);

  _f.resize(_vertices.size());
  _f.setZero();

  SparseVector<double> C(_vertices.size());

  C = _M*(2*sol - sol1 + dt*dt*B);

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
  	sol0.coeffRef(i) = fexacte(((_vertices[i]).GetVertice())[0],((_vertices[i]).GetVertice())[1],t);
  }
	return sol0;
}

double Calcul::FunctionInitial(const double x, const double y, double t) const
{
  double _x0 = 0;
  double _y0 = 0;
  if (x<=1)//(((x-_x0)*(x-_x0)+(y-_y0)*(y-_y0))<1)
  {
    return 1.;
  }
  else
  {
    return 0.;
  }
}


SparseVector<double> Calcul::SecondMember(const double t)
{
  SparseVector<double> B(_vertices.size());
  for (int u = 0; u < _vertices.size(); u++)
  {
    B.coeffRef(u) = f(_vertices[u].GetVertice()[0],_vertices[u].GetVertice()[1],t);
  }
  return B;
}

SparseVector<double> Calcul::SolutionExacte(const double t)
{
  SparseVector<double> B(_vertices.size());
  for (int u = 0; u < _vertices.size(); u++)
  {
    B.coeffRef(u) = fexacte(_vertices[u].GetVertice()[0],_vertices[u].GetVertice()[1],t);
  }
  return B;
}

double Calcul::fexacte(const double x, const double y, const double t) const
{
  double _x0 = 0.;
  double _y0 = 0.;
  double pi = 3.1415;
  double a = 1.;
  double b = 1.;
  return cos((x*x-100.)/_c)*cos((y*y-100.)/_c)*exp(-t);
}

double Calcul::f(const double x, const double y, const double t) const
{
  double _x0 = 0.;
  double _y0 = 0.;
  double pi = 3.1415;
  double a = 1.;
  double b = 1.;
  return cos((x*x-100.)/_c)*cos((y*y-100.)/_c)*exp(-t)+(2*exp(-t)*(_c*sin((x*x - 100.)/_c)*cos((y*y - 100.)/_c) + cos((x*x - 100.)/_c)*(2.*(x*x + y*y)*cos((y*y - 100.)/_c) + _c*sin((y*y - 100.)/_c))))/(_c*_c);
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

void Calcul::SaveNorm(SparseVector<double> sol, SparseVector<double> sol_e, double t)
{
  ofstream norme;
  norme.open("Resultats/norme.dat", ios::app);
  norme.precision(7);

  norme << t << " " << (sol-sol_e).norm() << endl;
  norme.close();
}

#define _MODULE_CPP
#endif
