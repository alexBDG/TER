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
  _maillage = maillage;
  _Sommets = maillage->GetP();
  _X = maillage->GetX();
  _Y = maillage->GetY();
  nb_pts = _vertices.size();
  _lambda = 2.;//59000000000;
  _mu = 1.;//26000000000;
  _rho = 1.;//2700.;
  _C12 = _lambda;
  _C66 = _mu;
  _C11 = _C12 + 2*_C66;
  _C22 = _C11;
  system("mkdir -p ./Resultats");
}


void Calcul::BuildMK(double dt)
{
  _M.resize(2*nb_pts,2*nb_pts);
  _K.resize(2*nb_pts,2*nb_pts);

  _M.setZero();
  _K.setZero();

  cout << "Création de M et K : " << endl;
  for (int u = 0; u < _triangles.size(); u++)
  {
    cout.flush();
    cout << "Progression : " << (double)u/((double)_triangles.size()-1.)*100 << "% \r";

    //Dérivés dues au changement de variables
    double dxpx, dxpy, dypx, dypy;
    dxpx = (_Y(u,2)-_Y(u,0))/((_Y(u,2)-_Y(u,0))*(_X(u,1)-_X(u,0))+(_X(u,0)-_X(u,2))*(_Y(u,1)-_Y(u,0)));
    dxpy = (_X(u,0)-_X(u,2))/((_Y(u,2)-_Y(u,0))*(_X(u,1)-_X(u,0))+(_X(u,0)-_X(u,2))*(_Y(u,1)-_Y(u,0)));
    dypx = (_Y(u,1)-_Y(u,0))/((_X(u,2)-_X(u,0))*(_Y(u,1)-_Y(u,0))+(_X(u,0)-_X(u,1))*(_Y(u,2)-_Y(u,0)));
    dypy = (_X(u,0)-_X(u,1))/((_X(u,2)-_X(u,0))*(_Y(u,1)-_Y(u,0))+(_X(u,0)-_X(u,1))*(_Y(u,2)-_Y(u,0)));


    for (int i = 0 ; i < 3 ; i++)
    {
      for (int j = 0 ; j < 3 ; j++)
      {

//METHODE 3/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if(i==j)
        {
          _M.coeffRef(_Sommets(u,i),_Sommets(u,j)) += _tri_area(u)/6.;
          _M.coeffRef(_Sommets(u,i)+nb_pts,_Sommets(u,j)+nb_pts) += _tri_area(u)/6.;
          if(i==0)
          {
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)) += _tri_area(u)*((2*_mu+_lambda)*(pow(dxpx+dypx,2)+_mu*pow(dxpy+dypy,2)));
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)+nb_pts) += _tri_area(u)*(_lambda+_mu)*(dxpx+dypx)*(dxpy+dypy);
            _K.coeffRef(_Sommets(u,i)+nb_pts,_Sommets(u,j)+nb_pts) += _tri_area(u)*((2*_mu+_lambda)*(pow(dxpy+dypy,2)+_mu*pow(dxpx+dypx,2)));
            _K.coeffRef(_Sommets(u,i)+nb_pts,_Sommets(u,j)) += _tri_area(u)*(_lambda+_mu)*(dxpx+dypx)*(dxpy+dypy);
          }
          else if(i==1)
          {
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)) += _tri_area(u)*((2*_mu+_lambda)*dxpx*dxpx+_mu*dxpy*dxpy);
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)+nb_pts) += _tri_area(u)*(_lambda+_mu)*dxpx*dxpy;
            _K.coeffRef(_Sommets(u,i)+nb_pts,_Sommets(u,j)+nb_pts) += _tri_area(u)*((2*_mu+_lambda)*dxpy*dxpy+_mu*dxpx*dxpx);
            _K.coeffRef(_Sommets(u,i)+nb_pts,_Sommets(u,j)) += _tri_area(u)*(_lambda+_mu)*dxpx*dxpy;
          }
          else if(i==2)
          {
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)) += _tri_area(u)*((2*_mu+_lambda)*dypx*dypx+_mu*dypy*dypy);
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)+nb_pts) += _tri_area(u)*(_lambda+_mu)*dypx*dypy;
            _K.coeffRef(_Sommets(u,i)+nb_pts,_Sommets(u,j)+nb_pts) += _tri_area(u)*((2*_mu+_lambda)*dypy*dypy+_mu*dypx*dypx);
            _K.coeffRef(_Sommets(u,i)+nb_pts,_Sommets(u,j)) += _tri_area(u)*(_lambda+_mu)*dypx*dypy;
          }
        }

        else
        {
          _M.coeffRef(_Sommets(u,i),_Sommets(u,j)) += _tri_area(u)/12.;
          _M.coeffRef(_Sommets(u,i)+nb_pts,_Sommets(u,j)+nb_pts) += _tri_area(u)/12.;
          if((i==0)and(j==1))
          {
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)) += _tri_area(u)*((2*_mu+_lambda)*(-dxpx-dypx)*dxpx+_mu*(-dxpy-dypy)*dxpy);
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)+nb_pts) += 0.5*_tri_area(u)*(_lambda*(-dxpy-dypy)*dxpx+_mu*(-dxpx-dypx)*dxpy+_lambda*(-dxpx-dypx)*dxpy+_mu*(-dxpy-dypy)*dxpx);
            _K.coeffRef(_Sommets(u,i)+nb_pts,_Sommets(u,j)+nb_pts) += _tri_area(u)*((2*_mu+_lambda)*(-dxpy-dypy)*dxpy+_mu*(-dxpx-dypx)*dxpx);
            _K.coeffRef(_Sommets(u,i)+nb_pts,_Sommets(u,j)) += 0.5*_tri_area(u)*(_lambda*(-dxpx-dypx)*dxpy+_mu*(-dxpy-dypy)*dxpx+_lambda*(-dxpy-dypy)*dxpx+_mu*(-dxpx-dypx)*dxpy);
          }
          if((i==1)and(j==0))
          {
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)) += _tri_area(u)*((2*_mu+_lambda)*(-dxpx-dypx)*dxpx+_mu*(-dxpy-dypy)*dxpy);
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)+nb_pts) += 0.5*_tri_area(u)*(_lambda*(-dxpx-dypx)*dxpy+_mu*(-dxpy-dypy)*dxpx+_lambda*(-dxpy-dypy)*dxpx+_mu*(-dxpx-dypx)*dxpy);
            _K.coeffRef(_Sommets(u,i)+nb_pts,_Sommets(u,j)+nb_pts) += _tri_area(u)*((2*_mu+_lambda)*(-dxpy-dypy)*dxpy+_mu*(-dxpx-dypx)*dxpx);
            _K.coeffRef(_Sommets(u,i)+nb_pts,_Sommets(u,j)) += 0.5*_tri_area(u)*(_lambda*(-dxpy-dypy)*dxpx+_mu*(-dxpx-dypx)*dxpy+_lambda*(-dxpx-dypx)*dxpy+_mu*(-dxpy-dypy)*dxpx);
          }
          if((i==0)and(j==2))
          {
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)) += _tri_area(u)*((2*_mu+_lambda)*(-dxpx-dypx)*dypx+_mu*(-dxpy-dypy)*dypy);
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)+nb_pts) += 0.5*_tri_area(u)*(_lambda*(-dxpy-dypy)*dypx+_mu*(-dxpx-dypx)*dypy+_lambda*(-dxpx-dypx)*dypy+_mu*(-dxpy-dypy)*dypx);
            _K.coeffRef(_Sommets(u,i)+nb_pts,_Sommets(u,j)+nb_pts) += _tri_area(u)*((2*_mu+_lambda)*(-dxpy-dypy)*dypy+_mu*(-dxpx-dypx)*dypx);
            _K.coeffRef(_Sommets(u,i)+nb_pts,_Sommets(u,j)) += 0.5*_tri_area(u)*(_lambda*(-dxpx-dypx)*dypy+_mu*(-dxpy-dypy)*dypx+_lambda*(-dxpy-dypy)*dypx+_mu*(-dxpx-dypx)*dypy);
          }
          if((i==2)and(j==0))
          {
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)) += _tri_area(u)*((2*_mu+_lambda)*(-dxpx-dypx)*dypx+_mu*(-dxpy-dypy)*dypy);
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)+nb_pts) += 0.5*_tri_area(u)*(_lambda*(-dxpx-dypx)*dypy+_mu*(-dxpy-dypy)*dypx+_lambda*(-dxpy-dypy)*dypx+_mu*(-dxpx-dypx)*dypy);
            _K.coeffRef(_Sommets(u,i)+nb_pts,_Sommets(u,j)+nb_pts) += _tri_area(u)*((2*_mu+_lambda)*(-dxpy-dypy)*dypy+_mu*(-dxpx-dypx)*dypx);
            _K.coeffRef(_Sommets(u,i)+nb_pts,_Sommets(u,j)) += 0.5*_tri_area(u)*(_lambda*(-dxpy-dypy)*dypx+_mu*(-dxpx-dypx)*dypy+_lambda*(-dxpx-dypx)*dypy+_mu*(-dxpy-dypy)*dypx);
          }
          if((i==1)and(j==2))
          {
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)) += _tri_area(u)*((2*_mu+_lambda)*dxpx*dypx+_mu*dxpy*dypy);
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)+nb_pts) += 0.5*_tri_area(u)*(_lambda*dxpy*dypx+_mu*dxpx*dypy+_lambda*dxpx*dypy+_mu*dxpy*dypx);
            _K.coeffRef(_Sommets(u,i)+nb_pts,_Sommets(u,j)+nb_pts) += _tri_area(u)*((2*_mu+_lambda)*dxpy*dypy+_mu*dxpx*dypx);
            _K.coeffRef(_Sommets(u,i)+nb_pts,_Sommets(u,j)) += 0.5*_tri_area(u)*(_lambda*dxpx*dypy+_mu*dxpy*dypx+_lambda*dxpy*dypx+_mu*dxpx*dypy);
          }
          if((i==2)and(j==1))
          {
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)) += _tri_area(u)*((2*_mu+_lambda)*dxpx*dypx+_mu*dxpy*dypy);
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)+nb_pts) += 0.5*_tri_area(u)*(_lambda*dxpx*dypy+_mu*dxpy*dypx+_lambda*dxpy*dypx+_mu*dxpx*dypy);
            _K.coeffRef(_Sommets(u,i)+nb_pts,_Sommets(u,j)+nb_pts) += _tri_area(u)*((2*_mu+_lambda)*dxpy*dypy+_mu*dxpx*dypx);
            _K.coeffRef(_Sommets(u,i)+nb_pts,_Sommets(u,j)) += 0.5*_tri_area(u)*(_lambda*dxpy*dypx+_mu*dxpx*dypy+_lambda*dxpx*dypy+_mu*dxpy*dypx);
          }
        }
      }
    }
  }
  _solver.compute(_M+dt*dt*(1./_rho)*_K);
}


void Calcul::BuildF(double t,double dt,SparseVector<double> & sol,SparseVector<double> & sol1)
{
  VectorXd dX0(2*nb_pts);

//  SparseVector<double> B(2*nb_pts);
//  B = SecondMember(t);

  _f.resize(2*nb_pts);
  _f.setZero();

  SparseVector<double> C(2*nb_pts);

  C = _M*(2*sol - sol1);// + dt*dt*(1./_rho)*B);

  VectorXd dC(2*nb_pts);
  dC = MatrixXd(C);

  dX0 = _solver.solve(dC);

  SparseVector<double> spX0;
  spX0 = dX0.sparseView();

  _f = spX0;
}


SparseVector<double> Calcul::InitialCondition(const double t)
{
	SparseVector<double> sol0(2*nb_pts);
  sol0.setZero();
  //CERCLE
  for (int i = 0; i < nb_pts; i++)
  {
    double X = ((_vertices[i]).GetVertice())[0];
    double Y = ((_vertices[i]).GetVertice())[1];
    sol0.coeffRef(i) = FunctionInitial(X,Y,t)*X/sqrt(X*X+Y*Y);
    sol0.coeffRef(i+nb_pts) = FunctionInitial(X,Y,t)*Y/sqrt(X*X+Y*Y);
  }/*
  //BANDE
  for (int i = 0; i < nb_pts; i++)
  {
    double X = ((_vertices[i]).GetVertice())[0];
    double Y = ((_vertices[i]).GetVertice())[1];
    sol0.coeffRef(i) = FunctionInitial(X,Y,t);
  }*/
	return sol0;
}

double Calcul::FunctionInitial(const double x, const double y, double t) const
{
  //CERCLE
  double _x0 = 0.;
  double _y0 = 0.;
  if (((x-_x0)*(x-_x0)+(y-_y0)*(y-_y0))<1)
  {
    return 1.;
  }
  else
  {
    return 0.;
  }/*
  //BANDE
  if (x<=-9.)
  {
    return 1.;
  }
  else
  {
    return 0.;
  }*/
}


SparseVector<double> Calcul::SecondMember(const double t)
{
  SparseVector<double> B(2*nb_pts);
  for (int u = 0; u < nb_pts; u++)
  {
    B.coeffRef(u) = f(_vertices[u].GetVertice()[0],_vertices[u].GetVertice()[1],t);
    B.coeffRef(u+nb_pts) = f(_vertices[u].GetVertice()[0],_vertices[u].GetVertice()[1],t);
  }
  return B;
}

SparseVector<double> Calcul::SolutionExacte(const double t)
{
  SparseVector<double> B(2*nb_pts);
  for (int u = 0; u < nb_pts; u++)
  {
    B.coeffRef(u) = fexacte(_vertices[u].GetVertice()[0],_vertices[u].GetVertice()[1],t);
    B.coeffRef(u+nb_pts) = fexacte(_vertices[u].GetVertice()[0],_vertices[u].GetVertice()[1],t);
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

  assert((sol.size() == 2*nb_pts) && "The size of the solution vector is not the same than the number of 2 * _vertices !");

	ofstream solution;
	solution.open(name_file, ios::out);
	solution.precision(7);

  solution << "# vtk DataFile Version 3.0 " << endl;
  solution << "2D Unstructured Grid" << endl;
  solution << "ASCII" << endl;
  solution << "DATASET UNSTRUCTURED_GRID" << endl;

  solution << "POINTS " << nb_pts << " float " << endl;
  for (int i = 0 ; i < nb_pts ; ++i)
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

  solution << "POINT_DATA " << nb_pts << endl;
  solution << "SCALARS sol float" << endl;
  solution << "LOOKUP_TABLE default" << endl;
  for (int i = 0 ; i < nb_pts ; ++i)
  {
    solution << float(sqrt(pow(sol.coeffRef(i),2)+pow(sol.coeffRef(i+nb_pts),2))) << endl;
  }
  solution << endl;

  solution << "VECTORS u float" << endl;
  for (int i = 0 ; i < nb_pts ; ++i)
  {
    solution << float(sol.coeffRef(i)) << " " << float(sol.coeffRef(i+nb_pts)) << " 10" << endl;
  }
  solution << endl;

	solution.close();
}

#define _MODULE_CPP
#endif
