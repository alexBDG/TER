#ifndef _MODULE_CPP

#include "Module.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <chrono>

using namespace Eigen;
using namespace std;

Calcul::Calcul(Mesh2D* maillage, DataFile* data_file) :
_vertices(maillage->GetVertices()),_tri_area(maillage->GetTrianglesArea()),
_triangles(maillage->GetTriangles()),_tri_center(maillage->GetTrianglesCenter()),
_edges(maillage->GetEdges()), _edges_bc(maillage->GetEdges_bc()),
_initial_condition(data_file-> Get_initial_condition_choice()),
_rho(data_file->Get_param_rho()), _lambda(data_file->Get_param_lambda()),
_mu(data_file->Get_param_mu()), _mesh_name(data_file->Get_mesh_name()),
_s1(data_file->Get_sensor1()), _s2(data_file->Get_sensor2()),
_coeff_bc(data_file->Get_param()),_file_name(data_file->Get_file_name()),
_x0(data_file->Get_param_x0()), _y0(data_file->Get_param_y0()),
_r(data_file->Get_param_r()), _c(data_file->Get_param_c()),
_if_exc_u(data_file->Get_bool_excitation_u()),_if_exc_l(data_file->Get_bool_excitation_l()),
_exc(data_file->Get_param_excitation()), _eps(data_file->Get_param_eps()),
_results(data_file->Get_results())
{
  _maillage = maillage;
  _Sommets = maillage->GetP();
  _X = maillage->GetX();
  _Y = maillage->GetY();
  _nb_pts = _vertices.size();

  /*int size_temp;
  for (int k = 0; k < 4; k++)
  {
    size_temp = _edges_bc[k].size();
    for (int i = 0; i < size_temp; i++) _edges_bc[k].push_back(_edges_bc[k][i] + _nb_pts);
  }*/

  /*
  cout << "----------------------------------" << endl;
  cout << "_coeff_bc =" << endl;
  for (int i=0; i<4; i++)
  {
    for (int j = 0; j < _coeff_bc[i].size(); j++) cout << _coeff_bc[i][j] << " ";
    cout << endl;
  }
  cout << "----------------------------------" << endl;
  */

  //Creation du repertoire contenant les solutions
  system(("mkdir -p ./" + _results).c_str());

  //Creation du fichier contenant les valeurs reçues par le capteur
  if ((_s1)or(_s2))
  {
    ofstream sensor;
    sensor.open((_results + "/sensor.dat").c_str(), ios::out);
    sensor.close();

    if ((_s1)and(_s2))
    {
      _sensor.setZero(2,2);
      _sensor << data_file->Get_param_s1_x(), data_file->Get_param_s1_y(),
              data_file->Get_param_s2_x(), data_file->Get_param_s2_y();
    }
    else if ((_s1)and(!_s2))
    {
      _sensor.setZero(1,2);
      _sensor << data_file->Get_param_s1_x(), data_file->Get_param_s1_y();
    }
    else if ((!_s1)and(_s2))
    {
      _sensor.setZero(1,2);
      _sensor << data_file->Get_param_s2_x(), data_file->Get_param_s2_y();
    }
  }
  else
  {
    _sensor.resize(0,0);
  }
}


void Calcul::BuildMK(double dt)
{
  _M.resize(2*_nb_pts,2*_nb_pts);
  _K.resize(2*_nb_pts,2*_nb_pts);

  _M.setZero();
  _K.setZero();

  cout << "Creation of M and K" << endl;
  for (int u = 0; u < _triangles.size(); u++)
  {
    cout.flush();
    cout << "Progression : " << (double)u/((double)_triangles.size())*100 << "%      \r";

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
          _M.coeffRef(_Sommets(u,i)+_nb_pts,_Sommets(u,j)+_nb_pts) += _tri_area(u)/6.;
          if(i==0)
          {
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)) += _tri_area(u)*((2*_mu+_lambda)*(pow(dxpx+dypx,2)+_mu*pow(dxpy+dypy,2)));
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)+_nb_pts) += _tri_area(u)*(_lambda+_mu)*(dxpx+dypx)*(dxpy+dypy);
            _K.coeffRef(_Sommets(u,i)+_nb_pts,_Sommets(u,j)+_nb_pts) += _tri_area(u)*((2*_mu+_lambda)*(pow(dxpy+dypy,2)+_mu*pow(dxpx+dypx,2)));
            _K.coeffRef(_Sommets(u,i)+_nb_pts,_Sommets(u,j)) += _tri_area(u)*(_lambda+_mu)*(dxpx+dypx)*(dxpy+dypy);
          }
          else if(i==1)
          {
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)) += _tri_area(u)*((2*_mu+_lambda)*dxpx*dxpx+_mu*dxpy*dxpy);
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)+_nb_pts) += _tri_area(u)*(_lambda+_mu)*dxpx*dxpy;
            _K.coeffRef(_Sommets(u,i)+_nb_pts,_Sommets(u,j)+_nb_pts) += _tri_area(u)*((2*_mu+_lambda)*dxpy*dxpy+_mu*dxpx*dxpx);
            _K.coeffRef(_Sommets(u,i)+_nb_pts,_Sommets(u,j)) += _tri_area(u)*(_lambda+_mu)*dxpx*dxpy;
          }
          else if(i==2)
          {
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)) += _tri_area(u)*((2*_mu+_lambda)*dypx*dypx+_mu*dypy*dypy);
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)+_nb_pts) += _tri_area(u)*(_lambda+_mu)*dypx*dypy;
            _K.coeffRef(_Sommets(u,i)+_nb_pts,_Sommets(u,j)+_nb_pts) += _tri_area(u)*((2*_mu+_lambda)*dypy*dypy+_mu*dypx*dypx);
            _K.coeffRef(_Sommets(u,i)+_nb_pts,_Sommets(u,j)) += _tri_area(u)*(_lambda+_mu)*dypx*dypy;
          }
        }

        else
        {
          _M.coeffRef(_Sommets(u,i),_Sommets(u,j)) += _tri_area(u)/12.;
          _M.coeffRef(_Sommets(u,i)+_nb_pts,_Sommets(u,j)+_nb_pts) += _tri_area(u)/12.;
          if((i==0)and(j==1))
          {
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)) += _tri_area(u)*((2*_mu+_lambda)*(-dxpx-dypx)*dxpx+_mu*(-dxpy-dypy)*dxpy);
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)+_nb_pts) += 0.5*_tri_area(u)*(_lambda*(-dxpy-dypy)*dxpx+_mu*(-dxpx-dypx)*dxpy+_lambda*(-dxpx-dypx)*dxpy+_mu*(-dxpy-dypy)*dxpx);
            _K.coeffRef(_Sommets(u,i)+_nb_pts,_Sommets(u,j)+_nb_pts) += _tri_area(u)*((2*_mu+_lambda)*(-dxpy-dypy)*dxpy+_mu*(-dxpx-dypx)*dxpx);
            _K.coeffRef(_Sommets(u,i)+_nb_pts,_Sommets(u,j)) += 0.5*_tri_area(u)*(_lambda*(-dxpx-dypx)*dxpy+_mu*(-dxpy-dypy)*dxpx+_lambda*(-dxpy-dypy)*dxpx+_mu*(-dxpx-dypx)*dxpy);
          }
          if((i==1)and(j==0))
          {
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)) += _tri_area(u)*((2*_mu+_lambda)*(-dxpx-dypx)*dxpx+_mu*(-dxpy-dypy)*dxpy);
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)+_nb_pts) += 0.5*_tri_area(u)*(_lambda*(-dxpx-dypx)*dxpy+_mu*(-dxpy-dypy)*dxpx+_lambda*(-dxpy-dypy)*dxpx+_mu*(-dxpx-dypx)*dxpy);
            _K.coeffRef(_Sommets(u,i)+_nb_pts,_Sommets(u,j)+_nb_pts) += _tri_area(u)*((2*_mu+_lambda)*(-dxpy-dypy)*dxpy+_mu*(-dxpx-dypx)*dxpx);
            _K.coeffRef(_Sommets(u,i)+_nb_pts,_Sommets(u,j)) += 0.5*_tri_area(u)*(_lambda*(-dxpy-dypy)*dxpx+_mu*(-dxpx-dypx)*dxpy+_lambda*(-dxpx-dypx)*dxpy+_mu*(-dxpy-dypy)*dxpx);
          }
          if((i==0)and(j==2))
          {
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)) += _tri_area(u)*((2*_mu+_lambda)*(-dxpx-dypx)*dypx+_mu*(-dxpy-dypy)*dypy);
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)+_nb_pts) += 0.5*_tri_area(u)*(_lambda*(-dxpy-dypy)*dypx+_mu*(-dxpx-dypx)*dypy+_lambda*(-dxpx-dypx)*dypy+_mu*(-dxpy-dypy)*dypx);
            _K.coeffRef(_Sommets(u,i)+_nb_pts,_Sommets(u,j)+_nb_pts) += _tri_area(u)*((2*_mu+_lambda)*(-dxpy-dypy)*dypy+_mu*(-dxpx-dypx)*dypx);
            _K.coeffRef(_Sommets(u,i)+_nb_pts,_Sommets(u,j)) += 0.5*_tri_area(u)*(_lambda*(-dxpx-dypx)*dypy+_mu*(-dxpy-dypy)*dypx+_lambda*(-dxpy-dypy)*dypx+_mu*(-dxpx-dypx)*dypy);
          }
          if((i==2)and(j==0))
          {
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)) += _tri_area(u)*((2*_mu+_lambda)*(-dxpx-dypx)*dypx+_mu*(-dxpy-dypy)*dypy);
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)+_nb_pts) += 0.5*_tri_area(u)*(_lambda*(-dxpx-dypx)*dypy+_mu*(-dxpy-dypy)*dypx+_lambda*(-dxpy-dypy)*dypx+_mu*(-dxpx-dypx)*dypy);
            _K.coeffRef(_Sommets(u,i)+_nb_pts,_Sommets(u,j)+_nb_pts) += _tri_area(u)*((2*_mu+_lambda)*(-dxpy-dypy)*dypy+_mu*(-dxpx-dypx)*dypx);
            _K.coeffRef(_Sommets(u,i)+_nb_pts,_Sommets(u,j)) += 0.5*_tri_area(u)*(_lambda*(-dxpy-dypy)*dypx+_mu*(-dxpx-dypx)*dypy+_lambda*(-dxpx-dypx)*dypy+_mu*(-dxpy-dypy)*dypx);
          }
          if((i==1)and(j==2))
          {
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)) += _tri_area(u)*((2*_mu+_lambda)*dxpx*dypx+_mu*dxpy*dypy);
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)+_nb_pts) += 0.5*_tri_area(u)*(_lambda*dxpy*dypx+_mu*dxpx*dypy+_lambda*dxpx*dypy+_mu*dxpy*dypx);
            _K.coeffRef(_Sommets(u,i)+_nb_pts,_Sommets(u,j)+_nb_pts) += _tri_area(u)*((2*_mu+_lambda)*dxpy*dypy+_mu*dxpx*dypx);
            _K.coeffRef(_Sommets(u,i)+_nb_pts,_Sommets(u,j)) += 0.5*_tri_area(u)*(_lambda*dxpx*dypy+_mu*dxpy*dypx+_lambda*dxpy*dypx+_mu*dxpx*dypy);
          }
          if((i==2)and(j==1))
          {
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)) += _tri_area(u)*((2*_mu+_lambda)*dxpx*dypx+_mu*dxpy*dypy);
            _K.coeffRef(_Sommets(u,i),_Sommets(u,j)+_nb_pts) += 0.5*_tri_area(u)*(_lambda*dxpx*dypy+_mu*dxpy*dypx+_lambda*dxpy*dypx+_mu*dxpx*dypy);
            _K.coeffRef(_Sommets(u,i)+_nb_pts,_Sommets(u,j)+_nb_pts) += _tri_area(u)*((2*_mu+_lambda)*dxpy*dypy+_mu*dxpx*dypx);
            _K.coeffRef(_Sommets(u,i)+_nb_pts,_Sommets(u,j)) += 0.5*_tri_area(u)*(_lambda*dxpy*dypx+_mu*dxpx*dypy+_lambda*dxpx*dypy+_mu*dxpy*dypx);
          }
        }
      }
    }
  }
  SparseMatrix<double, ColMajor> A(_M+(1./_rho)*dt*dt*_K);
  BoundaryConditionA(A);
  cout << "Building of the linear solver" << endl;
  _solver.analyzePattern(A);
  _solver.factorize(A);
}

void Calcul::BoundaryConditionA(SparseMatrix<double, ColMajor> A)
{
  SparseVector<double> V(2*_nb_pts);
  V.setZero();
  A = A.transpose();

  cout << "-------------------------------------------------" << endl;
  cout << "Modification of A depending on boundary condition" << endl;
  for (int k = 0; k< 4; k++)
  {
    cout.flush();
    cout << "Progression : " << (double)k/(double)4*100 << "%      \r";
    if(_coeff_bc[k][0] > 1)
    {
      for (int i = 0; i < _edges_bc[k].size(); i++)
      {
        A.col(_edges_bc[k][i]) = V;
        A.coeffRef(_edges_bc[k][i],_edges_bc[k][i]) = 1;

        A.col(_edges_bc[k][i] + _nb_pts) = V;
        A.coeffRef(_edges_bc[k][i] + _nb_pts,_edges_bc[k][i] + _nb_pts) = 1;
      }
    }
  }
  A = A.transpose();
  cout << "-------------------------------------------------" << endl;
}


void Calcul::BoundaryConditionB(SparseVector<double>& B, const double t)
{
  //  cout << " taille des edges : " << _edges.size() << endl;
  for (int k = 0; k< 4; k++)
  {
    if(_coeff_bc[k][0] > 30)
    {
      for (int i = 0; i < _edges_bc[k].size(); i++)
      {
        {
          if ((k==1)or(k==3))
          {
            if (_coeff_bc[k][4] > 0.)
            {
              B.coeffRef(_edges_bc[k][i]) = f(((_vertices[_edges_bc[k][i]]).GetVertice())[0],(_vertices[_edges_bc[k][i]].GetVertice())[1],t, _coeff_bc[k][1],_coeff_bc[k][2], _coeff_bc[k][3]);
              B.coeffRef(_edges_bc[k][i] + _nb_pts) = 0;
            }
            else
            {
              B.coeffRef(_edges_bc[k][i]) = 0;
              B.coeffRef(_edges_bc[k][i] + _nb_pts) = f(((_vertices[_edges_bc[k][i]]).GetVertice())[0],(_vertices[_edges_bc[k][i]].GetVertice())[1],t, _coeff_bc[k][1],_coeff_bc[k][2], _coeff_bc[k][3]);
            }
          }
          else
          {
            if (_coeff_bc[k][4] > 0.)
            {
              B.coeffRef(_edges_bc[k][i]) = 0;
              B.coeffRef(_edges_bc[k][i] + _nb_pts) = f(((_vertices[_edges_bc[k][i]]).GetVertice())[0],(_vertices[_edges_bc[k][i]].GetVertice())[1],t, _coeff_bc[k][1],_coeff_bc[k][2], _coeff_bc[k][3]);
            }
            else
            {
              B.coeffRef(_edges_bc[k][i]) = f(((_vertices[_edges_bc[k][i]]).GetVertice())[0],(_vertices[_edges_bc[k][i]].GetVertice())[1],t, _coeff_bc[k][1],_coeff_bc[k][2], _coeff_bc[k][3]);
              B.coeffRef(_edges_bc[k][i] + _nb_pts) = 0;
            }
          }
        }
      }
    }
    else if (_coeff_bc[k][0] < 1)
    {}
    else
    {
      for (int i = 0; i < _edges_bc[k].size(); i++)
      {
        if ((k==1)or(k==3))
        {
          if (_coeff_bc[k][4] > 0.)
          {
            B.coeffRef(_edges_bc[k][i]) = f(((_vertices[_edges_bc[k][i]]).GetVertice())[0],(_vertices[_edges_bc[k][i]].GetVertice())[1],t, _coeff_bc[k][1]);
            B.coeffRef(_edges_bc[k][i] + _nb_pts) = 0;
          }
          else
          {
            B.coeffRef(_edges_bc[k][i]) = 0;
            B.coeffRef(_edges_bc[k][i] + _nb_pts) = f(((_vertices[_edges_bc[k][i]]).GetVertice())[0],(_vertices[_edges_bc[k][i]].GetVertice())[1],t, _coeff_bc[k][1]);
          }
        }
        else
        {
          if (_coeff_bc[k][4] > 0.)
          {
            B.coeffRef(_edges_bc[k][i]) = 0;
            B.coeffRef(_edges_bc[k][i] + _nb_pts) = f(((_vertices[_edges_bc[k][i]]).GetVertice())[0],(_vertices[_edges_bc[k][i]].GetVertice())[1],t, _coeff_bc[k][1]);
          }
          else
          {
            B.coeffRef(_edges_bc[k][i]) = f(((_vertices[_edges_bc[k][i]]).GetVertice())[0],(_vertices[_edges_bc[k][i]].GetVertice())[1],t, _coeff_bc[k][1]);
            B.coeffRef(_edges_bc[k][i] + _nb_pts) = 0;
          }
        }
      }
    }
  }
}


void Calcul::BuildF(double t,double dt,SparseVector<double> & sol,SparseVector<double> & sol1)
{
  VectorXd dX0(2*_nb_pts);
  SparseVector<double> B(2*_nb_pts);


  _f.resize(2*_nb_pts);
  _f.setZero();
  B.setZero();
  //B = SecondMember(t);

  SparseVector<double> C(2*_nb_pts);

  C = _M*(2*sol - sol1); //+ dt*dt*(1./_rho)*B;

  BoundaryConditionB(C, t);

  VectorXd dC(2*_nb_pts);
  dC = MatrixXd(C);

  dX0 = _solver.solve(dC);

  SparseVector<double> spX0;
  spX0 = dX0.sparseView();

  _f = spX0;

  if((_if_exc_u)and(t<2*_exc[0][4])) excitation(_f, t, 2);
  if((_if_exc_l)and(t<2*_exc[1][4])) excitation(_f, t, 3);
}


SparseVector<double> Calcul::InitialCondition(const double t)
{
	SparseVector<double> sol0(2*_nb_pts);
  sol0.setZero();
  for (int i = 0; i < _nb_pts; i++)
  {
    double X = ((_vertices[i]).GetVertice())[0];
    double Y = ((_vertices[i]).GetVertice())[1];
    sol0.coeffRef(i) = InitialFunction_x(X,Y,t);
    sol0.coeffRef(i+_nb_pts) = InitialFunction_y(X,Y,t);
  }
	return sol0;
}

double Calcul::excitation(Eigen::SparseVector<double>& sol, double const t, int const k)
{
  vector<double> xy;
  xy.resize(2);
  for (int i = 0; i < _edges_bc[k].size(); i++)
  {
    xy[0] = ((_vertices[_edges_bc[k][i]]).GetVertice())[0];
    xy[1] = ((_vertices[_edges_bc[k][i]]).GetVertice())[1];



    if ((xy[k-2]>_exc[k-2][0])and(xy[k-2]<_exc[k-2][1]))
    {
      if (_exc[k-2][5] > 0.)
      {
        sol.coeffRef(_edges_bc[k][i]) = f(xy[0],xy[1],t, _exc[k-2][2],_exc[k-2][3], _exc[k-2][4]);
        sol.coeffRef(_edges_bc[k][i] + _nb_pts) = 0;
      }
      else
      {
        sol.coeffRef(_edges_bc[k][i]) = 0;
        sol.coeffRef(_edges_bc[k][i] + _nb_pts) = f(xy[0],xy[1],t, _exc[k-2][2],_exc[k-2][3], _exc[k-2][4]);
      }
    }
  }
}

double Calcul::InitialFunction_x(const double x, const double y, double t) const
{
  if (_initial_condition == "center")
  {
    if ((x-_x0)*(x-_x0)+(y-_y0)*(y-_y0) < _r) return x/sqrt(x*x+y*y);
    return 0;
  }
  else if (_initial_condition == "rectangular")
  {
    if ((x > _y0 - _r)and(x < _y0 + _r)) return 1;
    return 0;
  }
  else
  {
    return 0;
  }
}

double Calcul::InitialFunction_y(const double x, const double y, double t) const
{
  if (_initial_condition == "center")
  {
    if ((x-_x0)*(x-_x0)+(y-_y0)*(y-_y0) < _r) return y/sqrt(x*x+y*y);
    return 0;
  }
  else
  {
    return 0;
  }
}


SparseVector<double> Calcul::SecondMember(const double t)
{
  SparseVector<double> B(2*_nb_pts);
  for (int u = 0; u < _nb_pts; u++)
  {
    B.coeffRef(u) = 0.; //f(_vertices[u].GetVertice()[0],_vertices[u].GetVertice()[1],t);
    B.coeffRef(u+_nb_pts) = 0.; //f(_vertices[u].GetVertice()[0],_vertices[u].GetVertice()[1],t);
  }
  return B;
}

SparseVector<double> Calcul::SolutionExacte(const double t)
{
  SparseVector<double> B(2*_nb_pts);
  for (int u = 0; u < _nb_pts; u++)
  {
    B.coeffRef(u) = fexacte(_vertices[u].GetVertice()[0],_vertices[u].GetVertice()[1],t);
    B.coeffRef(u+_nb_pts) = fexacte(_vertices[u].GetVertice()[0],_vertices[u].GetVertice()[1],t);
  }
  return B;
}

double Calcul::fexacte(const double x, const double y, const double t) const
{
  double pi = 3.14159265358979323846;
  return cos((x*x-100.)/_c)*cos((y*y-100.)/_c)*exp(-t);
}

double Calcul::f(const double x, const double y, const double t, const double A, const double w, const double tm) const
{
  double pi = 3.14159265358979323846;
//  return cos((x*x-100.)/_c)*cos((y*y-100.)/_c)*exp(-t)+(2*exp(-t)*(_c*sin((x*x - 100.)/_c)*cos((y*y - 100.)/_c) + cos((x*x - 100.)/_c)*(2.*(x*x + y*y)*cos((y*y - 100.)/_c) + _c*sin((y*y - 100.)/_c))))/(_c*_c);
  if (t<=2*tm)
  {
    return A*cos(w*pi*t)*exp(-(t-tm)*(t-tm));
  }
  else
  {
    return 0.;
  }
}

// Sauvegarde la solution
void Calcul::SaveSol(SparseVector<double> sol, int n, const std::string nom)
{
	string name_file = _results + "/" + nom + std::to_string(n) + ".vtk";

  assert((sol.size() == 2*_nb_pts) && "The size of the solution vector is not the same than the number of 2 * _vertices !");

	ofstream solution;
	solution.open(name_file, ios::out);
	solution.precision(7);

  solution << "# vtk DataFile Version 3.0 " << endl;
  solution << "2D Unstructured Grid" << endl;
  solution << "ASCII" << endl;
  solution << "DATASET UNSTRUCTURED_GRID" << endl;

  solution << "POINTS " << _nb_pts << " float " << endl;
  for (int i = 0 ; i < _nb_pts ; ++i)
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

  solution << "POINT_DATA " << _nb_pts << endl;
  solution << "SCALARS sol float" << endl;
  solution << "LOOKUP_TABLE default" << endl;
  for (int i = 0 ; i < _nb_pts ; ++i)
  {
    solution << float(sqrt(pow(sol.coeffRef(i),2)+pow(sol.coeffRef(i+_nb_pts),2))) << endl;
  }
  solution << endl;

  solution << "VECTORS u float" << endl;
  for (int i = 0 ; i < _nb_pts ; ++i)
  {
    solution << float(sol.coeffRef(i)) << " " << float(sol.coeffRef(i+_nb_pts)) << " 10" << endl;
  }
  solution << endl;

	solution.close();
}

void Calcul::SaveSensor(SparseVector<double> sol, double t)
{
  if ((_s1)or(_s2))
  {
    assert((sol.size() == 2*_nb_pts) && "The size of the solution vector is not the same than the number of 2 * _vertices !");
    ofstream sensor;
    sensor.open(_results + "/sensor.dat", ios::app);
    sensor.precision(7);

    Matrix<double, Dynamic, Dynamic> val;
    double x, y;
    val.resize(_sensor.rows(),_sensor.cols());
    val.setZero();

    for (int i = 0 ; i < _nb_pts ; ++i)
    {
      x = ((_vertices[i]).GetVertice())[0];
      y = ((_vertices[i]).GetVertice())[1];
      for (int j = 0 ; j < val.rows(); j++)
      {
        if ((x<_sensor(j,0)+_eps)and(x>_sensor(j,0)-_eps)and(y<_sensor(j,1)+_eps)and(y>_sensor(j,1)-_eps))
        {
          val(j,0) += sol.coeffRef(i);
          val(j,1) += sol.coeffRef(i+_nb_pts);
        }
      }
    }
    sensor << t;
    for (int j = 0 ; j < val.rows(); j++)
    {
      sensor << " " << val(j,0) << " " << val(j,1);
    }
    sensor << endl;
    sensor.close();
  }
}
#define _MODULE_CPP
#endif
