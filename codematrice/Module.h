#ifndef _MODULE_H

#include "Mesh2D.h"
#include "Sparse"
#include "Dense"
#include <SparseCholesky>
#include <fstream>

class Calcul
{
  private:
    Mesh2D* _maillage;
		const std::vector<Vertice>& _vertices;
		const std::vector<Triangle>& _triangles;
		const Eigen::Matrix<double, Eigen::Dynamic, 2>& _tri_center;
  	const Eigen::VectorXd & _tri_area;

    Eigen::Matrix<double, Eigen::Dynamic, 3> _Sommets;
    Eigen::Matrix<double, Eigen::Dynamic, 3> _X;
    Eigen::Matrix<double, Eigen::Dynamic, 3> _Y;

    double _c, _v;

    Eigen::SparseVector<double> _f;

    Eigen::SparseMatrix<double> _M, _K, _M1, _K1;

    Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > _solver;

  public:
    Calcul(Mesh2D* maillage);
    void BuildMK(double dt);

    double Phi(const double x, const double y, int i, int som);
    double DxPhi(const double x, const double y, int i, int som);
    double DyPhi(const double x, const double y, int i, int som);

    Eigen::SparseVector<double> InitialCondition(const double t);
    double FunctionInitial(const double x, const double y, double t) const;
    const Eigen::SparseVector<double> & GetF() const {return _f;};

    void BuildF(double t,double dt,Eigen::SparseVector<double> & sol,Eigen::SparseVector<double> & sol1);

    Eigen::SparseVector<double> SecondMember(const double t);
    double f(const double x, const double y, const double t) const;
    Eigen::SparseVector<double> SolutionExacte(const double t);
    double fexacte(const double x, const double y, const double t) const;

    void SaveSol(Eigen::SparseVector<double> sol, int n, const std::string nom);
    void SaveNorm(Eigen::SparseVector<double> sol, Eigen::SparseVector<double> sol_e, double t);

};

#define _MODULE_H
#endif
