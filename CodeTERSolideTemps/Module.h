#ifndef _MODULE_H

#include "Mesh2D.h"
#include "Sparse"
#include "Dense"
#include <SparseLU>
#include <fstream>

class Calcul
{
  private:
    Mesh2D* _maillage;
		const std::vector<Vertice>& _vertices;
		const std::vector<Triangle>& _triangles;
    const std::vector<Edge>& _edges;
		const Eigen::Matrix<double, Eigen::Dynamic, 2>& _tri_center;
  	const Eigen::VectorXd & _tri_area;

    Eigen::Matrix<double, Eigen::Dynamic, 3> _Sommets;
    Eigen::Matrix<double, Eigen::Dynamic, 3> _X;
    Eigen::Matrix<double, Eigen::Dynamic, 3> _Y;

    double nb_pts;
    double _lambda, _mu, _rho, _vL, _vT;
    double _c, _a, _C11, _C22, _C66, _C12;

    Eigen::SparseVector<double> _f;

    Eigen::SparseMatrix<double, Eigen::RowMajor> _M, _K;

    Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>, Eigen::COLAMDOrdering<int> > _solver;

  public:
    Calcul(Mesh2D* maillage);
    void BuildMK(double dt);

    void BoundaryConditionA(Eigen::SparseMatrix<double, Eigen::ColMajor> A);

    void BoundaryConditionB(Eigen::SparseVector<double>& B, const double t);

    Eigen::SparseVector<double> InitialCondition(const double t);
    double FunctionInitial(const double x, const double y, double t) const;
    const Eigen::SparseVector<double> & GetF() const {return _f;};

    void BuildF(double t,double dt,Eigen::SparseVector<double> & sol,Eigen::SparseVector<double> & sol1);

    Eigen::SparseVector<double> SecondMember(const double t);
    double f(const double x, const double y, const double t) const;
    double f1(const double x, const double y, const double t) const;
    double f2(const double x, const double y, const double t) const;
    Eigen::SparseVector<double> SolutionExacte(const double t);
    double fexacte1(const double x, const double y, const double t) const;
    double fexacte2(const double x, const double y, const double t) const;

    void SaveSol(Eigen::SparseVector<double> sol, int n, const std::string nom);
    void SaveNorm(Eigen::SparseVector<double> sol, Eigen::SparseVector<double> sol_e, double t);

    void SaveSensor(Eigen::SparseVector<double> sol, double t);

};

#define _MODULE_H
#endif
