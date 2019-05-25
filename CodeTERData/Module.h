#ifndef _MODULE_H

#include "Mesh2D.h"
#include "DataFile.h"
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
    const std::vector<std::vector<double>>& _edges_bc;
		const Eigen::Matrix<double, Eigen::Dynamic, 2>& _tri_center;
  	const Eigen::VectorXd & _tri_area;
    const std::vector<std::vector<double>> _coeff_bc;
    const std::vector<std::vector<double>> _exc;

    Eigen::Matrix<double, Eigen::Dynamic, 3> _Sommets;
    Eigen::Matrix<double, Eigen::Dynamic, 3> _X;
    Eigen::Matrix<double, Eigen::Dynamic, 3> _Y;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> _sensor;

    std::string _initial_condition;
    std::string _mesh_name;
    std::string _file_name;
    std::string _results;
    bool _s1, _s2, _if_exc_u, _if_exc_l, _if_exc;
    double _lambda, _mu, _rho, _c;
    double _x0, _y0, _r;
    double s1_x, s1_y, s2_x, s2_y, _eps;
    double _nb_pts;


    Eigen::SparseVector<double> _f;

    Eigen::SparseMatrix<double, Eigen::RowMajor> _M, _K;

    Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>, Eigen::COLAMDOrdering<int> > _solver;

  public:
    Calcul(Mesh2D* maillage, DataFile* data_file);
    void BuildMK(double dt);

    void BoundaryConditionA(Eigen::SparseMatrix<double, Eigen::ColMajor> A);

    void BoundaryConditionB(Eigen::SparseVector<double>& B, const double t);

    double excitation(Eigen::SparseVector<double>& sol, double const t, int const k);

    Eigen::SparseVector<double> InitialCondition(const double t);
    double InitialFunction_x(const double x, const double y, double t) const;
    double InitialFunction_y(const double x, const double y, double t) const;
    const Eigen::SparseVector<double> & GetF() const {return _f;};

    void BuildF(double t,double dt,Eigen::SparseVector<double> & sol,Eigen::SparseVector<double> & sol1);

    Eigen::SparseVector<double> SecondMember(const double t);
    double f(const double x, const double y, const double t, const double A, const double w, const double tm) const;
    double f(const double x, const double y, const double t, const double A) const {return A;}
    Eigen::SparseVector<double> SolutionExacte(const double t);
    double fexacte(const double x, const double y, const double t) const;

    void SaveSol(Eigen::SparseVector<double> sol, int n, const std::string nom);

    void SaveSensor(Eigen::SparseVector<double> sol, double t);

};

#define _MODULE_H
#endif
