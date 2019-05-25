#ifndef _TIME_SCHEME_CPP

#include "TimeScheme.h"
#include <iostream>

using namespace Eigen;
using namespace std;

// Constructeur par défaut (ne pas oublier de mettre votre pointeur à 0 !!)
TimeScheme::TimeScheme() : _mod(0)
{}

// Destructeur (car on a des fonctions virtuelles)
TimeScheme::~TimeScheme()
{}

// Initialisation de vos différentes variables
void TimeScheme::Initialize(double dt, Calcul* mod)
{
  _dt = dt;
  _t = 0.;
  _mod = mod;
  _sol0 = mod->InitialCondition(_t);
  _mod->SaveSol(_sol0, 0, "solution_");
//  _mod->SaveSol(_mod->SolutionExacte(_t), 0, "FonctionExacte_");
  _t += dt;
  _sol1 = mod->InitialCondition(_t);
  _mod->SaveSol(_sol1, 1, "solution_");
//  _mod->SaveSol(_mod->SolutionExacte(_t), 1, "FonctionExacte_");
  _sol = _sol0;
  mod->BuildMK(dt);
}

void TimeScheme::SaveSolution(int n)
{
  _mod->SaveSol(_sol, n, "solution_");
//  _mod->SaveSol(_mod->SolutionExacte(_t), n, "FonctionExacte_");
}

void TimeScheme::SaveSensor(int n)
{
  _mod->SaveSensor(_sol, n*_dt);
}

// Renvoie _sol (pratique pour vérifier la résolution)
const SparseVector<double> & TimeScheme::GetIterateSolution() const
{
  return _sol;
}

// Euler Explicite
void EulerScheme::Advance()
{
  _t += _dt;
  _mod->BuildF(_t, _dt, _sol, _sol1);
  _sol1 = _sol;
  _sol = _mod->GetF();
}

#define _TIME_SCHEME_CPP
#endif
