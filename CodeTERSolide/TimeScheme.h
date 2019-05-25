#ifndef _TIME_SCHEME_H

#include "Module.h"
#include "Sparse"
#include "Dense"
class TimeScheme
{
  protected:
    // Pas de temps
    double _dt;
    // Temps en cours
    double _t;
    // Vecteur initial et vecteur solution
    Eigen::SparseVector<double> _sol0, _sol1, _sol;
    // Pointeur vers la classe Advection
    Calcul* _mod;

  public:
    // Constructeur par défaut
    TimeScheme();
    // Destructeur par défaut - Si la classe ne contient pas de destructeur par défaut
    // alors le compilateur en génère un implicitement.
    virtual ~TimeScheme();
    // Initialisation de vos différentes variables
    void Initialize(double dt, Calcul* mod);
    // Enregistre la solution un fichier
    void SaveSolution(int n);
    // Une étape du schéma en temps
    virtual void Advance() = 0;
    // Permet de récupérer _sol
    const Eigen::SparseVector<double> & GetIterateSolution() const;
};

class EulerScheme : public TimeScheme
{
  public:
    // Une étape du schéma en temps
    void Advance();
};

#define _TIME_SCHEME_H
#endif
