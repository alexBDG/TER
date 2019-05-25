#ifndef _MAIN_CC

#include <iostream>
#include <fstream>
#include <chrono>
#include "TimeScheme.h"

using namespace std;
using namespace Eigen;

int main(int argc, char** argv)
{

  // ---------------------------- Résolution  ----------------------------------
  Mesh2D* maillage = new Mesh2D();
  maillage->ReadMesh("gmsh/plaquemifine.mesh");
  Calcul* mod = new Calcul(maillage);
  TimeScheme* time_scheme = NULL;

  // Démarrage du chrono
  auto start = chrono::high_resolution_clock::now();

  time_scheme = new EulerScheme();

  double dt = 0.001*5;///6320;
  double nb_iterations = 10000;

  time_scheme->Initialize(dt,mod);

  cout << "-------------------------------------------------" << endl;
  cout << "Calcul des solutions dans le temps : " << endl;

  for (int n = 2; n <= nb_iterations; n++) // Boucle en temps
  {
    cout.flush();
    cout << "Progression : " << (double)n/(double)nb_iterations*100 << "% \r";
    time_scheme->Advance();
    if (n%10 == 0)
    {
      time_scheme->SaveSolution(n);
    }
  }

  cout << "\n-------------------------------------------------" << endl;

  // Fin du chrono
  auto finish = chrono::high_resolution_clock::now();
  double t = chrono::duration_cast<chrono::milliseconds>(finish-start).count();
  // Affichage du résultat
  cout << "Cela a pris "<< t/1000. << " seconds" << endl;

  delete time_scheme;
  delete mod;
  delete maillage;

  return 0;
}

#define _MAIN_CC
#endif
