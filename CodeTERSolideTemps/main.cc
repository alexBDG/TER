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
  maillage->ReadMesh("gmsh/plaque.mesh");
  Calcul* mod = new Calcul(maillage);
  TimeScheme* time_scheme = NULL;

  // Démarrage du chrono
  auto start = chrono::high_resolution_clock::now();

  time_scheme = new EulerScheme();

  double dt = 0.001;//0.00000001;
  double nb_iterations = 200;//30000;

  time_scheme->Initialize(dt,mod);

  cout << "-------------------------------------------------" << endl;
  cout << "Calcul des solutions dans le temps : " << endl;

  for (int n = 2; n <= nb_iterations; n++) // Boucle en temps
  {
    cout.flush();
    cout << "Progression : " << (double)n/(double)nb_iterations*100 << "%         \r";
    time_scheme->Advance();
//    time_scheme->SaveSensor(n);
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
  double temps(t/1000.),h(0),min(0),sec(0);
  cout << "Cela a pris ";
  while (temps - 3600 >= 0)
  {
    h++;
    temps -= 3600;
  }
  if (h > 1) cout << h << " heures ";
  if (h == 1) cout << h << " heure ";
  while (temps - 60 >= 0)
  {
    min++;
    temps -= 60;
  }
  if (min > 0) cout << min << " minutes et ";
  if (min == 1) cout << h << " minute et ";
  if (temps > 1)
  {
    cout << floor(temps) << " seconds" << endl;
  }
  else
  {
    cout << floor(temps) << " second" << endl;
  }

  delete time_scheme;
  delete mod;
  delete maillage;

  return 0;
}

#define _MAIN_CC
#endif
