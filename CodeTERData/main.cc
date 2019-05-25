#ifndef _MAIN_CC

#include <iostream>
#include <fstream>
#include <chrono>
#include "TimeScheme.h"
#include "DataFile.h"

using namespace std;
using namespace Eigen;

int main(int argc, char** argv)
{

  if (argc < 2)
  {
    cout << "Please, enter the name of your data file." << endl;
    abort();
  }
  const string data_file_name = argv[1];

  // ----------------------- Fichier de données --------------------------------
  DataFile* data_file = new DataFile(data_file_name);
  data_file->ReadDataFile();
  // ---------------------------------------------------------------------------

  // ------------------Définition du nombre d'itérations------------------------
  int nb_iterations = int(ceil((data_file->Get_tfinal()-data_file->Get_t0())/data_file->Get_dt()));
  data_file->Adapt_dt(data_file->Get_tfinal() / nb_iterations);
  // ---------------------------------------------------------------------------

  // -----------Définition de la fréquance d'ecriture des solutions-------------
  int it_saved = data_file->Get_its();
  // ---------------------------------------------------------------------------

  // --------------------Définition de la plaque a utiliser---------------------
  string mesh_name = data_file->Get_mesh_name();
  // ---------------------------------------------------------------------------

  // ------------------------------ Résolution----------------------------------
  Mesh2D* maillage = new Mesh2D();
  maillage->ReadMesh(mesh_name);
  Calcul* mod = new Calcul(maillage, data_file);
  TimeScheme* time_scheme = NULL;

  cout << "-------------------------------------------------" << endl;
  cout << "Search u such that : " << endl;
  cout << "             rho*dtt u - C*dxx u = f             " << endl;
  cout << "-------------------------------------------------" << endl;

  // Démarrage du chrono
  auto start = chrono::high_resolution_clock::now();

  time_scheme = new EulerScheme();
  time_scheme->Initialize(data_file->Get_dt(),mod);

  cout << "-------------------------------------------------" << endl;
  cout << "Solutions calculated" << endl;

  for (int n = 2; n <= nb_iterations; n++) // Boucle en temps
  {
    cout.flush();
    cout << "Progression : " << (double)n/(double)nb_iterations*100 << "%      \r";
    time_scheme->Advance();
    time_scheme->SaveSensor(n);
    if (n%it_saved == 0)
    {
      time_scheme->SaveSolution(n);
    }
  }

  cout << "-------------------------------------------------" << endl;

  // Fin du chrono
  auto finish = chrono::high_resolution_clock::now();
  double t = chrono::duration_cast<chrono::milliseconds>(finish-start).count();

  // Affichage du temps de calcul du résultat
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
