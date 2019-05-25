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
  maillage->ReadMesh("gmsh/plaquenulle.mesh");
  Calcul* mod = new Calcul(maillage);
  TimeScheme* time_scheme = NULL;

  // Démarrage du chrono
  auto start = chrono::high_resolution_clock::now();

  time_scheme = new EulerScheme();

  double dt = 0.01;
  double nb_iterations = 1000;

  
  int compt=0;
  double tsom,t_adv,t_save;
  
  auto start_init = chrono::high_resolution_clock::now();
  time_scheme->Initialize(dt,mod);
  auto finish_init = chrono::high_resolution_clock::now();
  double t_init = chrono::duration_cast<chrono::milliseconds>(finish_init-start_init).count();
  cout << "\n" << t_init/1000. << " pour l'initialisation" << endl;

      
  cout << "-------------------------------------------------" << endl;
  cout << "Calcul des solutions dans le temps : " << endl;

  for (int n = 2; n <= nb_iterations; n++) // Boucle en temps
  {
    cout.flush();
    cout << "Progression : " << (double)n/(double)nb_iterations*100 << "%       \r";

    auto start_adv = chrono::high_resolution_clock::now();
    time_scheme->Advance();
    auto finish_adv = chrono::high_resolution_clock::now();
    t_adv = chrono::duration_cast<chrono::milliseconds>(finish_adv-start_adv).count();
    //cout << "\n" << t_adv/1000. << " pour advance" << endl; 

    if (n%10 == 0)
    {
      if(compt<10){
	auto start_save = chrono::high_resolution_clock::now();
	time_scheme->SaveSolution(n);
	auto finish_save = chrono::high_resolution_clock::now();
        t_save = chrono::duration_cast<chrono::milliseconds>(finish_save-start_save).count();
	//cout << "\n" << t_save/1000. << " pour l'écriture des fichiers" << endl;
	tsom += t_save;
      }
      else{
		time_scheme->SaveSolution(n);
      }
      compt+=1;
    }
  }

  cout << "\n-------------------------------------------------" << endl;

  // Fin du chrono
  auto finish = chrono::high_resolution_clock::now();
  double t = chrono::duration_cast<chrono::milliseconds>(finish-start).count();

  // Affichage du résultat
  cout << "\n  " << t_init/1000 << "pour l'init";
  cout << "\n  " << tsom/10000 << "pour l'ecriture des fichiers " << endl;
  cout << "\n  " << t_adv/1000 << "pour un pas de temps";
  double temps(t/1000.),h(0),min(0),sec(0);
  cout << "\nCela a pris ";
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
