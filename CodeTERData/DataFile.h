#ifndef _DATA_FILE_H

#include <string>
#include <vector>
#include <iostream>
// Définition de la classe

class DataFile
{
private:
  double _t0, _tfinal, _dt;
  double _x0, _y0, _r, _w, _rho, _lambda, _mu, _c;
  double s1_x, s1_y, s2_x, s2_y, _eps;
  int _it_saved;

  // Premiere element de chaque vecteur  : contient un nombre associé aux conditions aux limites du bords en question
  // Deuxieme element de chaque vecteur  : Amplitude ou constante A
  // Troisieme element de chaque vecteur : pulsation w
  // quatrieme element de chaque vecteur : temps pour lequel l'amplitude du signal est maximale
  // cinquieme ou troisieme element      : nombre correspondant à l'orientation du déplacement
  std::vector<std::vector<double>> _coeff;
  std::vector<std::vector<double>> _exc;

  std::string _file_name, _mesh_name, _results;
  std::string _initial_condition_choice, _sigma;
  std::string _bc_left, _bc_right, _bc_up, _bc_down;
  std::string _sensor1, _sensor2, _excitation_l, _excitation_u;

  bool _if_excitation_l, _if_excitation_u;
  bool _if_mesh_name;
  bool _if_t0, _if_tfinal, _if_dt;
  bool _if_rho, _if_lambda, _if_mu, _if_c, _if_r;
  bool _if_initial_condition_choice;
  bool _if_bc_up, _if_bc_down, _if_bc_left, _if_bc_right;
  bool _if_sensor1, _if_sensor2, _if_eps;
  bool _if_results, _if_iteration_saved;

public: // Méthodes et opérateurs de la classe
  DataFile(std::string file_name);
  void ReadDataFile();
  void Adapt_dt(double dt){_dt = dt;};
  double Get_t0() const {return _t0;};
  double Get_tfinal() const {return _tfinal;};
  double Get_dt() const {return _dt;};
  double Get_param_x0() const { return _x0;};
  double Get_param_y0() const { return _y0;};
  double Get_param_r() const { return _r;};
  double Get_param_c() const { return _c;};
  double Get_param_rho() const { return _rho;};
  double Get_param_lambda() const { return _lambda;};
  double Get_param_mu() const { return _mu;};
  double Get_param_s1_x() const {return s1_x;};
  double Get_param_s1_y() const {return s1_y;};
  double Get_param_s2_x() const {return s2_x;};
  double Get_param_s2_y() const {return s2_y;};
  double Get_param_eps() const {return _eps;};
  int Get_its() const {return _it_saved;};
  bool Get_sensor1() const {return _if_sensor1;};
  bool Get_sensor2() const {return _if_sensor2;};
  bool Get_bool_excitation_l() const {return _if_excitation_l;};
  bool Get_bool_excitation_u() const {return _if_excitation_u;};
  std::string Get_bc(int edge) const;
  std::string Get_excitation_l_type() const {return _excitation_l;};
  std::string Get_excitation_u_type() const {return _excitation_u;};
  std::string Get_mesh_name() const {return _mesh_name;};
  std::string Get_initial_condition_choice() const {return _initial_condition_choice;};
  std::string Get_file_name() const {return _file_name;};
  std::string Get_results() const {return _results;};
  std::vector<std::vector<double>> Get_param() const {return _coeff;};
  std::vector<std::vector<double>> Get_param_excitation() const {return _exc;};
};

#define _DATA_FILE_H
#endif
