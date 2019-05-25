#ifndef _DATA_FILE_CPP

#include "DataFile.h"
#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;

DataFile::DataFile(std::string file_name)
: _file_name(file_name),
_if_results(false), _if_mesh_name(false),
_if_t0(false), _if_tfinal(false), _if_dt(false),
_if_rho(false), _if_lambda(false), _if_mu(false),
_if_c(false), _if_initial_condition_choice(false),
_if_bc_up(false), _if_bc_down(false), _if_bc_left(false),
_if_bc_right(false), _if_sensor1(false), _if_sensor2(false)
{}


void DataFile::ReadDataFile()
{
  ifstream data_file(_file_name.data());
  if (!data_file.is_open())
  {
    cout << "Unable to open file " << _file_name << endl;
    abort();
  }
  else
  {
    cout << "-------------------------------------------------" << endl;
    cout << "Reading data file " << _file_name << endl;
  }

  string file_line;
  _coeff.resize(4);
  _exc.resize(2);
  while (!data_file.eof())
  {
    getline(data_file, file_line);
    if (file_line.find("mesh") != std::string::npos)
    {
      data_file >> _mesh_name; _if_mesh_name = true;
    }

    if (file_line.find("t0") != std::string::npos)
    {
      data_file >> _t0; _if_t0 = true;
    }

    if (file_line.find("tfinal") != std::string::npos)
    {
      data_file >> _tfinal; _if_tfinal = true;
    }

    if (file_line.find("dt") != std::string::npos)
    {
      data_file >> _dt; _if_dt = true;
    }

    if (file_line.find("rho") != std::string::npos)
    {
      data_file >> _rho; _if_rho = true;
    }

    if (file_line.find("lambda") != std::string::npos)
    {
      data_file >> _lambda; _if_lambda = true;
    }

    if (file_line.find("mu") != std::string::npos)
    {
      data_file >> _mu; _if_mu = true;
    }

    if (file_line.find("celerity") != std::string::npos)
    {
      data_file >> _c; _if_c = true;
    }

    if (file_line.find("iteration_saved") != std::string::npos)
    {
      data_file >> _it_saved; _if_iteration_saved = true;
    }

    if (file_line.find("initial_condition") != std::string::npos)
    {
      data_file >> _initial_condition_choice; _if_initial_condition_choice = true;
      if (_initial_condition_choice == "center")
      {
        data_file >> _x0 >> _y0 >> _r;
      }
      else if (_initial_condition_choice == "rectangular")
      {
        data_file >> _y0 >> _r;
        _x0 = 0;
      }
      else if (_initial_condition_choice == "no")
      {
        _x0 = 0;
        _y0 = 0;
        _r = 0;
      }
      else
      {
        _if_initial_condition_choice = false;
      }
    }

    if (file_line.find("boundary_condition_down") != std::string::npos)
    {
      data_file >> _bc_down; _if_bc_down = true;
      if (_bc_down == "dirichlet_h")
      {
        _coeff[0].resize(3);
        _coeff[0][0] = 11;
        data_file >> _coeff[0][1] >> _sigma;
        if (_sigma == "c")
        {
          _coeff[0][2] = -1.; //cisaillement
        }
        else
        {
          _coeff[0][2] = 1.; //compression traction
        }
      }
      else if (_bc_down == "dirichlet_nh")
      {
        _coeff[0].resize(5);
        _coeff[0][0] = 31;
        data_file >> _coeff[0][1] >> _coeff[0][2] >> _coeff[0][3] >> _sigma;
        if (_sigma == "c")
        {
          _coeff[0][4] = -1.;
        }
        else
        {
          _coeff[0][4] = 1.;
        }
      }
      else if (_bc_down == "neumann")
      {
        _coeff[0].resize(1);
        _coeff[0][0] = 0;
      }
      else
      {
        _if_bc_down = false;
      }
    }

    if (file_line.find("boundary_condition_right") != std::string::npos)
    {
      data_file >> _bc_right; _if_bc_right = true;
      if (_bc_right == "dirichlet_h")
      {
        _coeff[1].resize(3);
        _coeff[1][0] = 11;
        data_file >> _coeff[1][1] >> _sigma;
        if (_sigma == "c")
        {
          _coeff[1][2] = -1.;
        }
        else
        {
          _coeff[1][2] = 1.;
        }
      }
      else if (_bc_right == "dirichlet_nh")
      {
        _coeff[1].resize(5);
        _coeff[1][0] = 31;
        data_file >> _coeff[1][1] >> _coeff[1][2] >> _coeff[1][3] >> _sigma;
        if (_sigma == "c")
        {
          _coeff[1][4] = -1.;
        }
        else
        {
          _coeff[1][4] = 1.;
        }
      }
      else if (_bc_right == "neumann")
      {
        _coeff[1].resize(1);
        _coeff[1][0] = 0;
      }
      else
      {
        _if_bc_right = false;
      }
    }

    if (file_line.find("boundary_condition_up") != std::string::npos)
    {
      data_file >> _bc_up; _if_bc_up = true;
      if (_bc_up == "dirichlet_h")
      {
        _coeff[2].resize(3);
        _coeff[2][0] = 11;
        data_file >> _coeff[2][1] >> _sigma;
        if (_sigma == "c")
        {
          _coeff[2][2] = -1.;
        }
        else
        {
          _coeff[2][2] = 1.;
        }
      }
      else if (_bc_up == "dirichlet_nh")
      {
        _coeff[2].resize(5);
        _coeff[2][0] = 31;
        data_file >> _coeff[2][1] >> _coeff[2][2] >> _coeff[2][3] >> _sigma;
        if (_sigma == "c")
        {
          _coeff[2][4] = -1.;
        }
        else
        {
          _coeff[2][4] = 1.;
        }
      }
      else if (_bc_up == "neumann")
      {
        _coeff[2].resize(1);
        _coeff[2][0] = 0;
      }
      else
      {
        _if_bc_up = false;
      }
    }

    if (file_line.find("boundary_condition_left") != std::string::npos)
    {
      data_file >> _bc_left; _if_bc_left = true;
      if (_bc_left == "dirichlet_h")
      {
        _coeff[3].resize(3);
        _coeff[3][0] = 11;
        data_file >> _coeff[3][1] >> _sigma;
        if (_sigma == "c")
        {
          _coeff[3][2] = -1.;
        }
        else
        {
          _coeff[3][2] = 1.;
        }
      }
      else if (_bc_left == "dirichlet_nh")
      {
        _coeff[3].resize(5);
        _coeff[3][0] = 31;
        data_file >> _coeff[3][1] >> _coeff[3][2] >> _coeff[3][3] >> _sigma;
        if (_sigma == "c")
        {
          _coeff[3][4] = -1.;
        }
        else
        {
          _coeff[3][4] = 1.;
        }
      }
      else if (_bc_left == "neumann")
      {
        _coeff[3].resize(1);
        _coeff[3][0] = 0;
      }
      else
      {
        _if_bc_left = false;
      }
    }

    if (file_line.find("excitation_left") != std::string::npos)
    {
      data_file >> _excitation_l;
      if (_excitation_l == "yes")
      {
        _if_excitation_l = true;

        _exc[1].resize(6);
        data_file >> _exc[1][0] >> _exc[1][1] >> _exc[1][2] >> _exc[1][3] >> _exc[1][4] >> _excitation_l;
        if (_excitation_l == "tc")
        {
          _exc[1][5] = 1.;
        }
        else
        {
          _exc[1][5] = -1;
        }
      }
    }

    if (file_line.find("excitation_up") != std::string::npos)
    {
      data_file >> _excitation_u;
      if (_excitation_u == "yes")
      {
        _if_excitation_u = true;

        _exc[0].resize(6);
        data_file >> _exc[0][0] >> _exc[0][1] >> _exc[0][2] >> _exc[0][3] >> _exc[0][4] >> _excitation_u;
        if (_excitation_u == "tc")
        {
          _exc[0][5] = 1.;
        }
        else
        {
          _exc[0][5] = -1;
        }
      }
    }

    if (file_line.find("state_of_sensor1") != std::string::npos)
    {
      data_file >> _sensor1; _if_sensor1 = true;
      if (_sensor1=="on")
      {
        data_file >> s1_x >> s1_y;
      }
      else
      {
        _if_sensor1 = false;
      }
    }

    if (file_line.find("state_of_sensor2") != std::string::npos)
    {
      data_file >> _sensor2; _if_sensor2 = true;
      if (_sensor2=="on")
      {
        data_file >> s2_x >> s2_y;
      }
      else
      {
        _if_sensor2 = false;
      }
    }

    if (file_line.find("precision") != std::string::npos)
    {
      data_file >> _eps; _if_eps = true;
    }

    if (file_line.find("results") != std::string::npos)
    {
      data_file >> _results; _if_results = true;
    }
  }

  if (!_if_mesh_name)
  {
    cout << "-------------------------------------------------" << endl;
    cout << "Do not forget to give the mesh name in the data file !" << endl;
    abort();
  }
  if (!_if_t0)
  {
    cout << "-------------------------------------------------" << endl;
    cout << "Beware - The default value (0.) is used for t0." << endl;
    _t0 = 0.;
  }
  if (!_if_tfinal)
  {
    cout << "-------------------------------------------------" << endl;
    cout << "Beware - The default value (0.1) is used for tfinal." << endl;
    _tfinal = 0.1;
  }
  if (!_if_dt)
  {
    cout << "-------------------------------------------------" << endl;
    cout << "Beware - The default value (0.001) is used for dt." << endl;
    _dt = 0.001;
  }
  if (!_if_rho)
  {
    cout << "-------------------------------------------------" << endl;
    cout << "Beware - The default value (1.) is used for rho." << endl;
    _rho = 1.;
  }
  if (!_if_lambda)
  {
    cout << "-------------------------------------------------" << endl;
    cout << "Beware - The default value (1.) is used for lambda." << endl;
    _lambda = 1.;
  }
  if (!_if_mu)
  {
    cout << "-------------------------------------------------" << endl;
    cout << "Beware - The default value (1.) is used for mu." << endl;
    _mu = 1.;
  }
  if (!_if_c)
  {
    cout << "-------------------------------------------------" << endl;
    cout << "Beware - The default value (100.) is used for the celerity." << endl;
    _c = 1.;
  }
  if (!_if_iteration_saved)
  {
    cout << "-------------------------------------------------" << endl;
    cout << "Beware - The default value (10) is used for the frequency of recording solutions." << endl;
    _it_saved = 10;
  }
  if (!_if_initial_condition_choice)
  {
    cout << "-------------------------------------------------" << endl;
    cout << "Beware - The default value (no initial condition) is used." << endl;
    _initial_condition_choice = "no";
    _x0 = 0;
    _y0 = 0;
    _r = 0;
  }
  if (!_if_bc_left)
  {
    cout << "-------------------------------------------------" << endl;
    cout << "Beware - The default Neumann boundary condition is used for the left boundary condition.";
    _bc_left == "neumann";
    _coeff[3].resize(1);
    _coeff[3][0] = 0.;
  }
  if (!_if_bc_right)
  {
    cout << "-------------------------------------------------" << endl;
    cout << "Beware - The default Neumann boundary condition is used for the right boundary condition.";
    _bc_right == "neumann";
    _coeff[1].resize(1);
    _coeff[1][0] = 0.;
  }
  if (!_if_bc_up)
  {
    cout << "-------------------------------------------------" << endl;
    cout << "Beware - The default Neumann boundary condition is used for the up boundary condition.";
    _bc_up == "neumann";
    _coeff[2].resize(1);
    _coeff[2][0] = 0.;
  }
  if (!_if_bc_down)
  {
    cout << "-------------------------------------------------" << endl;
    cout << "Beware - The default Neumann boundary condition is used for the down boundary condition.";
    _bc_down == "neumann";
    _coeff[0].resize(1);
    _coeff[0][0] = 0.;
  }

  if ((_if_excitation_l)and(_bc_left == "dirichlet_nh"))
  {
    cout << "-------------------------------------------------" << endl;
    cout << "Beware - You have chosen two different boundary conditions on the left border !";
    cout << "-------------------------------------------------" << endl;
    abort();
  }

  if ((_if_excitation_u)and(_bc_up == "dirichlet_nh"))
  {
    cout << "-------------------------------------------------" << endl;
    cout << "Beware - You have chosen two different boundary conditions on the up border !";
    cout << "-------------------------------------------------" << endl;
    abort();
  }

  if ((!_if_sensor1)and(_sensor1 != "off"))
  {
    cout << "-------------------------------------------------" << endl;
    cout << "Beware - The default state of sensor1 (off) is used.";
  }

  if ((!_if_sensor2)and(_sensor2 != "off"))
  {
    cout << "-------------------------------------------------" << endl;
    cout << "Beware - The default state of sensor2 (off) is used.";
  }

  if ((!_if_eps)and((_if_sensor1)or(_if_sensor2)))
  {
    cout << "-------------------------------------------------" << endl;
    cout << "Beware - The default value (0.5) is used for the precision of the sensors." << endl;
    _eps = 0.;
  }

  if (!_if_results)
  {
    cout << "-------------------------------------------------" << endl;
    cout << "Beware - The default results folder name (results) is used." << endl;
    _results = "results";
  }
  cout << "-------------------------------------------------" << endl;
}

// 1 : Down ; 2 : Right ; 3 = Up ; 4 = Left ;
string DataFile::Get_bc(int edge) const
{
  if (edge == 1)
  {
    return _bc_down;
  }
  else if (edge == 2)
  {
    return _bc_right;
  }
  else if (edge == 3)
  {
    return _bc_up;
  }
  else if (edge == 4)
  {
    return _bc_left;
  }
  else
  {
    return "";
  }
}

#define _DATA_FILE_CPP
#endif
