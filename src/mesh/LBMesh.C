//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html


#include "LBMesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/periodic_boundaries.h"
#include "libmesh/periodic_boundary_base.h"
#include "libmesh/unstructured_mesh.h"
#include "libmesh/node.h"


registerMooseObject("MooseApp", LBMesh);

InputParameters
LBMesh::validParams()
{
    InputParameters params = GeneratedMesh::validParams();
    // LBM cannot be used with single element meshes
    params.addRequiredParam<unsigned int>("nx", "Number of nodes in the X direction");
    params.addRequiredParam<unsigned int>("ny", "Number of nodes in the Y direction");
    params.addParam<unsigned int>("nz", 1, "Number of nodes in the Z direction");
    params.addParam<bool>("load_mesh_from_file", false, "Load mesh from file");
    params.addParam<std::string>("mesh_file", "", "Mesh file name");
    params.setParameters<bool>("allow_renumbering", false); // prevents renumbering of nodes

    params.addClassDescription(
      "Create structured mesh for Latice Boltzmann simulations.");

    return params;
}

LBMesh::LBMesh(const InputParameters & parameters)
    : GeneratedMesh(parameters),
    _lbnx(getParam<unsigned int>("nx")),
    _lbny(getParam<unsigned int>("ny")),
    _lbnz(getParam<unsigned int>("nz")),
    _load_mesh_from_file(getParam<bool>("load_mesh_from_file")),
    _mesh_file(getParam<std::string>("mesh_file"))
{   
    bool _gauss_lobatto_grid = getParam<bool>("gauss_lobatto_grid");
    Real _bias_x = getParam<Real>("bias_x");
    Real _bias_y = getParam<Real>("bias_y");
    Real _bias_z = getParam<Real>("bias_z");

    if (_gauss_lobatto_grid || _bias_x != 1.0 || _bias_y != 1.0 || _bias_z != 1.0)
        mooseError("Cannot apply Gauss-Lobatto mesh grading or biasing on Lattice Boltzmann Mesh.");
}

void
LBMesh::buildMesh()
{
  MooseEnum elem_type_enum = getParam<MooseEnum>("elem_type");

  if (!isParamValid("elem_type"))
  {
    switch (_dim)
    {
      case 1:
        mooseError("1D Lattice Boltzmann simulations are not supported.");
        break;
      case 2:
        elem_type_enum = "QUAD4";
        break;
      case 3:
        elem_type_enum = "HEX8";
        break;
    }
  }

  ElemType elem_type = Utility::string_to_enum<ElemType>(elem_type_enum);

  switch (_dim)
  {
    case 2:
      MeshTools::Generation::build_square(dynamic_cast<UnstructuredMesh &>(getMesh()),
                                          _lbnx - 1,
                                          _lbny - 1,
                                          _xmin,
                                          _xmax,
                                          _ymin,
                                          _ymax,
                                          elem_type,
                                          _gauss_lobatto_grid);
      break;
    case 3:
      MeshTools::Generation::build_cube(dynamic_cast<UnstructuredMesh &>(getMesh()),
                                        _lbnx - 1,
                                        _lbny - 1,
                                        _lbnz - 1,
                                        _xmin,
                                        _xmax,
                                        _ymin,
                                        _ymax,
                                        _zmin,
                                        _zmax,
                                        elem_type,
                                        _gauss_lobatto_grid);
      break;
  }
}


torch::Tensor
LBMesh::loadMeshFromFile()
{
  std::ifstream file(_mesh_file);
  if (!file.is_open())
  {
    mooseError("Cannot open file: " + _mesh_file);
  }
  std::vector<std::vector<unsigned int>> matrixData;
  std::string line;
    
  // Read each line from the file
  while (std::getline(file, line)) 
  {
      std::istringstream iss(line);
      unsigned int num;
      std::vector<unsigned int> row;
      
      // Parse each number in the line
      while (iss >> num) 
          row.push_back(num);

      matrixData.push_back(row);
  }
  file.close();
  
  torch::Tensor mesh = torch::ones({_lbnz, _lbny, _lbnx}, torch::kInt16);
  // reshape 
  for (auto i = 0; i < _lbnz; i ++)
  {
    for(auto j = 0; j < _lbny; j ++)
    {
      for (auto k = 0; k < _lbnx; k ++)
      {
        mesh.index_put_({i, j, k}, matrixData[i * _lbny + j][k]);
        // std::cout<<data[i * _lbnx *_lbny + k * _lbny + j]<<" ";
      }
    }
  }
  return mesh;
}

unsigned int
LBMesh::getNx() const
{
    return _lbnx;
}
unsigned int
LBMesh::getNy() const
{
    return _lbny;
}
unsigned int
LBMesh::getNz() const
{
    return _lbnz;
}
MooseEnum
LBMesh::getDim() const
{
    return _dim;
}
bool
LBMesh::load_mesh() const
{
    return _load_mesh_from_file;
}