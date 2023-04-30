#pragma once

#include "MooseMesh.h"
#include "libmesh/edge_edge2.h"

class DelPhiMesh : public MooseMesh
{
public:
  static InputParameters validParams();

  DelPhiMesh(const InputParameters & parameters);
  virtual ~DelPhiMesh() {}

  virtual std::unique_ptr<MooseMesh> safeClone() const override;
  virtual void buildMesh() override {}

  void build1DMesh(unsigned int subdomain_id,
                   unsigned int nodeset_id,
                   unsigned int bc_id_in,
                   unsigned int bc_id_out,
                   unsigned int n_elems,
                   Real x_min,
                   Real x_max);
};
