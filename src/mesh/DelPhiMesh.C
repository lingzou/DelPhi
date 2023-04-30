#include "DelPhiMesh.h"

registerMooseObject("delphiApp", DelPhiMesh);

InputParameters
DelPhiMesh::validParams()
{
  InputParameters params = MooseMesh::validParams();
  params.set<MooseEnum>("dim") = "3";
  return params;
}

DelPhiMesh::DelPhiMesh(const InputParameters & parameters)
  : MooseMesh(parameters)
{
}

std::unique_ptr<MooseMesh>
DelPhiMesh::safeClone() const
{
  return libmesh_make_unique<DelPhiMesh>(*this);
}

void
DelPhiMesh::build1DMesh(unsigned int subdomain_id,
                     unsigned int nodeset_id,
                     unsigned int bc_id_in,
                     unsigned int bc_id_out,
                     unsigned int n_elems,
                     Real x_min,
                     Real x_max)
{
  std::vector<unsigned int> node_ids;
  std::vector<unsigned int> elem_ids;

  setSubdomainName(subdomain_id, Moose::stringify(subdomain_id));

  // points
  Real delta_t = (x_max - x_min) / n_elems;
  Point p(0, 0, 0); // always starts from origin
  for (unsigned int i = 0; i <= n_elems; i++)
  {
    const Node * nd = getMesh().add_point(p);
    getMesh().boundary_info->add_node(nd, nodeset_id);
    node_ids.push_back(nd->id());
    p(0) += delta_t;
  }

  // elems
  for (unsigned int i = 0; i < n_elems; i++)
  {
    Elem * elem = getMesh().add_elem(new Edge2); // currently only support Edge2 elements
    elem->set_id(i);
    elem->subdomain_id() = subdomain_id;
    elem->set_node(0) = getMesh().node_ptr(node_ids[i]);
    elem->set_node(1) = getMesh().node_ptr(node_ids[i + 1]);
    elem_ids.push_back(elem->id());

    // BCs
    if (i == 0)
    {
      getMesh().boundary_info->add_side(elem, 0, bc_id_in);
    }
    if (i == (n_elems - 1))
    {
      getMesh().boundary_info->add_side(elem, 1, bc_id_out);
    }
  }
}
