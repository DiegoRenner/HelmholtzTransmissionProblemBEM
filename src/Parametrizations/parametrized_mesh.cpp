/**
 * \file parametrized_mesh.cpp
 * \brief This file defines the class for representing a mesh comprised
 *        of panels represented by parametrized curves.
 *
 *  This File is a part of the 2D-Parametric BEM package
 */

#include "parametrized_mesh.hpp"

#include <limits>

#include <Eigen/Dense>

namespace parametricbem2d {

ParametrizedMesh::ParametrizedMesh(PanelVector panels) : panels_(panels) {
  // std::cout << "ParametrizedMesh constructor called!" << std::endl;
  unsigned int N = getNumPanels();
  // Determining the split value by looping over the panels in the mesh. Non
  // zero split indicates two distinct boundaries in the mesh object. Used when
  // dealing with an annular domain
  for (unsigned int i = 0; i < N; ++i) {
    if ((panels[0]->operator()(-1.) - panels[i]->operator()(1.)).norm() <
        10*std::numeric_limits<double>::epsilon()) {
      // std::cout << "Break in continuity at position "<< i << std::endl;
      split_ = (i + 1) % N;
    }
  }
  // std::cout << "split : " << split_ << std::endl;
}

PanelVector ParametrizedMesh::getPanels() const {
  // Returning a copy of the stored panels
  return panels_;
}

unsigned int ParametrizedMesh::getNumPanels() const {
  // Returning the size of the PanelVector panels_
  return panels_.size();
}

Eigen::Vector2d ParametrizedMesh::getVertex(unsigned int i) const {
  assert(i < getNumPanels()); // Asserting requested index is within limits
  return panels_[i]->operator()(-1);
}

} // namespace parametricbem2d
