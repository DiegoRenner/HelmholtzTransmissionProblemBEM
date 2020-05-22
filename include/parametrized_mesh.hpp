/**
 * \file parametrized_mesh.hpp
 * \brief This file declares a class for representing a mesh comprising
 *        of panels in the form of parametrized curves.
 *
 *  This File is a part of the 2D-Parametric BEM package
 */

#ifndef PARAMETRIZEDMESHHPP
#define PARAMETRIZEDMESHHPP

#include "abstract_parametrized_curve.hpp"

/**
 * \class ParametrizedMesh
 * \brief This class represents a mesh which is comprised of panels in the
 *        form of parametrized curves. It is used in assembly of Galerkin
 *        matrices using the parametric BEM approach. It stores the panels
 *        using PanelVector and enforces additional constraints which
 *        require the end point of a panel to be the starting point of the
 *        next, such that the panels form a curved polygon.
 */
class ParametrizedMesh {
public:
  /**
   * Constructor using a PanelVector object which contains the component
   * panels of the mesh in the form of parametrized curves.
   */
  ParametrizedMesh(PanelVector panels);

  /**
   * This function is used for retrieving a PanelVector containing all the
   * parametrized curve panels in the parametrized mesh.
   *
   * @return A PanelVector containing all the parametrized panels in the mesh
   */
  PanelVector getPanels() const;

  /**
   * This function is used for getting the number of panels in the
   * parametrized mesh
   *
   * @return number of panels in the mesh
   */
  unsigned int getNumPanels() const;

  /**
   * This function is used for getting the ith vertex in the parametrized mesh
   *
   * @return ith vertex in the mesh as Eigen::Vector2d
   */
  Eigen::Vector2d getVertex(unsigned int i) const;

  /**
   * This function is used for getting the split value for the mesh. If split
   * is non zero, it indicates the position where the second boundary begins
   * in the mesh object; the domain is annular. A zero value indicates there is
   * only one boundary in the mesh object
   *
   * @return The position in the mesh where the second boundary begins
   */
  unsigned getSplit() const { return split_; }

  //void addPanels(const PanelVector& panels) {
  //  panels_.insert()
  //}

private:
  /**
   * Private const field for the PanelVector of the mesh
   */
  const PanelVector panels_;
  /**
   * Private unsigned field used to distinguish one boundary from another in the
   * mesh (annular domain). Indicates the starting position of second boundary.
   */
  unsigned split_;
}; // class ParametrizedMesh

#endif // PARAMETRIZEDMESHHPP
