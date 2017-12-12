anisotropy_code = """

class Conductivity : public Expression
{
public:

  // Create expression with 6 components
  Conductivity() : Expression(6) {}

  // Function for evaluating expression on each cell
  void eval(Array<double>& values, const Array<double>& x, const ufc::cell& cell) const
  {
    const uint D = cell.topological_dimension;
    const uint cell_index = cell.index;
    values[0] = (*c00)[cell_index];
    values[1] = (*c01)[cell_index];
    values[2] = (*c02)[cell_index];
    values[3] = (*c11)[cell_index];
    values[4] = (*c12)[cell_index];
    values[5] = (*c22)[cell_index];
  }

  // The data stored in mesh functions
  std::shared_ptr<MeshFunction<double> > c00;
  std::shared_ptr<MeshFunction<double> > c01;
  std::shared_ptr<MeshFunction<double> > c02;
  std::shared_ptr<MeshFunction<double> > c11;
  std::shared_ptr<MeshFunction<double> > c12;
  std::shared_ptr<MeshFunction<double> > c22;

};
"""


homogeneous_code = """

class Conductivity : public Expression
{
public:

  // Create expression with 1 components
  Conductivity() : Expression(1) {}

  // Function for evaluating expression on each cell
  void eval(Array<double>& values, const Array<double>& x, const ufc::cell& cell) const
  {
    const uint D = cell.topological_dimension;
    const uint cell_index = cell.index;
    values[0] = (*c00)[cell_index];
  }

  // The data stored in mesh functions
  std::shared_ptr<MeshFunction<double> > c00;
};
"""
