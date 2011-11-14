
function-factory-list(NAME) is:

  <ParameterList name=NAME>
    function-specification
  </ParameterList>

  * The parameter list name string NAME is arbitrary and meaningful only to the
    parent parameter list.
  
  * This list is given as input to the Amanzi::FunctionFactory::Create
    method which instantiates a new Amanzi::Function object.

function-specification is one of the following parameter lists.
New function types can added easily.

  * Constant function: f(x) = a, for all x

  <ParameterList name="function-constant">
    <Parameter name="value" type="double" value=double />
  </ParameterList>
  
  * Tabular function:  Given values x_i, y_i, i=0, ... n-1, f(x) is defined
    piecewise: f(x) = x_0, x <= x_0; (x, f(x)) is the linear interpolant
    between (x_{i-1},y_{i-1}) and {x_i,y_i) for x in [x_{i-1},x_i], i=1,...,n-1;
    f(x) = x_{n-1}, x > x_{n-1}.

  <ParameterList name="function-tabular">
    <Parameter name="x values" type="Array double" value=double-array />
    <Parameter name="y values" type="Array double" value=double-array />
  </ParameterList>
  
  * Smooth step function:  A C2 function f(x) such that f(x) = y0 for x < x0,
    f(x) = y1 for x > x1, and monotonically increasing for x in [x1, x2].

  <ParameterList name="function-smooth-step">
    <Parameter name="x0" type="double" value=double />
    <Parameter name="y0" type="double" value=double />
    <Parameter name="x1" type="double" value=double />
    <Parameter name="y1" type="double" value=double />
  </ParameterList>

  * Polynomial function: \[ f(x) = \sum_{j=0}^n c_j (x - x_0)^{p_j} \]
    - coefficients $c_j, j = 0, \dots n$
    - integer exponents $p_j, j = 0, \dots n$
    - reference point $x_0$

  <ParameterList name="function-polynomial">
    <Parameter name="coefficients" type="Array double" value=double-array />
    <Parameter name="exponents" type="Array int" value=int-array />
    <Parameter name="reference point" type="double" value=double />
  </ParameterList>
  
  * Multi-variable linear function: \[ f(x) = y_0 + \sum_{j=0}^{n-1} g_j (x_j - x_{0,j}) \]
    - constant term y_0 given by "y0"
    - linear term coefficients $(g_0, g_1, \dots, g_{n-1})$ given by "gradient"
    - reference point $x_0$ given by "x0".  If specified it must have the same
      number of values as "gradient".  Defaults to zero.

  <ParameterList name="function-linear">
    <Parameter name="y0" type="double" value=double />
    <Parameter name="gradient" type="Array double" value=double-array />
    <Parameter name="x0" type="Array double" value=double-array />
  </ParameterList>
  
  * Separable function: $f(x_0, x_1,\dots,x_{n-1}) = f_1(x_0) f_2(x_1,\dots,x_{n-1})$
    - $f_1$ is defined by the "function1" sublist, and $f_2$ by the "function2" sublist.

  <ParameterList name="function-separable">
    <ParameterList name="function1">
      function-specification
    </ParameterList>
    <ParameterList name="function2">
      function-specification
    </ParameterList>
  </ParameterList>
