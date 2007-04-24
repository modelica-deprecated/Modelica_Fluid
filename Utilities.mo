package Utilities 
  "Utility models to construct fluid components (should not be used directly) " 
  extends Modelica.Icons.Library;
  
  function checkBoundary "Check whether boundary definition is correct" 
    extends Modelica.Icons.Function;
    input String mediumName;
    input String substanceNames[:] "Names of substances";
    input Boolean singleState;
    input Boolean define_p;
    input Real X_boundary[:];
    input String modelName = "??? boundary ???";
  protected 
    Integer nX = size(X_boundary,1);
    String X_str;
  algorithm 
    assert(not singleState or singleState and define_p, "
Wrong value of parameter define_p (= false) in model \""   + modelName + "\":
The selected medium \"" + mediumName + "\" has Medium.singleState=true.
Therefore, an boundary density cannot be defined and
define_p = true is required.
");
    
    for i in 1:nX loop
      assert(X_boundary[i] >= 0.0, "
Wrong boundary mass fractions in medium \""
  + mediumName + "\" in model \"" + modelName + "\":
The boundary value X_boundary("
                              + String(i) + ") = " + String(
        X_boundary[i]) + "
is negative. It must be positive.
");
    end for;
    
    if nX > 0 and abs(sum(X_boundary) - 1.0) > 1.e-10 then
       X_str :="";
       for i in 1:nX loop
          X_str :=X_str + "   X_boundary[" + String(i) + "] = " + String(X_boundary[
          i]) + " \"" + substanceNames[i] + "\"\n";
       end for;
       Modelica.Utilities.Streams.error(
          "The boundary mass fractions in medium \"" + mediumName + "\" in model \"" + modelName + "\"\n" +
          "do not sum up to 1. Instead, sum(X_boundary) = " + String(sum(X_boundary)) + ":\n"
          + X_str);
    end if;
  end checkBoundary;
  
  function regRoot 
    "Anti-symmetric square root approximation with finite derivative in the origin" 
    extends Modelica.Icons.Function;
    input Real x;
    input Real delta=0.01 
      "Range of significant deviation from sqrt(abs(x))*sgn(x)";
    output Real y;
    annotation(derivative(zeroDerivative=delta)=Utilities.regRoot_der,
      Documentation(info="<html>
This function approximates sqrt(abs(x))*sgn(x), such that the derivative is finite and smooth in x=0. 
</p>
<p>
<table border=1 cellspacing=0 cellpadding=2> 
<tr><th>Function</th><th>Approximation</th><th>Range</th></tr>
<tr><td>y = regRoot(x)</td><td>y ~= sqrt(abs(x))*sgn(x)</td><td>abs(x) &gt;&gt delta</td></tr>
<tr><td>y = regRoot(x)</td><td>y ~= x/sqrt(delta)</td><td>abs(x) &lt;&lt  delta</td></tr>
</table>
<p>
With the default value of delta=0.01, the difference between sqrt(x) and regRoot(x) is 16% around x=0.01, 0.25% around x=0.1 and 0.0025% around x=1.
</p> 
</html>",
        revisions="<html>
<ul>
<li><i>15 Mar 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Created. </li>
</ul>
</html>"));
  algorithm 
    y := x/(x*x+delta*delta)^0.25;
    
  annotation (Documentation(info="<html>
This function approximates sqrt(x)*sign(x), such that the derivative is finite and smooth in x=0. 
</p>
<p>
<table border=1 cellspacing=0 cellpadding=2> 
<tr><th>Function</th><th>Approximation</th><th>Range</th></tr>
<tr><td>y = sqrtReg(x)</td><td>y ~= sqrt(abs(x))*sign(x)</td><td>abs(x) &gt;&gt delta</td></tr>
<tr><td>y = sqrtReg(x)</td><td>y ~= x/delta</td><td>abs(x) &lt;&lt  delta</td></tr>
</table>
<p>
With the default value of delta=0.01, the difference between sqrt(x) and sqrtReg(x) is 0.5% around x=0.1 and 0.005% around x=1.
</p> 
</html>",
        revisions="<html>
<ul>
<li><i>15 Mar 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Created. </li>
</ul>
</html>"));
  end regRoot;
  
  function regRoot_der "Derivative of regRoot" 
    extends Modelica.Icons.Function;
    input Real x;
    input Real delta=0.01 "Range of significant deviation from sqrt(x)";
    input Real dx "Derivative of x";
    output Real dy;
  algorithm 
    dy := dx*0.5*(x*x+2*delta*delta)/((x*x+delta*delta)^1.25);
    
  annotation (Documentation(info="<html>
</html>",
        revisions="<html>
<ul>
<li><i>15 Mar 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Created. </li>
</ul>
</html>"));
  end regRoot_der;
  
  function regSquare 
    "Anti-symmetric square approximation with non-zero derivative in the origin" 
    extends Modelica.Icons.Function;
    input Real x;
    input Real delta=0.01 "Range of significant deviation from x^2*sgn(x)";
    output Real y;
    annotation(Documentation(info="<html>
This function approximates x^2*sgn(x), such that the derivative is non-zero in x=0. 
</p>
<p>
<table border=1 cellspacing=0 cellpadding=2> 
<tr><th>Function</th><th>Approximation</th><th>Range</th></tr>
<tr><td>y = regSquare(x)</td><td>y ~= x^2*sgn(x)</td><td>abs(x) &gt;&gt delta</td></tr>
<tr><td>y = regSquare(x)</td><td>y ~= x*delta</td><td>abs(x) &lt;&lt  delta</td></tr>
</table>
<p>
With the default value of delta=0.01, the difference between x^2 and regSquare(x) is 41% around x=0.01, 0.4% around x=0.1 and 0.005% around x=1.
</p> 
</p> 
</html>",
        revisions="<html>
<ul>
<li><i>15 Mar 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Created. </li>
</ul>
</html>"));
  algorithm 
    y := x*sqrt(x*x+delta*delta);
    
  annotation (Documentation(info="<html>
This function approximates sqrt(x)*sign(x), such that the derivative is finite and smooth in x=0. 
</p>
<p>
<table border=1 cellspacing=0 cellpadding=2> 
<tr><th>Function</th><th>Approximation</th><th>Range</th></tr>
<tr><td>y = sqrtReg(x)</td><td>y ~= sqrt(abs(x))*sign(x)</td><td>abs(x) &gt;&gt delta</td></tr>
<tr><td>y = sqrtReg(x)</td><td>y ~= x/delta</td><td>abs(x) &lt;&lt  delta</td></tr>
</table>
<p>
With the default value of delta=0.01, the difference between sqrt(x) and sqrtReg(x) is 0.5% around x=0.1 and 0.005% around x=1.
</p> 
</html>",
        revisions="<html>
<ul>
<li><i>15 Mar 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Created. </li>
</ul>
</html>"));
  end regSquare;
  
  function regPow 
    "Anti-symmetric power approximation with non-zero derivative in the origin" 
    extends Modelica.Icons.Function;
    input Real x;
    input Real a;
    input Real delta=0.01 "Range of significant deviation from x^a*sgn(x)";
    output Real y;
    annotation(Documentation(info="<html>
This function approximates abs(x)^a*sign(x), such that the derivative is positive, finite and smooth in x=0. 
</p>
<p>
<table border=1 cellspacing=0 cellpadding=2> 
<tr><th>Function</th><th>Approximation</th><th>Range</th></tr>
<tr><td>y = regPow(x)</td><td>y ~= abs(x)^a*sgn(x)</td><td>abs(x) &gt;&gt delta</td></tr>
<tr><td>y = regPow(x)</td><td>y ~= x*delta^(a-1)</td><td>abs(x) &lt;&lt  delta</td></tr>
</table>
</html>",
        revisions="<html>
<ul>
<li><i>15 Mar 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Created. </li>
</ul>
</html>"));
  algorithm 
    y := x*(x*x+delta*delta)^((a-1)/2);
    
  annotation (Documentation(info="<html>
This function approximates sqrt(x)*sign(x), such that the derivative is finite and smooth in x=0. 
</p>
<p>
<table border=1 cellspacing=0 cellpadding=2> 
<tr><th>Function</th><th>Approximation</th><th>Range</th></tr>
<tr><td>y = sqrtReg(x)</td><td>y ~= sqrt(abs(x))*sign(x)</td><td>abs(x) &gt;&gt delta</td></tr>
<tr><td>y = sqrtReg(x)</td><td>y ~= x/delta</td><td>abs(x) &lt;&lt  delta</td></tr>
</table>
<p>
With the default value of delta=0.01, the difference between sqrt(x) and sqrtReg(x) is 0.5% around x=0.1 and 0.005% around x=1.
</p> 
</html>",
        revisions="<html>
<ul>
<li><i>15 Mar 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Created. </li>
</ul>
</html>"));
  end regPow;
  
  function regRoot2 
    "Anti-symmetric approximation of square root with discontinuous factor so that the first derivative is finite and continuous" 
    
    extends Modelica.Icons.Function;
    input Real x "abszissa value";
    input Real x_small(min=0)=0.01 
      "approximation of function for |x| <= x_small";
    input Real k1(min=0)=1 "y = if x>=0 then sqrt(k1*x) else -sqrt(k2*|x|)";
    input Real k2(min=0)=1 "y = if x>=0 then sqrt(k1*x) else -sqrt(k2*|x|)";
    input Boolean use_yd0 = false "= true, if yd0 shall be used";
    input Real yd0(min=0)=1 "Desired derivative at x=0: dy/dx = yd0";
    output Real y "ordinate value";
    annotation(smoothOrder=2, Documentation(info="<html>
<p>
Approximates the function
</p>
<pre>
   y = <b>if</b> x &ge; 0 <b>then</b> <b>sqrt</b>(k1*x) <b>else</b> -<b>sqrt</b>(k2*<b>abs</b>(x)), with k1, k2 > 0
</pre>
<p>
in such a way that within the region -x_small &le; x &le; x_small, 
the function is described by two polynomials of third order
(one in the region -x_small .. 0 and one within the region 0 .. x_small)
such that 
</p>
<ul>
<li> The derivative at x=0 is finite. </li>
<li> The overall function is continuous with a
     continuous first derivative everywhere.</li>
<li> If parameter use_yd0 = <b>false</b>, the two polynomials
     are constructed such that the second derivatives at x=0
     are identical. If use_yd0 = <b>true</b>, the derivative
     at x=0 is explicitly provided via the additional argument
     yd0. If necessary, the derivative yd0 is automatically 
     reduced in order that the polynomials are strict monotonically  
     increasing <i>[Fritsch and Carlson, 1980]</i>.</li>
</ul>
<p>
Typical screenshots for two different configurations
are shown below. The first one with k1=k2=1:
</p>
<p>
<img src=\"../Images/Components/regRoot2_a.png\">
</p>
<p>
and the second one with k1=1 and k2=3:
</p>
<p>
<img src=\"../Images/Components/regRoot2_b.png\">
</p>
 
<p>
The (smooth) derivative of the function with
k1=1, k2=3 is shown in the next figure:
<p>
<img src=\"../Images/Components/regRoot2_c.png\">
</p>
 
<p>
<b>Literature</b>
</p>
 
<dl>
<dt> Fritsch F.N. and Carlson R.E. (1980):</dt>
<dd> <b>Monotone piecewise cubic interpolation</b>.
     SIAM J. Numerc. Anal., Vol. 17, No. 2, April 1980, pp. 238-246</dd>
</dl>
</html>", revisions="<html>
<ul>
<li><i>Nov., 2005</i>
    by <a href=\"mailto:Martin.Otter@DLR.de\">Martin Otter</a>:<br>
    Designed and implementated.</li>
</ul>
</html>"));
  protected 
    encapsulated function regRoot2_utility 
      "Interpolating with two 3-order polynomials with a prescribed derivative at x=0" 
      import Modelica_Fluid.Utilities.evaluatePoly3_derivativeAtZero;
       input Real x;
       input Real x1 "approximation of function abs(x) < x1";
       input Real k1 "y = if x>=0 then sqrt(k1*x) else -sqrt(k2*|x|); k1 >= k2";
       input Real k2 "y = if x>=0 then sqrt(k1*x) else -sqrt(k2*|x|))";
       input Boolean use_yd0 "= true, if yd0 shall be used";
       input Real yd0(min=0) "Desired derivative at x=0: dy/dx = yd0";
       output Real y;
       annotation(smoothOrder=2);
    protected 
       Real x2;
       Real xsqrt1;
       Real xsqrt2;
       Real y1;
       Real y2;
       Real y1d;
       Real y2d;
       Real w;
       Real y0d;
       Real w1;
       Real w2;
    algorithm 
       x2 :=-x1*(k2/k1);
       //x2 :=-x1;
       if x <= x2 then
          y := -sqrt(k2*abs(x));
       else
          y1 :=sqrt(k1*x1);
          y2 :=-sqrt(k2*abs(x2));
          y1d :=sqrt(k1/x1)/2;
          y2d :=sqrt(k2/abs(x2))/2;
        
          if use_yd0 then
             y0d :=yd0;
          else
             /* Determine derivative, such that first and second derivative
              of left and right polynomial are identical at x=0:
           _
           Basic equations:
              y_right = a1*(x/x1) + a2*(x/x1)^2 + a3*(x/x1)^3
              y_left  = b1*(x/x2) + b2*(x/x2)^2 + b3*(x/x2)^3
              yd_right*x1 = a1 + 2*a2*(x/x1) + 3*a3*(x/x1)^2
              yd_left *x2 = b1 + 2*b2*(x/x2) + 3*b3*(x/x2)^2
              ydd_right*x1^2 = 2*a2 + 6*a3*(x/x1)
              ydd_left *x2^2 = 2*b2 + 6*b3*(x/x2)
           _
           Conditions (6 equations for 6 unknowns):
                     y1 = a1 + a2 + a3
                     y2 = b1 + b2 + b3
                 y1d*x1 = a1 + 2*a2 + 3*a3
                 y2d*x2 = b1 + 2*b2 + 3*b3
                    y0d = a1/x1 = b1/x2
                   y0dd = 2*a2/x1^2 = 2*b2/x2^2
           _
           Derived equations:
              b1 = a1*x2/x1
              b2 = a2*(x2/x1)^2
              b3 = y2 - b1 - b2
                 = y2 - a1*(x2/x1) - a2*(x2/x1)^2
              a3 = y1 - a1 - a2
           _
           Remaining equations
              y1d*x1 = a1 + 2*a2 + 3*(y1 - a1 - a2)
                     = 3*y1 - 2*a1 - a2
              y2d*x2 = a1*(x2/x1) + 2*a2*(x2/x1)^2 +
                       3*(y2 - a1*(x2/x1) - a2*(x2/x1)^2)
                     = 3*y2 - 2*a1*(x2/x1) - a2*(x2/x1)^2
              y0d    = a1/x1
           _
           Solving these equations results in y0d below
           (note, the denominator "(1-w)" is always non-zero, because w is negative) 
           */
             w :=x2/x1;
             y0d := ( (3*y2 - x2*y2d)/w - (3*y1 - x1*y1d)*w) /(2*x1*(1 - w));
          end if;
        
          /* Modify derivative y0d, such that the polynomial is 
           monotonically increasing. A sufficient condition is
             0 <= y0d <= sqrt(8.75*k_i/|x_i|)
        */
          w1 :=sqrt(8.75*k1/x1);
          w2 :=sqrt(8.75*k2/abs(x2));
          y0d :=min(y0d, 0.9*min(w1, w2));
        
          /* Perform interpolation in scaled polynomial:
           y_new = y/y1
           x_new = x/x1
        */
          y := y1*(if x >= 0 then evaluatePoly3_derivativeAtZero(x/x1,1,1,y1d*x1/y1,y0d*x1/y1) else 
                                  evaluatePoly3_derivativeAtZero(x/x1,x2/x1,y2/y1,y2d*x1/y1,y0d*x1/y1));
       end if;
    end regRoot2_utility;
  algorithm 
    y := smooth(2,if x >= x_small then sqrt(k1*x) else 
                  if x <= -x_small then -sqrt(k2*abs(x)) else 
                  if k1 >= k2 then regRoot2_utility(x,x_small,k1,k2,use_yd0,yd0) else 
                                  -regRoot2_utility(-x,x_small,k2,k1,use_yd0,yd0));
  end regRoot2;
  
  function regSquare2 
    "Anti-symmetric approximation of square with discontinuous factor so that the first derivative is non-zero and is continuous" 
    extends Modelica.Icons.Function;
    input Real x "abszissa value";
    input Real x_small(min=0)=0.01 
      "approximation of function for |x| <= x_small";
    input Real k1(min=0)=1 "y = (if x>=0 then k1 else k2)*x*|x|";
    input Real k2(min=0)=1 "y = (if x>=0 then k1 else k2)*x*|x|";
    input Boolean use_yd0 = false "= true, if yd0 shall be used";
    input Real yd0(min=0)=1 "Desired derivative at x=0: dy/dx = yd0";
    output Real y "ordinate value";
    annotation(smoothOrder=2, Documentation(info="<html>
<p>
Approximates the function
</p>
<pre>
   y = <b>if</b> x &ge; 0 <b>then</b> k1*x*x <b>else</b> -k2*x*x, with k1, k2 > 0
</pre>
<p>
in such a way that within the region -x_small &le; x &le; x_small, 
the function is described by two polynomials of third order
(one in the region -x_small .. 0 and one within the region 0 .. x_small)
such that
</p>
 
<ul>
<li> The derivative at x=0 is non-zero (in order that the
     inverse of the function does not have an infinite derivative). </li>
<li> The overall function is continuous with a
     continuous first derivative everywhere.</li>
<li> If parameter use_yd0 = <b>false</b>, the two polynomials
     are constructed such that the second derivatives at x=0
     are identical. If use_yd0 = <b>true</b>, the derivative
     at x=0 is explicitly provided via the additional argument
     yd0. If necessary, the derivative yd0 is automatically 
     reduced in order that the polynomials are strict monotonically  
     increasing <i>[Fritsch and Carlson, 1980]</i>.</li>
</ul>
</ul>
<p>
A typical screenshot for k1=1, k2=3 is shown in the next figure:
</p>
<p>
<img src=\"../Images/Components/regSquare2_b.png\">
</p>
 
<p>
The (smooth, non-zero) derivative of the function with
k1=1, k2=3 is shown in the next figure:
</p>
 
<p>
<img src=\"../Images/Components/regSquare2_c.png\">
</p>
 
<p>
<b>Literature</b>
</p>
 
<dl>
<dt> Fritsch F.N. and Carlson R.E. (1980):</dt>
<dd> <b>Monotone piecewise cubic interpolation</b>.
     SIAM J. Numerc. Anal., Vol. 17, No. 2, April 1980, pp. 238-246</dd>
</dl>
</html>", revisions="<html>
<ul>
<li><i>Nov., 2005</i>
    by <a href=\"mailto:Martin.Otter@DLR.de\">Martin Otter</a>:<br>
    Designed and implementated.</li>
</ul>
</html>"));
  protected 
    encapsulated function regSquare2_utility 
      "Interpolating with two 3-order polynomials with a prescribed derivative at x=0" 
      import Modelica_Fluid.Utilities.evaluatePoly3_derivativeAtZero;
       input Real x;
       input Real x1 "approximation of function abs(x) < x1";
       input Real k1 "y = (if x>=0 then k1 else -k2)*x*|x|; k1 >= k2";
       input Real k2 "y = (if x>=0 then k1 else -k2)*x*|x|";
       input Boolean use_yd0 = false "= true, if yd0 shall be used";
       input Real yd0(min=0)=1 "Desired derivative at x=0: dy/dx = yd0";
       output Real y;
       annotation(smoothOrder=2);
    protected 
       Real x2;
       Real y1;
       Real y2;
       Real y1d;
       Real y2d;
       Real w;
       Real w1;
       Real w2;
       Real y0d;
    algorithm 
       // x2 :=-x1*(k2/k1)^2;
       x2 := -x1;
       if x <= x2 then
          y := -k2*x^2;
       else
           y1 := k1*x1^2;
           y2 :=-k2*x2^2;
          y1d := k1*2*x1;
          y2d :=-k2*2*x2;
          if use_yd0 then
             y0d :=yd0;
          else
             /* Determine derivative, such that first and second derivative
              of left and right polynomial are identical at x=0:
              see derivation in function regRoot2
           */
             w :=x2/x1;
             y0d := ( (3*y2 - x2*y2d)/w - (3*y1 - x1*y1d)*w) /(2*x1*(1 - w));
          end if;
        
          /* Modify derivative y0d, such that the polynomial is 
           monotonically increasing. A sufficient condition is
             0 <= y0d <= sqrt(5)*k_i*|x_i|
        */
          w1 :=sqrt(5)*k1*x1;
          w2 :=sqrt(5)*k2*abs(x2);
          y0d :=min(y0d, 0.9*min(w1, w2));
        
          y := if x >= 0 then evaluatePoly3_derivativeAtZero(x,x1,y1,y1d,y0d) else 
                              evaluatePoly3_derivativeAtZero(x,x2,y2,y2d,y0d);
       end if;
    end regSquare2_utility;
  algorithm 
    y := smooth(2,if x >= x_small then k1*x^2 else 
                  if x <= -x_small then -k2*x^2 else 
                  if k1 >= k2 then regSquare2_utility(x,x_small,k1,k2,use_yd0,yd0) else 
                                  -regSquare2_utility(-x,x_small,k2,k1,use_yd0,yd0));
  end regSquare2;
  
  function evaluatePoly3_derivativeAtZero 
    "Evaluate polynomial of order 3 that passes the origin with a predefined derivative" 
    extends Modelica.Icons.Function;
    input Real x "Value for which polynomial shall be evaluated";
    input Real x1 "Abscissa value";
    input Real y1 "y1=f(x1)";
    input Real y1d "First derivative at y1";
    input Real y0d "First derivative at f(x=0)";
    output Real y;
    annotation(smoothOrder=3, Documentation(info="<html>
 
</html>"));
  protected 
    Real a1;
    Real a2;
    Real a3;
    Real xx;
  algorithm 
    a1 := x1*y0d;
    a2 := 3*y1 - x1*y1d - 2*a1;
    a3 := y1 - a2 - a1;
    xx := x/x1;
    y  := xx*(a1 + xx*(a2 + xx*a3));
  end evaluatePoly3_derivativeAtZero;
  
  function ReynoldsNumber_m_flow 
    "Return Reynolds number as a function of mass flow rate m_flow" 
    extends Modelica.Icons.Function;
      input SI.MassFlowRate m_flow "Mass flow rate";
    input SI.DynamicViscosity eta "Dynamic viscosity of medium";
    input SI.Diameter diameter "Diameter of pipe/orifice";
    output SI.ReynoldsNumber Re "Reynolds number";
  algorithm 
    Re :=abs(m_flow)*(4/Modelica.Constants.pi)/(diameter*eta);
    annotation (Documentation(info="<html>
 
</html>"));
  end ReynoldsNumber_m_flow;
  
  annotation (Documentation(info="<html>
 
</html>"));
end Utilities;
