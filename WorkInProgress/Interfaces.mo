package Interfaces 
  partial model PartialTwoPortTransport 
    "Partial element transporting fluid between two ports without storing mass or energy" 
    import Modelica.SIunits.*;
    import Modelica.Constants.*;
    import Modelica_Fluid;
    replaceable package Medium = PackageMedium extends 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component"  annotation (
        choicesAllMatching =                                                                            true);
    parameter Boolean allowFlowReversal = true 
      "Flow reversal at the ports is allowed by the equations"  annotation(Dialog(tab="Advanced"));
    
    Modelica_Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium = 
          Medium, m_flow(min=if allowFlowReversal then -inf else 0)) 
      annotation (extent=[-120, -10; -100, 10]);
    Modelica_Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium = 
          Medium, m_flow(max=if allowFlowReversal then +inf else 0)) 
      annotation (extent=[120, -10; 100, 10]);
    Medium.BaseProperties medium_a "Medium properties in port_a";
    Medium.BaseProperties medium_b "Medium properties in port_b";
    annotation (
      Coordsys(grid=[1, 1], component=[20, 20]),
      Diagram,
      Documentation(info="<html>
<p>
This component transports fluid between its two ports, without
storing mass or energy. Reversal and zero mass flow rate is taken
care of, for details see definition of built-in operator semiLinear().
<p>
When using this partial component, an equation for the momentum
balance has to be added by specifying a relationship
between the pressure drop <tt>dp</tt> and the mass flow rate <tt>m_flow</tt>.
</p>
</html>"));
  equation 
    // Properties in the ports
    port_a.p   = medium_a.p;
    port_a.h   = medium_a.h;
    port_a.Xi = medium_a.Xi;
    port_b.p   = medium_b.p;
    port_b.h   = medium_b.h;
    port_b.Xi = medium_b.Xi;
    
    /* Handle reverse and zero flow */
    port_a.H_flow   = semiLinear(port_a.m_flow, port_a.h,  port_b.h);
    port_a.mXi_flow = semiLinear(port_a.m_flow, port_a.Xi, port_b.Xi);
    
    /* Energy, mass and substance mass balance */
    port_a.H_flow + port_b.H_flow = 0;
    port_a.m_flow + port_b.m_flow = 0;
    port_a.mXi_flow + port_b.mXi_flow = zeros(Medium.nXi);
  end PartialTwoPortTransport;
  
  partial model PartialPressureLoss 
    "Generic pressure loss with constant turbulent loss factors" 
    import SI = Modelica.SIunits;
    import Modelica_Fluid.Utilities.regRoot2;
    import Modelica_Fluid.Utilities.regSquare2;
    
    parameter Modelica_Fluid.WorkInProgress.Utilities.PressureLossFactors 
      lossFactors "Loss factors for both flow directions";
    parameter Boolean from_dp=true 
      " = true, use m_flow = f(dp) else dp = f(m_flow)" 
      annotation (Evaluate=true, Dialog(tab="Advanced"));
    parameter Boolean use_Re = true 
      "= true, if turbulent region is defined by Re, otherwise by p_small or m_flow_small"
      annotation(Evaluate=true, Dialog(tab="Advanced"));
    parameter SI.AbsolutePressure dp_small = 1 
      "Turbulent flow if |dp| >= dp_small" 
      annotation(Dialog(tab="Advanced", enable=not use_Re and from_dp));
    parameter SI.MassFlowRate m_flow_small = 0.001 
      "Turbulent flow if |m_flow| >= m_flow_small" 
      annotation(Dialog(tab="Advanced", enable=not use_Re and not from_dp));
    SI.MassFlowRate m_flow(start=0) "Mass flow rate from port_a to port_b";
    SI.Pressure dp(start=0) "Pressure loss due to pressure loss component";
    input SI.Density d_a "Density at port_a"     annotation(Hide=true);
    input SI.Density d_b "Density at port_b"     annotation(Hide=true);
    input SI.DynamicViscosity eta_a 
      "Dynamic viscosity at port_a, if use_Re = true (otherwise not used)";
    input SI.DynamicViscosity eta_b 
      "Dynamic viscosity at port_b, if use_Re = true (otherwise not used)";
    final SI.ReynoldsNumber Re = noEvent(abs(m_flow))*(4/pi)/(lossFactors.D_Re*eta) if use_Re 
      "Reynolds number at smallest pipe diameter";
    SI.AbsolutePressure dp_turbulent 
      "The turbulent region is: |dp| >= dp_turbulent";
    SI.MassFlowRate m_flow_turbulent 
      "The turbulent region is: |m_flow| >= m_flow_turbulent";
    annotation (
      Diagram,
      Icon,
      Documentation(info="<html>
</html>"));
    
  protected 
    constant Real pi=Modelica.Constants.pi annotation(Hide=true);
    parameter SI.Diameter LD_a = lossFactors.D_a 
                                                annotation(Hide=true);
    parameter SI.Diameter LD_b = lossFactors.D_b 
                                                annotation(Hide=true);
    parameter Real k0=2*lossFactors.c0/(pi*lossFactors.D_Re^3);
    parameter Real k1=if lossFactors.zeta1_at_a then 8*lossFactors.zeta1/(pi*LD_a^2)^2 else 
                                                     8*lossFactors.zeta1/(pi*LD_b^2)^2;
    parameter Real k2=if lossFactors.zeta2_at_a then 8*lossFactors.zeta2/(pi*LD_a^2)^2 else 
                                                     8*lossFactors.zeta2/(pi*LD_b^2)^2;
    parameter Boolean withLaminar = lossFactors.zetaLaminarKnown annotation(Evaluate=true,Hide=true);
    Real yd0 
      "Derivative of dp=dp(m_flow) or m_flow=m_flow(dp) at zero, if use_Re and lossFactors.zetaLaminarKnown";
    SI.DynamicViscosity eta = if use_Re then (eta_a + eta_b)/2 else 0;
    
  equation 
  /*
Turbulent region:
   Re = m_flow*(4/pi)/(D_Re*eta)  
   dp = 0.5*zeta*d*v*|v|
      = 0.5*zeta*d*1/(d*A)^2 * m_flow * |m_flow|
      = 0.5*zeta/A^2 *1/d * m_flow * |m_flow|
      = k/d * m_flow * |m_flow| 
   k  = 0.5*zeta/A^2
      = 0.5*zeta/(pi*(D/2)^2)^2
      = 8*zeta/(pi*D^2)^2
   m_flow_turbulent = (pi/4)*D_Re*eta*Re_turbulent
   dp_turbulent     =  k/d *(D_Re*eta*pi/4)^2 * Re_turbulent^2
_
   The start of the turbulent region is computed with mean values
   of dynamic viscosity eta and density rho. Otherwise, one has
   to introduce different "delta" values for both flow directions.
   In order to simplify the approach, only one delta is used.  
_
Laminar region:
   dp = 0.5*zeta/(A^2*d) * m_flow * |m_flow|
      = 0.5 * c0/(|m_flow|*(4/pi)/(D_Re*eta)) / ((pi*(D_Re/2)^2)^2*d) * m_flow*|m_flow|
      = 0.5 * c0*(pi/4)*(D_Re*eta) * 16/(pi^2*D_Re^4*d) * m_flow*|m_flow|
      = 2*c0/(pi*D_Re^3) * eta/d * m_flow
      = k0 * eta/d * m_flow
   k0 = 2*c0/(pi*D_Re^3)
_
   In order that the derivative of dp=f(m_flow) is continuous 
   at m_flow=0, the mean values of eta and d are used in the
   laminar region: eta/d = (eta_a + eta_b)/(d_a + d_b)
   If lossFactors.zetaLaminarKnown = false then eta_a and eta_b are potentially zero
   (because dummy values) and therefore the division is only performed
   if zetaLaminarKnown = true.
*/
     if from_dp then
        dp_turbulent = if use_Re then (k1+k2)/(d_a+d_b)*(eta*lossFactors.D_Re*pi/4)^2
                                     *lossFactors.Re_turbulent^2 else dp_small;
        m_flow_turbulent = regRoot2(dp_turbulent, dp_turbulent, d_a/k1, d_b/k2) 
        "for information purposes";
        yd0 = if use_Re and withLaminar then (d_a+d_b)/(k0*(eta_a+eta_b)) else 0;
        m_flow = regRoot2(dp, dp_turbulent, d_a/k1, d_b/k2, use_Re and withLaminar, yd0);
     else
        m_flow_turbulent = if use_Re then (pi/4)*lossFactors.D_Re*eta*lossFactors.Re_turbulent else m_flow_small;
        dp_turbulent = regSquare2(m_flow_turbulent, m_flow_turbulent, k1/d_a, k2/d_b) 
        "for information purposes";
        yd0 = if use_Re and withLaminar then k0*(eta_a+eta_b)/(d_a+d_b) else 0;
        dp = regSquare2(m_flow, m_flow_turbulent, k1/d_a, k2/d_b, use_Re and withLaminar, yd0);
     end if;
    
  end PartialPressureLoss;
  
  model PressureLossWithoutIcon 
    "Generic pressure loss component with constant turbulent loss factors and without an icon" 
    extends Modelica_Fluid.WorkInProgress.Interfaces.PartialTwoPortTransport;
    extends Modelica_Fluid.WorkInProgress.Interfaces.PartialPressureLoss(
       m_flow = port_a.m_flow,
       dp = port_a.p - port_b.p,
       d_a = medium_a.d,
       d_b = medium_b.d,
       eta_a = if use_Re then Medium.dynamicViscosity(medium_a.state) else 0,
       eta_b = if use_Re then Medium.dynamicViscosity(medium_b.state) else 0);
    annotation (
      Diagram,
      Icon,
      Documentation(info="<html>
<p>
This model computes the pressure loss of a pipe
segment (orifice, bending etc.) with a minimum amount of data
provided via parameter <b>lossFactors</b>.
If available, data should be provided for <b>both flow directions</b>,
i.e., flow from port_a to port_b and from port_b to port_a, 
as well as for the <b>laminar</b> and the <b>turbulent</b> region.
It is also an option to provide the loss factor <b>only</b> for the
<b>turbulent</b> region for a flow from port_a to port_b.
</p>
<p>
The following equations are used:
</p>
<pre>   &Delta;p = 0.5*&zeta;*&rho;*v*|v|
      = 0.5*&zeta;/A^2 * (1/&rho;) * m_flow*|m_flow|
        Re = |v|*D*&rho;/&eta;
</pre>
<table border=1 cellspacing=0 cellpadding=2>
<tr><td><b>flow type</b></td>
    <td><b>&zeta;</b> = </td>
    <td><b>flow region</b></td></tr>
<tr><td>turbulent</td>
    <td><b>zeta1</b> = const.</td>
    <td>Re &ge;  Re_turbulent, v &ge; 0</td></tr>
<tr><td></td>
    <td><b>zeta2</b> = const.</td>
    <td>Re &ge; Re_turbulent, v &lt; 0</td></tr>
<tr><td>laminar</td>
    <td><b>c0</b>/Re</td>
    <td>both flow directions, Re small; c0 = const.</td></tr>
</table>
<p>
where
</p>
<ul>
<li> &Delta;p is the pressure drop: &Delta;p = port_a.p - port_b.p</li>
<li> v is the mean velocity.</li>
<li> &rho; is the density.</li>
<li> &zeta; is the loss factor that depends on the geometry of
     the pipe. In the turbulent flow regime, it is assumed that
     &zeta; is constant and is given by \"zeta1\" and
     \"zeta2\" depending on the flow direction.<br>
     When the Reynolds number Re is below \"Re_turbulent\", the
     flow is laminar for small flow velocities. For higher 
     velocities there is a transition region from 
     laminar to turbulent flow. The loss factor for
     laminar flow at small velocities is defined by the often occuring
     approximation c0/Re. If c0 is different for the two
     flow directions, the mean value has to be used 
     (c0 = (c0_ab + c0_ba)/2).<li>
<li> The equation \"&Delta;p = 0.5*&zeta;*&rho;*v*|v|\" is either with
     respect to port_a or to port_b, depending on the definition
     of the particular loss factor &zeta; (in some references loss
     factors are defined with respect to port_a, in other references
     with respect to port_b).</li>
 
<li> Re = |v|*D_Re*&rho;/&eta; = |m_flow|*D_Re/(A_Re*&eta;) 
     is the Reynolds number at the smallest cross
     section area. This is often at port_a or at port_b, but can
     also be between the two ports. In the record, the diameter
     D_Re of this smallest cross section area has to be provided, as
     well, as Re_turbulent, the absolute value of the 
     Reynolds number at which
     the turbulent flow starts. If Re_turbulent is different for
     the two flow directions, use the smaller value as Re_turbulent.</li>
<li> D is the diameter of the pipe. If the pipe has not a 
     circular cross section, D = 4*A/P, where A is the cross section
     area and P is the wetted perimeter.</li>
<li> A is the cross section area with A = &pi;(D/2)^2.
<li> &eta; is the dynamic viscosity.</li>
</ul>
<p>
The laminar and the transition region is usually of
not much technical interest because the operating point is
mostly in the turbulent regime. For simplification and for
numercial reasons, this whole region is described by two
polynomials of third order, one polynomial for m_flow &ge; 0 
and one for m_flow &lt; 0. The polynomials start at 
Re = |m_flow|*4/(&pi;*D_Re*&eta;), where D_Re is the
smallest diameter between port_a and port_b.
The common derivative
of the two polynomials at Re = 0 is
computed from the equation \"c0/Re\". Note, the pressure drop
equation above in the laminar region is always defined
with respect to the smallest diameter D_Re.
</p>
<p>
If no data for c0 is available, the derivative at Re = 0 is computed in such
a way, that the second derivatives of the two polynomials
are identical at Re = 0. The polynomials are constructed, such that
they smoothly touch the characteristic curves in the turbulent
regions. The whole characteristic is therefore <b>continuous</b>
and has a <b>finite</b>, <b>continuous first derivative everywhere</b>.
In some cases, the constructed polynomials would \"vibrate\". This is 
avoided by reducing the derivative at Re=0 in such a way that
the polynomials are guaranteed to be monotonically increasing.
The used sufficient criteria for monotonicity follows from:
</p>
 
<dl>
<dt> Fritsch F.N. and Carlson R.E. (1980):</dt>
<dd> <b>Monotone piecewise cubic interpolation</b>.
     SIAM J. Numerc. Anal., Vol. 17, No. 2, April 1980, pp. 238-246</dd>
</dl>
</html>"));
  end PressureLossWithoutIcon;
  
  
  
  partial model PartialTwoPortTransportWithDz 
    "Partial element transporting fluid between two ports without storing mass or energy" 
    import Modelica.SIunits.*;
    import Modelica.Constants.*;
    replaceable package Medium = 
        Modelica.Media.Interfaces.PartialTwoPhaseMedium                          extends 
      Modelica.Media.Interfaces.PartialTwoPhaseMedium "Medium in the component"
                                                                         annotation (
        choicesAllMatching =                                                                            true);
    parameter Boolean allowFlowReversal = true 
      "Flow reversal at the ports is allowed by the equations";
    
    Modelica_Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium = 
          Medium, m_flow(min=if allowFlowReversal then -inf else 0)) 
      annotation (extent=[-120, -10; -100, 10]);
    Modelica_Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium = 
          Medium, m_flow(max=if allowFlowReversal then +inf else 0)) 
      annotation (extent=[120, -10; 100, 10]);
    Medium.BaseProperties medium_a "Medium properties in port_a";
    Medium.BaseProperties medium_b "Medium properties in port_b";
    Medium.MassFlowRate m_flow 
      "Mass flow rate from port_a to port_b (m_flow > 0 is design flow direction)";
    Pressure dp(start=0) "Pressure difference between port_a and port_b";
    Real dz_in=0 "dz=hb-ha difference of height";
    constant Modelica.SIunits.Acceleration g=Modelica.Constants.g_n;
    parameter Medium.AbsolutePressure p_ambient=101325 "ambient pressure";
    Medium.SaturationProperties sat 
      "State vector to compute saturation properties";
    Medium.Density rho_l=Medium.bubbleDensity(sat) "density in liquid phase";
    Medium.Density rho_v=Medium.dewDensity(sat) "density in liquid phase";
    Integer help;
    Boolean liquid( start=true);
    annotation (
      Coordsys(grid=[1, 1], component=[20, 20]),
      Diagram,
      Documentation(info="<html>
<p>
This component transports fluid between its two ports, without
storing mass or energy. Reversal and zero mass flow rate is taken
care of, for details see definition of built-in operator semiLinear().
<p>
When using this partial component, an equation for the momentum
balance has to be added by specifying a relationship
between the pressure drop <tt>dp</tt> and the mass flow rate <tt>m_flow</tt>.
</p>
</html>"));
  equation 
    // Properties in the ports
    port_a.p   = medium_a.p;
    port_a.h   = medium_a.h;
    port_a.Xi = medium_a.Xi;
    port_b.p   = medium_b.p;
    port_b.h   = medium_b.h;
    port_b.Xi = medium_b.Xi;
    
  if m_flow > 0 then
    if pre(liquid) then
      if medium_a.d < 1*rho_v+0*rho_l then
        liquid = false;
        help = 1;
      else
        liquid = true;
        help = 2;
      end if;
    else
      if medium_a.d > 1*rho_l+0*rho_v then
        liquid = true;
        help = 3;
      else
        liquid = false;
        help = 4;
      end if;
    end if;
  else
    if pre(liquid) then
      if medium_b.d < 1*rho_v+0*rho_l then
        liquid = false;
        help = 5;
      else
        liquid = true;
        help = 6;
      end if;
    else
      if medium_b.d > 1*rho_l+0*rho_v then
        liquid = true;
        help = 7;
      else
        liquid = false;
        help = 8;
      end if;
    end if;
  end if;
    
    sat.psat = p_ambient;
    sat.Tsat = Medium.saturationTemperature(p_ambient);
    
    /* Handle reverse and zero flow */
    port_a.H_flow   = semiLinear(port_a.m_flow, port_a.h,  port_b.h);
    port_a.mXi_flow = semiLinear(port_a.m_flow, port_a.Xi, port_b.Xi);
    
    /* Energy, mass and substance mass balance */
    port_a.H_flow + port_b.H_flow = 0;
    port_a.m_flow + port_b.m_flow = 0;
    port_a.mXi_flow + port_b.mXi_flow = zeros(Medium.nXi);
    
    // Design direction of mass flow rate
    m_flow = port_a.m_flow;
    
    // Pressure difference between ports
    
    dp = port_a.p - port_b.p - dz_in*g*(if pre(liquid) then rho_l else rho_v);
    
  end PartialTwoPortTransportWithDz;
end Interfaces;
