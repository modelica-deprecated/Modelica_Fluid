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
    import Modelica_Fluid.WorkInProgress.Utilities.regRoot2;
    import Modelica_Fluid.WorkInProgress.Utilities.regSquare2;
    
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
    SI.AbsolutePressure dp_turbulent 
      "The turbulent region is: |dp| >= dp_turbulent (from_dp=true)";
    SI.MassFlowRate m_flow_turbulent 
      "The turbulent region is: |m_flow| >= m_flow_turbulent (from_dp=false)";
    
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
        m_flow_turbulent = 0 
        "m_flow_turbulent = regRoot2(dp_turbulent, dp_turbulent, d_a/k1, d_b/k2), but not needed";
        yd0 = if use_Re and withLaminar then (d_a+d_b)/(k0*(eta_a+eta_b)) else 0;
        m_flow = regRoot2(dp, dp_turbulent, d_a/k1, d_b/k2, withLaminar, yd0);
     else
        dp_turbulent = 0 
        "regSquare2(m_flow_turbulent, m_flow_turbulent, k1/d_a, k2/d_b), but not needed";
        m_flow_turbulent = if use_Re then (pi/4)*lossFactors.D_Re*eta*lossFactors.Re_turbulent else m_flow_small;
        yd0 = if use_Re and withLaminar then k0*(eta_a+eta_b)/(d_a+d_b) else 0;
        dp = regSquare2(m_flow, m_flow_turbulent, k1/d_a, k2/d_b, withLaminar, yd0);
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
       eta_a = if use_Re then Medium.dynamicViscosity(medium_a) else 0,
       eta_b = if use_Re then Medium.dynamicViscosity(medium_b) else 0);
    annotation (
      Diagram,
      Icon,
      Documentation(info="<html>
</html>"));
  end PressureLossWithoutIcon;
end Interfaces;
