package Components "Basic components for fluid models" 
  extends Modelica.Icons.Library;
  
  annotation (preferedView="info",
              Documentation(info="<html>
<p>
This package will contain basic component models of the fluid library.
It is currently empty as all components are being evaluated together 
with examples first (see sub-package Examples). 
</p>
</html>"));
  
model ShortPipe 
    
    import SI = Modelica.SIunits;
    
  extends FiniteVolume.Components.ShortVolume;
  public 
  parameter Boolean linearPressureDrop=true;
  parameter SI.AbsolutePressure dp_nominal(min=1.e-10) = 1 
      "Nominal pressure drop";
    
  parameter SI.MassFlowRate m_dot_nominal = 1E-3 
      "Nominal mass flow rate at nominal pressure drop";
equation 
    
  /*
  LongPipes.Components.PipeFriction friction[pipe.n](
    each from_dp=false, 
    each dp_nominal=500/pipe.n, 
    each roughness=1, 
    each diameter=30, 
    each length=length/pipe.n);
*/
    
  // Simple linear pressure drop in each segment
  dp_a/dp_nominal = if linearPressureDrop then m_dot_a/m_dot_nominal else abs(
    m_dot_a)*m_dot_a/m_dot_nominal^2;
  dp_b/dp_nominal = if linearPressureDrop then -m_dot_b/m_dot_nominal else abs(
    -m_dot_b)*(-m_dot_b)/m_dot_nominal^2;
    
end ShortPipe;

model LongPipe 
  replaceable package Medium = 
      Modelica_Media.Interfaces.PartialMedium               annotation (
            choicesAllMatching=true);
    
  Modelica_Fluid.Interfaces.FluidPort_a port_a(redeclare model Medium = 
          Medium) 
              annotation (extent=[-120, -10; -100, 10]);
  Modelica_Fluid.Interfaces.FluidPort_b port_b(redeclare model Medium = 
          Medium) 
              annotation (extent=[120, -10; 100, 10]);
    
    import SI = Modelica.SIunits;
    
    import Modelica.SIunits.Conversions.*;
    
  parameter Integer n=10;
    
  parameter SI.Length L=1 "Length of pipe";
  parameter SI.AbsolutePressure dp_nominal(min=1.e-10) = 1 
      "|frictionType = ConstantLaminar or ConstantTurbulent| Nominal pressure drop";
    
  parameter SI.MassFlowRate m_dot_nominal = 1E-3 
      "Nominal mass flow rate at nominal pressure drop";
    
  parameter SI.Area A_a=1;
  parameter SI.Area A_b=A_a;
    
  parameter SI.Length Z_a=0;
  parameter SI.Length Z_b=Z_a;
    
  parameter SI.ThermalConductivity k=0 "Thermal conductivity";
    
  parameter Real viscosityFactor1=0;
  parameter Real viscosityFactor2=1;
    
    parameter Boolean dynamicMomentumBalance=false 
      "If false, der(m_dot) is neglected in momentum balance" 
                                                 annotation(Evaluate=true,
      Dialog(tab="Level of Detail"));
    parameter Boolean includeKineticTerm=false 
      "If false, d*v^2 is neglected in momentum balance" 
                                             annotation(Evaluate=true,
      Dialog(tab="Level of Detail"));
    parameter Boolean includeThermalConductance=false 
      "If false, thermal conductance is neglected"  annotation(Evaluate=true, Dialog(tab=
          "Level of Detail"));
  //  parameter Boolean constantTemperature=false;
  //  parameter SI.Temperature T0;
    
  FiniteVolume.Components.SimpleShortPipe shortPipe[n](
      redeclare package Medium = Medium, 
      each L=L/n, 
      each dp_nominal=dp_nominal/n, 
      each A_a=A_a, 
      each Z_a=Z_a, 
      each m_dot_nominal=m_dot_nominal, 
      each dynamicMomentumBalance=dynamicMomentumBalance, 
      each k=k, 
      each includeKineticTerm=includeKineticTerm, 
      each includeThermalConductance=includeThermalConductance, 
      each viscosityFactor1=viscosityFactor1, 
      each viscosityFactor2=viscosityFactor2) 
                      annotation (extent=[-10, -10; 10, 10]);
equation 
    
  connect(port_a, shortPipe[1].port_a);
    
  connect(port_b, shortPipe[n].port_b);
    
  for i in 1:n - 1 loop
      
    connect(shortPipe[i].port_b, shortPipe[i + 1].port_a);
      
  end for;
    
annotation (Icon(
    Rectangle(extent=[-100, 60; 100, -60], style(color=0, fillColor=8)),
    Rectangle(extent=[-100, 34; 100, -36], style(
        color=69,
        gradient=2,
        fillColor=69)),
    Text(
      extent=[-120, 130; 116, 64],
      string="%name",
      style(gradient=2, fillColor=69))));
    
end LongPipe;
end Components;
