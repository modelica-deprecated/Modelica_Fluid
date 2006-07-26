package Boilers 
  extends Modelica_Fluid.Icons.VariantLibrary;
  model EquilibriumDrumBoiler 
    "Simple Evaporator with two states, see Astroem, Bell: Drum-boiler dynamics, Automatica 36, 2000, pp.363-378" 
    import Modelica.SIunits.Conversions.*;
    import Modelica.Constants;
    import Modelica_Fluid.Types;
    import Modelica_Fluid.Types.FlowDirection;
    import Modelica_Fluid.Types.FlowDirectionWithGlobalDefault;
  replaceable package Medium = 
      Modelica.Media.Interfaces.PartialTwoPhaseMedium 
    extends Modelica.Media.Interfaces.PartialTwoPhaseMedium "Medium model" 
      annotation (choicesAllMatching=true);
  parameter SI.Mass m_D "mass of surrounding drum metal";
  parameter Medium.SpecificHeatCapacity cp_D 
      "specific heat capacity of drum metal";
  parameter SI.Volume V_t "total volume inside drum";
  parameter Types.Init.Temp initType=
            Types.Init.NoInit "Initialization option" 
    annotation(Dialog(tab = "Initialization"));
  parameter Medium.AbsolutePressure p_start = Medium.p_default 
      "Start value of pressure" 
    annotation(Dialog(tab = "Initialization"));
  parameter SI.Volume V_l_start = V_t/2 
      "Start value of liquid volumeStart value of volume" 
    annotation(Dialog(tab = "Initialization"));
    
  parameter FlowDirectionWithGlobalDefault.Temp flowDirection=
            FlowDirectionWithGlobalDefault.UseGlobalFluidOption 
      "Unidirectional (port_a -> port_b) or bidirectional flow component" 
     annotation(Dialog(tab="Advanced"));
    
  Interfaces.FluidPort_a feedwater(redeclare package Medium = Medium,
                     m_flow(min=if allowFlowReversal then -Constants.inf else 0)) 
    annotation (extent=[-110,-10; -90,10]);
  Interfaces.FluidPort_b steam(redeclare package Medium = Medium,
                     m_flow(max=if allowFlowReversal then +Constants.inf else 0)) 
    annotation (extent=[110,-10; 90,10]);
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort 
    annotation (extent=[-10,-110; 10,-90]);
  Modelica.Blocks.Interfaces.RealOutput V(
    redeclare type SignalType = SI.Volume) "liquid volume" 
    annotation (extent=[30, 100; 50, 120], rotation=90);
    
  Medium.SaturationProperties sat 
      "State vector to compute saturation properties";
  Medium.AbsolutePressure p(start=p_start, stateSelect=StateSelect.prefer) 
      "pressure inside drum boiler";
  Medium.Temperature T "temperature inside drum boiler";
  SI.Volume V_v "volume of vapour phase";
  SI.Volume V_l(start=V_l_start, stateSelect=StateSelect.prefer) 
      "volumes of liquid phase";
  Medium.SpecificEnthalpy h_v=Medium.dewEnthalpy(sat) 
      "specific enthalpy of vapour";
  Medium.SpecificEnthalpy h_l=Medium.bubbleEnthalpy(sat) 
      "specific enthalpy of liquid";
  Medium.Density rho_v=Medium.dewDensity(sat) "density in vapour phase";
  Medium.Density rho_l=Medium.bubbleDensity(sat) "density in liquid phase";
  SI.Mass m "total mass of drum boiler";
  SI.Energy U "internal energy";
  Medium.Temperature T_D=heatPort.T "temperature of drum";
  SI.HeatFlowRate q_F=heatPort.Q_flow "heat flow rate from furnace";
  Medium.SpecificEnthalpy h_W=feedwater.h "feed water enthalpy";
  Medium.SpecificEnthalpy h_S=steam.h "steam enthalpy";
  SI.MassFlowRate qm_W=feedwater.m_flow "feed water mass flow rate";
  SI.MassFlowRate qm_S=steam.m_flow "steam mass flow rate";
  /*outer Modelica_Fluid.Components.FluidOptions fluidOptions 
    "Global default options";*/
  protected 
  parameter Boolean allowFlowReversal=
     flowDirection == FlowDirection.Bidirectional 
      "= false, if flow only from port_a to port_b, otherwise reversing flow allowed"
     annotation(Evaluate=true, Hide=true);
  equation 
  // balance equations  
  m = rho_v*V_v + rho_l*V_l + m_D "Total mass";
  U = rho_v*V_v*h_v + rho_l*V_l*h_l - p*V_t + m_D*cp_D*T_D "Total energy";
  der(m) = qm_W + qm_S "Mass balance";
  der(U) = q_F + qm_W*h_W + qm_S*h_S "Energy balance";
  V_t = V_l + V_v;
    
  // Properties of saturated liquid and steam
  sat.psat = p;
  sat.Tsat = T;
  sat.Tsat = Medium.saturationTemperature(p);
    
  // ideal heat transfer between metal and water
  T_D = T;
    
  // boundary conditions at the ports
  feedwater.p = p;
  feedwater.H_flow = semiLinear(feedwater.m_flow, feedwater.h, h_l);
  steam.p = p;
  steam.H_flow = semiLinear(steam.m_flow, steam.h, h_v);
    
  // liquid volume 
  V = V_l;
    
  // Check that two-phase equilibrium is actually possible
  assert(p<Medium.fluidConstants[1].criticalPressure-10000,
         "Evaporator model requires subcritical pressure");
  initial equation 
  // Initial conditions
  if initType == Types.Init.NoInit then
    // no initial equations
  elseif initType == Types.Init.InitialValues then
   p = p_start;
   V_l = V_l_start;
  elseif initType == Types.Init.SteadyState then
    der(p) = 0;
    der(V_l) = 0;
  elseif initType == Types.Init.SteadyStateHydraulic then
    der(p) = 0;
    V_l = V_l_start;
  else
    assert(false, "Unsupported initialization option");
  end if;
    
  annotation (
    Coordsys(grid=[1, 1], component=[20, 20]),
    Diagram,
    Icon(
      Rectangle(extent=[-100, 59; 100, -61], style(
          color=0,
          gradient=2,
          fillColor=8)),
      Rectangle(extent=[-100, 34; 100, -36], style(
          color=69,
          gradient=2,
          fillColor=69)),
      Ellipse(extent=[18, 0; 48, -29], style(pattern=0, fillColor=7)),
      Ellipse(extent=[-1, 29; 29, 0], style(pattern=0, fillColor=7)),
      Ellipse(extent=[48, 34; 78, 5], style(pattern=0, fillColor=7)),
      Ellipse(extent=[-31, 1; -1, -28], style(pattern=0, fillColor=7)),
      Ellipse(extent=[47, 14; 77, -15], style(pattern=0, fillColor=7)),
      Ellipse(extent=[-72, 25; -42, -4], style(pattern=0, fillColor=7)),
      Ellipse(extent=[71, 0; 101, -29], style(pattern=0, fillColor=7)),
      Ellipse(extent=[74, 14; 104, -15], style(pattern=0, fillColor=7)),
      Ellipse(extent=[71, 29; 101, 0], style(pattern=0, fillColor=7)),
      Text(
        extent=[-139,111; 144,57],
        string="%name",
        style(gradient=2, fillColor=69)),
      Line(points=[0, -60; 0, -100], style(color=42)),
      Line(points=[40, 99; 40, 60])),
    Documentation(revisions="<html>
<ul>
<li><i>2 Nov 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
     Initialization options fixed</li>
<li><i>6 Sep 2005</i><br>
    Model by Ruediger Franke modified after the 45th Design Meeting</li>
</ul>
</html>",
        info="<html>
Model of a simple evaporator with two states. The model assumes two-phase equilibrium inside the component; saturated steam goes out of the steam outlet.
<p>
References: Astroem, Bell: Drum-boiler dynamics, Automatica 36, 2000, pp.363-378
</html>"));
  equation 
    
  end EquilibriumDrumBoiler;
end Boilers;
