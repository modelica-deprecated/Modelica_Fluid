within FluidSandbox;
package HeatTransfer "Component models depicting heat transfer" 
  extends Icons.VariantLibrary;
  model EvaporatingVessel 
    "Valve for water/steam flows with linear pressure drop" 
    
    extends Interfaces.PartialComponent;
    extends FluidInterface.PartialTwoSidedVolume(h_a=h_l, h_b=h_v);
    
    import Modelica_Fluid.Types;
    
    // Drum parameters
    parameter Modelica.SIunits.Mass m_D "mass of surrounding drum metal";
    parameter Medium.SpecificHeatCapacity cp_D 
      "specific heat capacity of drum metal";
    parameter Modelica.SIunits.Volume V_t "total volume inside drum";
    
    // Drum boiler variables  
    Medium.AbsolutePressure p(start=p_start, stateSelect=StateSelect.prefer) 
      "pressure inside drum boiler";
    Modelica.SIunits.Volume V_l(start=V_l_start, stateSelect=StateSelect.prefer) 
      "volumes of liquid phase";
    Modelica.SIunits.Volume V_v "volume of vapour phase";
    Medium.Temperature T "temperature inside drum boiler";
    Medium.SaturationProperties sat 
      "State vector to compute saturation properties";
    Medium.SpecificEnthalpy h_v=Medium.dewEnthalpy(sat) 
      "specific enthalpy of vapour";
    Medium.SpecificEnthalpy h_l=Medium.bubbleEnthalpy(sat) 
      "specific enthalpy of liquid";
    Medium.Density rho_v=Medium.dewDensity(sat) "density in vapour phase";
    Medium.Density rho_l=Medium.bubbleDensity(sat) "density in liquid phase";
    Modelica.SIunits.Mass m "total mass of drum boiler";
    Modelica.SIunits.Energy U "internal energy";
    Medium.Temperature T_D=heatPort.T "temperature of drum";
    Modelica.SIunits.HeatFlowRate q_F=heatPort.Q_flow 
      "heat flow rate from furnace";
    
    // initialization
    parameter Types.Init.Temp initType=Types.Init.NoInit 
      "Initialization option" 
      annotation(Dialog(tab = "Initialization"));
    parameter Modelica.SIunits.Volume V_l_start=V_t/2 
      "Start value of liquid volumeStart value of volume" 
      annotation(Dialog(tab = "Initialization"));
    parameter Medium.AbsolutePressure p_start=Medium.p_default 
      "Start value of pressure" 
      annotation(Dialog(tab = "Initialization"));
    
    // output connector  
    Modelica.Blocks.Interfaces.RealOutput V(redeclare type SignalType = 
          Modelica.SIunits.Volume) "liquid volume" 
      annotation (extent=[30, 100; 50, 120], rotation=90);
    
  equation 
  // balance equations  
    m = rho_v*V_v + rho_l*V_l + m_D "Total mass";
    U = rho_v*V_v*h_v + rho_l*V_l*h_l - p*V_t + m_D*cp_D*T_D "Total energy";
    der(m) = m_flow_net "Mass balance";
    der(U) = q_F + H_flow_net "Energy balance";
    V_t = V_l + V_v;
    zeros(Medium.nXi) = mXi_flow_net 
      "May not be used for media with nXi<>0, included for equation count reasons";
    
  // Properties of saturated liquid and steam
    sat.psat = p;
    sat.Tsat = Medium.saturationTemperature(p);
    
  // ideal heat transfer between metal and water
    T_D = T;
    
  // liquid volume 
    V = V_l;
    
  // component properties
  // assuming that always steam is leaving the drum boiler
    medium.p = p;
    medium.h = h_v;
    
  // Check that two-phase equilibrium is actually possible
    assert(p < Medium.fluidConstants[1].criticalPressure - 10000,
      "Evaporator model requires subcritical pressure");
    
    assert(Medium.nXi==0, "Usable for two phase media only.");
    
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
    
    annotation (Documentation(info=
            "Taken from PowerFluid library by Ruediger Franke."), Icon(
        Rectangle(extent=[-100,59; 100,-61], style(
            color=0,
            gradient=2,
            fillColor=8)),
        Rectangle(extent=[-100,34; 100,-36], style(
            color=69,
            gradient=2,
            fillColor=69)),
        Ellipse(extent=[-72,25; -42,-4], style(pattern=0, fillColor=7)),
        Ellipse(extent=[-31,1; -1,-28], style(pattern=0, fillColor=7)),
        Ellipse(extent=[-1,29; 29,0], style(pattern=0, fillColor=7)),
        Ellipse(extent=[18,0; 48,-29], style(pattern=0, fillColor=7)),
        Ellipse(extent=[47,14; 77,-15], style(pattern=0, fillColor=7)),
        Ellipse(extent=[48,34; 78,5], style(pattern=0, fillColor=7)),
        Ellipse(extent=[71,29; 101,0], style(pattern=0, fillColor=7)),
        Ellipse(extent=[71,0; 101,-29], style(pattern=0, fillColor=7)),
        Ellipse(extent=[74,14; 104,-15], style(pattern=0, fillColor=7))));
  end EvaporatingVessel;
end HeatTransfer;
