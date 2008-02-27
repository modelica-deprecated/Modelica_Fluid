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
  
  model HeatedPipe "Model for heating or cooling a fluid" 
    
    extends Interfaces.PartialComponent;
    extends FluidInterface.PartialPortsAndMediumOnlyBA;
    
    import Modelica_Fluid.Types;
    
    parameter Types.Init.Temp initType=
              Types.Init.NoInit "Initialization option" 
      annotation(Evaluate=true, Dialog(tab = "Initialization"));
    parameter Medium.AbsolutePressure p_start = Medium.p_default 
      "Start value of pressure" 
      annotation(Dialog(tab = "Initialization"));
    parameter Boolean use_T_start = true 
      "= true, use T_start, otherwise h_start" 
      annotation(Dialog(tab = "Initialization"), Evaluate=true);
    parameter Medium.Temperature T_start=
      if use_T_start then Medium.T_default else Medium.temperature_phX(p_start,h_start,X_start) 
      "Start value of temperature" 
      annotation(Dialog(tab = "Initialization", enable = use_T_start));
    parameter Medium.SpecificEnthalpy h_start=
      if use_T_start then Medium.specificEnthalpy_pTX(p_start, T_start, X_start) else Medium.h_default 
      "Start value of specific enthalpy" 
      annotation(Dialog(tab = "Initialization", enable = not use_T_start));
    parameter Medium.MassFraction X_start[Medium.nX] = Medium.X_default 
      "Start value of mass fractions m_i/m" 
      annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
    
    import SI = Modelica.SIunits;
    import Modelica.SIunits.Conversions.*;
    
    parameter SI.Volume V_lumped=1 "Volume";
    
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b heatPort 
       annotation (extent=[-10,-110; 10,-90]);
    
    annotation (Icon(
        Rectangle(extent=[-100,60; 100,-60], style(
            color=0,
            rgbcolor={0,0,0},
            gradient=2,
            fillColor=1,
            rgbfillColor={255,0,0})),
        Text(
          extent=[-150,140; 150,80],
          string="%name",
          style(gradient=2, fillColor=69)),
        Line(points=[0,-60; 0,-90], style(color=1, rgbcolor={255,0,0}))), Diagram);
    Volumes.Volume lumpedVolume(
      initType=initType,
      p_start=p_start,
      use_T_start=use_T_start,
      T_start=T_start,
      h_start=h_start,
      X_start=X_start,
      redeclare package FluidInterface = FluidInterface,
      redeclare package Medium = Medium,
      V=V_lumped,
      n_a=1,
      n_b=1,
      provide_p=provide_p_b,
      provide_T=provide_T_b)                     annotation (extent=[-10,10; 10,
          -10]);
    PressureLosses.WallFriction pipe(
      provide_p_b=false,
      provide_T_b=false,
      redeclare package FluidInterface = FluidInterface,
      redeclare package Medium = Medium,
      redeclare package WallFriction = WallFriction,
      length=length,
      roughness=roughness,
      diameter=diameter,
      provide_p_a=provide_p_a,
      provide_T_a=provide_T_a,
      provide_m_flow_ab=provide_m_flow_a) 
              annotation (extent=[-66,-10; -46,10]);
    
    replaceable package WallFriction = 
        PressureLosses.WallFrictionCorrelations.Laminar 
      extends PressureLosses.WallFrictionCorrelations.PartialWallFriction 
      "Characteristic of wall friction"  annotation(choicesAllMatching=true);
    
    parameter SI.Length length "Length of pipe";
    parameter SI.Length roughness=2.5e-5 
      "Absolute roughness of pipe (default = smooth steel pipe)";
    parameter SI.Diameter diameter "Inner (hydraulic) diameter of pipe";
    
    // sensors
    parameter Boolean provide_p_a=false "Provide pressure at port a?" annotation(Evaluate=true, Dialog(descriptionLabel=true, tab="Sensors"));
    parameter Boolean provide_T_a=false "Provide temperature at port a?" annotation(Evaluate=true, Dialog(descriptionLabel=true, tab="Sensors"));
    parameter Boolean provide_m_flow_a=false 
      "Provide mass flow rate at port a?"                                        annotation(Evaluate=true, Dialog(descriptionLabel=true, tab="Sensors"));
    Modelica.Blocks.Interfaces.RealOutput m_flow_a(redeclare type SignalType = 
          SI.MassFlowRate) if provide_m_flow_a 
      annotation (extent=[-120,70; -100,90], rotation=180);
    Modelica.Blocks.Interfaces.RealOutput p_a(redeclare type SignalType = 
          SI.Pressure) if provide_p_a 
      annotation (extent=[-120,50; -100,70], rotation=180);
    Modelica.Blocks.Interfaces.RealOutput T_a(redeclare type SignalType = 
          SI.Temperature) if provide_T_a 
      annotation (extent=[-120,30; -100,50], rotation=180);
    
    parameter Boolean provide_p_b=false "Provide pressure at port b?" annotation(Evaluate=true, Dialog(descriptionLabel=true, tab="Sensors"));
    parameter Boolean provide_T_b=false "Provide temperature at port b?" annotation(Evaluate=true, Dialog(descriptionLabel=true, tab="Sensors"));
    Modelica.Blocks.Interfaces.RealOutput T_b(redeclare type SignalType = 
          SI.Temperature) if provide_T_b  annotation (extent=[100,30; 120,50]);
    Modelica.Blocks.Interfaces.RealOutput p_b(redeclare type SignalType = 
          SI.Pressure) if provide_p_b  annotation (extent=[100,50; 120,70]);
  equation 
    connect(port_a, pipe.port_a) annotation (points=[-102,0; -66,0], style(
      color=69,
      rgbcolor={0,127,255},
      gradient=3,
      fillColor=1,
      rgbfillColor={255,0,0}));
    connect(lumpedVolume.thermalPort, heatPort) annotation (points=[0,-9.8; 0,
          -100], style(
        color=42,
        rgbcolor={191,0,0},
        smooth=0));
    connect(pipe.port_b, lumpedVolume.port_a[1]) annotation (points=[-46,0; -10,0],
        style(
        color=69,
        rgbcolor={0,127,255},
        smooth=0));
    connect(lumpedVolume.port_b[1], port_b) annotation (points=[10,0; 100,0],
        style(
        color=69,
        rgbcolor={0,127,255},
        smooth=0));
    connect(pipe.m_flow_ab, m_flow_a) annotation (points=[-62,11; -62,80; -110,80],
        style(
        color=74,
        rgbcolor={0,0,127},
        smooth=0));
    connect(pipe.p_a, p_a) annotation (points=[-67,8; -68,8; -68,60; -110,60],
        style(
        color=74,
        rgbcolor={0,0,127},
        smooth=0));
    connect(pipe.T_a, T_a) annotation (points=[-67,5; -76,5; -76,40; -110,40],
        style(
        color=74,
        rgbcolor={0,0,127},
        smooth=0));
    connect(lumpedVolume.T, T_b) annotation (points=[-11,-5; -14,-5; -14,40; 110,
          40], style(
        color=74,
        rgbcolor={0,0,127},
        smooth=0));
    connect(lumpedVolume.p, p_b) annotation (points=[-11,-8; -20,-8; -20,60; 110,
          60], style(
        color=74,
        rgbcolor={0,0,127},
        smooth=0));
  end HeatedPipe;
end HeatTransfer;
