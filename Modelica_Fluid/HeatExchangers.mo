within Modelica_Fluid;
package HeatExchangers "Evaporators and condensor components"
  extends Modelica_Fluid.Icons.VariantLibrary;
  model EquilibriumDrumBoiler
    "Simple Evaporator with two states, see Astroem, Bell: Drum-boiler dynamics, Automatica 36, 2000, pp.363-378"
    import Modelica.SIunits.Conversions.*;
    import Modelica.Constants;
    import Modelica_Fluid.Types;
    import Modelica_Fluid.Types.FlowDirection;
    outer Modelica_Fluid.System system "System properties";
    replaceable package Medium = Modelica.Media.Water.StandardWater 
      constrainedby Modelica.Media.Interfaces.PartialTwoPhaseMedium
      "Medium model" 
        annotation (choicesAllMatching=true);
    parameter SI.Mass m_D "mass of surrounding drum metal";
    parameter Medium.SpecificHeatCapacity cp_D
      "specific heat capacity of drum metal";
    parameter SI.Volume V_t "total volume inside drum";
    parameter Types.Init initType=Types.Init.NoInit "Initialization option" 
    annotation(Dialog(tab = "Initialization"));
    parameter Medium.AbsolutePressure p_start=Medium.p_default
      "Start value of pressure" 
    annotation(Dialog(tab = "Initialization"));
    parameter SI.Volume V_l_start=V_t/2
      "Start value of liquid volumeStart value of volume" 
    annotation(Dialog(tab = "Initialization"));

    parameter FlowDirection flowDirection=system.flowDirection
      "Unidirectional (port_a -> port_b) or bidirectional flow component" 
     annotation(Dialog(tab="Advanced"));

    Modelica_Fluid.Interfaces.FluidPort_a feedwater(redeclare package Medium = 
          Medium, m_flow(min=if allowFlowReversal then -Constants.inf else 0)) 
    annotation (Placement(transformation(extent={{-110,-10},{-90,10}}, rotation=
             0)));
    Modelica_Fluid.Interfaces.FluidPort_b steam(redeclare package Medium = Medium,
        m_flow(max=if allowFlowReversal then +Constants.inf else 0)) 
    annotation (Placement(transformation(extent={{110,-10},{90,10}}, rotation=0)));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort 
    annotation (Placement(transformation(extent={{-10,-110},{10,-90}}, rotation=
             0)));
    Modelica.Blocks.Interfaces.RealOutput V "liquid volume" 
    annotation (Placement(transformation(
          origin={100,110},
          extent={{-10,-10},{10,10}},
          rotation=90)));

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
    Medium.SpecificEnthalpy h_W=inStream(feedwater.h_outflow)
      "Feed water enthalpy (specific enthalpy close to feedwater port when mass flows in to the boiler)";
    Medium.SpecificEnthalpy h_S=inStream(steam.h_outflow)
      "steam enthalpy (specific enthalpy close to steam port when mass flows in to the boiler)";
    SI.MassFlowRate qm_W=feedwater.m_flow "feed water mass flow rate";
    SI.MassFlowRate qm_S=steam.m_flow "steam mass flow rate";
  /*outer Modelica_Fluid.Components.FluidOptions fluidOptions 
    "Global default options";*/
  protected
    parameter Boolean allowFlowReversal=flowDirection == FlowDirection.Bidirectional
      "= false, if flow only from port_a to port_b, otherwise reversing flow allowed"
     annotation(Evaluate=true, Hide=true);
  equation
  // balance equations
    m = rho_v*V_v + rho_l*V_l + m_D "Total mass";
    U = rho_v*V_v*h_v + rho_l*V_l*h_l - p*V_t + m_D*cp_D*T_D "Total energy";
    der(m) = qm_W + qm_S "Mass balance";
    der(U) = q_F
              + feedwater.m_flow*actualStream(feedwater.h_outflow)
              + steam.m_flow*actualStream(steam.h_outflow) "Energy balance";
    V_t = V_l + V_v;

  // Properties of saturated liquid and steam
    sat.psat = p;
    sat.Tsat = T;
    sat.Tsat = Medium.saturationTemperature(p);

  // ideal heat transfer between metal and water
    T_D = T;

  // boundary conditions at the ports
    feedwater.p = p;
    feedwater.h_outflow = h_l;
  // feedwater.H_flow = semiLinear(feedwater.m_flow, feedwater.h, h_l);
    steam.p = p;
    steam.h_outflow = h_v;
  //steam.H_flow = semiLinear(steam.m_flow, steam.h, h_v);

  // liquid volume
    V = V_l;

  // Check that two-phase equilibrium is actually possible
    assert(p < Medium.fluidConstants[1].criticalPressure - 10000,
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
      Diagram(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics),
      Icon(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics={
          Rectangle(
            extent={{-100,60},{100,-61}},
            lineColor={0,0,0},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={192,192,192}),
          Rectangle(
            extent={{-100,34},{100,-36}},
            lineColor={0,0,0},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={0,127,255}),
          Ellipse(
            extent={{18,0},{48,-29}},
            lineColor={0,0,255},
            pattern=LinePattern.None,
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-1,29},{29,0}},
            lineColor={0,0,255},
            pattern=LinePattern.None,
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{48,34},{78,5}},
            lineColor={0,0,255},
            pattern=LinePattern.None,
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-31,1},{-1,-28}},
            lineColor={0,0,255},
            pattern=LinePattern.None,
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{47,14},{77,-15}},
            lineColor={0,0,255},
            pattern=LinePattern.None,
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-72,25},{-42,-4}},
            lineColor={0,0,255},
            pattern=LinePattern.None,
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{71,0},{101,-29}},
            lineColor={0,0,255},
            pattern=LinePattern.None,
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{74,14},{104,-15}},
            lineColor={0,0,255},
            pattern=LinePattern.None,
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{71,29},{101,0}},
            lineColor={0,0,255},
            pattern=LinePattern.None,
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-150,70},{150,110}},
            lineColor={0,0,0},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={0,127,255},
            textString="%name"),
          Line(points={{0,-60},{0,-100}}, color={191,0,0}),
          Line(points={{100,100},{100,60}}, color={0,0,127})}),
      Documentation(revisions="<html>
<ul>
<li><i>2 Nov 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
     Initialization options fixed</li>
<li><i>6 Sep 2005</i><br>
    Model by Ruediger Franke modified after the 45th Design Meeting</li>
</ul>
</html>",   info="<html>
Model of a simple evaporator with two states. The model assumes two-phase equilibrium inside the component; saturated steam goes out of the steam outlet.
<p>
References: Astroem, Bell: Drum-boiler dynamics, Automatica 36, 2000, pp.363-378
</html>"));
  equation

  end EquilibriumDrumBoiler;
  annotation (Documentation(info="<html>
 
</html>"));
  model BasicHX "Simple heat exchanger model"

    //General
    parameter Integer n(min=1) = 1 "Spatial segmentation";
    replaceable package Medium_1 = Modelica.Media.Water.StandardWater constrainedby
      Modelica.Media.Interfaces.PartialMedium "Fluid 1" 
                                                      annotation(choicesAllMatching, Dialog(tab="General",group="Fluid 1"));
    replaceable package Medium_2 = Modelica.Media.Water.StandardWater constrainedby
      Modelica.Media.Interfaces.PartialMedium "Fluid 2" 
                                                      annotation(choicesAllMatching,Dialog(tab="General", group="Fluid 2"));
    parameter SI.Area Ah_1 "Heat transfer area" annotation(Dialog(tab="General",group="Fluid 1"));
    parameter SI.Area Ah_2 "Heat transfer area" annotation(Dialog(tab="General",group="Fluid 2"));
    parameter SI.Area Ac_1 "Cross sectional area" annotation(Dialog(tab="General",group="Fluid 1"));
    parameter SI.Area Ac_2 "Cross sectional area" annotation(Dialog(tab="General",group="Fluid 2"));
    parameter SI.Length P_1 "Flow channel perimeter" annotation(Dialog(tab="General",group="Fluid 1"));
    parameter SI.Length P_2 "Flow channel perimeter" annotation(Dialog(tab="General",group="Fluid 2"));
    parameter SI.Length length(min=0) "Length of flow path for both fluids";
    parameter SI.Length s_wall(min=0) "Wall thickness";
    //Wall
    parameter SI.Density d_wall "Density of wall material" annotation(Dialog(tab="General", group="Solid material properties"));
    parameter SI.SpecificHeatCapacity c_wall
      "Specific heat capacity of wall material" annotation(Dialog(tab="General", group="Solid material properties"));
    final parameter SI.Mass m_wall=sum(wall.m) "Wall mass";
    parameter SI.ThermalConductivity k_wall
      "Thermal conductivity of wall material" 
      annotation (Dialog(group="Solid material properties"));

    //Initialization pipe 1
    parameter Types.Init initType=Types.Init.InitialValues
      "Initialization option" 
      annotation(Evaluate=true, Dialog(tab = "Initialization"));
    parameter SI.Temperature Twall_start "Start value of wall temperature" 
                                                                          annotation(Dialog(tab="Initialization", group="Wall"));
    parameter SI.Temperature dT "Start value for port_b.T - port_a.T" 
      annotation (Dialog(tab="Initialization", group="Wall"));
    parameter Boolean use_T_start=true "Use T_start if true, otherwise h_start"
      annotation(Evaluate=true, Dialog(tab = "Initialization"));
    parameter Medium_1.AbsolutePressure p_a_start1=Medium_1.p_default
      "Start value of pressure" 
      annotation(Dialog(tab = "Initialization", group = "Fluid 1"));
    parameter Medium_1.AbsolutePressure p_b_start1=Medium_1.p_default
      "Start value of pressure" 
      annotation(Dialog(tab = "Initialization", group = "Fluid 1"));
    parameter Medium_1.Temperature T_start_1=if use_T_start then Medium_1.
        T_default else Medium_1.temperature_phX(
          (p_a_start1 + p_b_start1)/2,
          h_start_1,
          X_start_1) "Start value of temperature" 
      annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Fluid 1", enable = use_T_start));
    parameter Medium_1.SpecificEnthalpy h_start_1=if use_T_start then Medium_1.specificEnthalpy_pTX(
          (p_a_start1 + p_b_start1)/2,
          T_start_1,
          X_start_1[1:Medium_1.nXi]) else Medium_1.h_default
      "Start value of specific enthalpy" 
      annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Fluid 1", enable = not use_T_start));
    parameter Medium_1.MassFraction X_start_1[Medium_1.nX]=Medium_1.X_default
      "Start value of mass fractions m_i/m" 
      annotation (Dialog(tab="Initialization", group = "Fluid 1", enable=(Medium_1.nXi > 0)));
    parameter Medium_1.MassFlowRate mflow_start_1
      "Start value of mass flow rate" annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Fluid 1"));
    //Initialization pipe 2

    parameter Medium_2.AbsolutePressure p_a_start2=Medium_2.p_default
      "Start value of pressure" 
      annotation(Dialog(tab = "Initialization", group = "Fluid 2"));
    parameter Medium_2.AbsolutePressure p_b_start2=Medium_2.p_default
      "Start value of pressure" 
      annotation(Dialog(tab = "Initialization", group = "Fluid 2"));
    parameter Medium_2.Temperature T_start_2=if use_T_start then Medium_2.
        T_default else Medium_2.temperature_phX(
          (p_a_start2 + p_b_start2)/2,
          h_start_2,
          X_start_2) "Start value of temperature" 
      annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Fluid 2", enable = use_T_start));
    parameter Medium_2.SpecificEnthalpy h_start_2=if use_T_start then Medium_2.specificEnthalpy_pTX(
          (p_a_start2 + p_b_start2)/2,
          T_start_2,
          X_start_2[1:Medium_2.nXi]) else Medium_2.h_default
      "Start value of specific enthalpy" 
      annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Fluid 2", enable = not use_T_start));
    parameter Medium_2.MassFraction X_start_2[Medium_2.nX]=Medium_2.X_default
      "Start value of mass fractions m_i/m" 
      annotation (Dialog(tab="Initialization", group = "Fluid 2", enable=Medium_2.nXi>0));
    parameter Medium_2.MassFlowRate mflow_start_2
      "Start value of mass flow rate"    annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Fluid 2"));
    //Advanced
    parameter Boolean static=false
      "= true, use quasistatic mass and energy balances" 
                             annotation(Evaluate=true, Dialog(tab="General", group="Model options"));

    //Pressure drop and heat transfer
    replaceable package WallFriction = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent
      constrainedby
      Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.PartialWallFriction
      "Characteristic of wall friction"                                                            annotation(choicesAllMatching, Dialog(tab="General", group="Pressure drop"));
    parameter SI.Length roughness_1=2.5e-5
      "Absolute roughness of pipe (default = smooth steel pipe)" annotation(Dialog(tab="General", group="Fluid 1"));
    parameter SI.Length roughness_2=2.5e-5
      "Absolute roughness of pipe (default = smooth steel pipe)" annotation(Dialog(tab="General", group="Fluid 2"));
    parameter SI.DynamicViscosity eta_nominal_M1=0.01
      "Nominal dynamic viscosity (e.g. eta_liquidWater = 1e-3, eta_air = 1.8e-5)"
                                                                                             annotation(Dialog(tab="General", group="Fluid 1"));
    parameter SI.DynamicViscosity eta_nominal_M2=0.01
      "Nominal dynamic viscosity (e.g. eta_liquidWater = 1e-3, eta_air = 1.8e-5)"
                                                                                         annotation(Dialog(tab="General", group="Fluid 2"));
    parameter Boolean use_eta_nominal=false
      "= true, if eta_nominal is used, otherwise computed from medium" annotation(Evaluate=true, Dialog(tab="General", group="Pressure drop"));
    replaceable model HeatTransfer_1 = 
        Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PipeHT_constAlpha 
      constrainedby
      Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PartialPipeHeatTransfer
      "Heat transfer model" annotation(choicesAllMatching, Dialog(tab="General", group="Fluid 1"));
    replaceable model HeatTransfer_2 = 
        Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PipeHT_constAlpha 
      constrainedby
      Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PartialPipeHeatTransfer
      "Heat transfer model" annotation(choicesAllMatching, Dialog(tab="General", group="Fluid 2"));
    //Display variables
    SI.HeatFlowRate Q_flow_1 "Total heat flow rate of pipe 1";
    SI.HeatFlowRate Q_flow_2 "Total heat flow rate of pipe 2";

    Modelica_Fluid.Pipes.DistributedPipe pipe_1(
      redeclare package Medium = Medium_1,
      isCircular=false,
      diameter=0,
      n=n,
      static=static,
      length=length,
      area_h=Ah_1,
      redeclare HeatTransfer_1 heatTransfer,
      initType=initType,
      use_T_start=use_T_start,
      T_start=T_start_1,
      h_start=h_start_1,
      X_start=X_start_1,
      mflow_start=mflow_start_1,
      perimeter=P_1,
      area=Ac_1,
      redeclare package WallFriction = WallFriction,
      roughness=roughness_1,
      use_eta_nominal=use_eta_nominal,
      eta_nominal=eta_nominal_M1) 
                               annotation (Placement(transformation(extent={{
              -40,-60},{20,0}}, rotation=0)));

    Modelica_Fluid.Pipes.DistributedPipe pipe_2(
      redeclare package Medium = Medium_2,
      n=n,
      static=static,
      length=length,
      isCircular=false,
      diameter=0,
      redeclare HeatTransfer_2 heatTransfer,
      use_T_start=use_T_start,
      T_start=T_start_2,
      h_start=h_start_2,
      X_start=X_start_2,
      initType=initType,
      mflow_start=mflow_start_2,
      perimeter=P_2,
      area=Ac_2,
      area_h=Ah_2,
      p_a_start=p_a_start1,
      p_b_start=p_b_start2,
      redeclare package WallFriction = WallFriction,
      roughness=roughness_2,
      use_eta_nominal=use_eta_nominal,
      eta_nominal=eta_nominal_M2,
      show_Re=false) 
                annotation (Placement(transformation(extent={{-40,88},{20,28}},
            rotation=0)));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{
              -100,-100},{100,100}}),
                        graphics),
                         Icon(coordinateSystem(preserveAspectRatio=false,
            extent={{-100,-100},{100,100}}), graphics={
          Rectangle(
            extent={{-100,-26},{100,-30}},
            lineColor={0,0,0},
            fillColor={95,95,95},
            fillPattern=FillPattern.Forward),
          Rectangle(
            extent={{-100,30},{100,26}},
            lineColor={0,0,0},
            fillColor={95,95,95},
            fillPattern=FillPattern.Forward),
          Rectangle(
            extent={{-100,60},{100,30}},
            lineColor={0,0,0},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={0,63,125}),
          Rectangle(
            extent={{-100,-30},{100,-60}},
            lineColor={0,0,0},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={0,63,125}),
          Rectangle(
            extent={{-100,26},{100,-26}},
            lineColor={0,0,0},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={0,128,255}),
          Text(
            extent={{-100,-60},{100,-100}},
            lineColor={0,0,255},
            textString="%name")}),
      Documentation(info="<html>
Simple model of a heat exchanger consisting of two pipes and one wall in between. For both fluids geometry parameters, such as heat transfer area and cross section as well as heat transfer and pressure drop correlations may be chosen. The flow scheme be cocurrent or counterflow, defined by the respective flow directions of the fluids entering the component.
</html>"));
    Modelica_Fluid.Interfaces.FluidPort_b port_b1(redeclare package Medium = 
          Medium_1) annotation (Placement(transformation(extent={{100,-12},{120,
              8}}, rotation=0)));
    Modelica_Fluid.Interfaces.FluidPort_a port_a1(redeclare package Medium = 
          Medium_1) annotation (Placement(transformation(extent={{-120,-12},{
              -100,8}}, rotation=0)));
    Modelica_Fluid.Interfaces.FluidPort_a port_a2(redeclare package Medium = 
          Medium_2) annotation (Placement(transformation(extent={{-120,36},{
              -100,56}}, rotation=0)));
    Modelica_Fluid.Interfaces.FluidPort_b port_b2(redeclare package Medium = 
          Medium_2) annotation (Placement(transformation(extent={{100,-56},{120,
              -36}}, rotation=0)));

    Modelica_Fluid.Thermal.WallConstProps wall(
      n=n,
      d_wall=d_wall,
      c_wall=c_wall,
      T_start=Twall_start,
      k_wall=k_wall,
      dT=dT,
      s=s_wall,
      area_h=(Ah_1 + Ah_2)/2,
      initType=initType) 
      annotation (Placement(transformation(extent={{-28,-14},{10,44}}, rotation=
             0)));

  equation
    Q_flow_1 = sum(pipe_1.heatTransfer.Q_flow);
    Q_flow_2 = sum(pipe_2.heatTransfer.Q_flow);
    connect(pipe_2.port_b, port_b2) annotation (Line(
        points={{20,58},{60,58},{60,-46},{110,-46}},
        color={0,127,255},
        thickness=0.5));
    connect(pipe_1.port_b, port_b1) annotation (Line(
        points={{20,-30},{42,-30},{42,-2},{110,-2}},
        color={0,127,255},
        thickness=0.5));
    connect(pipe_1.port_a, port_a1) annotation (Line(
        points={{-40,-30},{-75.3,-30},{-75.3,-2},{-110,-2}},
        color={0,127,255},
        thickness=0.5));
    connect(pipe_2.port_a, port_a2) annotation (Line(
        points={{-40,58},{-76,58},{-76,46},{-110,46}},
        color={0,127,255},
        thickness=0.5));
    connect(pipe_2.thermalPort, wall.thermalPort_a) annotation (Line(points={{-10,
            41.8},{-10,29.5},{-9,29.5}},     color={191,0,0}));
    connect(wall.thermalPort_b, pipe_1.thermalPort) annotation (Line(points={{-9,0.5},
            {-9,-7.75},{-10,-7.75},{-10,-13.8}},         color={191,0,0}));
  end BasicHX;
end HeatExchangers;
