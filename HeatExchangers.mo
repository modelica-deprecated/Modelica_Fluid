within Modelica_Fluid;
package HeatExchangers "Evaporators and condensor components"
  extends Modelica_Fluid.Icons.VariantLibrary;
  model EquilibriumDrumBoiler
    "Simple Evaporator with two states, see Astroem, Bell: Drum-boiler dynamics, Automatica 36, 2000, pp.363-378"
    extends Modelica_Fluid.Interfaces.PartialTwoPort(
      final port_a_exposesState=true,
      final port_b_exposesState=true,
      redeclare replaceable package Medium = 
          Modelica.Media.Water.StandardWater 
          constrainedby Modelica.Media.Interfaces.PartialTwoPhaseMedium);
    import Modelica.SIunits.Conversions.*;
    import Modelica.Constants;
    import Modelica_Fluid.Types;

    parameter SI.Mass m_D "mass of surrounding drum metal";
    parameter Medium.SpecificHeatCapacity cp_D
      "specific heat capacity of drum metal";
    parameter SI.Volume V_t "total volume inside drum";
    parameter Types.Init initType=Types.Init.NoInit "Initialization option" 
    annotation(Dialog(tab = "Initialization"));
    parameter Medium.AbsolutePressure p_start=system.p_start
      "Start value of pressure" 
    annotation(Dialog(tab = "Initialization"));
    parameter SI.Volume V_l_start=V_t/2
      "Start value of liquid volumeStart value of volume" 
    annotation(Dialog(tab = "Initialization"));

    parameter Boolean allowFlowReversal = system.allowFlowReversal
      "allow flow reversal, false restricts to design direction (port_a -> port_b)"
      annotation(Dialog(tab="Assumptions"), Evaluate=true);

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
    Medium.SpecificEnthalpy h_W=inStream(port_a.h_outflow)
      "Feed water enthalpy (specific enthalpy close to feedwater port when mass flows in to the boiler)";
    Medium.SpecificEnthalpy h_S=inStream(port_b.h_outflow)
      "steam enthalpy (specific enthalpy close to steam port when mass flows in to the boiler)";
    SI.MassFlowRate qm_W=port_a.m_flow "feed water mass flow rate";
    SI.MassFlowRate qm_S=port_b.m_flow "steam mass flow rate";
  /*outer Modelica_Fluid.Components.FluidOptions fluidOptions 
    "Global default options";*/
  equation
  // balance equations
    m = rho_v*V_v + rho_l*V_l + m_D "Total mass";
    U = rho_v*V_v*h_v + rho_l*V_l*h_l - p*V_t + m_D*cp_D*T_D "Total energy";
    der(m) = qm_W + qm_S "Mass balance";
    der(U) = q_F
              + port_a.m_flow*actualStream(port_a.h_outflow)
              + port_b.m_flow*actualStream(port_b.h_outflow) "Energy balance";
    V_t = V_l + V_v;

  // Properties of saturated liquid and steam
    sat.psat = p;
    sat.Tsat = T;
    sat.Tsat = Medium.saturationTemperature(p);

  // ideal heat transfer between metal and water
    T_D = T;

  // boundary conditions at the ports
    port_a.p = p;
    port_a.h_outflow = h_l;
    port_b.p = p;
    port_b.h_outflow = h_v;

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
            lineColor={0,0,255},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={0,127,255},
            textString="%name"),
          Line(points={{0,-61},{0,-100}}, color={191,0,0}),
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
    outer Modelica_Fluid.System system "System properties";
    //General
    parameter Integer nNodes(min=1) = 1 "Spatial segmentation";
    replaceable package Medium_1 = Modelica.Media.Water.StandardWater constrainedby
      Modelica.Media.Interfaces.PartialMedium "Fluid 1" 
                                                      annotation(choicesAllMatching, Dialog(tab="General",group="Fluid 1"));
    replaceable package Medium_2 = Modelica.Media.Water.StandardWater constrainedby
      Modelica.Media.Interfaces.PartialMedium "Fluid 2" 
                                                      annotation(choicesAllMatching,Dialog(tab="General", group="Fluid 2"));
    parameter SI.Area crossArea_1 "Cross sectional area" annotation(Dialog(tab="General",group="Fluid 1"));
    parameter SI.Area crossArea_2 "Cross sectional area" annotation(Dialog(tab="General",group="Fluid 2"));
    parameter SI.Length perimeter_1 "Flow channel perimeter" annotation(Dialog(tab="General",group="Fluid 1"));
    parameter SI.Length perimeter_2 "Flow channel perimeter" annotation(Dialog(tab="General",group="Fluid 2"));
    parameter SI.Length length(min=0) "Length of flow path for both fluids";
    parameter SI.Length s_wall(min=0) "Wall thickness";
    // Heat transfer
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

    parameter SI.Area area_h_1 "Heat transfer area" annotation(Dialog(tab="General",group="Fluid 1"));
    parameter SI.Area area_h_2 "Heat transfer area" annotation(Dialog(tab="General",group="Fluid 2"));
   //Wall
    parameter SI.Density d_wall "Density of wall material" annotation(Dialog(tab="General", group="Solid material properties"));
    parameter SI.SpecificHeatCapacity c_wall
      "Specific heat capacity of wall material" annotation(Dialog(tab="General", group="Solid material properties"));
    final parameter SI.Area area_h=(area_h_1 + area_h_2)/2 "Heat transfer area";
    final parameter SI.Mass m_wall=d_wall*area_h*s_wall "Wall mass";
    parameter SI.ThermalConductivity k_wall
      "Thermal conductivity of wall material" 
      annotation (Dialog(group="Solid material properties"));

    // Assumptions
    parameter Boolean allowFlowReversal = system.allowFlowReversal
      "allow flow reversal, false restricts to design direction (port_a -> port_b)"
      annotation(Dialog(tab="Assumptions"), Evaluate=true);
    parameter Modelica_Fluid.Types.Dynamics dynamicsType=system.dynamicsType
      "Dynamics option" 
      annotation(Evaluate=true, Dialog(tab = "Assumptions"));

    //Initialization pipe 1
    parameter Types.Init initType=Types.Init.InitialValues
      "Initialization option" 
      annotation(Evaluate=true, Dialog(tab = "Initialization"));
    parameter SI.Temperature Twall_start "Start value of wall temperature" 
                                                                          annotation(Dialog(tab="Initialization", group="Wall"));
    parameter SI.Temperature dT "Start value for pipe_1.T - pipe_2.T" 
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
    parameter Medium_1.MassFlowRate m_flow_start_1 = system.m_flow_start
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
    parameter Medium_2.MassFlowRate m_flow_start_2 = system.m_flow_start
      "Start value of mass flow rate"    annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Fluid 2"));

    //Pressure drop and heat transfer
    replaceable package WallFriction_1 = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent
      constrainedby
      Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.PartialWallFriction
      "Characteristic of wall friction"                                                            annotation(choicesAllMatching, Dialog(tab="General", group="Fluid 1"));
    replaceable package WallFriction_2 = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent
      constrainedby
      Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.PartialWallFriction
      "Characteristic of wall friction"                                                            annotation(choicesAllMatching, Dialog(tab="General", group="Fluid 2"));
    parameter SI.Length roughness_1=2.5e-5
      "Absolute roughness of pipe (default = smooth steel pipe)" annotation(Dialog(tab="General", group="Fluid 1"));
    parameter SI.Length roughness_2=2.5e-5
      "Absolute roughness of pipe (default = smooth steel pipe)" annotation(Dialog(tab="General", group="Fluid 2"));
    parameter SI.DynamicViscosity eta_nominal_1=Medium_1.dynamicViscosity(Medium_1.setState_pTX(Medium_1.p_default, Medium_1.T_default, Medium_1.X_default))
      "Nominal dynamic viscosity (e.g. eta_liquidWater = 1e-3, eta_air = 1.8e-5)"
                                                                                             annotation(Dialog(tab="Advanced", group="Fluid 1", enable=use_eta_nominal));
    parameter SI.DynamicViscosity eta_nominal_2=Medium_2.dynamicViscosity(Medium_2.setState_pTX(Medium_2.p_default, Medium_2.T_default, Medium_2.X_default))
      "Nominal dynamic viscosity (e.g. eta_liquidWater = 1e-3, eta_air = 1.8e-5)"
                                                                                         annotation(Dialog(tab="Advanced", group="Fluid 2", enable=use_eta_nominal));
    parameter SI.Density d_nominal_1 = Medium_1.density_pTX(Medium_1.p_default, Medium_1.T_default, Medium_1.X_default)
      "Nominal density (e.g. d_liquidWater = 995, d_air = 1.2)" 
       annotation(Dialog(tab="Advanced", group="Fluid 1", enable=use_eta_nominal));
    parameter SI.Density d_nominal_2 = Medium_1.density_pTX(Medium_2.p_default, Medium_2.T_default, Medium_2.X_default)
      "Nominal density (e.g. d_liquidWater = 995, d_air = 1.2)" 
       annotation(Dialog(tab="Advanced", group="Fluid 2", enable=use_eta_nominal));
    parameter Boolean use_eta_nominal=false
      "= true, if eta_ and d_nominal are used, otherwise computed from media" annotation(Evaluate=true, Dialog(tab="Advanced", group="Pressure loss"));
    //Display variables
    SI.HeatFlowRate Q_flow_1 "Total heat flow rate of pipe 1";
    SI.HeatFlowRate Q_flow_2 "Total heat flow rate of pipe 2";

    Modelica_Fluid.Pipes.DistributedPipe_Old pipe_1(
      redeclare package Medium = Medium_1,
      isCircular=false,
      diameter=0,
      nNodes=nNodes,
      allowFlowReversal=allowFlowReversal,
      dynamicsType=dynamicsType,
      length=length,
      redeclare model HeatTransfer = HeatTransfer_1(area=area_h_1),
      initType=initType,
      use_T_start=use_T_start,
      T_start=T_start_1,
      h_start=h_start_1,
      X_start=X_start_1,
      m_flow_start=m_flow_start_1,
      perimeter=perimeter_1,
      crossArea=crossArea_1,
      redeclare package WallFriction = WallFriction_1,
      roughness=roughness_1,
      use_eta_nominal=use_eta_nominal,
      eta_nominal=eta_nominal_1,
      d_nominal=d_nominal_1)   annotation (Placement(transformation(extent={{-40,-80},
              {20,-20}},        rotation=0)));

    Modelica_Fluid.Pipes.DistributedPipe_Old pipe_2(
      redeclare package Medium = Medium_2,
      nNodes=nNodes,
      allowFlowReversal=allowFlowReversal,
      dynamicsType=dynamicsType,
      length=length,
      isCircular=false,
      diameter=0,
      redeclare model HeatTransfer = HeatTransfer_2(area=area_h_2),
      use_T_start=use_T_start,
      T_start=T_start_2,
      h_start=h_start_2,
      X_start=X_start_2,
      initType=initType,
      m_flow_start=m_flow_start_2,
      perimeter=perimeter_2,
      crossArea=crossArea_2,
      p_a_start=p_a_start1,
      p_b_start=p_b_start2,
      redeclare package WallFriction = WallFriction_2,
      roughness=roughness_2,
      use_eta_nominal=use_eta_nominal,
      eta_nominal=eta_nominal_2,
      d_nominal=d_nominal_2,
      show_Re=false) 
                annotation (Placement(transformation(extent={{20,88},{-40,28}},
            rotation=0)));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
              -100},{100,100}},
          grid={1,1}),  graphics),
                         Icon(coordinateSystem(preserveAspectRatio=false,
            extent={{-100,-100},{100,100}},
          grid={1,1}), graphics={
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
            extent={{-150,110},{150,70}},
            lineColor={0,0,255},
            textString="%name"),
          Line(
            points={{30,-85},{-60,-85}},
            color={0,128,255},
            smooth=Smooth.None),
          Polygon(
            points={{20,-70},{60,-85},{20,-100},{20,-70}},
            lineColor={0,128,255},
            smooth=Smooth.None,
            fillColor={0,128,255},
            fillPattern=FillPattern.Solid),
          Line(
            points={{30,77},{-60,77}},
            color={0,128,255},
            smooth=Smooth.None),
          Polygon(
            points={{-50,92},{-90,77},{-50,62},{-50,92}},
            lineColor={0,128,255},
            smooth=Smooth.None,
            fillColor={0,128,255},
            fillPattern=FillPattern.Solid)}),
      Documentation(info="<html>
Simple model of a heat exchanger consisting of two pipes and one wall in between. 
For both fluids geometry parameters, such as heat transfer area and cross section as well as heat transfer and pressure drop correlations may be chosen. 
The flow scheme may be concurrent or counterflow, defined by the respective flow directions of the fluids entering the component.
The design flow direction with positive m_flow variables is counterflow.
</html>"));
    Modelica_Fluid.Interfaces.FluidPort_b port_b1(redeclare package Medium = 
          Medium_1) annotation (Placement(transformation(extent={{100,-12},{120,
              8}}, rotation=0)));
    Modelica_Fluid.Interfaces.FluidPort_a port_a1(redeclare package Medium = 
          Medium_1) annotation (Placement(transformation(extent={{-120,-12},{
              -100,8}}, rotation=0)));
    Modelica_Fluid.Interfaces.FluidPort_b port_b2(redeclare package Medium = 
          Medium_2) annotation (Placement(transformation(extent={{-120,36},{
              -100,56}}, rotation=0)));
    Modelica_Fluid.Interfaces.FluidPort_a port_a2(redeclare package Medium = 
          Medium_2) annotation (Placement(transformation(extent={{100,-56},{120,
              -36}}, rotation=0)));

    Modelica.Thermal.HeatTransfer.Components.ThermalConductor[nNodes] wall1(G=2*
          k_wall/s_wall*area_h/nNodes*ones(nNodes), dT(each start=0.5*dT))           annotation (
        Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-10,-10})));
    Modelica.Thermal.HeatTransfer.Components.ThermalConductor[nNodes] wall2(G=2*
          k_wall/s_wall*area_h/nNodes*ones(nNodes), dT(each start=0.5*dT)) 
                                       annotation (
        Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-10,20})));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor[nNodes] wall(
      C=c_wall*m_wall/nNodes*ones(nNodes),
      T(start=Twall_start*ones(nNodes), each fixed=(initType <> Types.Init.SteadyState)),
      der_T(start=zeros(nNodes), each fixed=(initType == Types.Init.SteadyState))) if dynamicsType <> Types.Dynamics.SteadyState 
      annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-40,5})));

  equation
    Q_flow_1 = sum(pipe_1.heatTransfer.Q_flow);
    Q_flow_2 = sum(pipe_2.heatTransfer.Q_flow);
    connect(pipe_2.port_b, port_b2) annotation (Line(
        points={{-40,58},{-76,58},{-76,46},{-110,46}},
        color={0,127,255},
        thickness=0.5));
    connect(pipe_1.port_b, port_b1) annotation (Line(
        points={{20,-50},{42,-50},{42,-2},{110,-2}},
        color={0,127,255},
        thickness=0.5));
    connect(pipe_1.port_a, port_a1) annotation (Line(
        points={{-40,-50},{-75.3,-50},{-75.3,-2},{-110,-2}},
        color={0,127,255},
        thickness=0.5));
    connect(pipe_2.port_a, port_a2) annotation (Line(
        points={{20,58},{65,58},{65,-46},{110,-46}},
        color={0,127,255},
        thickness=0.5));
    connect(pipe_1.heatPorts, wall1.port_a) annotation (Line(
        points={{-9.7,-34.4},{-9.7,-27.2},{-10,-27.2},{-10,-20}},
        color={191,0,0},
        smooth=Smooth.None));
    connect(wall1.port_b, wall2.port_a) annotation (Line(
        points={{-10,0},{-10,10}},
        color={191,0,0},
        smooth=Smooth.None));
    connect(wall1.port_b, wall.port) annotation (Line(
        points={{-10,0},{-10,5},{-30,5}},
        color={191,0,0},
        smooth=Smooth.None));
    connect(wall2[1:nNodes].port_b, pipe_2.heatPorts[nNodes:-1:1]) annotation (Line(
        points={{-10,30},{-10,36.2},{-10,42.4},{-10.3,42.4}},
        color={191,0,0},
        smooth=Smooth.None));
  end BasicHX;
end HeatExchangers;
