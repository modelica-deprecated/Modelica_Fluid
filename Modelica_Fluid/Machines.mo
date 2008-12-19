within Modelica_Fluid;
package Machines
  "Devices for converting between energy held in a fluid and mechanical energy"
  extends Modelica_Fluid.Icons.VariantLibrary;
  model SweptVolume
    "varying cylindric volume depending on the postition of the piston"
    extends Vessels.BaseClasses.PartialLumpedVolumePorts(Qs_flow=0);
    parameter SI.Area pistonCrossArea "cross sectional area of pistion";
    parameter SI.Volume clearance "remaining volume at zero piston stroke";

    annotation (Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
              -100},{100,100}}),
                        graphics),
                         Icon(coordinateSystem(preserveAspectRatio=true,
            extent={{-100,-100},{100,100}}), graphics={
          Rectangle(
            extent={{-44,36},{44,-90}},
            lineColor={0,0,255},
            pattern=LinePattern.None,
            lineThickness=1,
            fillColor={170,213,255},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{-44,60},{-40,60},{-40,-32},{-44,-32},{-44,60}},
            lineColor={95,95,95},
            smooth=Smooth.None,
            fillColor={135,135,135},
            fillPattern=FillPattern.Backward),
          Polygon(
            points={{40,60},{44,60},{44,-34},{40,-34},{40,60}},
            lineColor={95,95,95},
            smooth=Smooth.None,
            fillColor={135,135,135},
            fillPattern=FillPattern.Backward),
          Rectangle(
            extent={{-40,40},{40,30}},
            lineColor={95,95,95},
            fillColor={135,135,135},
            fillPattern=FillPattern.Forward),
          Rectangle(
            extent={{-6,92},{6,40}},
            lineColor={95,95,95},
            fillColor={135,135,135},
            fillPattern=FillPattern.Forward),
          Polygon(
            points={{-40,-86},{40,-86},{40,70},{44,70},{44,-90},{-44,-90},{-44,
                70},{-40,70},{-40,-86}},
            lineColor={95,95,95},
            smooth=Smooth.None,
            fillColor={135,135,135},
            fillPattern=FillPattern.Backward)}),
      Documentation(info="<html>
<p> Mixing volume with varying size. The size of the volume is given by:</p>
<ul>
  <li>cross sectional piston area</li>
  <li>piston stroke given by the flange position s</li>
  <li>clearance (volume at flang position = 0)</li>
</ul> 
 
<p> The flange position has to be equal or greater than zero. Otherwise the simulation stops. The force of the flange results from the pressure difference between medium and ambient pressure and the cross sectional piston area. For using the component, a top level instance of the ambient model with the inner attribute is needed.</p>
<p> The pressure at both fluid ports equals the medium pressure in the volume. No suction nor discharge valve is included in the model.</p>
<p>The thermal port is directly connected to the medium. The temperature of the thermal port equals the medium temperature. The heat capacity of the cylinder and the piston are not includes in the model.</p>
</html>",
        revisions="<html>
<ul>
<li><i>29 Oct 2007</i>
    by Carsten Heinrich:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
    Modelica.Mechanics.Translational.Interfaces.Flange_b flange
      "translation flange for piston" annotation (Placement(transformation(
            extent={{-10,90},{10,110}},   rotation=0)));

  equation
    assert(flange.s >= 0, "Piston stroke (given by flange.s) must not be smaller than zero!");

    // volume size
    fluidVolume = clearance + flange.s * pistonCrossArea;

    flange.f = (medium.p - system.p_ambient) * pistonCrossArea;

    // energy balances
    Ws_flow = medium.p * pistonCrossArea * (-der(flange.s));

    // definition of ports pressure
    ports_p_static = medium.p;
  end SweptVolume;

  model Pump "Centrifugal pump with mechanical connector for the shaft"
    extends Modelica_Fluid.Machines.BaseClasses.PartialPump;
    SI.Angle phi "Shaft angle";
    SI.AngularVelocity omega "Shaft angular velocity";
    Modelica.Mechanics.Rotational.Interfaces.Flange_a shaft 
    annotation (Placement(transformation(extent={{80,4},{110,32}}, rotation=0),
          iconTransformation(extent={{80,-34},{110,-6}})));
  equation
    phi = shaft.phi;
    omega = der(phi);
    N = Modelica.SIunits.Conversions.to_rpm(omega);
    W_single = omega*shaft.tau;
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
              100}}), graphics={Rectangle(
            extent={{60,-12},{84,-26}},
            lineColor={0,0,0},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={95,95,95})}),
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}}),
            graphics),
    Documentation(info="<HTML>
<p>This model describes a centrifugal pump (or a group of <tt>nParallel</tt> pumps) with a mechanical rotational connector for the shaft, to be used when the pump drive has to be modelled explicitly. In the case of <tt>nParallel</tt> pumps, the mechanical connector is relative to a single pump.
<p>The model extends <tt>PartialPump</tt>
 </HTML>",
       revisions="<html>
<ul>
<li><i>31 Oct 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
  end Pump;

  model ControlledPump
    "Centrifugal pump with ideally controlled mass flow rate"
    import Modelica.SIunits.Conversions.NonSIunits.AngularVelocity_rpm;
    extends Modelica_Fluid.Machines.BaseClasses.PartialPump(
      N_nominal=1500,
      N(start=N_nominal),
      redeclare replaceable function flowCharacteristic = 
          Modelica_Fluid.Machines.BaseClasses.PumpCharacteristics.quadraticFlow
          ( q_nominal={0, q_op, 1.5*q_op},
            head_nominal={2*head_op, head_op, 0}));

    // nominal values
    parameter Medium.AbsolutePressure p_a_nominal
      "Nominal inlet pressure for predefined pump characteristics";
    parameter Medium.AbsolutePressure p_b_nominal
      "Nominal outlet pressure, fixed if not control_m_flow and not use_p_set";
    parameter Medium.MassFlowRate m_flow_nominal
      "Nominal mass flow rate, fixed if control_m_flow and not use_m_flow_set";

    // what to control
    parameter Boolean control_m_flow = true
      "= false to control outlet pressure port_b.p instead of m_flow" 
      annotation(Evaluate = true);
    parameter Boolean use_m_flow_set = false
      "= true to use input signal m_flow_set instead of m_flow_nominal" 
      annotation (Dialog(enable = control_m_flow));
    parameter Boolean use_p_set = false
      "= true to use input signal p_set instead of p_b_nominal" 
      annotation (Dialog(enable = not control_m_flow));

    // exemplary characteristics
    final parameter SI.VolumeFlowRate q_op = m_flow_nominal/d_nominal
      "operational volume flow rate according to nominal values";
    final parameter SI.Height head_op = (p_b_nominal-p_a_nominal)/(d_nominal*g)
      "operational pump head according to nominal values";

    Modelica.Blocks.Interfaces.RealInput m_flow_set if use_m_flow_set
      "Prescribed mass flow rate" 
      annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=-90,
          origin={-50,82})));
    Modelica.Blocks.Interfaces.RealInput p_set if use_p_set
      "Prescribed outlet pressure" 
      annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=-90,
          origin={50,82})));

  protected
    Modelica.Blocks.Interfaces.RealInput m_flow_set_internal
      "Needed to connect to conditional connector";
    Modelica.Blocks.Interfaces.RealInput p_set_internal
      "Needed to connect to conditional connector";
  equation
    // Ideal control
    if control_m_flow then
      m_flow = m_flow_set_internal;
    else
      dp = p_set_internal - port_a.p;
    end if;

    // Internal connector value when use_m_flow_set = false
    if not use_m_flow_set then
      m_flow_set_internal = m_flow_nominal;
    end if;
    if not use_p_set then
      p_set_internal = p_b_nominal;
    end if;
    connect(m_flow_set, m_flow_set_internal);
    connect(p_set, p_set_internal);

    annotation (defaultComponentName="pump",
      Icon(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},{100,
              100}}), graphics={Text(
            visible=use_p_set,
            extent={{82,108},{176,92}},
            textString="p_set"), Text(
            visible=use_m_flow_set,
            extent={{-20,108},{170,92}},
            textString="m_flow_set")}),
      Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},{
              100,100}}), graphics),
      Documentation(info="<HTML>
<p>
This model describes a centrifugal pump (or a group of <tt>nParallel</tt> pumps) 
with ideally controlled mass flow rate or pressure.
</p>
<p>
Nominal values are used to predefine an exemplary pump characteristics and to define the operation of the pump. 
The input connectors <tt>m_flow_set</tt> or <tt>p_set</tt> can optionally be enabled to provide time varying set points. 
</p>
<p>
Use this model if the pump characteristics is of secondary interest. 
The actual characteristics can be configured later on for the appropriate rotational speed N. 
Then the model can be replaced with a Pump with rotational shaft or with a PrescribedPump.
</p>
</HTML>",
        revisions="<html>
<ul>
<li><i>15 Dec 2008</i>
    by Ruediger Franke</a>:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
  end ControlledPump;

  model PrescribedPump "Centrifugal pump with ideally controlled speed"
    extends Modelica_Fluid.Machines.BaseClasses.PartialPump;
    parameter Boolean use_N_input = false
      "Get the rotational speed from the input connector";
    parameter Modelica.SIunits.Conversions.NonSIunits.AngularVelocity_rpm
      N_const =                                                                     N_nominal
      "Constant rotational speed" annotation(Dialog(enable = not use_N_input));
    Modelica.Blocks.Interfaces.RealInput N_in(unit="1/min") if use_N_input
      "Prescribed rotational speed" 
      annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=-90,
          origin={0,100}), iconTransformation(
          extent={{-20,-20},{20,20}},
          rotation=-90,
          origin={0,100})));
    annotation (defaultComponentName="pump",
      Icon(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},{100,
              100}}), graphics={Text(
            visible=use_N_input,
            extent={{14,98},{178,82}},
            textString="N_in [rpm]")}),
      Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},{
              100,100}}), graphics),
      Documentation(info="<HTML>
<p>This model describes a centrifugal pump (or a group of <tt>nParallel</tt> pumps) with prescribed speed, either fixed or provided by an external signal.
<p>The model extends <tt>PartialPump</tt>
<p>If the <tt>N_in</tt> input connector is wired, it provides rotational speed of the pumps (rpm); otherwise, a constant rotational speed equal to <tt>n_const</tt> (which can be different from <tt>N_nominal</tt>) is assumed.</p>
</HTML>",
        revisions="<html>
<ul>
<li><i>31 Oct 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));

  protected
    Modelica.Blocks.Interfaces.RealInput N_in_internal(unit="1/min")
      "Needed to connect to conditional connector";
  equation
    // Connect statement active only if use_p_in = true
    connect(N_in, N_in_internal);
    // Internal connector value when use_p_in = false
    if not use_N_input then
      N_in_internal = N_const;
    end if;
    // Set N with a lower limit to avoid singularities at zero speed
    N = max(N_in_internal,1e-3) "Rotational speed";

  end PrescribedPump;

  model PrescribedPumpNPSH
    "Centrifugal pump with ideally controlled speed and NPSHa computation"
    extends Modelica_Fluid.Machines.PrescribedPump(
                 redeclare replaceable package Medium = 
      Modelica.Media.Water.WaterIF97_ph constrainedby
        Modelica.Media.Interfaces.PartialTwoPhaseMedium);
    SI.Length NPSHa "Net Positive Suction Head available";
    Medium.AbsolutePressure pv "Saturation pressure of inlet liquid";

  equation
    // NPSHa computation
    pv = Medium.saturationPressure(Tin);
    NPSHa = (port_a.p-pv)/(d*Modelica.Constants.g_n);

    // Check for cavitation
    assert(port_a.p >= pv, "Cavitation occurs at the inlet");
    annotation (defaultComponentName="pump",
  Documentation(info="<html>Same as the PrescribedPump model, with added computation of Net Positive Suction Head available. Requires a two-phase medium model.
</html>", revisions="<html>
<ul>
<li><i>30 Jul 2007</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Model added to the Fluid library</li>
</ul>
 
</html>"), Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}), graphics));
  end PrescribedPumpNPSH;

  package BaseClasses "Base classes for Turbomachinery components"
    extends Modelica_Fluid.Icons.BaseClassLibrary;

  partial model PartialPump "Base model for centrifugal pumps"
      import Modelica.SIunits.Conversions.NonSIunits.*;
      import Modelica.Constants;
    extends Modelica_Fluid.Interfaces.PartialTwoPort(
      port_a_exposesState = (V > 0),
      port_b_exposesState = (V > 0),
      port_a(
        p(start=p_a_start),
        m_flow(start = m_flow_start,
               min = if allowFlowReversal and not checkValve then -Constants.inf else 0)),
      port_b(
        p(start=p_b_start),
        m_flow(start = -m_flow_start,
               max = if allowFlowReversal and not checkValve then +Constants.inf else 0)));

    parameter Integer nParallel(min=1) = 1 "Number of pumps in parallel" 
      annotation(Dialog(group="Characteristics"));
    replaceable function flowCharacteristic = 
        PumpCharacteristics.baseFlow
        "Head vs. q_flow characteristic at nominal speed and density" 
      annotation(Dialog(group="Characteristics"), choicesAllMatching=true);
    parameter AngularVelocity_rpm N_nominal
        "Nominal rotational speed for flow characteristic" 
      annotation(Dialog(group="Characteristics"));
    parameter Medium.Density d_nominal = Medium.density_pTX(Medium.p_default, Medium.T_default, Medium.X_default)
        "Nominal fluid density" 
      annotation(Dialog(group="Characteristics"));
    parameter Boolean usePowerCharacteristic = false
        "Use powerCharacteristic (vs. efficiencyCharacteristic)" 
       annotation(Evaluate=true,Dialog(group="Characteristics"));
    replaceable function efficiencyCharacteristic = 
      PumpCharacteristics.constantEfficiency(eta_nominal = 0.8) constrainedby
        PumpCharacteristics.baseEfficiency
        "Efficiency vs. q_flow at nominal speed and density" 
      annotation(Dialog(group="Characteristics",enable = not usePowerCharacteristic),
                 choicesAllMatching=true);
    replaceable function powerCharacteristic = 
          PumpCharacteristics.quadraticPower (
         q_nominal={0,0,0},W_nominal={0,0,0})
        "Power consumption vs. q_flow at nominal speed and density" 
      annotation(Dialog(group="Characteristics", enable = usePowerCharacteristic),
                 choicesAllMatching=true);

    // Assumptions
    parameter Types.Dynamics energyDynamics=system.energyDynamics
        "Formulation of energy balance" 
      annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));
    parameter Boolean checkValve=false "= true to prevent reverse flow" 
      annotation(Dialog(tab="Assumptions"));
    parameter SI.Volume V = 0
        "Fluid volume inside the pump; consider if checkValve" 
      annotation(Dialog(tab="Assumptions"));

    // Initialization
    parameter Medium.AbsolutePressure p_a_start=system.p_start
        "Guess value for inlet pressure" 
      annotation(Dialog(tab="Initialization"));
    parameter Medium.AbsolutePressure p_b_start=2*p_a_start
        "Guess value for outlet pressure" 
      annotation(Dialog(tab="Initialization"));
    parameter Boolean use_T_start = true
        "Use T_start if true, otherwise h_start" 
      annotation(Dialog(tab = "Initialization"), Evaluate = true);
    parameter Medium.Temperature T_start=
      if use_T_start then system.T_start else Medium.temperature_phX(p_a_start,h_start,X_start)
        "Guess value for temperature" 
      annotation(Dialog(tab = "Initialization", enable = use_T_start));
    parameter Medium.SpecificEnthalpy h_start=
      if use_T_start then Medium.specificEnthalpy_pTX(p_a_start, T_start, X_start) else Medium.h_default
        "Guess value for specific enthalpy" 
      annotation(Dialog(tab = "Initialization", enable = not use_T_start));
    parameter Medium.MassFraction X_start[Medium.nX] = Medium.X_default
        "Guess value for mass fractions m_i/m" 
      annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
    parameter SI.MassFlowRate m_flow_start = 1
        "Guess value for mass flow rate (total)" 
      annotation(Dialog(tab="Initialization"));
    final parameter SI.Acceleration g=system.g;

    SI.Pressure dp = port_b.p - port_a.p "Pressure increase";
    SI.Height head = dp/(d*g) "Pump head";
    Medium.Density d "Liquid density at the inlet port_a";
    Medium.SpecificEnthalpy h(start=h_start)
        "Enthalpy of the liquid stored in the pump if M>0";
    Medium.Temperature Tin "Liquid inlet temperature";
    SI.MassFlowRate m_flow = port_a.m_flow "Mass flow rate (total)";
    SI.MassFlowRate m_flow_single = m_flow/nParallel
        "Mass flow rate (single pump)";
    SI.VolumeFlowRate q_flow = m_flow/d "Volume flow rate (total)";
    SI.VolumeFlowRate q_flow_single = q_flow/nParallel
        "Volume flow rate (single pump)";
    AngularVelocity_rpm N "Shaft rotational speed";
    SI.Power W_single "Power Consumption (single pump)";
    SI.Power W_tot = W_single*nParallel "Power Consumption (total)";
    constant SI.Power W_eps=1e-8
        "Small coefficient to avoid numerical singularities in efficiency computations";
    Real eta "Global Efficiency";
    Real s(start = m_flow_start)
        "Curvilinear abscissa for the flow curve in parametric form (either mass flow rate or head)";
    Medium.ThermodynamicState port_a_state_inflow
        "Medium state close to inlet for inflowing mass flow";
    protected
    constant SI.Height unitHead = 1;
    constant SI.MassFlowRate unitMassFlowRate = 1;
    input SI.HeatFlowRate Qs_flow = 0
        "Heat flow, e.g. to consider a housing in an extending class";
  equation
    // Flow equations
    if noEvent(s > 0 or (not checkValve)) then
      // Flow characteristics when check valve is open or with no check valve
      head = (N/N_nominal)^2*flowCharacteristic(q_flow_single*(N_nominal/N));
      q_flow_single = s*unitMassFlowRate/d;
    else
      // Flow characteristics when check valve is closed
      head = (N/N_nominal)^2*flowCharacteristic(0) - s*unitHead;
      q_flow_single = 0;
    end if;

    // Power consumption
    if usePowerCharacteristic then
      W_single = (N/N_nominal)^3*(d/d_nominal)*powerCharacteristic(q_flow_single*(N_nominal/N))
          "Power consumption (single pump)";
      eta = (dp*q_flow_single)/(W_single + W_eps) "Hydraulic efficiency";
    else
      eta = efficiencyCharacteristic(q_flow_single*(N_nominal/N));
      W_single = dp*q_flow_single/eta;
    end if;

    // Medium states close to the ports when mass flows in to the respective port
    // The inlet inflow state is used also in case of flow reversal, to avoid
    // discontinuities.
    port_a_state_inflow = Medium.setState_phX(port_a.p, inStream(port_a.h_outflow), inStream(port_a.Xi_outflow));
    // port_b_state_inflow = Medium.setState_phX(port_b.p, inStream(port_b.h_outflow), inStream(port_b.Xi_outflow));

    // Inflow density and temperature at the inlet port
    d = Medium.density(port_a_state_inflow);
    Tin = Medium.temperature(port_a_state_inflow);

    // Mass balances
    port_a.m_flow + port_b.m_flow = 0 "Mass balance";
    port_a.Xi_outflow  = inStream(port_b.Xi_outflow);
    port_b.Xi_outflow = inStream(port_a.Xi_outflow);
    port_a.C_outflow = inStream(port_b.C_outflow);
    port_b.C_outflow = inStream(port_a.C_outflow);

    // Energy balances
    if energyDynamics <> Types.Dynamics.SteadyState and V > 0 then
       // Dynamic energy balance
       // mass variations and p/d are neglected
       nParallel*V*d_nominal*der(h) = port_a.m_flow*actualStream(port_a.h_outflow) +
                     port_b.m_flow*actualStream(port_b.h_outflow) +
                     W_tot + Qs_flow;
       port_b.h_outflow = h;
       port_a.h_outflow = h;
    else
      /* In the following two equations the extra term W_tot/m_flow is
       present where a 0/0 term appears if m_flow = 0. In order to avoid
       numerical problems, this term is analytically simplified:
         W_tot/m_flow = nParallel*W_single/(d*q_flow_single*nParallel)
                      = nParallel*dp*q_flow_single/(eta*d*q_flow_single*nParallel)
                      = dp/(eta*d)
    */
      port_b.h_outflow  = inStream(port_a.h_outflow) + dp/(d*eta);
      port_a.h_outflow  = inStream(port_b.h_outflow) + dp/(d*eta);
      h = Medium.h_default
          "Unused (set to an arbitrary value within the medium region)";
      assert(abs(Qs_flow) < Modelica.Constants.small, "Specify V > 0 for using Qs_flow.");
    end if;

  initial equation
    if V > 0 then
      if energyDynamics == Types.Dynamics.FixedInitial then
        h = h_start;
      elseif energyDynamics == Types.Dynamics.SteadyStateInitial then
        der(h) = 0;
      end if;
    end if;

    annotation (
      Icon(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},{100,
                100}}), graphics={
            Polygon(
              points={{-48,-60},{-72,-100},{72,-100},{48,-60},{-48,-60}},
              lineColor={0,0,255},
              pattern=LinePattern.None,
              fillColor={0,0,191},
              fillPattern=FillPattern.Solid),
            Ellipse(
              extent={{-80,80},{80,-80}},
              lineColor={0,0,0},
              fillPattern=FillPattern.Sphere),
            Polygon(
              points={{-28,30},{-28,-30},{50,-2},{-28,30}},
              lineColor={0,0,0},
              pattern=LinePattern.None,
              fillPattern=FillPattern.HorizontalCylinder,
              fillColor={255,255,255}),
            Text(
              extent={{-148,-110},{152,-150}},
              textString="%name",
              lineColor={0,0,255}),
            Line(
              points={{-80,0},{-100,0}},
              color={0,128,255},
              smooth=Smooth.None),
            Line(
              points={{100,0},{80,0}},
              color={0,128,255},
              smooth=Smooth.None)}),
      Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},
                {100,100}}),
              graphics),
      Documentation(info="<HTML>
<p>This is the base model for the <tt>Pump</tt> and <tt>
PumpMech</tt> pump models.
<p>The model describes a centrifugal pump, or a group of <tt>nParallel</tt> identical pumps. The pump model is based on the theory of kinematic similarity: the pump characteristics are given for nominal operating conditions (rotational speed and fluid density), and then adapted to actual operating condition, according to the similarity equations. 
<p><b>Modelling options</b></p>
<p> The nominal hydraulic characteristic (head vs. volume flow rate) is given by the the replaceable function <tt>flowCharacteristic</tt>. 
<p> The pump energy balance can be specified in two alternative ways:
<ul>
<li><tt>usePowerCharacteristic = false</tt> (default option): the replaceable function <tt>efficiencyCharacteristic</tt> (efficiency vs. volume flow rate in nominal conditions) is used to determine the efficiency, and then the power consumption. The default is a constant efficiency of 0.8.
<li><tt>usePowerCharacteristic = true</tt>: the replaceable function <tt>powerCharacteristic</tt> (power consumption vs. volume flow rate in nominal conditions) is used to determine the power consumption, and then the efficiency.
</ul>
<p>
Several functions are provided in the package <tt>PumpCharacteristics</tt> to specify the characteristics as a function of some operating points at nominal conditions.
<p>Depending on the value of the <tt>checkValve</tt> parameter, the model either supports reverse flow conditions, or includes a built-in check valve to avoid flow reversal.
<p>It is possible to take into account the heat capacity of the fluid inside the pump by specifying its volume <tt>V</tt>; 
this is necessary to avoid singularities in the computation of the outlet enthalpy in case of zero flow rate. 
If zero flow rate conditions are always avoided, this dynamic effect can be neglected by leaving the default value <tt>V = 0</tt>, thus avoiding a fast state variable in the model.
</HTML>",
        revisions="<html>
<ul>
<li><i>Dec 2008</i>
    by R6uuml;diger Franke:<br>
       Changed fluid mass M to volume V(*d_nominal) and added Qs_flow</li>
<li><i>31 Oct 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
  equation

  end PartialPump;

  package PumpCharacteristics "Functions for pump characteristics"
      import NonSI = Modelica.SIunits.Conversions.NonSIunits;

    partial function baseFlow "Base class for pump flow characteristics"
      extends Modelica.Icons.Function;
      input SI.VolumeFlowRate q_flow "Volumetric flow rate";
      output SI.Height head "Pump head";
    end baseFlow;

    partial function basePower
        "Base class for pump power consumption characteristics"
      extends Modelica.Icons.Function;
      input SI.VolumeFlowRate q_flow "Volumetric flow rate";
      output SI.Power consumption "Power consumption";
    end basePower;

    partial function baseEfficiency "Base class for efficiency characteristics"
      extends Modelica.Icons.Function;
      input SI.VolumeFlowRate q_flow "Volumetric flow rate";
      output Real eta "Efficiency";
    end baseEfficiency;

    function linearFlow "Linear flow characteristic"
      extends baseFlow;
      input SI.VolumeFlowRate q_nominal[2]
          "Volume flow rate for two operating points (single pump)" 
                                                                  annotation(Dialog);
      input SI.Height head_nominal[2] "Pump head for two operating points" annotation(Dialog);
      /* Linear system to determine the coefficients:
  head_nominal[1] = c[1] + q_nominal[1]*c[2];
  head_nominal[2] = c[1] + q_nominal[2]*c[2];
  */
      protected
      Real c[2] = Modelica.Math.Matrices.solve([ones(2),q_nominal],head_nominal)
          "Coefficients of linear head curve";
    algorithm
      // Flow equation: head = q*c[1] + c[2];
      head := c[1] + q_flow*c[2];
    end linearFlow;

    function quadraticFlow "Quadratic flow characteristic"
      extends baseFlow;
      input SI.VolumeFlowRate q_nominal[3]
          "Volume flow rate for three operating points (single pump)" 
                                                                    annotation(Dialog);
      input SI.Height head_nominal[3] "Pump head for three operating points" annotation(Dialog);
      protected
      Real q_nominal2[3] = {q_nominal[1]^2,q_nominal[2]^2, q_nominal[3]^2}
          "Squared nominal flow rates";
      /* Linear system to determine the coefficients:
  head_nominal[1] = c[1] + q_nominal[1]*c[2] + q_nominal[1]^2*c[3];
  head_nominal[2] = c[1] + q_nominal[2]*c[2] + q_nominal[2]^2*c[3];
  head_nominal[3] = c[1] + q_nominal[3]*c[2] + q_nominal[3]^2*c[3];
  */
      Real c[3] = Modelica.Math.Matrices.solve([ones(3), q_nominal, q_nominal2],head_nominal)
          "Coefficients of quadratic head curve";
    algorithm
      // Flow equation: head  = c[1] + q_flow*c[2] + q_flow^2*c[3];
      head := c[1] + q_flow*c[2] + q_flow^2*c[3];
    end quadraticFlow;

    function polynomialFlow "Polynomial flow characteristic"
      extends baseFlow;
      input SI.VolumeFlowRate q_nominal[:]
          "Volume flow rate for N operating points (single pump)" 
                                                                annotation(Dialog);
      input SI.Height head_nominal[:] "Pump head for N operating points" annotation(Dialog);
      protected
      Integer N = size(q_nominal,1) "Number of nominal operating points";
      Real q_nominal_pow[N,N] = {{q_nominal[j]^(i-1) for j in 1:N} for i in 1:N}
          "Rows: different operating points; columns: increasing powers";
      /* Linear system to determine the coefficients (example N=3):
  head_nominal[1] = c[1] + q_nominal[1]*c[2] + q_nominal[1]^2*c[3];
  head_nominal[2] = c[1] + q_nominal[2]*c[2] + q_nominal[2]^2*c[3];
  head_nominal[3] = c[1] + q_nominal[3]*c[2] + q_nominal[3]^2*c[3];
  */
      Real c[N] = Modelica.Math.Matrices.solve(q_nominal_pow,head_nominal)
          "Coefficients of polynomial head curve";
    algorithm
      // Flow equation (example N=3): head  = c[1] + q_flow*c[2] + q_flow^2*c[3];
      // Note: the implementation is numerically efficient only for low values of Na
      head := sum(q_flow^(i-1)*c[i] for i in 1:N);
    end polynomialFlow;

    function constantEfficiency "Constant efficiency characteristic"
       extends baseEfficiency;
       input Real eta_nominal "Nominal efficiency" annotation(Dialog);
    algorithm
      eta := eta_nominal;
    end constantEfficiency;

    function linearPower "Linear power consumption characteristic"
      extends basePower;
      input SI.VolumeFlowRate q_nominal[2]
          "Volume flow rate for two operating points (single pump)" annotation(Dialog);
      input SI.Power W_nominal[2] "Power consumption for two operating points"   annotation(Dialog);
      /* Linear system to determine the coefficients:
  W_nominal[1] = c[1] + q_nominal[1]*c[2];
  W_nominal[2] = c[1] + q_nominal[2]*c[2];
  */
      protected
      Real c[2] = Modelica.Math.Matrices.solve([ones(3),q_nominal],W_nominal)
          "Coefficients of linear power consumption curve";
    algorithm
      consumption := c[1] + q_flow*c[2];
    end linearPower;

    function quadraticPower "Quadratic power consumption characteristic"
      extends basePower;
      input SI.VolumeFlowRate q_nominal[3]
          "Volume flow rate for three operating points (single pump)" 
                                                                    annotation(Dialog);
      input SI.Power W_nominal[3]
          "Power consumption for three operating points"                         annotation(Dialog);
      protected
      Real q_nominal2[3] = {q_nominal[1]^2,q_nominal[2]^2, q_nominal[3]^2}
          "Squared nominal flow rates";
      /* Linear system to determine the coefficients:
  W_nominal[1] = c[1] + q_nominal[1]*c[2] + q_nominal[1]^2*c[3];
  W_nominal[2] = c[1] + q_nominal[2]*c[2] + q_nominal[2]^2*c[3];
  W_nominal[3] = c[1] + q_nominal[3]*c[2] + q_nominal[3]^2*c[3];
  */
      Real c[3] = Modelica.Math.Matrices.solve([ones(3),q_nominal,q_nominal2],W_nominal)
          "Coefficients of quadratic power consumption curve";
    algorithm
      consumption := c[1] + q_flow*c[2] + q_flow^2*c[3];
    end quadraticPower;

  end PumpCharacteristics;
  end BaseClasses;
  annotation (Documentation(info="<html>
 
</html>"));
end Machines;
