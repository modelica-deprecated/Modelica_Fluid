within Modelica_Fluid;
package Pumps "Pump components"
  extends Modelica_Fluid.Icons.VariantLibrary;
  model Pump "Centrifugal pump with ideally controlled speed"
    extends Modelica_Fluid.Pumps.BaseClasses.PartialPump;
    parameter Boolean use_N_input = false
      "Get the rotational speed from the input connector";
    parameter Modelica.SIunits.Conversions.NonSIunits.AngularVelocity_rpm
      N_const =                                                                     N_nom
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
    annotation (
      Icon(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},{100,
              100}}), graphics={Text(
            visible=use_N_input,
            extent={{14,98},{178,82}},
            textString="N_in [rpm]")}),
      Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},{
              100,100}}), graphics),
      Documentation(info="<HTML>
<p>This model describes a centrifugal pump (or a group of <tt>Np</tt> pumps in parallel) with controlled speed, either fixed or provided by an external signal.
<p>The model extends <tt>PartialPump</tt>
<p>If the <tt>N_in</tt> input connector is wired, it provides rotational speed of the pumps (rpm); otherwise, a constant rotational speed equal to <tt>n_const</tt> (which can be different from <tt>N_nom</tt>) is assumed.</p>
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
    // Connect statement active only if usePressureInput = true
    connect(N_in, N_in_internal);
    // Internal connector value when usePressureInput = false
    if not use_N_input then
      N_in_internal = N_const;
    end if;
    // Set N with a lower limit to avoid singularities at zero speed
    N = max(N_in_internal,1e-3) "Rotational speed";

  end Pump;

  model PumpShaft "Centrifugal pump with mechanical connector for the shaft"
    extends Modelica_Fluid.Pumps.BaseClasses.PartialPump;
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
<p>This model describes a centrifugal pump (or a group of <tt>Np</tt> pumps in parallel) with a mechanical rotational connector for the shaft, to be used when the pump drive has to be modelled explicitly. In the case of <tt>Np</tt> pumps in parallel, the mechanical connector is relative to a single pump.
<p>The model extends <tt>PartialPump</tt>
 </HTML>",
       revisions="<html>
<ul>
<li><i>31 Oct 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
  end PumpShaft;

  model PumpNPSH
    "Centrifugal pump with ideally controlled speed and NPSHa computation"
    extends Pump(redeclare replaceable package Medium = 
      Modelica.Media.Water.WaterIF97_ph constrainedby
        Modelica.Media.Interfaces.PartialTwoPhaseMedium);
    SI.Length NPSHa "Net Positive Suction Head available";
    Medium.AbsolutePressure pv "Saturation pressure of inlet liquid";

  equation
    // NPSHa computation
    pv = Medium.saturationPressure(Tin);
    NPSHa = (inlet.p-pv)/(d*Modelica.Constants.g_n);

    // Check for cavitation
    assert(inlet.p >= pv, "Cavitation occurs at the inlet");
    annotation (Documentation(info="<html>Same as the Pump model, with added computation of Net Positive Suction Head available. Requires a two-phase medium model.
</html>", revisions="<html>
<ul>
<li><i>30 Jul 2007</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Model added to the Fluid library</li>
</ul>
 
</html>"), Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}), graphics));
  end PumpNPSH;

  package BaseClasses "Base classes for Turbomachinery components"
    extends Modelica_Fluid.Icons.BaseClassLibrary;

  partial model PartialPump "Base model for centrifugal pumps"
      import Modelica.SIunits.Conversions.NonSIunits.*;
      import Modelica.Constants;
    replaceable package Medium = Modelica.Media.Interfaces.PartialMedium
        "Medium model" 
                     annotation(choicesAllMatching=true);
    replaceable function flowCharacteristic = 
        PumpCharacteristics.baseFlow
        "Head vs. q_flow characteristic at nominal speed and density" 
      annotation(Dialog(group="Characteristics"), choicesAllMatching=true);
    parameter Boolean usePowerCharacteristic = false
        "Use powerCharacteristic (vs. efficiencyCharacteristic)" 
       annotation(Evaluate=true,Dialog(group="Characteristics"));
    replaceable function powerCharacteristic = 
          PumpCharacteristics.quadraticPower (
         q_nom={0,0,0},W_nom={0,0,0})
        "Power consumption vs. q_flow at nominal speed and density" 
      annotation(Dialog(group="Characteristics", enable = usePowerCharacteristic),
                 choicesAllMatching=true);
    replaceable function efficiencyCharacteristic = 
      PumpCharacteristics.constantEfficiency(eta_nom = 0.8) constrainedby
        PumpCharacteristics.baseEfficiency
        "Efficiency vs. q_flow at nominal speed and density" 
      annotation(Dialog(group="Characteristics",enable = not usePowerCharacteristic),
                 choicesAllMatching=true);
    parameter AngularVelocity_rpm N_nom = 1500 "Nominal rotational speed" 
      annotation(Dialog(group="Characteristics"));
    parameter Medium.Density d_nom = Medium.density_pTX(Medium.p_default, Medium.T_default, Medium.X_default)
        "Nominal fluid density" 
      annotation(Dialog(group="Characteristics"));
    parameter Integer Np(min=1) = 1 "Number of pumps in parallel";
    parameter SI.Mass M = 0 "Fluid mass inside the pump";
    parameter Boolean checkValve=false "Reverse flow stopped";
    parameter Boolean allowFlowReversal = system.allowFlowReversal
        "allow flow reversal, false restricts to design direction (port_a -> port_b)"
      annotation(Dialog(tab="Assumptions"), Evaluate=true);
    parameter Medium.AbsolutePressure pin_start
        "Guess value for inlet pressure" 
      annotation(Dialog(tab="Initialization"));
    parameter Medium.AbsolutePressure pout_start
        "Guess value for outlet pressure" 
      annotation(Dialog(tab="Initialization"));
    parameter Boolean use_T_start = true
        "Use T_start if true, otherwise h_start" 
      annotation(Dialog(tab = "Initialization"), Evaluate = true);
    parameter Medium.Temperature T_start=
      if use_T_start then Medium.T_default else Medium.temperature_phX(pin_start,h_start,X_start)
        "Guess value for temperature" 
      annotation(Dialog(tab = "Initialization", enable = use_T_start));
    parameter Medium.SpecificEnthalpy h_start=
      if use_T_start then Medium.specificEnthalpy_pTX(pin_start, T_start, X_start) else Medium.h_default
        "Guess value for specific enthalpy" 
      annotation(Dialog(tab = "Initialization", enable = not use_T_start));
    parameter Medium.MassFraction X_start[Medium.nX] = Medium.X_default
        "Guess value for mass fractions m_i/m" 
      annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
    parameter SI.MassFlowRate m_flow_start = 0
        "Guess value for mass flow rate (total)" 
      annotation(Dialog(tab="Initialization"));
    outer Modelica_Fluid.System system "System properties";
    parameter SI.Acceleration g=system.g;
    parameter Types.Init initType=
              Types.Init.NoInit "Initialization option" 
      annotation(Evaluate=true, Dialog(tab = "Initialization"));
    Modelica_Fluid.Interfaces.FluidPort_a inlet(
        redeclare package Medium = Medium,
        p(start=pin_start),
        m_flow(start = m_flow_start,
               min = if allowFlowReversal and not checkValve then -Constants.inf else 0)) 
    annotation (Placement(transformation(extent={{
              -90,-10},{-70,10}}), iconTransformation(extent={{-90,-10},{-70,10}})));
    Modelica_Fluid.Interfaces.FluidPort_b outlet(
                                  redeclare package Medium = Medium,
        p(start=pout_start),
        m_flow(start = -m_flow_start,
               max = if allowFlowReversal and not checkValve then +Constants.inf else 0)) 
    annotation(Placement(transformation(extent={{
              70,-10},{90,10}}), iconTransformation(extent={{70,-10},{90,10}})));
    SI.Pressure dp = outlet.p - inlet.p "Pressure increase";
    SI.Height head = dp/(d*g) "Pump head";
    Medium.Density d "Liquid density at the inlet";
    Medium.SpecificEnthalpy h(start=h_start)
        "Enthalpy of the liquid stored in the pump if M>0";
    Medium.Temperature Tin "Liquid inlet temperature";
    SI.MassFlowRate m_flow = inlet.m_flow "Mass flow rate (total)";
    SI.MassFlowRate m_flow_single = m_flow/Np "Mass flow rate (single pump)";
    SI.VolumeFlowRate q_flow = m_flow/d "Volume flow rate (total)";
    SI.VolumeFlowRate q_flow_single = q_flow/Np
        "Volume flow rate (single pump)";
    AngularVelocity_rpm N "Shaft rotational speed";
    SI.Power W_single "Power Consumption (single pump)";
    SI.Power W_tot = W_single*Np "Power Consumption (total)";
    constant SI.Power W_eps=1e-8
        "Small coefficient to avoid numerical singularities in efficiency computations";
    Real eta "Global Efficiency";
    Real s(start = m_flow_start)
        "Curvilinear abscissa for the flow curve in parametric form (either mass flow rate or head)";
    Medium.ThermodynamicState inlet_state_inflow
        "Medium state close to inlet for inflowing mass flow";
    protected
    constant SI.Height unitHead = 1;
    constant SI.MassFlowRate unitMassFlowRate = 1;
  equation
    // Flow equations
    if noEvent(s > 0 or (not checkValve)) then
      // Flow characteristics when check valve is open or with no check valve
      head = (N/N_nom)^2*flowCharacteristic(q_flow_single*(N_nom/N));
      q_flow_single = s*unitMassFlowRate/d;
    else
      // Flow characteristics when check valve is closed
      head = (N/N_nom)^2*flowCharacteristic(0) - s*unitHead;
      q_flow_single = 0;
    end if;

    // Power consumption
    if usePowerCharacteristic then
      W_single = (N/N_nom)^3*(d/d_nom)*powerCharacteristic(q_flow_single*(N_nom/N))
          "Power consumption (single pump)";
      eta = (dp*q_flow_single)/(W_single + W_eps) "Hydraulic efficiency";
    else
      eta = efficiencyCharacteristic(q_flow_single*(N_nom/N));
      W_single = dp*q_flow_single/eta;
    end if;

    // Medium states close to the ports when mass flows in to the respective port
    // The inlet inflow state is used also in case of flow reversal, to avoid
    // discontinuities.
    inlet_state_inflow = Medium.setState_phX(inlet.p, inlet.h_outflow, inlet.Xi_outflow);
    // outlet_state_inflow = Medium.setState_phX(outlet.p, outlet.h_outflow, outlet.Xi_outflow);

    // Inflow density and temperature at the inlet port
    d = Medium.density(inlet_state_inflow);
    Tin = Medium.temperature(inlet_state_inflow);

    // Mass balances
    inlet.m_flow + outlet.m_flow = 0 "Mass balance";
    inlet.Xi_outflow  = inStream(outlet.Xi_outflow);
    outlet.Xi_outflow = inStream(inlet.Xi_outflow);
    inlet.C_outflow = inStream(outlet.C_outflow);
    outlet.C_outflow = inStream(inlet.C_outflow);

    // Energy balances
    if M > 0 then
       // Dynamic energy balance
       // mass variations and p/d are neglected
       Np*M*der(h) = inlet.m_flow*actualStream(inlet.h_outflow) +
                     outlet.m_flow*actualStream(outlet.h_outflow) +
                     W_tot;
       outlet.h_outflow  = h;
       inlet.h_outflow   = h;
    else
      /* In the following two equations the extra term W_tot/m_flow is
       present where a 0/0 term appears if m_flow = 0. In order to avoid
       numerical problems, this term is analytically simplified:
         W_tot/m_flow = Np*W_single/(d*q_flow_single*Np)
                      = Np*dp*q_flow_single/(eta*d*q_flow_single*Np)
                      = dp/(eta*d)
    */
      outlet.h_outflow  = inStream(inlet.h_outflow) + dp/(d*eta);
      inlet.h_outflow  = inStream(outlet.h_outflow) + dp/(d*eta);
      h = Medium.h_default
          "Unused (set to an arbitrary value within the medium region)";
    end if;

  initial equation
    if initType == Types.Init.NoInit or not M > 0 then
    // no initial equations
    elseif initType == Types.Init.InitialValues then
      h = h_start;
    elseif initType == Types.Init.SteadyState then
      der(h) = 0;
    else
      assert(false, "Unsupported initialization option");
    end if;
    annotation (
      Icon(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},{100,
                100}}), graphics={
            Polygon(
              points={{-40,-64},{-60,-100},{60,-100},{40,-64},{-40,-64}},
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
              lineColor={0,0,255})}),
      Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},{
                100,100}}),
              graphics),
      Documentation(info="<HTML>
<p>This is the base model for the <tt>Pump</tt> and <tt>
PumpMech</tt> pump models.
<p>The model describes a centrifugal pump, or a group of <tt>Np</tt> identical pumps in parallel. The pump model is based on the theory of kinematic similarity: the pump characteristics are given for nominal operating conditions (rotational speed and fluid density), and then adapted to actual operating condition, according to the similarity equations. 
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
<p>If the <tt>in_Np</tt> input connector is wired, it provides the number of pumps in parallel; otherwise,  <tt>Np_n</tt> parallel pumps are assumed.</p>
<p>It is possible to take into account the heat capacity of the fluid inside the pump by specifying its mass <tt>M</tt> at nominal conditions; this is necessary to avoid singularities in the computation of the outlet enthalpy in case of zero flow rate. If zero flow rate conditions are always avoided, this dynamic effect can be neglected by leaving the default value <tt>M = 0</tt>, thus avoiding a fast state variable in the model.
<p>If <tt>computeNPSHa = true</tt>, the available net positive suction head is also computed; this requires a two-phase medium model to provide the fluid saturation pressure.
</HTML>",
        revisions="<html>
<ul>
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
      input SI.VolumeFlowRate q_nom[2]
          "Volume flow rate for two operating points (single pump)" 
                                                                  annotation(Dialog);
      input SI.Height head_nom[2] "Pump head for two operating points" annotation(Dialog);
      /* Linear system to determine the coefficients:
  head_nom[1] = c[1] + q_nom[1]*c[2];
  head_nom[2] = c[1] + q_nom[2]*c[2];
  */
      protected
      Real c[2] = Modelica.Math.Matrices.solve([ones(2),q_nom],head_nom)
          "Coefficients of linear head curve";
    algorithm
      // Flow equation: head = q*c[1] + c[2];
      head := c[1] + q_flow*c[2];
    end linearFlow;

    function quadraticFlow "Quadratic flow characteristic"
      extends baseFlow;
      input SI.VolumeFlowRate q_nom[3]
          "Volume flow rate for three operating points (single pump)" 
                                                                    annotation(Dialog);
      input SI.Height head_nom[3] "Pump head for three operating points" annotation(Dialog);
      protected
      Real q_nom2[3] = {q_nom[1]^2,q_nom[2]^2, q_nom[3]^2}
          "Squared nominal flow rates";
      /* Linear system to determine the coefficients:
  head_nom[1] = c[1] + q_nom[1]*c[2] + q_nom[1]^2*c[3];
  head_nom[2] = c[1] + q_nom[2]*c[2] + q_nom[2]^2*c[3];
  head_nom[3] = c[1] + q_nom[3]*c[2] + q_nom[3]^2*c[3];
  */
      Real c[3] = Modelica.Math.Matrices.solve([ones(3), q_nom, q_nom2],head_nom)
          "Coefficients of quadratic head curve";
    algorithm
      // Flow equation: head  = c[1] + q_flow*c[2] + q_flow^2*c[3];
      head := c[1] + q_flow*c[2] + q_flow^2*c[3];
    end quadraticFlow;

    function polynomialFlow "Polynomial flow characteristic"
      extends baseFlow;
      input SI.VolumeFlowRate q_nom[:]
          "Volume flow rate for N operating points (single pump)" 
                                                                annotation(Dialog);
      input SI.Height head_nom[:] "Pump head for N operating points" annotation(Dialog);
      protected
      Integer N = size(q_nom,1) "Number of nominal operating points";
      Real q_nom_pow[N,N] = {{q_nom[j]^(i-1) for j in 1:N} for i in 1:N}
          "Rows: different operating points; columns: increasing powers";
      /* Linear system to determine the coefficients (example N=3):
  head_nom[1] = c[1] + q_nom[1]*c[2] + q_nom[1]^2*c[3];
  head_nom[2] = c[1] + q_nom[2]*c[2] + q_nom[2]^2*c[3];
  head_nom[3] = c[1] + q_nom[3]*c[2] + q_nom[3]^2*c[3];
  */
      Real c[N] = Modelica.Math.Matrices.solve(q_nom_pow,head_nom)
          "Coefficients of polynomial head curve";
    algorithm
      // Flow equation (example N=3): head  = c[1] + q_flow*c[2] + q_flow^2*c[3];
      // Note: the implementation is numerically efficient only for low values of Na
      head := sum(q_flow^(i-1)*c[i] for i in 1:N);
    end polynomialFlow;

    function constantEfficiency "Constant efficiency characteristic"
       extends baseEfficiency;
       input Real eta_nom "Nominal efficiency" annotation(Dialog);
    algorithm
      eta := eta_nom;
    end constantEfficiency;

    function linearPower "Linear power consumption characteristic"
      extends basePower;
      input SI.VolumeFlowRate q_nom[2]
          "Volume flow rate for two operating points (single pump)" annotation(Dialog);
      input SI.Power W_nom[2] "Power consumption for two operating points"   annotation(Dialog);
      /* Linear system to determine the coefficients:
  W_nom[1] = c[1] + q_nom[1]*c[2];
  W_nom[2] = c[1] + q_nom[2]*c[2];
  */
      protected
      Real c[2] = Modelica.Math.Matrices.solve([ones(3),q_nom],W_nom)
          "Coefficients of linear power consumption curve";
    algorithm
      consumption := c[1] + q_flow*c[2];
    end linearPower;

    function quadraticPower "Quadratic power consumption characteristic"
      extends basePower;
      input SI.VolumeFlowRate q_nom[3]
          "Volume flow rate for three operating points (single pump)" 
                                                                    annotation(Dialog);
      input SI.Power W_nom[3] "Power consumption for three operating points" annotation(Dialog);
      protected
      Real q_nom2[3] = {q_nom[1]^2,q_nom[2]^2, q_nom[3]^2}
          "Squared nominal flow rates";
      /* Linear system to determine the coefficients:
  W_nom[1] = c[1] + q_nom[1]*c[2] + q_nom[1]^2*c[3];
  W_nom[2] = c[1] + q_nom[2]*c[2] + q_nom[2]^2*c[3];
  W_nom[3] = c[1] + q_nom[3]*c[2] + q_nom[3]^2*c[3];
  */
      Real c[3] = Modelica.Math.Matrices.solve([ones(3),q_nom,q_nom2],W_nom)
          "Coefficients of quadratic power consumption curve";
    algorithm
      consumption := c[1] + q_flow*c[2] + q_flow^2*c[3];
    end quadraticPower;

  end PumpCharacteristics;
  end BaseClasses;
  annotation (Documentation(info="<html>
 
</html>"));
end Pumps;
