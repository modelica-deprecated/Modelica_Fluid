within Modelica_Fluid;
package Pumps "Pump components" 
  extends Modelica_Fluid.Icons.VariantLibrary;
  model Pump "Centrifugal pump with ideally controlled speed" 
    extends Modelica_Fluid.Pumps.BaseClasses.PartialPump;
    parameter Modelica.SIunits.Conversions.NonSIunits.AngularVelocity_rpm 
      N_const =                                                                     N_nom 
      "Constant rotational speed" annotation(Dialog(enable = not use_N_input));
    parameter Boolean use_N_input = false 
      "Get the rotational speed from the input connector";
    Modelica.Blocks.Interfaces.RealInput N_in(redeclare type SignalType = 
          Modelica.SIunits.Conversions.NonSIunits.AngularVelocity_rpm) if use_N_input 
      "Prescribed rotational speed" 
      annotation (extent=[-36,34; -16,54],   rotation=-90);
  protected 
    Modelica.Blocks.Interfaces.RealInput N_in_internal(
      redeclare type SignalType = 
          Modelica.SIunits.Conversions.NonSIunits.AngularVelocity_rpm) 
      "Needed to connect to conditional connector";
  equation 
    // Connect statement active only if usePressureInput = true
    connect(N_in, N_in_internal);
    // Internal connector value when usePressureInput = false
    if not use_N_input then
      N_in_internal = N_const;
    end if;
    N = N_in_internal "Rotational speed";
    annotation (
      Icon(
        Text(extent=[-58,58; -30,38], string="n")),
      Diagram,
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
  end Pump;
  
  model PumpShaft "Centrifugal pump with mechanical connector for the shaft" 
    extends Modelica_Fluid.Pumps.BaseClasses.PartialPump;
    SI.Angle phi "Shaft angle";
    SI.AngularVelocity omega "Shaft angular velocity";
    Modelica.Mechanics.Rotational.Interfaces.Flange_a shaft 
    annotation (extent=[80,4; 110,32]);
  equation 
    phi = shaft.phi;
    omega = der(phi);
    N = Modelica.SIunits.Conversions.to_rpm(omega);
    W_single = omega*shaft.tau;
  annotation (
    Icon(Rectangle(extent=[60,26; 84,12], style(
            color=10,
            rgbcolor={95,95,95},
            gradient=2,
            fillColor=10,
            rgbfillColor={95,95,95}))),
    Diagram,
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
      Modelica.Media.Water.WaterIF97_ph extends 
        Modelica.Media.Interfaces.PartialTwoPhaseMedium);
    SI.Length NPSHa "Net Positive Suction Head available";
    Medium.AbsolutePressure pv "Saturation pressure of inlet liquid";
    
  equation 
    // NPSHa computation
    pv = Medium.saturationPressure(Tin);
    NPSHa = (inlet.p-pv)/(d*Modelica.Constants.g_n);
    
    // Check for cavitation
    assert(inlet.p >= pv, "Cavitation occurs at the inlet (propably a valve is closed or a tank/reservoir is empty)");
    annotation (Documentation(info="<html>Same as the Pump model, with added computation of Net Positive Suction Head available. Requires a two-phase medium model.
</html>", revisions="<html>
<ul>
<li><i>30 Jul 2007</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Model added to the Fluid library</li>
</ul>

</html>"));
  end PumpNPSH;
  
  package BaseClasses "Base classes for Turbomachinery components" 
    extends Modelica_Fluid.Icons.BaseClassLibrary;
    
  partial model PartialPump "Base model for centrifugal pumps" 
      import Modelica.SIunits.Conversions.NonSIunits.*;
      import Modelica.Constants;
    replaceable package Medium = Modelica.Media.Interfaces.PartialMedium 
        "Medium model" 
                     annotation(choicesAllMatching=true);
      
  /*
  Medium.BaseProperties fluid(p(start=pin_start),h(start=h_start)) 
    "Fluid properties at the inlet";
*/
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
      PumpCharacteristics.constantEfficiency(eta_nom = 0.8) extends 
        PumpCharacteristics.baseEfficiency 
        "Efficiency vs. q_flow at nominal speed and density" 
      annotation(Dialog(group="Characteristics",enable = not usePowerCharacteristic),
                 choicesAllMatching=true);
    parameter AngularVelocity_rpm N_nom = 1500 "Nominal rotational speed" 
      annotation(Dialog(group="Characteristics"));
    parameter Medium.Density d_nom = 1000 "Nominal fluid density" 
      annotation(Dialog(group="Characteristics"));
    parameter Integer Np(min=1) = 1 "Number of pumps in parallel";
    parameter SI.Mass M = 0 "Fluid mass inside the pump";
    parameter Boolean checkValve=false "Reverse flow stopped";
    parameter Types.FlowDirection.Temp flowDirection=
                     Types.FlowDirection.Bidirectional 
        "Unidirectional (inlet -> outlet) or bidirectional flow component" 
       annotation(Dialog(tab="Advanced"));
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
    outer Modelica_Fluid.Ambient ambient "Ambient conditions";
    parameter SI.Acceleration g=ambient.g;
  //  parameter Choices.Init.Options.Temp initOpt=Choices.Init.Options.noInit 
  //    "Initialisation option";
    Modelica_Fluid.Interfaces.FluidPort_a inlet(
                                 redeclare package Medium = Medium,
        p(start=pin_start),
        m_flow(start = m_flow_start,
               min = if allowFlowReversal and not checkValve then -Constants.inf else 0)) 
    annotation (extent=[-100,-40; -60,0]);
    Modelica_Fluid.Interfaces.FluidPort_b outlet(
                                  redeclare package Medium = Medium,
        p(start=pout_start),
        m_flow(start = -m_flow_start,
               max = if allowFlowReversal and not checkValve then +Constants.inf else 0)) 
    annotation (extent=[40,12; 80,52]);
    SI.Pressure dp = outlet.p - inlet.p "Pressure increase";
    SI.Height head = dp/(d*g) "Pump head";
    Medium.Density d "Liquid density at the inlet";
    Medium.SpecificEnthalpy h_out(start=h_start) 
        "Enthalpy of the liquid flowing out of the pump";
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
  //  outer Modelica_Fluid.Components.FluidOptions fluidOptions 
  //    "Global default options";
      
    Medium.ThermodynamicState inlet_state_inflow 
        "Medium state close to inlet for inflowing mass flow";
    Medium.EnthalpyFlowRate inlet_H_flow;
    Medium.EnthalpyFlowRate outlet_H_flow;
    protected 
   parameter Boolean allowFlowReversal=
       flowDirection == Modelica_Fluid.Types.FlowDirection.Bidirectional 
        "= false, if flow only from port_a to port_b, otherwise reversing flow allowed"
       annotation(Evaluate=true, Hide=true);
    constant SI.Height unitHead = 1;
    constant SI.MassFlowRate unitMassFlowRate = 1;
  equation 
    // Flow equations
    if noEvent(s > 0 or (not checkValve)) then
      // Flow characteristics when check valve is open
      // q_flow_single = s;
      q_flow_single = s*unitMassFlowRate/d;
      // head = (N/N_nom)^2*flowCharacteristic(q_flow_single*N_nom/(noEvent(if abs(N) > 1e-6 then N else 1e-6)));
      head = noEvent((((if abs(N) > 1e-6 then N else 1e-6))/N_nom)^2*flowCharacteristic(q_flow_single*N_nom/((if abs(N) > 1e-6 then N else 1e-6))));
    else
      // Flow characteristics when check valve is closed
      head = (N/N_nom)^2*flowCharacteristic(0) - s*unitHead;
      q_flow_single = 0;
    end if;
      
    // Power consumption  
    if usePowerCharacteristic then
      W_single = (N/N_nom)^3*(d/d_nom)*powerCharacteristic(q_flow_single*N_nom/(noEvent(if abs(N) > 1e-6 then N else 1e-6))) 
          "Power consumption (single pump)";
      eta = (dp*q_flow_single)/(W_single + W_eps) "Hydraulic efficiency";
    else
      eta = efficiencyCharacteristic(q_flow_single*N_nom/(noEvent(if abs(N) > 1e-6 then N else 1e-10)));
      W_single = dp*q_flow/eta;
    end if;
      
    // Medium states close to the ports when mass flows in to the respective port
    inlet_state_inflow = Medium.setState_phX(inlet.p, inlet.h_outflow, inlet.Xi_outflow);
    // outlet_state_inflow = Medium.setState_phX(outlet.p, outlet.h_outflow, outlet.Xi_outflow);
      
    // Inflow density and temperature at the inlet port
    d = Medium.density(inlet_state_inflow);
    Tin = Medium.temperature(inlet_state_inflow);
      
    // Mass and energy balances
    inlet.m_flow + outlet.m_flow = 0 "Mass balance";
    inlet.h_outflow   = h_out;
    outlet.h_outflow  = h_out;
    inlet.Xi_outflow  = inflow(outlet.Xi_outflow);
    outlet.Xi_outflow = inflow(inlet.Xi_outflow);
    inlet.C_outflow = inflow(outlet.C_outflow);
    outlet.C_outflow = inflow(inlet.C_outflow);
    inlet_H_flow=semiLinear(inlet.m_flow, inflow(inlet.h_outflow), h_out) 
        "Enthalpy flow at the inlet";
    outlet_H_flow=semiLinear(outlet.m_flow, inflow(outlet.h_outflow), h_out) 
        "Enthalpy flow at the outlet";
    if M > 0 then
      // M * der(h_out) = m_flow_single*(inlet.h - outlet.h) + W_single 
      (M*Np) * der(h_out) = inlet_H_flow + outlet_H_flow + W_tot 
          "Dynamic energy balance (density variations neglected)";
    else
      inlet_H_flow + outlet_H_flow + W_tot = 0 "Static energy balance";
    end if;
      
  /*
initial equation 
  if initOpt == Choices.Init.Options.noInit then
    // do nothing
  elseif initOpt == Choices.Init.Options.steadyState then
    if ThermalCapacity then
      der(h)=0;
    end if;
  else
    assert(false, "Unsupported initialisation option");
  end if;
*/
      
    annotation (
      Icon(
        Polygon(points=[-40,-64; -60,-100; 60,-100; 40,-64; -40,-64],
            style(pattern=0, fillColor=74)),
        Ellipse(extent=[-60,40; 60,-80],   style(gradient=3)),
        Polygon(points=[-30,12; -30,-48; 48,-20; -30,12],   style(
            pattern=0,
            gradient=2,
            fillColor=7)),
        Text(extent=[-100,-110; 100,-136], string="%name"),
        Text(extent=[-10,60; 18,40],  string="Np")),
      Diagram,
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
