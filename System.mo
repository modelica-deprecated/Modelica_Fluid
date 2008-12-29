within Modelica_Fluid;
model System
  "System properties and default values (ambient, flow direction, initialization)"

  replaceable package Medium = 
    Modelica.Media.Interfaces.PartialMedium "Default Medium model" 
      annotation (choicesAllMatching = true);
  parameter Medium.AbsolutePressure p_ambient=101325 "Default ambient pressure"
    annotation(Dialog(group="Environment"));
  parameter Medium.Temperature T_ambient=293.15 "Default ambient temperature" 
    annotation(Dialog(group="Environment"));
  parameter SI.Acceleration g=Modelica.Constants.g_n
    "Constant gravity acceleration" 
    annotation(Dialog(group="Environment"));

  // Assumptions
  parameter Boolean allowFlowReversal = true
    "= false to restrict to design flow direction (port_a -> port_b)" 
    annotation(Dialog(tab="Assumptions"), Evaluate=true);
  parameter Types.Dynamics energyDynamics=
    Types.Dynamics.DynamicFreeInitial "Default formulation of energy balances" 
    annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));
  parameter Types.Dynamics massDynamics=
    energyDynamics "Default formulation of mass balances" 
    annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));
  final parameter Types.Dynamics substanceDynamics=
    massDynamics "Default formulation of substance balances" 
    annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));
  final parameter Types.Dynamics traceDynamics=
    massDynamics "Default formulation of trace substance balances" 
    annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));
  parameter Types.Dynamics momentumDynamics=
    massDynamics
    "Default formulation of momentum balances, if options available" 
    annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));

  // Initialization
  parameter Medium.AbsolutePressure p_start = Medium.p_default
    "Default start value for pressures" 
    annotation(Dialog(tab = "Initialization"));
  parameter Medium.MassFlowRate m_flow_start = 0
    "Default start value for mass flow rates" 
    annotation(Dialog(tab = "Initialization"));
  parameter Medium.Temperature T_start = Medium.T_default
    "Default start value for temperatures" 
    annotation(Dialog(tab = "Initialization"));

  annotation (
    defaultComponentName="system",
    defaultComponentPrefixes="inner",
    missingInnerMessage="
Your model is using an outer \"system\" component but
an inner \"system\" component is not defined.
For simulation drag Modelica_Fluid.System into your model
to specify system properties. The default System setting
is used for the current simulation.
",  Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
            100}}), graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-150,150},{150,110}},
          lineColor={0,0,255},
          textString="%name"),
        Line(points={{-86,-30},{82,-30}}, color={0,0,0}),
        Line(points={{-82,-68},{-52,-30}}, color={0,0,0}),
        Line(points={{-48,-68},{-18,-30}}, color={0,0,0}),
        Line(points={{-14,-68},{16,-30}}, color={0,0,0}),
        Line(points={{22,-68},{52,-30}}, color={0,0,0}),
        Line(points={{74,84},{74,14}}, color={0,0,0}),
        Polygon(
          points={{60,14},{88,14},{74,-18},{60,14}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{16,20},{60,-18}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid,
          textString="g"),
        Text(
          extent={{-90,82},{74,50}},
          lineColor={0,0,0},
          textString="defaults"),
        Line(
          points={{-82,14},{-42,-20},{2,30}},
          color={0,0,0},
          thickness=0.5),
        Ellipse(
          extent={{-10,40},{12,18}},
          pattern=LinePattern.None,
          lineColor={0,0,0},
          fillColor={255,0,0},
          fillPattern=FillPattern.Solid)}),
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
            100,100}}),
            graphics),
    Documentation(info="<html>
<p>
 A system component is needed in each fluid model to provide system-wide settings, such as ambient conditions and overall modeling assumptions.
 The system settings are propagated to the fluid models using the inner/outer mechanism. 
</p>
<p>
 A model should never directly use system parameters. 
 Instead a local parameter should be declared, which uses the global setting as default. 
 The only exception currently made is the gravity system.g.
</p>
</html>"));

end System;