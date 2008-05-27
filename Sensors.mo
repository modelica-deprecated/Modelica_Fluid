within Modelica_Fluid;
package Sensors
  "Ideal sensor components (provide the variables in a fluid connector as signals)"
  extends Modelica_Fluid.Icons.VariantLibrary;

  annotation (preferedView="info", Documentation(info="<html>
<p align = justify>
Package <b>Sensors</b> consists of idealized sensor components that
provide variables of a medium model and/or fluid ports as
output signals. These signals can be, e.g., further processed
with components of the Modelica.Blocks library.
Also more realistic sensor models can be built, by further
processing (e.g., by attaching block Modelica.Blocks.FirstOrder to
model the time constant of the sensor).
 
</p>
 
<p align = justify>For the thermodynamic state variables temperature, specific entalpy, specific entropy and density the fluid library provides two different types of sensors: <b>one port</b> and <b>two port</b> sensors. </p>
 
<ul>
<li>
The <b>one port</b> sensors have the advantage of easily introducing them to and removing them from a complete model, due to no connections have to be broken. But if the sensor is placed in a three or more way connection without any explicit junction model, the resulting temperature signal can be something unexpected. <a href= \"Modelica_Fluid.Test.TestComponents.Sensors.TestTemperatureSensor\">Modelica_Fluid.Test.TestComponents.Sensors.TestTemperatureSensor </a> provides a test case, which shows the arising trouble. That means, that one port sensors should only be used in combination with explicit junction models (<a href= \"Modelica_Fluid.Junctions\">Modelica_Fluid.Junctions</a>) or in the case that the user is really sure to connect the sensor only between two other components! </li> 
 
<li> The <b>two port</b> sensors are the save way for getting the expecting results. If they are connected in series the output signal corresponds to the passing fluid. To provide a temperature signal of a volume, e.g. a tank, the two port sensors can be used only with one port connected to the volume.</li>
</ul>
 
</html>",
      revisions="<html>
<ul>
<li><i>31 Oct 2007</i>
    by Carsten Heinrich<br>
       updated sensor models, included one and two port sensors for thermodynamic state variables</li>
</ul>
</html>"));

  model Pressure "Ideal pressure sensor"
    extends Sensors.BaseClasses.PartialAbsoluteSensor;
    extends Modelica.Icons.RotationalSensor;
    Modelica.Blocks.Interfaces.RealOutput p(final quantity="Pressure",
                                            final unit="Pa",
                                            displayUnit="bar",
                                            min=0) "Pressure at port" 
      annotation (Placement(transformation(extent={{100,-10},{120,10}},
            rotation=0)));
    annotation (
    Diagram(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics={Line(points={{70,0},{100,0}}, color={0,0,127}), 
            Line(points={{0,-70},{0,-100}}, color={0,127,255})}),
    Icon(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics={
          Line(points={{70,0},{100,0}}, color={0,0,127}),
          Line(points={{0,-70},{0,-100}}, color={0,127,255}),
          Text(extent={{-150,80},{150,120}}, textString=
                                               "%name"),
          Text(
            extent={{212,-51},{52,-103}},
            lineColor={0,0,0},
            textString=
                 "p")}),
      Documentation(info="<HTML>
<p>
This component monitors the absolute pressure at its fluid port. The sensor is 
ideal, i.e., it does not influence the fluid.
</p>
</HTML>
"));
  equation
    p = port.p;
  end Pressure;

  model DensityOnePort "Ideal one port density sensor"
    extends Sensors.BaseClasses.PartialAbsoluteSensor;
    extends Modelica.Icons.RotationalSensor;
    Modelica.Blocks.Interfaces.RealOutput d(final quantity="Density",
                                            final unit="kg/m3",
                                            displayUnit="g/cm3",
                                            min=0) "Density in port medium" 
      annotation (Placement(transformation(extent={{100,-10},{120,10}},
            rotation=0)));

  annotation (defaultComponentName="density",
    Diagram(graphics={Line(points={{0,-70},{0,-100}}, color={0,0,127}), Line(
              points={{70,0},{100,0}}, color={0,0,127})}),
    Icon(graphics={
          Line(points={{0,-70},{0,-100}}, color={0,0,127}),
          Text(extent={{-150,80},{150,120}}, textString=
                                               "%name"),
          Text(
            extent={{170,-53},{10,-105}},
            lineColor={0,0,0},
            textString=
                 "d"),
          Line(points={{70,0},{100,0}}, color={0,0,127})}),
    Documentation(info="<HTML>
<p>
This component monitors the density of the medium in the flow
between fluid ports. The sensor is 
ideal, i.e., it does not influence the fluid.
</p>
<p>If using the one port sensor please read the <a href = Modelica_Fluid.Sensors>Information</a>  first.</p>
 
</HTML>
"));
  equation
    d = Medium.density(Medium.setState_phX(port.p, inflow(port.h_outflow), inflow(port.Xi_outflow)));
  end DensityOnePort;

  model DensityTwoPort "Ideal two port density sensor"
    extends Sensors.BaseClasses.PartialFlowSensor;
    extends Modelica.Icons.RotationalSensor;
    Modelica.Blocks.Interfaces.RealOutput d(final quantity="Density",
                                            final unit="kg/m3",
                                            displayUnit="g/cm3",
                                            min=0)
      "Density of the passing fluid" 
      annotation (Placement(transformation(
          origin={0,-110},
          extent={{-10,-10},{10,10}},
          rotation=270)));
    parameter Medium.MassFlowRate m_flow_small(min=0) = 1e-4
      "For bi-directional flow, density is regularized in the region |m_flow| < m_flow_small (m_flow_small > 0 required)"
      annotation(Dialog(tab="Advanced"));
  annotation (defaultComponentName="density",
    Diagram(graphics={
          Line(points={{-100,0},{-70,0}}, color={0,128,255}),
          Line(points={{70,0},{100,0}}, color={0,128,255}),
          Line(points={{0,-70},{0,-100}}, color={0,0,127})}),
    Icon(graphics={
          Text(extent={{-150,80},{150,120}}, textString=
                                               "%name"),
          Text(
            extent={{168,-71},{8,-123}},
            lineColor={0,0,0},
            textString=
                 "d"),
          Line(points={{0,-70},{0,-100}}, color={0,0,127}),
          Line(points={{-100,0},{-70,0}}, color={0,128,255}),
          Line(points={{70,0},{100,0}}, color={0,128,255})}),
    Documentation(info="<HTML>
<p>
This component monitors the volume flow rate flowing from port_a to port_b. 
The sensor is ideal, i.e. it does not influence the fluid.
</p>
</HTML>
"));
  protected
    Medium.Density d_a_inflow "Density of inflowing fluid at port_a";
    Medium.Density d_b_inflow
      "Density of inflowing fluid at port_b or d_a_inflow, if uni-directional flow";
  equation
    if allowFlowReversal then
       d_a_inflow = Medium.density(Medium.setState_phX(port_b.p, port_b.h_outflow, port_b.Xi_outflow));
       d_b_inflow = Medium.density(Medium.setState_phX(port_a.p, port_a.h_outflow, port_a.Xi_outflow));
       d = Modelica_Fluid.Utilities.regStep(port_a.m_flow, d_a_inflow, d_b_inflow, m_flow_small);
    else
       d = Medium.density(Medium.setState_phX(port_b.p, port_b.h_outflow, port_b.Xi_outflow));
       d_a_inflow = d;
       d_b_inflow = d;
    end if;
  end DensityTwoPort;

  model TemperatureOnePort "Ideal one port temperature sensor"
      extends Sensors.BaseClasses.PartialAbsoluteSensor;

    Modelica.Blocks.Interfaces.RealOutput T(final quantity="ThermodynamicTemperature",
                                            final unit = "K", displayUnit = "degC", min=0)
      "Temperature in port medium" 
      annotation (Placement(transformation(extent={{60,-10},{80,10}}, rotation=
              0)));

  annotation (defaultComponentName="temperature",
    Diagram(graphics={
          Line(points={{0,-70},{0,-100}}, color={0,0,127}),
          Ellipse(
            extent={{-20,-98},{20,-60}},
            lineColor={0,0,0},
            lineThickness=0.5,
            fillColor={191,0,0},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-12,40},{12,-68}},
            lineColor={191,0,0},
            fillColor={191,0,0},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{-12,40},{-12,80},{-10,86},{-6,88},{0,90},{6,88},{10,86},{
                12,80},{12,40},{-12,40}},
            lineColor={0,0,0},
            lineThickness=0.5),
          Line(
            points={{-12,40},{-12,-64}},
            color={0,0,0},
            thickness=0.5),
          Line(
            points={{12,40},{12,-64}},
            color={0,0,0},
            thickness=0.5),
          Line(points={{-40,-20},{-12,-20}}, color={0,0,0}),
          Line(points={{-40,20},{-12,20}}, color={0,0,0}),
          Line(points={{-40,60},{-12,60}}, color={0,0,0}),
          Line(points={{12,0},{60,0}}, color={0,0,127})}),
      Icon(graphics={
          Ellipse(
            extent={{-20,-88},{20,-50}},
            lineColor={0,0,0},
            lineThickness=0.5,
            fillColor={191,0,0},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-12,50},{12,-58}},
            lineColor={191,0,0},
            fillColor={191,0,0},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{-12,50},{-12,90},{-10,96},{-6,98},{0,100},{6,98},{10,96},{
                12,90},{12,50},{-12,50}},
            lineColor={0,0,0},
            lineThickness=0.5),
          Line(
            points={{-12,50},{-12,-54}},
            color={0,0,0},
            thickness=0.5),
          Line(
            points={{12,50},{12,-54}},
            color={0,0,0},
            thickness=0.5),
          Line(points={{-40,-10},{-12,-10}}, color={0,0,0}),
          Line(points={{-40,30},{-12,30}}, color={0,0,0}),
          Line(points={{-40,70},{-12,70}}, color={0,0,0}),
          Text(
            extent={{120,-40},{0,-90}},
            lineColor={0,0,0},
            textString=
                 "T"),
          Text(extent={{-150,110},{150,150}}, textString=
                                                  "%name"),
          Line(points={{12,0},{60,0}}, color={0,0,127})}),
      Documentation(info="<HTML>
<p>
This component monitors the temperature of the medium in the flow
between fluid ports. The sensor is 
ideal, i.e., it does not influence the fluid.
</p>
</HTML>
"));
  equation
    T = Medium.temperature(Medium.setState_phX(port.p, inflow(port.h_outflow), inflow(port.Xi_outflow)));
  end TemperatureOnePort;

  model TemperatureTwoPort "Ideal two port temperature sensor"
    extends Sensors.BaseClasses.PartialFlowSensor;

    Modelica.Blocks.Interfaces.RealOutput T( final quantity="ThermodynamicTemperature",
                                             final unit="K",
                                             min = 0,
                                             displayUnit="degC")
      "Temperature of the passing fluid" 
      annotation (Placement(transformation(
          origin={0,-110},
          extent={{-10,-10},{10,10}},
          rotation=270)));
    parameter Medium.MassFlowRate m_flow_small(min=0) = 1e-4
      "For bi-directional flow, temperature is regularized in the region |m_flow| < m_flow_small (m_flow_small > 0 required)"
      annotation(Dialog(tab="Advanced"));

  annotation (defaultComponentName="temperature",
    Diagram(graphics),
    Icon(graphics={
          Text(extent={{-150,110},{150,150}}, textString=
                                                "%name"),
          Line(points={{0,-70},{0,-100}}, color={0,0,127}),
          Line(points={{-92,0},{100,0}}, color={0,128,255}),
          Ellipse(
            extent={{-20,-88},{20,-50}},
            lineColor={0,0,0},
            lineThickness=0.5,
            fillColor={191,0,0},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-12,50},{12,-58}},
            lineColor={191,0,0},
            fillColor={191,0,0},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{-12,50},{-12,90},{-10,96},{-6,98},{0,100},{6,98},{10,96},{
                12,90},{12,50},{-12,50}},
            lineColor={0,0,0},
            lineThickness=0.5),
          Line(
            points={{-12,50},{-12,-54}},
            color={0,0,0},
            thickness=0.5),
          Line(
            points={{12,50},{12,-54}},
            color={0,0,0},
            thickness=0.5),
          Line(points={{-40,-10},{-12,-10}}, color={0,0,0}),
          Line(points={{-40,30},{-12,30}}, color={0,0,0}),
          Line(points={{-40,70},{-12,70}}, color={0,0,0}),
          Text(
            extent={{120,-40},{0,-90}},
            lineColor={0,0,0},
            textString=
                 "T")}),
    Documentation(info="<HTML>
<p>
This component monitors the temperature of the passing fluid. 
The sensor is ideal, i.e. it does not influence the fluid.
</p>
</HTML>
"));
  protected
    Medium.Temperature T_a_inflow "Temperature of inflowing fluid at port_a";
    Medium.Temperature T_b_inflow
      "Temperature of inflowing fluid at port_b or T_a_inflow, if uni-directional flow";
  equation
    if allowFlowReversal then
       T_a_inflow = Medium.temperature(Medium.setState_phX(port_b.p, port_b.h_outflow, port_b.Xi_outflow));
       T_b_inflow = Medium.temperature(Medium.setState_phX(port_a.p, port_a.h_outflow, port_a.Xi_outflow));
       T = Modelica_Fluid.Utilities.regStep(port_a.m_flow, T_a_inflow, T_b_inflow, m_flow_small);
    else
       T = Medium.temperature(Medium.setState_phX(port_b.p, port_b.h_outflow, port_b.Xi_outflow));
       T_a_inflow = T;
       T_b_inflow = T;
    end if;
  end TemperatureTwoPort;

  model SpecificEnthalpyOnePort "Ideal one port specific enthalphy sensor"
    extends Sensors.BaseClasses.PartialAbsoluteSensor;
    extends Modelica.Icons.RotationalSensor;
    Modelica.Blocks.Interfaces.RealOutput h_out(final quantity="SpecificEnergy",
                                                final unit="J/kg")
      "Specific enthalpy in port medium" 
      annotation (Placement(transformation(extent={{100,-10},{120,10}},
            rotation=0)));

  annotation (defaultComponentName="specificEnthalpy",
    Diagram(graphics={Line(points={{0,-70},{0,-100}}, color={0,0,127}), Line(
              points={{70,0},{100,0}}, color={0,0,127})}),
    Icon(graphics={
          Line(points={{0,-70},{0,-100}}, color={0,0,127}),
          Text(extent={{-150,80},{150,120}}, textString=
                                               "%name"),
          Text(
            extent={{212,-51},{52,-103}},
            lineColor={0,0,0},
            textString=
                 "h"),
          Line(points={{70,0},{100,0}}, color={0,0,127})}),
    Documentation(info="<HTML>
<p>
This component monitors the specific enthalphy of the medium in the flow
between fluid ports. The sensor is ideal, i.e., it does not influence the fluid.
</p>
</HTML>
"));
  equation
    h_out = inflow(port.h_outflow);
  end SpecificEnthalpyOnePort;

  model SpecificEnthalpyTwoPort
    "Ideal two port sensor for the specific enthalpy"
    extends Sensors.BaseClasses.PartialFlowSensor;
    extends Modelica.Icons.RotationalSensor;
    Modelica.Blocks.Interfaces.RealOutput h_out(final quantity="SpecificEnergy",
                                                final unit="J/kg")
      "Specific enthalpy of the passing fluid" 
      annotation (Placement(transformation(
          origin={0,-110},
          extent={{-10,-10},{10,10}},
          rotation=270)));

    parameter Medium.MassFlowRate m_flow_small(min=0) = 1e-4
      "For bi-directional flow, specific enthalpy is regularized in the region |m_flow| < m_flow_small (m_flow_small > 0 required)"
      annotation(Dialog(tab="Advanced"));

  annotation (defaultComponentName="specificEnthalpy",
    Diagram(graphics={
          Line(points={{-100,0},{-70,0}}, color={0,128,255}),
          Line(points={{70,0},{100,0}}, color={0,128,255}),
          Line(points={{0,-70},{0,-100}}, color={0,0,127})}),
    Icon(graphics={
          Text(extent={{-150,80},{150,120}}, textString=
                                               "%name"),
          Text(
            extent={{168,-71},{8,-123}},
            lineColor={0,0,0},
            textString=
                 "h"),
          Line(points={{0,-70},{0,-100}}, color={0,0,127}),
          Line(points={{-100,0},{-70,0}}, color={0,128,255}),
          Line(points={{70,0},{100,0}}, color={0,128,255})}),
    Documentation(info="<HTML>
<p>
This component monitors the specific enthalpy of passing fluid. 
The sensor is ideal, i.e. it does not influence the fluid.
</p>
</HTML>
"));
  equation
    if allowFlowReversal then
       h_out = Modelica_Fluid.Utilities.regStep(port_a.m_flow, port_b.h_outflow, port_a.h_outflow, m_flow_small);
    else
       h_out = port_b.h_outflow;
    end if;
  end SpecificEnthalpyTwoPort;

  model SpecificEntropyOnePort "Ideal one port specific entropy sensor"
    extends Sensors.BaseClasses.PartialAbsoluteSensor;
    extends Modelica.Icons.RotationalSensor;
    Modelica.Blocks.Interfaces.RealOutput s(final quantity="SpecificEntropy",
                                            final unit="J/(kg.K)")
      "Specific entropy in port medium" 
      annotation (Placement(transformation(extent={{100,-10},{120,10}},
            rotation=0)));

  annotation (defaultComponentName="specificEntropy",
    Diagram(graphics={Line(points={{0,-70},{0,-100}}, color={0,0,127}), Line(
              points={{70,0},{100,0}}, color={0,0,127})}),
    Icon(graphics={
          Line(points={{0,-70},{0,-100}}, color={0,0,127}),
          Text(extent={{-150,80},{150,120}}, textString=
                                               "%name"),
          Text(
            extent={{170,-55},{10,-107}},
            lineColor={0,0,0},
            textString=
                 "s"),
          Line(points={{70,0},{100,0}}, color={0,0,127})}),
    Documentation(info="<HTML>
<p>
This component monitors the specific enthalphy of the medium in the flow
between fluid ports. The sensor is ideal, i.e., it does not influence the fluid.
</p>
</HTML>
"));
  equation
    s = Medium.specificEntropy(Medium.setState_phX(port.p, inflow(port.h_outflow), inflow(port.Xi_outflow)));
  end SpecificEntropyOnePort;

  model SpecificEntropyTwoPort "Ideal two port sensor for the specific entropy"
    extends Sensors.BaseClasses.PartialFlowSensor;
    extends Modelica.Icons.RotationalSensor;
    Modelica.Blocks.Interfaces.RealOutput s(final quantity="SpecificEntropy",
                                            final unit="J/(kg.K)")
      "Specific entropy of the passing fluid" 
      annotation (Placement(transformation(
          origin={0,-110},
          extent={{-10,-10},{10,10}},
          rotation=270)));
    parameter Medium.MassFlowRate m_flow_small(min=0) = 1e-4
      "For bi-directional flow, specific entropy is regularized in the region |m_flow| < m_flow_small (m_flow_small > 0 required)"
      annotation(Dialog(tab="Advanced"));

  annotation (defaultComponentName="specificEntropy",
    Diagram(graphics={
          Line(points={{-100,0},{-70,0}}, color={0,128,255}),
          Line(points={{70,0},{100,0}}, color={0,128,255}),
          Line(points={{0,-70},{0,-100}}, color={0,0,127})}),
    Icon(graphics={
          Text(extent={{-150,80},{150,120}}, textString=
                                               "%name"),
          Text(
            extent={{168,-71},{8,-123}},
            lineColor={0,0,0},
            textString=
                 "s"),
          Line(points={{0,-70},{0,-100}}, color={0,0,127}),
          Line(points={{-100,0},{-70,0}}, color={0,128,255}),
          Line(points={{70,0},{100,0}}, color={0,128,255})}),
    Documentation(info="<HTML>
<p>
This component monitors the specific entropy of the passing fluid. 
The sensor is ideal, i.e. it does not influence the fluid.
</p>
</HTML>
"));
  protected
    Medium.SpecificEntropy s_a_inflow
      "Specific entropy of inflowing fluid at port_a";
    Medium.SpecificEntropy s_b_inflow
      "Specific entropy of inflowing fluid at port_b or s_a_inflow, if uni-directional flow";
  equation
    if allowFlowReversal then
       s_a_inflow = Medium.specificEntropy(Medium.setState_phX(port_b.p, port_b.h_outflow, port_b.Xi_outflow));
       s_b_inflow = Medium.specificEntropy(Medium.setState_phX(port_a.p, port_a.h_outflow, port_a.Xi_outflow));
       s = Modelica_Fluid.Utilities.regStep(port_a.m_flow, s_a_inflow, s_b_inflow, m_flow_small);
    else
       s = Medium.specificEntropy(Medium.setState_phX(port_b.p, port_b.h_outflow, port_b.Xi_outflow));
       s_a_inflow = s;
       s_b_inflow = s;
    end if;
  end SpecificEntropyTwoPort;

  model MassFlowRate "Ideal sensor for mass flow rate"
    extends Sensors.BaseClasses.PartialFlowSensor;
    extends Modelica.Icons.RotationalSensor;
    Modelica.Blocks.Interfaces.RealOutput m_flow(quantity="MassFlowRate",
                                                 final unit="kg/s")
      "Mass flow rate from port_a to port_b" annotation (Placement(
          transformation(
          origin={0,-110},
          extent={{-10,-10},{10,10}},
          rotation=270)));

  annotation (
    Diagram(graphics={
          Line(points={{-100,0},{-70,0}}, color={0,128,255}),
          Line(points={{70,0},{100,0}}, color={0,128,255}),
          Line(points={{0,-70},{0,-100}}, color={0,0,127})}),
    Icon(graphics={
          Text(extent={{-150,80},{150,120}}, textString=
                                               "%name"),
          Line(points={{70,0},{100,0}}, color={0,128,255}),
          Text(
            extent={{178,-81},{18,-133}},
            lineColor={0,0,0},
            textString=
                 "m_flow"),
          Line(points={{0,-70},{0,-100}}, color={0,0,127}),
          Line(points={{-100,0},{-70,0}}, color={0,128,255})}),
    Documentation(info="<HTML>
<p>
This component monitors the mass flow rate flowing from port_a to port_b. 
The sensor is ideal, i.e., it does not influence the fluid.
</p>
</HTML>
"));
  equation
    m_flow = port_a.m_flow;
  end MassFlowRate;

  model VolumeFlowRate "Ideal sensor for volume flow rate"
    extends Sensors.BaseClasses.PartialFlowSensor;
    extends Modelica.Icons.RotationalSensor;
    Modelica.Blocks.Interfaces.RealOutput V_flow(final quantity="VolumeFlowRate",
                                                 final unit="m3/s")
      "Volume flow rate from port_a to port_b" 
      annotation (Placement(transformation(
          origin={0,-110},
          extent={{-10,-10},{10,10}},
          rotation=270)));
    parameter Medium.MassFlowRate m_flow_small(min=0) = 1e-4
      "For bi-directional flow, density is regularized in the region |m_flow| < m_flow_small (m_flow_small > 0 required)"
      annotation(Dialog(tab="Advanced"));

  annotation (
    Diagram(graphics={
          Line(points={{-100,0},{-70,0}}, color={0,128,255}),
          Line(points={{70,0},{100,0}}, color={0,128,255}),
          Line(points={{0,-70},{0,-100}}, color={0,0,127})}),
    Icon(graphics={
          Text(extent={{-150,80},{150,120}}, textString=
                                               "%name"),
          Text(
            extent={{188,-71},{28,-123}},
            lineColor={0,0,0},
            textString=
                 "V_flow"),
          Line(points={{0,-70},{0,-100}}, color={0,0,127}),
          Line(points={{-100,0},{-70,0}}, color={0,128,255}),
          Line(points={{70,0},{100,0}}, color={0,128,255})}),
    Documentation(info="<HTML>
<p>
This component monitors the volume flow rate flowing from port_a to port_b. 
The sensor is ideal, i.e. it does not influence the fluid.
</p>
</HTML>
"));
  protected
    Medium.Density d_a_inflow "Density of inflowing fluid at port_a";
    Medium.Density d_b_inflow
      "Density of inflowing fluid at port_b or d_a_inflow, if uni-directional flow";
    Medium.Density d "Density of the passing fluid";
  equation
    if allowFlowReversal then
       d_a_inflow = Medium.density(Medium.setState_phX(port_b.p, port_b.h_outflow, port_b.Xi_outflow));
       d_b_inflow = Medium.density(Medium.setState_phX(port_a.p, port_a.h_outflow, port_a.Xi_outflow));
       d = Modelica_Fluid.Utilities.regStep(port_a.m_flow, d_a_inflow, d_b_inflow, m_flow_small);
    else
       d = Medium.density(Medium.setState_phX(port_b.p, port_b.h_outflow, port_b.Xi_outflow));
       d_a_inflow = d;
       d_b_inflow = d;
    end if;
    V_flow = port_a.m_flow/d;
  end VolumeFlowRate;

  model RelativePressure "Ideal relative pressure sensor"
    extends Modelica.Icons.TranslationalSensor;
    replaceable package Medium = 
      Modelica.Media.Interfaces.PartialMedium "Medium in the sensor"  annotation (
        choicesAllMatching = true);

    Modelica_Fluid.Interfaces.FluidPort_a port_a(m_flow(min=0),
                                  redeclare package Medium = Medium) 
      annotation (Placement(transformation(extent={{-110,-10},{-90,10}},
            rotation=0)));
    Modelica_Fluid.Interfaces.FluidPort_b port_b(m_flow(min=0),
                                  redeclare package Medium = Medium) 
      annotation (Placement(transformation(extent={{110,-12},{90,8}}, rotation=
              0)));

    Modelica.Blocks.Interfaces.RealOutput p_rel(final quantity="Pressure",
                                                final unit="Pa",
                                                displayUnit="bar")
      "Relative pressure signal" annotation (Placement(transformation(
          origin={0,-90},
          extent={{10,-10},{-10,10}},
          rotation=90)));
    annotation (
      Icon(graphics={
          Line(points={{-100,0},{-70,0}}, color={0,127,255}),
          Line(points={{70,0},{100,0}}, color={0,127,255}),
          Line(points={{0,-30},{0,-80}}, color={0,0,127}),
          Text(extent={{-150,40},{150,80}}, textString=
                                              "%name"),
          Text(
            extent={{92,-62},{34,-122}},
            lineColor={0,0,0},
            textString=
                 "p_rel")}),
      Diagram(graphics={
          Line(points={{-100,0},{-70,0}}, color={0,127,255}),
          Line(points={{70,0},{100,0}}, color={0,127,255}),
          Line(points={{0,-30},{0,-80}}, color={0,0,127}),
          Text(
            extent={{64,-74},{32,-102}},
            lineColor={0,0,0},
            textString=
                 "p_rel")}),
      Documentation(info="<HTML>
<p>
The relative pressure \"port_a.p - port_b.p\" is determined between
the two ports of this component and is provided as output signal. The
sensor should be connected in parallel with other equipment, no flow
through the sensor is allowed.
</p>
</HTML>
"));
  equation
    // Zero flow equations for connectors
    port_a.m_flow = 0;
    port_b.m_flow = 0;

    // No contribution of specific quantities
    port_a.h_outflow = 0;
    port_b.h_outflow = 0;
    port_a.Xi_outflow = zeros(Medium.nXi);
    port_b.Xi_outflow = zeros(Medium.nXi);
    port_a.C_outflow  = zeros(Medium.nC);
    port_b.C_outflow  = zeros(Medium.nC);

    // Relative pressure
    p_rel = port_a.p - port_b.p;
  end RelativePressure;

  model RelativeTemperature "Ideal relative temperature sensor"
    extends Modelica.Icons.TranslationalSensor;
    replaceable package Medium = 
      Modelica.Media.Interfaces.PartialMedium "Medium in the sensor"  annotation (
        choicesAllMatching = true);
    Modelica_Fluid.Interfaces.FluidPort_a port_a(m_flow(min=0),
                                  redeclare package Medium = Medium) 
      annotation (Placement(transformation(extent={{-110,-10},{-90,10}},
            rotation=0)));
    Modelica_Fluid.Interfaces.FluidPort_b port_b(m_flow(min=0),
                                  redeclare package Medium = Medium) 
      annotation (Placement(transformation(extent={{110,-10},{90,10}}, rotation
            =0)));

    Modelica.Blocks.Interfaces.RealOutput T_rel(final quantity="ThermodynamicTemperature",
                                                final unit = "K", displayUnit = "degC", min=0)
      "Relative temperature signal"                                                                               annotation (Placement(
          transformation(
          origin={0,-90},
          extent={{10,-10},{-10,10}},
          rotation=90)));
    annotation (
      Icon(graphics={
          Line(points={{-100,0},{-70,0}}, color={0,127,255}),
          Line(points={{70,0},{100,0}}, color={0,127,255}),
          Line(points={{0,-30},{0,-80}}, color={0,0,127}),
          Text(extent={{-150,40},{150,80}}, textString=
                                              "%name"),
          Text(
            extent={{92,-62},{34,-122}},
            lineColor={0,0,0},
            textString=
                 "p_rel")}),
      Diagram(graphics={
          Line(points={{-100,0},{-70,0}}, color={0,127,255}),
          Line(points={{70,0},{100,0}}, color={0,127,255}),
          Line(points={{0,-30},{0,-80}}, color={0,0,127}),
          Text(
            extent={{64,-74},{32,-102}},
            lineColor={0,0,0},
            textString=
                 "T_rel")}),
      Documentation(info="<HTML>
<p>
The relative temperature \"port_a.T - port_b.T\" is determined between
the two ports of this component and is provided as output signal. The
sensor should be connected in parallel with other equipment, no flow
through the sensor is allowed.
</p>
</HTML>
"));
  equation
    // Zero flow equations for connectors
    port_a.m_flow = 0;
    port_b.m_flow = 0;

    // No contribution of specific quantities
    port_a.h_outflow = 0;
    port_b.h_outflow = 0;
    port_a.Xi_outflow = zeros(Medium.nXi);
    port_b.Xi_outflow = zeros(Medium.nXi);
    port_a.C_outflow  = zeros(Medium.nC);
    port_b.C_outflow  = zeros(Medium.nC);

    // Relative temperature
    T_rel = Medium.temperature(Medium.setState_phX(port_a.p, inflow(port_a.h_outflow), inflow(port_a.Xi_outflow))) -
            Medium.temperature(Medium.setState_phX(port_b.p, inflow(port_b.h_outflow), inflow(port_b.Xi_outflow)));
  end RelativeTemperature;

/*
  model SpecificInternalEnergy "Ideal specific internal energy sensor" 
    extends BaseClasses.Sensors.PartialSensor;
    extends Modelica.Icons.RotationalSensor;
    Medium.BaseProperties medium;
    Modelica.Blocks.Interfaces.RealOutput u(unit = "J/kg") 
      "Specific internal energy in port medium" annotation (extent=[100,-10; 120,10]);
    
  annotation (
    Diagram(
        Line(points=[70,0; 100,0], style(rgbcolor={0,0,127})),
        Line(points=[0,-70; 0,-100], style(color=69))),
    Icon(
        Line(points=[70,0; 100,0], style(rgbcolor={0,0,127})),
        Line(points=[0,-70; 0,-100], style(color=69)),
        Text(extent=[-126,160; 138,98], string="%name"),
        Text(
          extent=[212,-51; 52,-103],
          style(color=0),
          string="u")),
    Documentation(info="<HTML>
<p>
This component monitors the specific internal energy of the medium in the fluid port. The sensor is 
ideal, i.e., it does not influence the fluid.
</p>
</HTML>
"));
  equation 
    port.p   = medium.p;
    port.h   = medium.h;
    port.Xi = medium.Xi;
    u = medium.u;
  end SpecificInternalEnergy;
  
  model RelDensity "Ideal relative density sensor" 
    extends Interfaces.PartialRelativeSensor;
    extends Modelica.Icons.TranslationalSensor;
    Medium.BaseProperties medium_a;
    Medium.BaseProperties medium_b;
    Modelica.Blocks.Interfaces.RealOutput d_rel(redeclare type SignalType = 
          SI.Density) "Relative density signal" annotation (extent=[-10, -80; 10, -100], rotation=90);
    annotation (
      Icon(
        Line(points=[-100, 0; -70, 0], style(color=69)),
        Line(points=[70, 0; 100, 0], style(color=69)),
        Line(points=[0, -30; 0, -80], style(rgbcolor={0,0,127})),
        Text(extent=[-140, 94; 144, 34], string="%name"),
        Text(
          extent=[92, -62; 34, -122],
          string="d_rel",
          style(color=0))),
      Diagram(
        Line(points=[-100, 0; -70, 0], style(color=69)),
        Line(points=[70, 0; 100, 0], style(color=69)),
        Line(points=[0, -30; 0, -80], style(rgbcolor={0,0,127})),
        Text(
          extent=[64, -74; 32, -102],
          string="d_rel",
          style(color=0))),
      Documentation(info="<HTML>
<p>
The relative density \"d(port_a) - d(port_b)\" is determined between
the two ports of this component and is provided as output signal.
</p>
</HTML>
"));
  equation 
    port_a.p   = medium_a.p;
    port_a.h   = medium_a.h;
    port_a.Xi = medium_a.Xi;
    port_b.p   = medium_b.p;
    port_b.h   = medium_b.h;
    port_b.Xi = medium_b.Xi;
    
    d_rel = medium_a.d - medium_b.d;
  end RelDensity;
  
  model RelTemperature "Ideal relative temperature sensor" 
    extends Interfaces.PartialRelativeSensor;
    extends Modelica.Icons.TranslationalSensor;
    Medium.BaseProperties medium_a;
    Medium.BaseProperties medium_b;
    Modelica.Blocks.Interfaces.RealOutput T_rel(redeclare type SignalType = 
          SI.Temperature) "Relative temperature signal" annotation (extent=[-10, -80; 10, -100], rotation=90);
    annotation (
      Icon(
        Line(points=[-100, 0; -70, 0], style(color=69)),
        Line(points=[70, 0; 100, 0], style(color=69)),
        Line(points=[0, -30; 0, -80], style(rgbcolor={0,0,127})),
        Text(extent=[-140, 94; 144, 34], string="%name"),
        Text(
          extent=[92, -62; 34, -122],
          string="T_rel",
          style(color=0))),
      Diagram(
        Line(points=[-100, 0; -70, 0], style(color=69)),
        Line(points=[70, 0; 100, 0], style(color=69)),
        Line(points=[0, -30; 0, -80], style(rgbcolor={0,0,127})),
        Text(
          extent=[64, -74; 32, -102],
          string="T_rel",
          style(color=0))),
      Documentation(info="<HTML>
<p>
The relative temperature \"T(port_a) - T(port_b)\" is determined between
the two ports of this component and is provided as output signal.
</p>
</HTML>
"));
  equation 
    port_a.p   = medium_a.p;
    port_a.h   = medium_a.h;
    port_a.Xi = medium_a.Xi;
    port_b.p   = medium_b.p;
    port_b.h   = medium_b.h;
    port_b.Xi = medium_b.Xi;
    
    T_rel = medium_a.T - medium_b.T;
  end RelTemperature;
  
  model RelSpecificEnthalpy "Ideal relative specific enthalpy sensor" 
    extends Interfaces.PartialRelativeSensor;
    extends Modelica.Icons.TranslationalSensor;
    Modelica.Blocks.Interfaces.RealOutput h_rel(redeclare type SignalType = 
          SI.SpecificEnthalpy) "Relative specific enthalpy signal" annotation (extent=[-10, -80; 10, -100], rotation=90);
    annotation (
      Icon(
        Line(points=[-100, 0; -70, 0], style(color=69)),
        Line(points=[70, 0; 100, 0], style(color=69)),
        Line(points=[0, -30; 0, -80], style(rgbcolor={0,0,127})),
        Text(extent=[-140, 94; 144, 34], string="%name"),
        Text(
          extent=[92, -62; 34, -122],
          string="h_rel",
          style(color=0))),
      Diagram(
        Line(points=[-100, 0; -70, 0], style(color=69)),
        Line(points=[70, 0; 100, 0], style(color=69)),
        Line(points=[0, -30; 0, -80], style(rgbcolor={0,0,127})),
        Text(
          extent=[64, -74; 32, -102],
          string="h_rel",
          style(color=0))),
      Documentation(info="<HTML>
<p>
The relative specific enthalpy \"port_a.h - port_b.h\" is determined between
the two ports of this component and is provided as output signal.
</p>
</HTML>
"));
  equation 
    h_rel = port_a.h - port_b.h;
  end RelSpecificEnthalpy;
  
  model RelSpecificInternalEnergy 
    "Ideal relative specific internal energy sensor" 
    extends Interfaces.PartialRelativeSensor;
    extends Modelica.Icons.TranslationalSensor;
    Medium.BaseProperties medium_a;
    Medium.BaseProperties medium_b;
    Modelica.Blocks.Interfaces.RealOutput u_rel(redeclare type SignalType = 
          SI.SpecificEnergy) "Relative specific internal energy signal" annotation (extent=[-10, -80; 10, -100], rotation=90);
    annotation (
      Icon(
        Line(points=[-100, 0; -70, 0], style(color=69)),
        Line(points=[70, 0; 100, 0], style(color=69)),
        Line(points=[0, -30; 0, -80], style(rgbcolor={0,0,127})),
        Text(extent=[-140, 94; 144, 34], string="%name"),
        Text(
          extent=[92, -62; 34, -122],
          string="u_rel",
          style(color=0))),
      Diagram(
        Line(points=[-100, 0; -70, 0], style(color=69)),
        Line(points=[70, 0; 100, 0], style(color=69)),
        Line(points=[0, -30; 0, -80], style(rgbcolor={0,0,127})),
        Text(
          extent=[64, -74; 32, -102],
          string="u_rel",
          style(color=0))),
      Documentation(info="<HTML>
<p>
The relative specific internal energy \"u(port_a) - u(port_b)\" is determined between
the two ports of this component and is provided as output signal.
</p>
</HTML>
"));
  equation 
    port_a.p   = medium_a.p;
    port_a.h   = medium_a.h;
    port_a.Xi = medium_a.Xi;
    port_b.p   = medium_b.p;
    port_b.h   = medium_b.h;
    port_b.Xi = medium_b.Xi;
    
    u_rel = medium_a.u - medium_b.u;
  end RelSpecificInternalEnergy;
*/

  package BaseClasses
    extends Modelica_Fluid.Icons.BaseClassLibrary;

    partial model PartialAbsoluteSensor
      "Partial component to model a sensor that measures a potential variable"

      replaceable package Medium=Modelica.Media.Interfaces.PartialMedium
        "Medium in the sensor" 
        annotation(choicesAllMatching=true);

      Modelica_Fluid.Interfaces.FluidPort_a port(redeclare package Medium=Medium, m_flow(min=0)) 
        annotation (Placement(transformation(
            origin={0,-100},
            extent={{-10,-10},{10,10}},
            rotation=90)));

      annotation (Documentation(info="<html>
<p>
Partial component to model an <b>absolute sensor</b>. Can be used for pressure sensor models.
Use for other properties such as temperature or density is discouraged, because the enthalpy at the connector can have different meanings, depending on the connection topology. Use <tt>PartialFlowSensor</tt> instead.
as signal.
</p>
</html>"),
        Diagram(coordinateSystem(
            preserveAspectRatio=false,
            extent={{-100,-100},{100,100}},
            grid={1,1}), graphics),
        Icon(coordinateSystem(
            preserveAspectRatio=false,
            extent={{-100,-100},{100,100}},
            grid={1,1}), graphics));
    equation
      port.m_flow = 0;
      port.h_outflow = 0;
      port.Xi_outflow = zeros(Medium.nXi);
      port.C_outflow = zeros(Medium.nC);
    end PartialAbsoluteSensor;

    partial model PartialFlowSensor
      "Partial component to model sensors that measure flow properties"
      import Modelica.Constants;
      import Modelica_Fluid.Types.FlowDirection;

      replaceable package Medium=Modelica.Media.Interfaces.PartialMedium
        "Medium in the sensor" 
        annotation(choicesAllMatching=true);

      Modelica_Fluid.Interfaces.FluidPort_a port_a(
        redeclare package Medium=Medium,
        m_flow(min=if allowFlowReversal then -Constants.inf else 0.0)) 
        annotation (Placement(transformation(extent={{-110,-10},{-90,10}},
              rotation=0)));
      Modelica_Fluid.Interfaces.FluidPort_b port_b(
        redeclare package Medium=Medium,
        m_flow(max=if allowFlowReversal then +Constants.inf else 0.0)) 
        annotation (Placement(transformation(extent={{110,-10},{90,10}},
              rotation=0)));

      parameter FlowDirection.Temp flowDirection=FlowDirection.Bidirectional
        "Unidirectional (port_a -> port_b) or bidirectional flow component" 
         annotation(Dialog(tab="Advanced"));

      annotation (Documentation(info="<html>
<p>
Partial component to model a <b>sensor</b> that measures any intensive properties
of a flow, e.g., to get temperature or density in the flow
between fluid connectors.<br>
The model includes zero-volume balance equations. Sensor models inheriting from
this partial class should add a medium instance to calculate the measured property.
</p>
</html>"),
        Diagram(coordinateSystem(
            preserveAspectRatio=false,
            extent={{-100,-100},{100,100}},
            grid={1,1}), graphics));

    protected
      parameter Boolean allowFlowReversal=
         flowDirection == Modelica_Fluid.Types.FlowDirection.Bidirectional
        "= false, if flow only from port_a to port_b, otherwise reversing flow allowed"
         annotation(Evaluate=true, Hide=true);
    equation
      // mass balance
      0 = port_a.m_flow + port_b.m_flow;

      // momentum equation (no pressure loss)
      port_a.p = port_b.p;

      // isenthalpic state transformation (no storage and no loss of energy)
      port_a.h_outflow = inflow(port_b.h_outflow);
      port_b.h_outflow = inflow(port_a.h_outflow);

      port_a.Xi_outflow = inflow(port_b.Xi_outflow);
      port_b.Xi_outflow = inflow(port_a.Xi_outflow);

      port_a.C_outflow = inflow(port_b.C_outflow);
      port_b.C_outflow = inflow(port_a.C_outflow);
    end PartialFlowSensor;

  end BaseClasses;
end Sensors;
