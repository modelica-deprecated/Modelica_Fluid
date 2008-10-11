within Modelica_Fluid;
model System "System properties and default values"

  parameter Modelica.Media.Interfaces.PartialMedium.AbsolutePressure p_ambient=
                101325 "Default ambient pressure" 
      annotation(Dialog(group="Defaults"));
  parameter Modelica.Media.Interfaces.PartialMedium.Temperature T_ambient=
                293.15 "Default ambient temperature" 
      annotation(Dialog(group="Defaults"));
  parameter SI.Acceleration g=Modelica.Constants.g_n
    "Constant gravity acceleration"                                                  annotation(Dialog(group="Environment"));
  parameter Modelica_Fluid.Types.FlowDirection flowDirection=
      Modelica_Fluid.Types.FlowDirection.Bidirectional
    "Default for unidirectional (port_a -> port_b) or bidirectional flow" 
     annotation(Dialog(tab="Advanced"));
  parameter Types.Init initType=
            Types.Init.NoInit "Default initialization option" 
    annotation(Evaluate=true, Dialog(tab = "Initialization"));

  annotation (
    preferedView="info",
    defaultComponentName="system",
    defaultComponentPrefixes="inner",
    missingInnerMessage="Your model is using an outer \"system\" component. An inner \"system\" component is not defined. For simulation drag Modelica_Fluid.System into your model to specify system properties.",
    Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
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
 
</html>"));

/*
    Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape dummyShape 
    "Just temporarily, to force Dymola to open an animation window (only then animation setup is available for diagram animation)"
      annotation (extent=[-60,20; -40,40]);
*/
end System;
