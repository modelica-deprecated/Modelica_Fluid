within Modelica_Fluid;
model Ambient "Ambient field component"

  parameter Modelica.Media.Interfaces.PartialMedium.AbsolutePressure
    default_p_ambient = 101325 "Default ambient pressure" 
      annotation(Dialog(group="Defaults"));
  parameter Modelica.Media.Interfaces.PartialMedium.Temperature
    default_T_ambient = 293.15 "Default ambient temperature" 
      annotation(Dialog(group="Defaults"));
  parameter SI.Acceleration g=Modelica.Constants.g_n
    "Constant gravity acceleration"                                                  annotation(Dialog(group="Environment"));

  annotation (
    preferedView="info",
    defaultComponentName="ambient",
    defaultComponentPrefixes="inner",
    missingInnerMessage="Your model is using an outer \"ambient\" component. An inner \"ambient\" component is not defined. For simulation drag Modelica_Fluid.Components.Ambient into your model and specify ambient conditions.",
    Icon(graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-150,150},{150,110}},
          lineColor={0,0,255},
          textString=
               "%name"),
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
          textString=
               "g"),
        Text(
          extent={{-90,82},{74,50}},
          lineColor={0,0,0},
          textString=
                 "defaults")}),
    Diagram(graphics),
    Documentation(info="<html>
 
</html>"));

/*
    Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape dummyShape 
    "Just temporarily, to force Dymola to open an animation window (only then animation setup is available for diagram animation)"
      annotation (extent=[-60,20; -40,40]);
*/
end Ambient;
