within Modelica_Fluid.Test.TestComponents.Volumes;
model TestSweptVolume
  "Enclosed medium with fixed quantity in an adiabatic volume with varying size"
  import Modelica_Fluid;
  extends Modelica.Icons.Example;
  Modelica_Fluid.Volumes.SweptVolume sweptVolume(
    pistonCrossArea=0.05*0.05*Modelica.Constants.pi/4,
    clearance=0.05*0.05*Modelica.Constants.pi/4*0.03,
    redeclare package Medium = Modelica.Media.Air.DryAirNasa) 
    annotation (Placement(transformation(extent={{0,-10},{20,10}}, rotation=0)));
  Modelica.Mechanics.Translational.Sources.Position position 
    annotation (Placement(transformation(extent={{-20,-40},{0,-20}}, rotation=0)));
  Modelica.Blocks.Sources.Sine sine(
    phase=Modelica.Constants.pi/2,
    amplitude=0.1,
    offset=0.1) annotation (Placement(transformation(extent={{-80,-40},{-60,-20}},
          rotation=0)));
  inner Modelica_Fluid.System system  annotation (Placement(transformation(
          extent={{80,-40},{100,-20}}, rotation=0)));

  annotation (Diagram(graphics={Text(
          extent={{-100,80},{100,60}},
          lineColor={0,0,0},
          textString=
            "Enclosed medium with fixed quantity in an adiabatic volume with varying size")}),
    experiment(StopTime=10, Tolerance=1e-007),
    experimentSetupOutput);

equation
  connect(position.flange,   sweptVolume.flange) annotation (Line(
      points={{0,-30},{10,-30},{10,-10}},
      color={0,127,0},
      pattern=LinePattern.None,
      smooth=Smooth.None));
  connect(sine.y, position.s_ref) annotation (Line(
      points={{-59,-30},{-22,-30}},
      color={0,0,127},
      pattern=LinePattern.None,
      smooth=Smooth.None));

end TestSweptVolume;
