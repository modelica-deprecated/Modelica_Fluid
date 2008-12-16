within Modelica_Fluid.Test.TestComponents.Sensors;
model TestSpecificEntropy
  import Modelica_Fluid;
  annotation (
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
            100,100}}),
            graphics),
    experiment(Tolerance=1e-006),
    experimentSetupOutput);
  inner Modelica_Fluid.System system  annotation (Placement(transformation(
          extent={{-100,-100},{-80,-80}}, rotation=0)));
  Modelica_Fluid.Sources.Boundary_ph boundary_prescribed_1(
    use_h_in=true,
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    p=system.p_ambient) annotation (Placement(transformation(extent={{
            -40,10},{-20,30}}, rotation=0)));
  Modelica_Fluid.Sensors.SpecificEntropy specificEntropy(redeclare package
      Medium =         Modelica.Media.Water.StandardWater) 
    annotation (Placement(transformation(extent={{-10,20},{10,40}}, rotation=0)));
  Modelica_Fluid.Sources.Boundary_ph boundary_prescribed_2(
    use_h_in=true,
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    p=system.p_ambient) annotation (Placement(transformation(extent={{
            -40,-30},{-20,-10}}, rotation=0)));
  Modelica_Fluid.Sensors.SpecificEntropyTwoPort specificEntropy1(redeclare
      package Medium = Modelica.Media.Water.StandardWater) 
    annotation (Placement(transformation(extent={{-10,-30},{10,-10}}, rotation=
            0)));
  Modelica.Blocks.Sources.Sine sine1(amplitude=1600e3, offset=1800e3) 
                                    annotation (Placement(transformation(extent=
           {{-80,-10},{-60,10}}, rotation=0)));
equation
  connect(boundary_prescribed_1.ports[1], specificEntropy.port) 
                                                            annotation (Line(
      points={{-20,20},{0,20}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(sine1.y, boundary_prescribed_1.h_in) annotation (Line(
      points={{-59,0},{-52,0},{-52,20},{-42,20}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(sine1.y, boundary_prescribed_2.h_in) annotation (Line(
      points={{-59,0},{-52,0},{-52,-20},{-42,-20}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(boundary_prescribed_2.ports[1], specificEntropy1.port_a) 
                                                               annotation (Line(
      points={{-20,-20},{-10,-20}},
      color={0,127,255},
      smooth=Smooth.None));
end TestSpecificEntropy;
