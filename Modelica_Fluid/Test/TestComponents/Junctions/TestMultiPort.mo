within Modelica_Fluid.Test.TestComponents.Junctions;
model TestMultiPort
  import Modelica_Fluid;
  extends Modelica.Icons.Example;

  Modelica_Fluid.Junctions.MultiPort multiPort(redeclare package Medium = 
        Modelica.Media.Air.DryAirNasa, nPorts_b=2) 
                                            annotation (Placement(
        transformation(extent={{-28,-30},{-20,-10}},
                                                   rotation=0)));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,
            -100},{100,100}}),
                      graphics));
  Modelica_Fluid.Sources.FixedBoundary_pTX source2(
    T=278.15,
    p=5e5,
    redeclare package Medium = Modelica.Media.Air.DryAirNasa) 
    annotation (Placement(transformation(
        origin={80,-20},
        extent={{-10,-10},{10,10}},
        rotation=180)));
  Modelica_Fluid.Sources.FixedBoundary_pTX source3(
    T=283.15,
    p=2e5,
    redeclare package Medium = Modelica.Media.Air.DryAirNasa) 
    annotation (Placement(transformation(
        origin={10,70},
        extent={{-10,-10},{10,10}},
        rotation=270)));
  inner Modelica_Fluid.System system 
    annotation (Placement(transformation(extent={{-100,80},{-80,100}}, rotation=
           0)));
  Modelica_Fluid.Sources.PrescribedBoundary_pTX source1(          p=5e5,
      redeclare package Medium = Modelica.Media.Air.DryAirNasa,
    T=system.T_ambient,
    usePressureInput=true) 
    annotation (Placement(transformation(extent={{-60,-30},{-40,-10}}, rotation=
           0)));
  Modelica.Blocks.Sources.Ramp ramp(
    duration=1,
    height=-6.5e5,
    offset=7e5) annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=270,
        origin={-80,10})));
  Modelica_Fluid.PressureLosses.WallFrictionAndGravity pipe1(redeclare package
      Medium = 
        Modelica.Media.Air.DryAirNasa,
    length=1,
    diameter=0.1)                      annotation (Placement(transformation(
          extent={{40,-30},{60,-10}}, rotation=0)));
  Modelica_Fluid.PressureLosses.WallFrictionAndGravity pipe2(redeclare package
      Medium = 
        Modelica.Media.Air.DryAirNasa,
    length=1,
    diameter=0.1) 
    annotation (Placement(transformation(
        origin={10,24},
        extent={{-10,-10},{10,10}},
        rotation=90)));
equation
  connect(ramp.y, source1.p_in) annotation (Line(points={{-80,-1},{-80,-1},{-80,
          -14},{-62,-14}},
        color={0,0,127}));
  connect(pipe1.port_b, source2.ports[1]) 
                                      annotation (Line(points={{60,-20},{70,-20}},
        color={0,127,255}));
  connect(pipe2.port_b, source3.ports[1]) 
                                      annotation (Line(points={{10,34},{10,47},
          {10,60}},         color={0,127,255}));
  connect(multiPort.ports_b[1], pipe2.port_a) 
                                             annotation (Line(
      points={{-20,-18},{10,-18},{10,14}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(multiPort.ports_b[2], pipe1.port_a) 
                                             annotation (Line(
      points={{-20,-22},{35,-22},{35,-20},{40,-20}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(source1.ports[1], multiPort.port_a) annotation (Line(
      points={{-40,-20},{-28,-20}},
      color={0,127,255},
      smooth=Smooth.None));
end TestMultiPort;