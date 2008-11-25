within Modelica_Fluid.Test.TestComponents.Junctions;
model TestGeneircJunction
  import Modelica_Fluid;
  extends Modelica.Icons.Example;

  Modelica_Fluid.Junctions.GenericJunction junction(             redeclare
      package Medium = 
        Modelica.Media.Air.DryAirNasa, V=20e-6,
    initType=Modelica_Fluid.Types.Init.SteadyState,
    n_a=1,
    n_b=2,
    modelStructure=Modelica_Fluid.Types.ModelStructure.avb,
    dp_nominal=100000,
    m_flow_nominal=0.01,
    p_start=100000)                         annotation (Placement(
        transformation(extent={{10,-30},{30,-10}}, rotation=0)));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,
            -100},{100,100}}),
                      graphics));
  Modelica_Fluid.Sources.FixedBoundary_pTX source2(
    T=278.15,
    p=5e5,
    redeclare package Medium = Modelica.Media.Air.DryAirNasa) 
    annotation (Placement(transformation(
        origin={90,-20},
        extent={{-10,-10},{10,10}},
        rotation=180)));
  Modelica_Fluid.Sources.FixedBoundary_pTX source3(
    T=283.15,
    p=2e5,
    redeclare package Medium = Modelica.Media.Air.DryAirNasa) 
    annotation (Placement(transformation(
        origin={40,70},
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
  Modelica_Fluid.PressureLosses.WallFrictionAndGravity pipe(redeclare package
      Medium = 
        Modelica.Media.Air.DryAirNasa,
    length=1,
    diameter=0.1)                      annotation (Placement(transformation(
          extent={{-30,-30},{-10,-10}},
                                      rotation=0)));
  Modelica_Fluid.PressureLosses.WallFrictionAndGravity pipe1(redeclare package
      Medium = 
        Modelica.Media.Air.DryAirNasa,
    length=1,
    diameter=0.1)                      annotation (Placement(transformation(
          extent={{50,-30},{70,-10}}, rotation=0)));
  Modelica_Fluid.PressureLosses.WallFrictionAndGravity pipe2(redeclare package
      Medium = 
        Modelica.Media.Air.DryAirNasa,
    length=1,
    diameter=0.1) 
    annotation (Placement(transformation(
        origin={40,24},
        extent={{-10,-10},{10,10}},
        rotation=90)));
equation
  connect(ramp.y, source1.p_in) annotation (Line(points={{-80,-1},{-80,-1},{-80,
          -14},{-62,-14}},
        color={0,0,127}));
  connect(source1.port, pipe.port_a) annotation (Line(points={{-40,-20},{-30,
          -20}}, color={0,127,255}));
  connect(pipe1.port_b, source2.port) annotation (Line(points={{70,-20},{80,-20}},
        color={0,127,255}));
  connect(pipe2.port_b, source3.port) annotation (Line(points={{40,34},{40,47},
          {40,60}},         color={0,127,255}));
  connect(junction.ports_a[1], pipe.port_b) annotation (Line(
      points={{10,-20},{-10,-20}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(junction.ports_b[1], pipe2.port_a) annotation (Line(
      points={{30,-18},{40,-18},{40,14}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(junction.ports_b[2], pipe1.port_a) annotation (Line(
      points={{30,-22},{45,-22},{45,-20},{50,-20}},
      color={0,127,255},
      smooth=Smooth.None));
end TestGeneircJunction;
