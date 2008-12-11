within Modelica_Fluid.Test.TestComponents.Fittings;
model TestJunctionVolume
  extends Modelica.Icons.Example;

  Modelica_Fluid.Fittings.TJunctionVolume junction(             redeclare
      package Medium = 
        Modelica.Media.Air.DryAirNasa, V=20e-6,
    initType=Modelica_Fluid.Types.Init.SteadyState,
    p_start=100000)                         annotation (Placement(
        transformation(extent={{20,-30},{40,-10}}, rotation=0)));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
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
        origin={30,70},
        extent={{-10,-10},{10,10}},
        rotation=270)));
  inner Modelica_Fluid.System system 
    annotation (Placement(transformation(extent={{-100,80},{-80,100}}, rotation=
           0)));
  Modelica_Fluid.Sources.PrescribedBoundary_pTX source1(          p=5e5,
      redeclare package Medium = Modelica.Media.Air.DryAirNasa,
    T=system.T_ambient,
    usePressureInput=true) 
    annotation (Placement(transformation(extent={{-40,-30},{-20,-10}}, rotation=
           0)));
  Modelica.Blocks.Sources.Ramp ramp(
    duration=1,
    height=-6.5e5,
    offset=7e5) annotation (Placement(transformation(extent={{-90,-24},{-70,-4}},
          rotation=0)));
  Modelica_Fluid.Fittings.WallFrictionAndGravity pipe(      redeclare package
      Medium = 
        Modelica.Media.Air.DryAirNasa,
    length=1,
    diameter=0.1)                      annotation (Placement(transformation(
          extent={{-12,-30},{8,-10}}, rotation=0)));
  Modelica_Fluid.Fittings.WallFrictionAndGravity pipe1(      redeclare package
      Medium = 
        Modelica.Media.Air.DryAirNasa,
    length=1,
    diameter=0.1)                      annotation (Placement(transformation(
          extent={{50,-30},{70,-10}}, rotation=0)));
  Modelica_Fluid.Fittings.WallFrictionAndGravity pipe2(      redeclare package
      Medium = 
        Modelica.Media.Air.DryAirNasa,
    length=1,
    diameter=0.1) 
    annotation (Placement(transformation(
        origin={30,24},
        extent={{-10,-10},{10,10}},
        rotation=90)));
equation
  connect(ramp.y, source1.p_in) annotation (Line(points={{-69,-14},{-42,-14}},
        color={0,0,127}));
  connect(source1.ports[1], pipe.port_a) 
                                     annotation (Line(points={{-20,-20},{-12,
          -20}}, color={0,127,255}));
  connect(pipe.port_b, junction.port_1) 
    annotation (Line(points={{8,-20},{20,-20}}, color={0,127,255}));
  connect(pipe1.port_b, source2.ports[1]) 
                                      annotation (Line(points={{70,-20},{80,-20}},
        color={0,127,255}));
  connect(junction.port_2, pipe1.port_a) annotation (Line(points={{40,-20},{50,
          -20}}, color={0,127,255}));
  connect(pipe2.port_b, source3.ports[1]) 
                                      annotation (Line(points={{30,34},{30,47},
          {30,60},{30,60}}, color={0,127,255}));
  connect(pipe2.port_a, junction.port_3) 
    annotation (Line(points={{30,14},{30,-10}}, color={0,127,255}));
end TestJunctionVolume;
