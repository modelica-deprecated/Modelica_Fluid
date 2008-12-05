within Modelica_Fluid.Test.TestComponents.Sources;
model TestSources "Test model for models in source package"
  import Modelica_Fluid;
  replaceable package Medium=Modelica.Media.Air.MoistAir(extraPropertiesNames={"CO2"});
  Modelica_Fluid.Sources.PrescribedBoundary_pTX boundary(redeclare package
      Medium = Medium,
    useTraceInput=true,
    usePressureInput=true) 
    annotation (Placement(transformation(extent={{-50,40},{-30,60}})));
  Modelica_Fluid.Sources.PrescribedBoundary_phX boundary1(redeclare package
      Medium = Medium,
    usePressureInput=true,
    useCompositionInput=false,
    useTraceInput=true) 
    annotation (Placement(transformation(extent={{-50,0},{-30,20}})));
  Modelica_Fluid.Sources.PrescribedMassFlowRate_TX boundary2(redeclare package
      Medium = Medium,
    useTraceInput=true,
    m_flow=0.1) 
    annotation (Placement(transformation(extent={{-50,-40},{-30,-20}})));
  Modelica_Fluid.Sources.PrescribedMassFlowRate_hX boundary3(redeclare package
      Medium = Medium,
    useTraceInput=true,
    m_flow=0.1) 
    annotation (Placement(transformation(extent={{-50,-80},{-30,-60}})));
  Modelica_Fluid.Sources.FixedBoundary boundary4(redeclare package Medium = 
        Medium) 
    annotation (Placement(transformation(extent={{80,40},{60,60}})));
  Modelica_Fluid.Sources.FixedBoundary_pTX boundary5(redeclare package Medium
      = Medium) annotation (Placement(transformation(extent={{80,0},{60,20}})));
  Modelica_Fluid.Sources.FixedBoundary_phX boundary6(nPorts=2, redeclare
      package Medium = Medium) 
    annotation (Placement(transformation(extent={{80,-60},{60,-40}})));
  Modelica_Fluid.Pipes.StaticPipe pipe(
    length=1,
    diameter=0.25,
    redeclare package Medium = Medium) 
    annotation (Placement(transformation(extent={{-20,40},{0,60}})));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,
            -100},{100,100}}), graphics), experiment(StopTime=2));
  Modelica_Fluid.Pipes.StaticPipe pipe1(
    length=1,
    diameter=0.25,
    redeclare package Medium = Medium) 
    annotation (Placement(transformation(extent={{0,0},{20,20}})));
  Modelica_Fluid.Pipes.StaticPipe pipe2(
    length=1,
    diameter=0.25,
    redeclare package Medium = Medium) 
    annotation (Placement(transformation(extent={{0,-40},{20,-20}})));
  Modelica_Fluid.Pipes.StaticPipe pipe3(
    length=1,
    diameter=0.25,
    redeclare package Medium = Medium) 
    annotation (Placement(transformation(extent={{0,-80},{20,-60}})));
  inner Modelica_Fluid.System system 
                                   annotation (Placement(transformation(extent={{60,70},
            {80,90}},          rotation=0)));
  Modelica.Blocks.Sources.Ramp C(duration=1, height=1.519E-3)
    "substance concentration, raising to 1000 PPM CO2" 
    annotation (Placement(transformation(extent={{-100,40},{-80,60}})));
  Modelica.Blocks.Sources.Ramp P(
    duration=2,
    height=-100,
    offset=101325 + 50) "Pressure" 
    annotation (Placement(transformation(extent={{-100,0},{-80,20}})));
  Modelica_Fluid.Sensors.TraceSubstancesOnePort traceSubstance(redeclare
      package Medium = Medium) 
    annotation (Placement(transformation(extent={{0,70},{20,90}})));
  Modelica_Fluid.Sensors.TraceSubstancesTwoPort traceSubstance1(redeclare
      package Medium = Medium) 
    annotation (Placement(transformation(extent={{22,40},{42,60}})));
  Modelica_Fluid.Sensors.TraceSubstancesOnePort traceSubstance2(redeclare
      package Medium = Medium) 
    annotation (Placement(transformation(extent={{-22,14},{-2,34}})));
  Modelica_Fluid.Sensors.TraceSubstancesOnePort traceSubstance3(redeclare
      package Medium = Medium) 
    annotation (Placement(transformation(extent={{30,-20},{50,0}})));
  Modelica_Fluid.Sensors.TraceSubstancesOnePort traceSubstance4(redeclare
      package Medium = Medium) 
    annotation (Placement(transformation(extent={{-20,-64},{0,-44}})));
equation
  connect(boundary.ports[1], pipe.port_a) annotation (Line(
      points={{-30,50},{-20,50}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(boundary1.ports[1], pipe1.port_a) annotation (Line(
      points={{-30,10},{0,10}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(pipe1.port_b, boundary5.ports[1]) annotation (Line(
      points={{20,10},{60,10}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(boundary2.ports[1], pipe2.port_a) annotation (Line(
      points={{-30,-30},{0,-30}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(pipe2.port_b, boundary6.ports[1]) annotation (Line(
      points={{20,-30},{46,-30},{46,-48},{60,-48}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(boundary3.ports[1], pipe3.port_a) annotation (Line(
      points={{-30,-70},{0,-70}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(pipe3.port_b, boundary6.ports[2]) annotation (Line(
      points={{20,-70},{46,-70},{46,-52},{60,-52}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(C.y, boundary.C_in[1]) annotation (Line(
      points={{-79,50},{-70,50},{-70,42},{-52,42}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(P.y, boundary1.p_in) annotation (Line(
      points={{-79,10},{-72,10},{-72,18},{-52,18}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(P.y, boundary.p_in) annotation (Line(
      points={{-79,10},{-72,10},{-72,58},{-52,58}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(C.y, boundary1.C_in[1]) annotation (Line(
      points={{-79,50},{-70,50},{-70,2},{-52,2}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(pipe.port_b, traceSubstance1.port_a) annotation (Line(
      points={{0,50},{22,50}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(traceSubstance1.port_b, boundary4.ports[1]) annotation (Line(
      points={{42,50},{60,50}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(pipe3.port_a, traceSubstance4.port) annotation (Line(
      points={{0,-70},{-6,-70},{-6,-64},{-10,-64}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(traceSubstance.port, pipe.port_b) annotation (Line(
      points={{10,70},{6,70},{6,50},{0,50}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(pipe1.port_a, traceSubstance2.port) annotation (Line(
      points={{0,10},{-6,10},{-6,14},{-12,14}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(pipe2.port_b, traceSubstance3.port) annotation (Line(
      points={{20,-30},{31,-30},{31,-20},{40,-20}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(C.y, boundary2.C_in[1]) annotation (Line(
      points={{-79,50},{-64,50},{-64,-38},{-50,-38}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(C.y, boundary3.C_in[1]) annotation (Line(
      points={{-79,50},{-64,50},{-64,-78},{-50,-78}},
      color={0,0,127},
      smooth=Smooth.None));
end TestSources;
