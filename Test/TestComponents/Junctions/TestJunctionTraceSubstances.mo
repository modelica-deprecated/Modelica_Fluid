within Modelica_Fluid.Test.TestComponents.Junctions;
model TestJunctionTraceSubstances
  import Modelica_Fluid;
  extends Modelica.Icons.Example;
  package Medium=Modelica.Media.Air.MoistAir(extraPropertiesNames={"CO2"});
  Modelica_Fluid.Junctions.IdealTJunction junction(redeclare package Medium = 
        Medium)                             annotation (Placement(
        transformation(extent={{0,-30},{20,-10}},  rotation=0)));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
            -100},{100,100}}),
                      graphics={Text(
          extent={{-80,-58},{92,-120}},
          lineColor={255,0,0},
          textString=
              "Note: Multiport has been removed due to bug. See ticket #50.")}));
  Modelica_Fluid.Sources.FixedBoundary_pTX source2(
    T=278.15,
    redeclare package Medium = Medium,
    p=100000,
    X=Medium.X_default) 
    annotation (Placement(transformation(
        origin={90,-20},
        extent={{-10,-10},{10,10}},
        rotation=180)));
  Modelica_Fluid.Sources.FixedBoundary_pTX source3(
    T=283.15,
    redeclare package Medium = Medium,
    p=100000,
    X=Medium.X_default) 
    annotation (Placement(transformation(
        origin={-30,70},
        extent={{-10,-10},{10,10}},
        rotation=270)));
  inner Modelica_Fluid.System system 
    annotation (Placement(transformation(extent={{-100,80},{-80,100}}, rotation=
           0)));
  Modelica_Fluid.Sources.PrescribedBoundary_pTX source1(
    T=system.T_ambient,
    usePressureInput=true,
    redeclare package Medium = Medium,
    nPorts=2,
    useTraceInput=true,
    useCompositionInput=true,
    p=500000) 
    annotation (Placement(transformation(extent={{-68,-28},{-48,-8}},  rotation=
           0)));
  Modelica.Blocks.Sources.Ramp ramp(
    duration=1,
    height=-40,
    offset=1E5 + 20) 
                annotation (Placement(transformation(extent={{-100,24},{-80,44}},
          rotation=0)));
  Modelica_Fluid.PressureLosses.WallFrictionAndGravity pipe1(
    length=1,
    diameter=0.1,
    redeclare package Medium = Medium) annotation (Placement(transformation(
          extent={{-26,-30},{-6,-10}},rotation=0)));
  Modelica_Fluid.PressureLosses.WallFrictionAndGravity pipe2(
    length=1,
    diameter=0.1,
    redeclare package Medium = Medium) annotation (Placement(transformation(
          extent={{56,-30},{76,-10}}, rotation=0)));
  Modelica_Fluid.PressureLosses.WallFrictionAndGravity pipe3(
    length=1,
    diameter=0.1,
    redeclare package Medium = Medium) 
    annotation (Placement(transformation(
        origin={-30,30},
        extent={{-10,-10},{10,10}},
        rotation=90)));
  Modelica_Fluid.Junctions.TJunctionVolume junction1(
    redeclare package Medium = Medium,
    V=0.1,
    C_start={1E-4},
    dynamicsType=Modelica_Fluid.Types.Dynamics.SteadyStateMomentum,
    initType=Modelica_Fluid.Types.Init.InitialValues,
    X_start=Medium.X_default)               annotation (Placement(
        transformation(extent={{30,-30},{50,-10}}, rotation=0)));
  Modelica_Fluid.Sources.FixedBoundary_pTX source4(
    T=283.15,
    redeclare package Medium = Medium,
    p=100000,
    X=Medium.X_default) 
    annotation (Placement(transformation(
        origin={10,70},
        extent={{-10,-10},{10,10}},
        rotation=270)));
  Modelica_Fluid.PressureLosses.WallFrictionAndGravity pipe4(
    length=1,
    diameter=0.1,
    redeclare package Medium = Medium) 
    annotation (Placement(transformation(
        origin={10,30},
        extent={{-10,-10},{10,10}},
        rotation=90)));
  Modelica_Fluid.Sources.FixedBoundary_pTX source5(
    T=283.15,
    redeclare package Medium = Medium,
    p=100000,
    X=Medium.X_default) 
    annotation (Placement(transformation(
        origin={40,70},
        extent={{-10,-10},{10,10}},
        rotation=270)));
  Modelica_Fluid.PressureLosses.WallFrictionAndGravity pipe5(
    length=1,
    diameter=0.1,
    redeclare package Medium = Medium) 
    annotation (Placement(transformation(
        origin={40,30},
        extent={{-10,-10},{10,10}},
        rotation=90)));
  Modelica_Fluid.Sensors.TraceSubstancesOnePort traceSubstance2(redeclare
      package Medium = Medium)
    annotation (Placement(transformation(extent={{56,4},{76,24}})));
  Modelica_Fluid.Sensors.TraceSubstancesOnePort traceSubstance(redeclare
      package Medium = Medium)
    annotation (Placement(transformation(extent={{-70,16},{-50,36}})));
  Modelica.Blocks.Sources.Ramp C(duration=1, height=1.519E-3)
    "substance concentration, raising to 1000 PPM CO2" 
    annotation (Placement(transformation(extent={{-100,-60},{-80,-40}})));
  Modelica.Blocks.Sources.Constant X(k=0.02) 
    annotation (Placement(transformation(extent={{-100,-8},{-80,12}})));
  Modelica.Blocks.Sources.RealExpression X2(y=1 - X.y) "Concentration of X[2]"
    annotation (Placement(transformation(extent={{-100,-32},{-80,-12}})));
equation
  connect(ramp.y, source1.p_in) annotation (Line(points={{-79,34},{-74,34},{-74,
          -10},{-70,-10}},
        color={0,0,127}));
  connect(pipe1.port_b, junction.port_1) 
    annotation (Line(points={{-6,-20},{0,-20}}, color={0,127,255}));
  connect(pipe2.port_b, source2.ports[1]) 
                                      annotation (Line(points={{76,-20},{76,-20},
          {80,-20}},
        color={0,127,255}));
  connect(junction.port_2, junction1.port_1) annotation (Line(
      points={{20,-20},{30,-20}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(junction1.port_2, pipe2.port_a) annotation (Line(
      points={{50,-20},{56,-20}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(pipe3.port_b, source3.ports[1]) annotation (Line(
      points={{-30,40},{-30,60}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(pipe4.port_a, junction.port_3) annotation (Line(
      points={{10,20},{10,-10}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(source4.ports[1], pipe4.port_b) annotation (Line(
      points={{10,60},{10,40}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(source5.ports[1], pipe5.port_b) annotation (Line(
      points={{40,60},{40,40}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(pipe5.port_a, junction1.port_3) annotation (Line(
      points={{40,20},{40,-10}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(source1.ports[1], pipe3.port_a) annotation (Line(
      points={{-48,-16},{-30,-16},{-30,20}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(source1.ports[2], pipe1.port_a) annotation (Line(
      points={{-48,-20},{-26,-20}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(traceSubstance.port, pipe1.port_a) annotation (Line(
      points={{-60,16},{-44,16},{-44,-20},{-26,-20}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(C.y, source1.C_in[1]) annotation (Line(
      points={{-79,-50},{-76,-50},{-76,-26},{-70,-26}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(pipe2.port_a, traceSubstance2.port) annotation (Line(
      points={{56,-20},{54,-20},{54,4},{66,4}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(X.y, source1.X_in[1]) annotation (Line(
      points={{-79,2},{-76,2},{-76,-22},{-70,-22}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(X2.y, source1.X_in[2]) annotation (Line(
      points={{-79,-22},{-70,-22}},
      color={0,0,127},
      smooth=Smooth.None));
end TestJunctionTraceSubstances;
