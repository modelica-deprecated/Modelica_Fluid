within Modelica_Fluid.Examples;
package TraceSubstances "Library demonstrating the usage of trace substances"
  extends Modelica.Icons.Library;
  model RoomCO2 "Demonstrates a room volume with CO2 accumulation"
    extends Modelica.Icons.Example;
    package Medium=Modelica.Media.Air.MoistAir(extraPropertiesNames={"CO2"});
    Modelica.Blocks.Sources.Constant C(k=0.1*1.519E-3)
      "substance concentration, raising to 1000 PPM CO2" 
      annotation (Placement(transformation(extent={{-100,-28},{-80,-8}})));
    Sources.FixedBoundary boundary4(redeclare package Medium = Medium) 
      annotation (Placement(transformation(extent={{80,-20},{60,0}})));
    Sensors.TraceSubstancesOnePort traceSubstanceVolume(redeclare package
        Medium = Medium) 
      annotation (Placement(transformation(extent={{0,20},{20,40}})));
    inner System system              annotation (Placement(transformation(extent={{60,60},
              {80,80}},          rotation=0)));
    Sources.PrescribedMassFlowRate_TX boundary1(
      useTraceInput=true,
      m_flow=100/1.2/3600*5,
      redeclare package Medium = Medium,
      nPorts=2) 
      annotation (Placement(transformation(extent={{-60,-20},{-40,0}})));
    Volumes.Volume volume(
      C_start={1.519E-3},
      V=100,
      redeclare package Medium = Medium,
      initType=Modelica_Fluid.Types.Init.InitialValues,
      nPorts=2) annotation (Placement(transformation(extent={{-20,0},{0,20}})));
    PressureLosses.WallFrictionAndGravity pipeFriction(
      redeclare package Medium = Medium,
      length=1,
      show_Re=true,
      diameter=0.15) 
      annotation (Placement(transformation(extent={{20,-20},{40,0}})));
    Sensors.TraceSubstancesOnePort traceSubstanceSource(redeclare package
        Medium = Medium) 
      annotation (Placement(transformation(extent={{-40,20},{-20,40}})));
  equation
    connect(C.y, boundary1.C_in[1]) annotation (Line(
        points={{-79,-18},{-60,-18}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(pipeFriction.port_b, boundary4.ports[1]) annotation (Line(
        points={{40,-10},{60,-10}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(volume.ports[2], pipeFriction.port_a) annotation (Line(
        points={{-10,-2},{-10,0},{-6,0},{-6,-10},{20,-10}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(traceSubstanceVolume.port, pipeFriction.port_a) annotation (Line(
        points={{10,20},{10,-10},{20,-10}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(boundary1.ports[1], volume.ports[1]) annotation (Line(
        points={{-40,-8},{-10,-8},{-10,2}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(boundary1.ports[2], traceSubstanceSource.port) annotation (Line(
        points={{-40,-12},{-30,-12},{-30,20}},
        color={0,127,255},
        smooth=Smooth.None));
    annotation (
      Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
              100,100}}), graphics),
      experiment(StopTime=7200),
      Documentation(info="<html>
This example consists of a volume with a carbon dioxide concentration that corresponds to about 1000 PPM.
There is an air stream with a carbon dioxide concentration of about 100 PPM that causes 
an air exchange rate of about 5 air changes per hour. After 2 hours of simulation time, the 
volume's carbon dioxide concentration is close to the concentration of the supply air.
</html>"));
  end RoomCO2;

end TraceSubstances;
