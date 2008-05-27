within Modelica_Fluid.Test.TestComponents.Pipes;
model LumpedPipe
  import Modelica_Fluid;
  extends Modelica.Icons.Example;
  replaceable package Medium = 
      Modelica_Fluid.Media.Water.ConstantPropertyLiquidWater;
  //Modelica.Media.Water.StandardWater;

  Modelica_Fluid.Sources.FixedBoundary_pTX source(
    redeclare package Medium = Medium,
    p=5.0e5,
    T=300) annotation (Placement(transformation(extent={{-76,4},{-64,16}},
          rotation=0)));
  Modelica_Fluid.Pipes.LumpedPipe pipe1(
    redeclare package Medium = Medium,
    use_T_start=true,
    redeclare package WallFriction = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
    length=10,
    diameter=2.54e-2,
    initType=Modelica_Fluid.Types.Init.InitialValues,
    p_a_start=500000,
    p_b_start=100000,
    T_start=300)      annotation (Placement(transformation(extent={{-20,0},{0,
            20}}, rotation=0)));

  annotation (
    Diagram(graphics),
    experiment(StopTime=5),
    experimentSetupOutput(equdistant=false),
    Documentation(info="<html>
Simulation starts with both valves open. At t=1, valve 1 closes; at t=2 valve 2 closes, and the simulation fails.
</html>"));
  Modelica_Fluid.Sources.FixedBoundary_pTX sink(
    redeclare package Medium = Medium,
    p=100000,
    T=300)   annotation (Placement(transformation(extent={{56,4},{44,16}},
          rotation=0)));

  inner Ambient ambient annotation (Placement(transformation(extent={{-100,60},
            {-80,80}}, rotation=0)));
equation
  connect(source.port, pipe1.port_a) annotation (Line(points={{-64,10},{-20,10}}, 
        color={0,127,255}));
  connect(pipe1.port_b, sink.port) annotation (Line(
      points={{0,10},{44,10}},
      color={0,127,255},
      smooth=Smooth.None));
end LumpedPipe;
