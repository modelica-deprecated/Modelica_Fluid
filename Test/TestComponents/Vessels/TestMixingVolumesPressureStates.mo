within Modelica_Fluid.Test.TestComponents.Vessels;
model TestMixingVolumesPressureStates
  "Test case where in one of the mixing volumes a pressure state appears"
  import Modelica_Fluid;
  extends Modelica.Icons.Example;
  package Medium = Modelica.Media.Water.StandardWater;
  annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,
            -100},{100,100}}),
                      graphics),
                       experiment(StopTime=10));
  Modelica_Fluid.Vessels.Volume mixingVolume1(
    V=1e-3,
    redeclare package Medium = Medium,
    p_start=system.p_ambient,
    use_T_start=true,
    T_start=system.T_ambient,
    nPorts=2,
    portDiameters={0,0},
    neglectPortDiameters=true,
    energyDynamics=Modelica_Fluid.Types.Dynamics.SteadyStateInitial,
    massDynamics=Modelica_Fluid.Types.Dynamics.SteadyStateInitial) 
                 annotation (Placement(transformation(extent={{-30,38},{-10,58}},
          rotation=0)));

  Modelica_Fluid.Sources.MassFlowSource_h flowSource2(
    m_flow=1,
    h=2e5,
    redeclare package Medium = Medium) 
                   annotation (Placement(transformation(extent={{-100,30},{-80,
            50}}, rotation=0)));
  Modelica_Fluid.Vessels.Volume mixingVolume2(
    V=1e-3,
    p_start=system.p_ambient,
    use_T_start=false,
    h_start=1e5,
    redeclare package Medium = Medium,
    nPorts=2,
    portDiameters={0,0}) 
                 annotation (Placement(transformation(extent={{10,38},{30,58}},
          rotation=0)));
  Modelica_Fluid.Sensors.Temperature Tmix_in(
                                         redeclare package Medium = Medium) 
    annotation (Placement(transformation(extent={{-50,70},{-30,90}}, rotation=0)));
  Modelica_Fluid.Sensors.Temperature Tmix_out(
                                          redeclare package Medium = Medium) 
    annotation (Placement(transformation(extent={{30,68},{50,88}}, rotation=0)));
  Modelica_Fluid.Sources.Boundary_ph sink2(             redeclare package
      Medium = Medium,
    h=Medium.h_default,
    p=101325) 
    annotation (Placement(transformation(extent={{100,30},{80,50}}, rotation=0)));
  inner Modelica_Fluid.System system 
    annotation (Placement(transformation(extent={{-100,-100},{-80,-80}},
          rotation=0)));
  Modelica_Fluid.Pipes.BaseClasses.WallFriction.TestWallFrictionAndGravity
    simpleGenericOrifice2(
    redeclare package Medium = Medium,
    diameter=0.2,
    length=1,
    redeclare package WallFriction = 
        Modelica_Fluid.Pipes.BaseClasses.WallFriction.Laminar) 
                                annotation (Placement(transformation(extent={{
            50,30},{70,50}}, rotation=0)));
  Modelica_Fluid.Pipes.BaseClasses.WallFriction.TestWallFrictionAndGravity
    simpleGenericOrifice1(
    redeclare package Medium = Medium,
    diameter=0.2,
    length=1,
    redeclare package WallFriction = 
        Modelica_Fluid.Pipes.BaseClasses.WallFriction.Laminar) 
                                annotation (Placement(transformation(extent={{
            -70,30},{-50,50}}, rotation=0)));
equation
  connect(simpleGenericOrifice2.port_b,sink2. ports[1]) 
                                                    annotation (Line(
      points={{70,40},{80,40}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(flowSource2.ports[1],simpleGenericOrifice1.port_a) 
                                                          annotation (Line(
      points={{-80,40},{-70,40}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(simpleGenericOrifice1.port_b, Tmix_in.port) annotation (Line(
      points={{-50,40},{-40,40},{-40,70}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(Tmix_out.port, simpleGenericOrifice2.port_a) annotation (Line(
      points={{40,68},{40,40},{50,40}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(simpleGenericOrifice1.port_b, mixingVolume1.ports[1]) annotation (
      Line(
      points={{-50,40},{-20,40}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(mixingVolume1.ports[2], mixingVolume2.ports[2]) annotation (Line(
      points={{-20,36},{-20,34},{20,34},{20,36}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(mixingVolume2.ports[1], simpleGenericOrifice2.port_a) annotation (
      Line(
      points={{20,40},{50,40}},
      color={0,127,255},
      smooth=Smooth.None));
end TestMixingVolumesPressureStates;
