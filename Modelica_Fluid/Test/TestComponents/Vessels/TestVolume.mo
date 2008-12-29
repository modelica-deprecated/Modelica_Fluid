within Modelica_Fluid.Test.TestComponents.Vessels;
model TestVolume
  extends Modelica.Icons.Example;
  Modelica_Fluid.Vessels.Volume Volume(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    V=1,
    use_T_start=false,
    h_start=3e6,
    nPorts=2,
    use_portDiameters=false,
    energyDynamics=Modelica_Fluid.Types.Dynamics.SteadyState,
    massDynamics=Modelica_Fluid.Types.Dynamics.SteadyState) 
         annotation (Placement(transformation(extent={{-40,14},{-20,34}},
          rotation=0)));
  Modelica_Fluid.Sources.MassFlowSource_h FlowSource(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    m_flow=1,
    h=3e6) annotation (Placement(transformation(extent={{-82,0},{-62,20}},
          rotation=0)));
  Modelica_Fluid.Sources.Boundary_pT Sink( redeclare package Medium = 
        Modelica.Media.Water.StandardWater, p=101325,
    T=system.T_ambient) 
    annotation (Placement(transformation(extent={{60,0},{40,20}}, rotation=0)));
  Modelica_Fluid.Valves.ValveLinear Valve(   redeclare package Medium = 
        Modelica.Media.Water.StandardWater,
    dp_nominal=10000,
    m_flow_nominal=0.1)                                    annotation (Placement(
        transformation(extent={{2,0},{22,20}}, rotation=0)));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
            -100},{100,100}}),
                      graphics),
                       experiment(StopTime=5));
  Modelica.Blocks.Sources.Step Step1(
    startTime=1,
    height=-0.5,
    offset=1) annotation (Placement(transformation(extent={{-40,48},{-20,68}},
          rotation=0)));
  inner Modelica_Fluid.System system 
    annotation (Placement(transformation(extent={{-100,-100},{-80,-80}},
          rotation=0)));
equation
  connect(Valve.port_b, Sink.ports[1]) 
    annotation (Line(points={{22,10},{40,10}}, color={0,127,255}));
  connect(Step1.y, Valve.opening) annotation (Line(points={{-19,58},{12,58},{
          12,18}},
                color={0,0,127}));
  connect(FlowSource.ports[1], Volume.ports[1]) annotation (Line(
      points={{-62,10},{-30,10},{-30,16}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(Volume.ports[2], Valve.port_a) annotation (Line(
      points={{-30,12},{-30,10},{2,10}},
      color={0,127,255},
      smooth=Smooth.None));
end TestVolume;
