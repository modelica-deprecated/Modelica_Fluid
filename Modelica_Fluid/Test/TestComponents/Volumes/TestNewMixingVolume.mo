within Modelica_Fluid.Test.TestComponents.Volumes;
model TestNewMixingVolume
  extends Modelica.Icons.Example;
  Modelica_Fluid.Volumes.MixingVolume Volume(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    initType=Modelica_Fluid.Types.Init.SteadyState,
    V=1,
    use_T_start=false,
    h_start=3e6) 
         annotation (Placement(transformation(extent={{-42,0},{-22,20}},
          rotation=0)));
  Modelica_Fluid.Sources.PrescribedMassFlowRate_hX FlowSource(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    m_flow=1,
    h=3e6) annotation (Placement(transformation(extent={{-82,0},{-62,20}},
          rotation=0)));
  Modelica_Fluid.Sources.FixedBoundary_pTX Sink(
                                           redeclare package Medium = 
        Modelica.Media.Water.StandardWater, p=101325,
    T=system.T_ambient) 
    annotation (Placement(transformation(extent={{60,0},{40,20}}, rotation=0)));
  Modelica_Fluid.ControlValves.ValveLinear Valve(
                                             redeclare package Medium = 
        Modelica.Media.Water.StandardWater,
    dp_nominal=10000,
    m_flow_nominal=0.1)                                    annotation (Placement(
        transformation(extent={{2,0},{22,20}}, rotation=0)));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
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
  connect(FlowSource.port, Volume.port_a) annotation (Line(points={{-62,10},{
          -42.2,10}}, color={0,127,255}));
  connect(Volume.port_b, Valve.port_a) 
    annotation (Line(points={{-22,10},{2,10}}, color={0,127,255}));
  connect(Valve.port_b, Sink.port) 
    annotation (Line(points={{22,10},{40,10}}, color={0,127,255}));
  connect(Step1.y, Valve.stemPosition) annotation (Line(points={{-19,58},{12,58},{12,
          19}}, color={0,0,127}));
end TestNewMixingVolume;
