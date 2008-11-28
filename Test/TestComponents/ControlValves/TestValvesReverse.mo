within Modelica_Fluid.Test.TestComponents.ControlValves;
model TestValvesReverse "Test case for valves with reverse and zero flow"
  extends Modelica.Icons.Example;
  package Medium = Modelica.Media.Water.StandardWater;
  Modelica_Fluid.Sources.FixedBoundary_pTX SourceP1(
    T=system.T_ambient,
    redeclare package Medium = 
        Modelica.Media.Water.StandardWaterOnePhase,
    p=1000000) 
  annotation (Placement(transformation(extent={{-100,26},{-80,46}}, rotation=0)));
  Modelica_Fluid.Sources.FixedBoundary_pTX SourceP2(
    T=system.T_ambient,
    redeclare package Medium = 
        Modelica.Media.Water.StandardWaterOnePhase,
    p=800000) 
  annotation (Placement(transformation(extent={{-100,-50},{-80,-30}}, rotation=
            0)));
  Modelica_Fluid.Sources.FixedBoundary_pTX SinkP1(
    T=system.T_ambient,
    redeclare package Medium = 
        Modelica.Media.Water.StandardWaterOnePhase,
    p=100000) 
  annotation (Placement(transformation(extent={{82,-4},{62,16}}, rotation=0)));
  Modelica_Fluid.ControlValves.ValveIncompressible V1(
    dp_nominal=9e5,
    m_flow_nominal=1.5,
  redeclare package Medium = Medium,
    CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
    Av=0.1) annotation (Placement(transformation(extent={{-50,58},{-30,78}},
          rotation=0)));
  Modelica_Fluid.ControlValves.ValveIncompressible V2(
    dp_nominal=5e5,
    m_flow_nominal=1.2,
  redeclare package Medium = Medium,
    CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
    Av=0.1) annotation (Placement(transformation(extent={{-38,26},{-18,46}},
          rotation=0)));
  Modelica_Fluid.ControlValves.ValveIncompressible V3(
    dp_nominal=3e5,
    m_flow_nominal=1.1,
  redeclare package Medium = Medium,
    CvData=Modelica_Fluid.Types.CvTypes.OpPoint) 
            annotation (Placement(transformation(extent={{-38,-38},{-18,-18}},
          rotation=0)));
  Modelica_Fluid.ControlValves.ValveIncompressible V4(
    dp_nominal=8e5,
    m_flow_nominal=1.3,
  redeclare package Medium = Medium,
    CvData=Modelica_Fluid.Types.CvTypes.OpPoint) 
            annotation (Placement(transformation(extent={{-40,-78},{-20,-58}},
          rotation=0)));
  Modelica_Fluid.ControlValves.ValveIncompressible V5(
    dp_nominal=4e5,
    m_flow_nominal=2,
  redeclare package Medium = Medium,
    CvData=Modelica_Fluid.Types.CvTypes.OpPoint) 
            annotation (Placement(transformation(extent={{30,-4},{50,16}},
          rotation=0)));

annotation (
  Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
            100}}),
          graphics),
  experiment(StopTime=4, Tolerance=1e-006),
  Documentation(info=""));
  Modelica_Fluid.Sources.FixedBoundary_pTX SinkP2(
    T=system.T_ambient,
    redeclare package Medium = 
        Modelica.Media.Water.StandardWaterOnePhase,
    p=100000) 
  annotation (Placement(transformation(extent={{4,58},{-16,78}}, rotation=0)));
  Modelica_Fluid.Sources.FixedBoundary_pTX SinkP3(
    T=system.T_ambient,
    redeclare package Medium = 
        Modelica.Media.Water.StandardWaterOnePhase,
    p=100000) 
  annotation (Placement(transformation(extent={{26,-78},{6,-58}}, rotation=0)));
  Modelica.Blocks.Sources.Ramp CloseLoad(
    duration=1,
    offset=1,
    startTime=1,
    height=-0.99) annotation (Placement(transformation(extent={{8,26},{28,46}},
          rotation=0)));
  Modelica.Blocks.Sources.Ramp OpenRelief(
    duration=2,
    height=1,
    offset=0,
    startTime=1) 
              annotation (Placement(transformation(extent={{-92,70},{-72,90}},
          rotation=0)));
  Modelica.Blocks.Sources.Ramp CloseValves(
    duration=2,
    offset=1,
    startTime=1,
    height=-1) 
              annotation (Placement(transformation(extent={{-96,-12},{-76,8}},
          rotation=0)));

  inner Modelica_Fluid.System system 
                                   annotation (Placement(transformation(extent=
            {{60,68},{80,88}}, rotation=0)));
equation
  connect(V1.port_b, SinkP2.ports[1]) 
                                  annotation (Line(points={{-30,68},{-16,68}}));
  connect(V4.port_b, SinkP3.ports[1]) 
                                  annotation (Line(points={{-20,-68},{6,-68}}));
  connect(SourceP1.ports[1], V1.port_a) 
                                   annotation (Line(points={{-80,36},{-68,36},{
          -68,68},{-50,68}}, color={0,127,255}));
  connect(SourceP1.ports[1], V2.port_a) 
                                   annotation (Line(points={{-80,36},{-38,36}},
        color={0,127,255}));
  connect(V2.port_b, V5.port_a) 
                             annotation (Line(points={{-18,36},{5,36},{5,6},{30,
          6}}, color={0,127,255}));
  connect(V3.port_b, V5.port_a) 
                             annotation (Line(points={{-18,-28},{6,-28},{6,6},{
          30,6}}, color={0,127,255}));
  connect(SourceP2.ports[1], V4.port_a) 
                                   annotation (Line(points={{-80,-40},{-60,-40},
          {-60,-68},{-40,-68}}, color={0,127,255}));
  connect(SourceP2.ports[1], V3.port_a) 
                                   annotation (Line(points={{-80,-40},{-60,-40},
          {-60,-28},{-38,-28}}, color={0,127,255}));
  connect(OpenRelief.y, V1.stemPosition) annotation (Line(points={{-71,80},{-40,
          80},{-40,77}}, color={0,0,127}));
  connect(OpenRelief.y, V4.stemPosition) annotation (Line(points={{-71,80},{-64,
          80},{-64,-52},{-30,-52},{-30,-59}}, color={0,0,127}));
  connect(CloseValves.y, V2.stemPosition) annotation (Line(points={{-75,-2},{
          -46,-2},{-46,54},{-28,54},{-28,45}}, color={0,0,127}));
  connect(CloseValves.y, V3.stemPosition) annotation (Line(points={{-75,-2},{
          -28,-2},{-28,-19}}, color={0,0,127}));
  connect(CloseLoad.y, V5.stemPosition) annotation (Line(points={{29,36},{40,36},
          {40,15}}, color={0,0,127}));
  connect(V5.port_b, SinkP1.ports[1]) 
    annotation (Line(points={{50,6},{62,6}}, color={0,127,255}));
end TestValvesReverse;
