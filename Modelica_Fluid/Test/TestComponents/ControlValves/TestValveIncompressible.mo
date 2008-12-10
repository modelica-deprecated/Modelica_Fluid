within Modelica_Fluid.Test.TestComponents.ControlValves;
model TestValveIncompressible "Test case with different valve characteristics"
  extends Modelica.Icons.Example;
  package Medium = Modelica.Media.Water.StandardWater;
  Modelica_Fluid.Sources.FixedBoundary_pTX SourceP1(
                                               p=10e5,
  redeclare package Medium = Modelica.Media.Water.StandardWater,
    T=system.T_ambient) 
  annotation (Placement(transformation(extent={{-100,30},{-80,50}}, rotation=0)));
  Modelica_Fluid.Valves.ValveIncompressible V1(
    d_nominal=1000,
    dp_nominal=9e5,
    m_flow_nominal=1.5,
  redeclare package Medium = Modelica.Media.Water.StandardWater,
    CvData=Modelica_Fluid.Types.CvTypes.Cv,
    Cv=10)  annotation (Placement(transformation(extent={{-50,30},{-30,50}},
          rotation=0)));

annotation (
  Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
            100}}),
          graphics),
  experiment(StopTime=4, Tolerance=1e-006),
  Documentation(info=""));
  Modelica_Fluid.Sources.FixedBoundary_pTX SinkP1(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    T=system.T_ambient,
    p=100000) 
  annotation (Placement(transformation(extent={{22,30},{2,50}}, rotation=0)));
  Modelica.Blocks.Sources.Ramp Opening(
    duration=2,
    height=1,
    offset=0,
    startTime=1) 
              annotation (Placement(transformation(extent={{-92,74},{-72,94}},
          rotation=0)));
  Modelica_Fluid.Sources.FixedBoundary_pTX SourceP2(
                                               p=10e5,
  redeclare package Medium = Modelica.Media.Water.StandardWater,
    T=system.T_ambient) 
  annotation (Placement(transformation(extent={{-100,-10},{-80,10}}, rotation=0)));
  Modelica_Fluid.Valves.ValveIncompressible V2(
    d_nominal=1000,
    dp_nominal=9e5,
    m_flow_nominal=1.5,
  redeclare package Medium = Modelica.Media.Water.StandardWater,
    CvData=Modelica_Fluid.Types.CvTypes.Cv,
    Cv=10,
    redeclare function valveCharacteristic = 
        Modelica_Fluid.Valves.BaseClasses.ValveCharacteristics.equalPercentage)
            annotation (Placement(transformation(extent={{-50,-10},{-30,10}},
          rotation=0)));
  Modelica_Fluid.Sources.FixedBoundary_pTX SinkP2(
  redeclare package Medium = Modelica.Media.Water.StandardWater,
    T=system.T_ambient,
    p=100000) 
  annotation (Placement(transformation(extent={{22,-10},{2,10}}, rotation=0)));
  Modelica_Fluid.Sources.FixedBoundary_pTX SourceP3(
                                               p=10e5,
  redeclare package Medium = Modelica.Media.Water.StandardWater,
    T=system.T_ambient) 
  annotation (Placement(transformation(extent={{-100,-50},{-80,-30}}, rotation=
            0)));
  Modelica_Fluid.Valves.ValveIncompressible V3(
    d_nominal=1000,
    dp_nominal=9e5,
    m_flow_nominal=1.5,
  redeclare package Medium = Modelica.Media.Water.StandardWater,
    CvData=Modelica_Fluid.Types.CvTypes.Cv,
    Cv=10,
    redeclare function valveCharacteristic = 
        Modelica_Fluid.Valves.BaseClasses.ValveCharacteristics.equalPercentage
        (                                                                              rangeability=10)) 
            annotation (Placement(transformation(extent={{-50,-50},{-30,-30}},
          rotation=0)));
  Modelica_Fluid.Sources.FixedBoundary_pTX SinkP3(
                                             p=1e5,
  redeclare package Medium = Modelica.Media.Water.StandardWater,
    T=system.T_ambient) 
  annotation (Placement(transformation(extent={{22,-50},{2,-30}}, rotation=0)));

  inner Modelica_Fluid.System system 
                                   annotation (Placement(transformation(extent=
            {{58,72},{78,92}}, rotation=0)));
equation
  connect(V1.port_b, SinkP1.ports[1]) 
                                  annotation (Line(points={{-30,40},{2,40}}));
  connect(Opening.y, V1.opening) 
  annotation (Line(points={{-71,84},{-40,84},{-40,48}}, color={0,0,255}));
  connect(SourceP1.ports[1], V1.port_a) 
                                   annotation (Line(points={{-80,40},{-50,40}},
        color={0,127,255}));
  connect(Opening.y, V2.opening) annotation (Line(points={{-71,84},{-64,84},
          {-64,20},{-40,20},{-40,8}}, color={0,0,127}));
  connect(Opening.y, V3.opening) annotation (Line(points={{-71,84},{-64,84},
          {-64,-22},{-40,-22},{-40,-32}}, color={0,0,127}));
  connect(SourceP2.ports[1], V2.port_a) 
    annotation (Line(points={{-80,0},{-50,0}}, color={0,127,255}));
  connect(V2.port_b,SinkP2.ports[1]) 
                                  annotation (Line(points={{-30,0},{2,0}},
        color={0,127,255}));
  connect(SourceP3.ports[1], V3.port_a) 
                                    annotation (Line(points={{-80,-40},{-50,-40}},
        color={0,127,255}));
  connect(V3.port_b, SinkP3.ports[1]) 
                                  annotation (Line(points={{-30,-40},{2,-40}},
        color={0,127,255}));
end TestValveIncompressible;
