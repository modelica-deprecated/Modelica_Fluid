within Modelica_Fluid.Test.TestComponents.ControlValves;
model TestValveCoefficients
  "Test case to compare different choices of flow coefficient"
  extends Modelica.Icons.Example;
  Modelica_Fluid.Sources.FixedBoundary_pTX SourceP1(
  redeclare package Medium = Modelica.Media.Water.StandardWater,
    T=system.T_ambient,
    p=200000) 
  annotation (Placement(transformation(extent={{-94,18},{-74,38}},  rotation=0)));
  Modelica_Fluid.ControlValves.ValveIncompressible V1(
    d_nominal=1000,
  redeclare package Medium = Modelica.Media.Water.StandardWater,
    m_flow_nominal=1,
    CvData=Modelica_Fluid.Types.CvTypes.Av,
    Av=240e-6,
    dp_nominal=100000) 
            annotation (Placement(transformation(extent={{-44,18},{-24,38}},
          rotation=0)));
  Modelica_Fluid.Sources.FixedBoundary_pTX SinkP1(
  redeclare package Medium = Modelica.Media.Water.StandardWater,
    T=system.T_ambient,
    p=100000) 
  annotation (Placement(transformation(extent={{28,18},{8,38}}, rotation=0)));
  Modelica.Blocks.Sources.Constant Opening(k=1) 
              annotation (Placement(transformation(extent={{-96,62},{-76,82}},
          rotation=0)));
  inner Modelica_Fluid.System system 
                                   annotation (Placement(transformation(extent={{64,60},
            {84,80}},          rotation=0)));
  Modelica_Fluid.Sources.FixedBoundary_pTX SourceP2(
  redeclare package Medium = Modelica.Media.Water.StandardWater,
    T=system.T_ambient,
    p=200000) 
  annotation (Placement(transformation(extent={{-94,-18},{-74,2}},  rotation=0)));
  Modelica_Fluid.ControlValves.ValveIncompressible V2(
    d_nominal=1000,
  redeclare package Medium = Modelica.Media.Water.StandardWater,
    CvData=Modelica_Fluid.Types.CvTypes.Kv,
    m_flow_nominal=1,
    Kv=9,
    dp_nominal=100000) 
            annotation (Placement(transformation(extent={{-44,-18},{-24,2}},
          rotation=0)));
  Modelica_Fluid.Sources.FixedBoundary_pTX SinkP2(
  redeclare package Medium = Modelica.Media.Water.StandardWater,
    T=system.T_ambient,
    p=100000) 
  annotation (Placement(transformation(extent={{28,-18},{8,2}}, rotation=0)));
  Modelica_Fluid.Sources.FixedBoundary_pTX SourceP3(
  redeclare package Medium = Modelica.Media.Water.StandardWater,
    T=system.T_ambient,
    p=200000) 
  annotation (Placement(transformation(extent={{-94,-56},{-74,-36}},rotation=0)));
  Modelica_Fluid.ControlValves.ValveIncompressible V3(
    d_nominal=1000,
  redeclare package Medium = Modelica.Media.Water.StandardWater,
    CvData=Modelica_Fluid.Types.CvTypes.Cv,
    Cv=10,
    m_flow_nominal=1,
    dp_nominal=100000) 
            annotation (Placement(transformation(extent={{-44,-56},{-24,-36}},
          rotation=0)));
  Modelica_Fluid.Sources.FixedBoundary_pTX SinkP3(
  redeclare package Medium = Modelica.Media.Water.StandardWater,
    T=system.T_ambient,
    p=100000) 
  annotation (Placement(transformation(extent={{28,-56},{8,-36}},
                                                                rotation=0)));
  Modelica_Fluid.Sources.FixedBoundary_pTX SourceP4(
  redeclare package Medium = Modelica.Media.Water.StandardWater,
    T=system.T_ambient,
    p=200000) 
  annotation (Placement(transformation(extent={{-94,-88},{-74,-68}},rotation=0)));
  Modelica_Fluid.ControlValves.ValveIncompressible V4(
    d_nominal=1000,
  redeclare package Medium = Modelica.Media.Water.StandardWater,
    Cv=10,
    CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
    dp_nominal=100000,
    m_flow_nominal=2.4) 
            annotation (Placement(transformation(extent={{-44,-88},{-24,-68}},
          rotation=0)));
  Modelica_Fluid.Sources.FixedBoundary_pTX SinkP4(
  redeclare package Medium = Modelica.Media.Water.StandardWater,
    T=system.T_ambient,
    p=100000) 
  annotation (Placement(transformation(extent={{28,-88},{8,-68}},
                                                                rotation=0)));
equation

  connect(V1.port_b,SinkP1. port) annotation (Line(points={{-24,28},{8,28}}));
  connect(Opening.y,V1. stemPosition) 
  annotation (Line(points={{-75,72},{-34,72},{-34,37}}, color={0,0,255}));
  connect(SourceP1.port,V1. port_a) 
                                   annotation (Line(points={{-74,28},{-44,28}},
        color={0,127,255}));
  connect(V2.port_b,SinkP2. port) annotation (Line(points={{-24,-8},{8,-8}}));
  connect(SourceP2.port,V2. port_a) 
                                   annotation (Line(points={{-74,-8},{-44,-8}},
        color={0,127,255}));
  connect(V3.port_b,SinkP3. port) annotation (Line(points={{-24,-46},{8,-46}}));
  connect(SourceP3.port,V3. port_a) 
                                   annotation (Line(points={{-74,-46},{-44,-46}},
        color={0,127,255}));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,
            -100},{100,100}}), graphics));
  connect(V2.stemPosition, Opening.y) annotation (Line(
      points={{-34,1},{-34,8},{-62,8},{-62,72},{-75,72}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(V3.stemPosition, Opening.y) annotation (Line(
      points={{-34,-37},{-34,-26},{-62,-26},{-62,72},{-75,72}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(V4.port_b,SinkP4. port) annotation (Line(points={{-24,-78},{8,-78}}));
  connect(SourceP4.port,V4. port_a) 
                                   annotation (Line(points={{-74,-78},{-44,-78}},
        color={0,127,255}));
  connect(V4.stemPosition, Opening.y) annotation (Line(
      points={{-34,-69},{-34,-58},{-62,-58},{-62,72},{-75,72}},
      color={0,0,127},
      smooth=Smooth.None));
end TestValveCoefficients;
