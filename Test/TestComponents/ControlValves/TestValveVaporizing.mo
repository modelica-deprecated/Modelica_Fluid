within Modelica_Fluid.Test.TestComponents.ControlValves;
model TestValveVaporizing "Test case for vaporizing fluid valve"
  extends Modelica.Icons.Example;
  Modelica_Fluid.Sources.FixedBoundary_pTX SourceP1(
  redeclare package Medium = Modelica.Media.Water.StandardWater,
    p=200000,
    T(displayUnit="degC") = 373.15) 
  annotation (Placement(transformation(extent={{-76,-22},{-56,-2}}, rotation=0)));
  Modelica_Fluid.ControlValves.ValveVaporizing V1(
  redeclare package Medium = Modelica.Media.Water.StandardWater,
    m_flow_nominal=1,
    CvData=Modelica_Fluid.Types.CvTypes.Av,
    Av=240e-6,
    dp_nominal=100000,
    d_nominal=1000) 
            annotation (Placement(transformation(extent={{-26,-22},{-6,
            -2}},
          rotation=0)));
  Modelica_Fluid.Sources.PrescribedBoundary_pTX SinkP1(
  redeclare package Medium = Modelica.Media.Water.StandardWater,
    usePressureInput=true,
    p=100000,
    T=293.15) 
  annotation (Placement(transformation(extent={{46,-22},{26,-2}},
                                                                rotation=0)));
  Modelica.Blocks.Sources.Constant Opening(k=1) 
              annotation (Placement(transformation(extent={{-78,22},{
            -58,42}},
          rotation=0)));
  inner Modelica_Fluid.System system 
                                   annotation (Placement(transformation(extent={{64,60},
            {84,80}},          rotation=0)));
  Modelica.Blocks.Sources.Ramp downstreamPressure(
    offset=2e5,
    duration=0.6,
    startTime=0.2,
    height=-1.5e5) 
    annotation (Placement(transformation(extent={{14,20},{34,40}})));
equation

  connect(V1.port_b,SinkP1. port) annotation (Line(points={{-6,-12},{26,
          -12}}));
  connect(Opening.y,V1. stemPosition) 
  annotation (Line(points={{-57,32},{-16,32},{-16,-3}}, color={0,0,255}));
  connect(SourceP1.port,V1. port_a) 
                                   annotation (Line(points={{-56,-12},{
          -26,-12}},
        color={0,127,255}));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,
            -100},{100,100}}), graphics), experiment);
  connect(downstreamPressure.y, SinkP1.p_in) annotation (Line(
      points={{35,30},{60,30},{60,-6},{48,-6}},
      color={0,0,127},
      smooth=Smooth.None));
end TestValveVaporizing;
